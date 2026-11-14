// SPDX-License-Identifier: GPL-3.0-or-later
//
// dvbt-rx: DVB-T signalling from IQ samples.
//
// Reads complex float32 IQ (interleaved I/Q, little-endian) on stdin, centred on a
// DVB-T multiplex, and prints what the multiplex signals about itself (EN 300 744's
// transmission parameter signalling) as JSON lines on stdout. Reception is done by
// SDRangel's modemdvb library (GPL-3.0-or-later, Jon Beniston M7RCE), which this
// program links; it is therefore a separate program under the GPL, built by the
// operator (build.sh) and run by a ZetaSDR node as an external decoder.
//
//   dvbt-rx --rate <input samples per second> [--bandwidth auto|<Hz>] [--every <s>]
//
// With `--bandwidth auto` (the default) the channel widths in `widths` are tried in
// turn until one acquires. The input is first low-passed to the width being tried, so
// the neighbouring multiplexes a wider input carries do not reach the receiver. A line is printed when lock is gained or lost, when the
// signalled parameters change, and otherwise every `--every` seconds of signal.

#include <fftw3.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "dsp/fftengine.h"
#include "dvbt.h"
#include "dvbtofdm.h"
#include "resampler.h"

namespace {

class FftwEngine : public FFTEngine
{
public:
    ~FftwEngine() override { release(); }

    void configure(int n, bool inverse) override
    {
        release();
        m_in = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * n);
        m_out = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * n);
        m_plan = fftwf_plan_dft_1d(n, m_in, m_out, inverse ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_ESTIMATE);
    }

    void transform() override { fftwf_execute(m_plan); }
    Complex *in() override { return reinterpret_cast<Complex *>(m_in); }
    Complex *out() override { return reinterpret_cast<Complex *>(m_out); }

private:
    void release()
    {
        if (m_plan)
        {
            fftwf_destroy_plan(m_plan);
            fftwf_free(m_in);
            fftwf_free(m_out);
            m_plan = nullptr;
        }
    }

    fftwf_complex *m_in = nullptr;
    fftwf_complex *m_out = nullptr;
    fftwf_plan m_plan = nullptr;
};

// A low-pass that leaves only the wanted channel: the cut is at half the channel width, in the
// gap between this multiplex's outermost carriers and a neighbour's, and the transition (about
// 1 % of the input rate, 80 dB down beyond it) fits inside that gap. Without it a neighbour 20 dB
// or more above the wanted multiplex takes the receiver's timing. Filtering is done a block at a
// time in the frequency domain (overlap-save), which costs far less than the taps would directly.
class ChannelFilter
{
public:
    ChannelFilter(double rate, double cutoff)
    {
        m_time = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * m_size);
        m_freq = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * m_size);
        m_forward = fftwf_plan_dft_1d(m_size, m_time, m_freq, FFTW_FORWARD, FFTW_ESTIMATE);
        m_inverse = fftwf_plan_dft_1d(m_size, m_freq, m_time, FFTW_BACKWARD, FFTW_ESTIMATE);

        // A Kaiser-windowed sinc, then its spectrum (with the inverse transform's 1/size)
        const double pi = 3.14159265358979323846;
        const double beta = 7.857; // 80 dB
        const int centre = (m_taps - 1) / 2;
        const double width = 2.0 * cutoff / rate;
        std::vector<double> taps(m_taps);
        double sum = 0.0;

        for (int i = 0; i < m_taps; i++)
        {
            const double x = i - centre;
            const double sinc = x == 0.0 ? width : std::sin(pi * width * x) / (pi * x);
            const double r = x / centre;
            taps[i] = sinc * besselI0(beta * std::sqrt(1.0 - r * r)) / besselI0(beta);
            sum += taps[i];
        }

        for (int i = 0; i < m_size; i++)
        {
            m_time[i][0] = i < m_taps ? (float) (taps[i] / sum) : 0.0f;
            m_time[i][1] = 0.0f;
        }

        fftwf_execute(m_forward);
        m_response.resize(m_size);

        for (int i = 0; i < m_size; i++) {
            m_response[i] = dvbt::Cf(m_freq[i][0], m_freq[i][1]) / (float) m_size;
        }

        m_block.assign(m_size, dvbt::Cf(0.0f, 0.0f));
        m_fill = m_taps - 1;
    }

    ~ChannelFilter()
    {
        fftwf_destroy_plan(m_forward);
        fftwf_destroy_plan(m_inverse);
        fftwf_free(m_time);
        fftwf_free(m_freq);
    }

    ChannelFilter(const ChannelFilter&) = delete;
    ChannelFilter& operator=(const ChannelFilter&) = delete;

    //! Appends the filtered samples to out; they lag the input by up to a block
    void process(const dvbt::Cf *in, int count, std::vector<dvbt::Cf>& out)
    {
        const int overlap = m_taps - 1;

        while (count > 0)
        {
            const int n = std::min(count, m_size - m_fill);
            std::copy(in, in + n, m_block.begin() + m_fill);
            m_fill += n;
            in += n;
            count -= n;

            if (m_fill < m_size) {
                break;
            }

            memcpy(m_time, m_block.data(), sizeof(fftwf_complex) * m_size);
            fftwf_execute(m_forward);

            for (int i = 0; i < m_size; i++)
            {
                const dvbt::Cf v = dvbt::Cf(m_freq[i][0], m_freq[i][1]) * m_response[i];
                m_freq[i][0] = v.real();
                m_freq[i][1] = v.imag();
            }

            fftwf_execute(m_inverse);
            const dvbt::Cf *filtered = reinterpret_cast<const dvbt::Cf *>(m_time);
            out.insert(out.end(), filtered + overlap, filtered + m_size); // The first are wrapped round
            std::copy(m_block.end() - overlap, m_block.end(), m_block.begin());
            m_fill = overlap;
        }
    }

private:
    static double besselI0(double x)
    {
        double sum = 1.0;
        double term = 1.0;

        for (int k = 1; k < 40; k++)
        {
            term *= (x / (2.0 * k)) * (x / (2.0 * k));
            sum += term;
        }

        return sum;
    }

    static constexpr int m_taps = 1025;
    static constexpr int m_size = 8192;
    fftwf_complex *m_time;
    fftwf_complex *m_freq;
    fftwf_plan m_forward;
    fftwf_plan m_inverse;
    std::vector<dvbt::Cf> m_response;
    std::vector<dvbt::Cf> m_block;  //!< The last taps - 1 inputs, then the new ones
    int m_fill;
};

// The channel widths `--bandwidth auto` tries, Hz.
const int widths[] = {8000000, 7000000, 6000000, 5000000};
// Signal seconds a width is given to get past acquisition before the next is tried.
const double trySeconds = 0.5;
// Signal seconds between lines while nothing is locked.
const double searchEvery = 5.0;
// The library resamples between these ratios of input to receiver rate.
const double minRatio = 0.8;
const double maxRatio = 2.2;

bool verbose = false;
bool filtering = true;

// One channel width's receiver: the resampler to its sample rate and the demodulator.
struct Chain
{
    int m_bandwidth;
    double m_rate;
    bool m_resampling;
    dvbt::Resampler m_resampler;
    dvbt::OFDMDemod m_demod;
    std::vector<dvbt::Cf> m_resampled;
    std::unique_ptr<ChannelFilter> m_filter; //!< Null when the input is no wider than the channel
    std::vector<dvbt::Cf> m_filtered;

    Chain(int bandwidth, double inputRate) :
        m_bandwidth(bandwidth),
        m_rate(dvbt::sampleRate(bandwidth)),
        m_resampling(std::abs(inputRate / m_rate - 1.0) > 1e-6),
        m_demod(bandwidth)
    {
        if (m_resampling) {
            m_resampler.configure(inputRate, m_rate, 0.4165 * m_rate); // Half the occupied bandwidth
        }

        if (filtering && (bandwidth < 0.98 * inputRate)) {
            m_filter.reset(new ChannelFilter(inputRate, bandwidth / 2.0));
        }

        m_demod.setLogCallback([bandwidth](const char *message) {
            if (verbose) {
                fprintf(stderr, "dvbt-rx: %d Hz: %s\n", bandwidth, message);
            }
        });
    }

    void process(const dvbt::Cf *samples, int count)
    {
        if (m_filter)
        {
            m_filtered.clear();
            m_filter->process(samples, count, m_filtered);
            samples = m_filtered.data();
            count = (int) m_filtered.size();
        }

        if (m_resampling)
        {
            m_resampled.clear();
            m_resampler.process(samples, count, m_resampled);
            m_demod.process(m_resampled.data(), (int) m_resampled.size());
        }
        else
        {
            m_demod.process(samples, count);
        }
    }
};

// What a line reports apart from its measurements; a change prints a line at once.
struct Report
{
    bool m_locked = false;
    int m_bandwidth = 0;
    std::string m_parameters;

    bool operator!=(const Report& other) const
    {
        return (m_locked != other.m_locked) || (m_bandwidth != other.m_bandwidth) || (m_parameters != other.m_parameters);
    }
};

std::string quoted(const char *text) { return std::string("\"") + text + "\""; }

Report describe(const Chain& chain, const dvbt::OFDMDemod::Status& status)
{
    Report report;
    report.m_locked = (status.m_state == dvbt::OFDMDemod::Locked) && status.m_tpsValid;

    if (report.m_locked)
    {
        const dvbt::TPS& tps = status.m_tps;
        const bool hierarchical = tps.m_hierarchy != dvbt::Hierarchy::None;
        report.m_bandwidth = chain.m_bandwidth;
        report.m_parameters = "\"mode\":" + quoted(dvbt::toString(tps.m_mode))
            + ",\"guard\":" + quoted(dvbt::toString(tps.m_guard))
            + ",\"constellation\":" + quoted(dvbt::toString(tps.m_constellation))
            + ",\"code_rate\":" + quoted(dvbt::toString(tps.m_hpCodeRate))
            + ",\"hierarchy\":" + (hierarchical ? quoted(dvbt::toString(tps.m_hierarchy)) : std::string("null"))
            + ",\"lp_code_rate\":" + (hierarchical ? quoted(dvbt::toString(tps.m_lpCodeRate)) : std::string("null"))
            + ",\"cell_id\":" + (status.m_cellId >= 0 ? std::to_string(status.m_cellId) : std::string("null"));
    }

    return report;
}

// How uneven the channel is across the multiplex: the spread of its frequency response, dB,
// "null" until a symbol has been equalised. Multipath shows here before it shows anywhere else.
std::string rippleOf(const dvbt::OFDMDemod& demod)
{
    std::vector<dvbt::Cf> cells;
    std::vector<float> response;
    demod.getDisplayData(cells, 1, response, 64);

    if (response.size() < 2) {
        return "null";
    }

    const auto range = std::minmax_element(response.begin(), response.end());
    char text[32];
    snprintf(text, sizeof(text), "%.1f", *range.second - *range.first);
    return text;
}

void print(const Report& report, const dvbt::OFDMDemod::Status& status, const std::string& ripple)
{
    if (report.m_locked)
    {
        printf("{\"locked\":true,\"lock\":\"locked\",\"bandwidth_hz\":%d,%s,\"offset_hz\":%ld,\"clock_ppm\":%.1f,\"snr_db\":%.1f,\"ripple_db\":%s,\"resyncs\":%llu}\n",
            report.m_bandwidth, report.m_parameters.c_str(), std::lround(status.m_frequencyOffset), status.m_clockOffset, status.m_snr,
            ripple.c_str(), (unsigned long long) status.m_resyncs);
    }
    else
    {
        printf("{\"locked\":false,\"lock\":\"searching\",\"bandwidth_hz\":null,\"mode\":null,\"guard\":null,\"constellation\":null,\"code_rate\":null,"
            "\"hierarchy\":null,\"lp_code_rate\":null,\"cell_id\":null,\"offset_hz\":null,\"clock_ppm\":null,\"snr_db\":null,\"ripple_db\":null,\"resyncs\":null}\n");
    }

    fflush(stdout);
}

int usage()
{
    fprintf(stderr, "usage: dvbt-rx --rate <Hz> [--bandwidth auto|<Hz>] [--every <seconds>] [--no-filter] [--verbose] < iq.cf32\n");
    return 2;
}

} // namespace

FFTEngine *FFTEngine::create(const QString&, const QString&) { return new FftwEngine(); }

int main(int argc, char **argv)
{
    double rate = 0.0;
    double every = 2.0;
    int bandwidth = 0; // 0 is auto

    for (int i = 1; i < argc; i++)
    {
        const bool more = i + 1 < argc;

        if (!strcmp(argv[i], "--verbose")) {
            verbose = true;
        } else if (!strcmp(argv[i], "--no-filter")) {
            filtering = false;
        } else if (!strcmp(argv[i], "--rate") && more) {
            rate = atof(argv[++i]);
        } else if (!strcmp(argv[i], "--every") && more) {
            every = atof(argv[++i]);
        } else if (!strcmp(argv[i], "--bandwidth") && more) {
            i++;
            bandwidth = strcmp(argv[i], "auto") ? atoi(argv[i]) : 0;
        } else {
            return usage();
        }
    }

    if ((rate <= 0.0) || (every <= 0.0)) {
        return usage();
    }

    std::vector<int> candidates;

    for (int width : widths)
    {
        const double ratio = rate / dvbt::sampleRate(width);

        if (((bandwidth == 0) || (bandwidth == width)) && (ratio >= minRatio) && (ratio <= maxRatio)) {
            candidates.push_back(width);
        }
    }

    if (candidates.empty())
    {
        fprintf(stderr, "dvbt-rx: no channel width to try at %.0f samples/s (bandwidth %d)\n", rate, bandwidth);
        return 2;
    }

    size_t candidate = 0;
    std::unique_ptr<Chain> chain(new Chain(candidates[candidate], rate));
    std::vector<dvbt::Cf> input(32768);
    Report last;
    bool printed = false;
    double seconds = 0.0;       // Of signal read
    double printedAt = 0.0;
    double acquiringSince = 0.0;
    size_t count;

    while ((count = fread(input.data(), sizeof(dvbt::Cf), input.size(), stdin)) > 0)
    {
        chain->process(input.data(), (int) count);
        seconds += count / rate;
        dvbt::OFDMDemod::Status status = chain->m_demod.getStatus();

        if (status.m_state != dvbt::OFDMDemod::Acquire)
        {
            acquiringSince = seconds;
        }
        else if ((candidates.size() > 1) && (seconds - acquiringSince >= trySeconds))
        {
            candidate = (candidate + 1) % candidates.size();
            chain.reset(new Chain(candidates[candidate], rate));
            acquiringSince = seconds;
            status = chain->m_demod.getStatus();
        }

        const Report report = describe(*chain, status);
        const double period = report.m_locked ? every : searchEvery;

        if (!printed || (report != last) || (seconds - printedAt >= period))
        {
            // A first line that says nothing is locked waits for every width to have had its turn
            if (printed || report.m_locked || (seconds >= trySeconds * candidates.size() + 1.0))
            {
                print(report, status, report.m_locked ? rippleOf(chain->m_demod) : std::string("null"));
                last = report;
                printed = true;
                printedAt = seconds;
            }
        }
    }

    return 0;
}
