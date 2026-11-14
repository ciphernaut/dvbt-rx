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
// turn until one acquires. A line is printed when lock is gained or lost, when the
// signalled parameters change, and otherwise every `--every` seconds of signal.

#include <fftw3.h>

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

// One channel width's receiver: the resampler to its sample rate and the demodulator.
struct Chain
{
    int m_bandwidth;
    double m_rate;
    bool m_resampling;
    dvbt::Resampler m_resampler;
    dvbt::OFDMDemod m_demod;
    std::vector<dvbt::Cf> m_resampled;

    Chain(int bandwidth, double inputRate) :
        m_bandwidth(bandwidth),
        m_rate(dvbt::sampleRate(bandwidth)),
        m_resampling(std::abs(inputRate / m_rate - 1.0) > 1e-6),
        m_demod(bandwidth)
    {
        if (m_resampling) {
            m_resampler.configure(inputRate, m_rate, 0.4165 * m_rate); // Half the occupied bandwidth
        }

        m_demod.setLogCallback([bandwidth](const char *message) {
            if (verbose) {
                fprintf(stderr, "dvbt-rx: %d Hz: %s\n", bandwidth, message);
            }
        });
    }

    void process(const dvbt::Cf *samples, int count)
    {
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

void print(const Report& report, const dvbt::OFDMDemod::Status& status)
{
    if (report.m_locked)
    {
        printf("{\"locked\":true,\"lock\":\"locked\",\"bandwidth_hz\":%d,%s,\"offset_hz\":%ld,\"clock_ppm\":%.1f,\"snr_db\":%.1f}\n",
            report.m_bandwidth, report.m_parameters.c_str(), std::lround(status.m_frequencyOffset), status.m_clockOffset, status.m_snr);
    }
    else
    {
        printf("{\"locked\":false,\"lock\":\"searching\",\"bandwidth_hz\":null,\"mode\":null,\"guard\":null,\"constellation\":null,\"code_rate\":null,"
            "\"hierarchy\":null,\"lp_code_rate\":null,\"cell_id\":null,\"offset_hz\":null,\"clock_ppm\":null,\"snr_db\":null}\n");
    }

    fflush(stdout);
}

int usage()
{
    fprintf(stderr, "usage: dvbt-rx --rate <Hz> [--bandwidth auto|<Hz>] [--every <seconds>] [--verbose] < iq.cf32\n");
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
                print(report, status);
                last = report;
                printed = true;
                printedAt = seconds;
            }
        }
    }

    return 0;
}
