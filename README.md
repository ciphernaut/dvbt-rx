# dvbt-rx

DVB-T signalling from IQ samples, as JSON lines.

`dvbt-rx` reads complex float32 IQ on stdin, centred on a DVB-T multiplex, and prints what the
multiplex signals about itself (the transmission parameter signalling of EN 300 744): the channel
width it locked at, 2K or 8K carriers, the guard interval, the constellation, the code rate, the
hierarchy and the cell id when the multiplex sends one. It does not decode the transport stream:
no picture, sound or service names yet.

It is a small wrapper around the `modemdvb` library of [SDRangel](https://github.com/f4exb/sdrangel)
(Jon Beniston, M7RCE), which does the reception. It was written as the external decoder behind
the `dvbt_decode` mode row of [ZetaSDR](https://github.com/ciphernaut/ZetaSDR), which runs it as
a separate process; it works as well from a shell.

## Build

Needs `g++` (C++17), FFTW's single-precision library (Debian: `libfftw3-dev`) and a checkout of
SDRangel. Only the library's DVB-T reception sources are compiled; Qt is not needed (the few Qt
and SDRangel headers the library includes are replaced by the stand-ins under `shim/`).

```bash
git clone --depth 1 https://github.com/f4exb/sdrangel
./build.sh sdrangel            # writes ./dvbt-rx
sudo install dvbt-rx /usr/local/bin/
```

Tested against SDRangel commit `97964ee4` (2026-11-13) on x86-64. `modemdvb` is new and changes
quickly; if a later commit does not build, try that one.

## Use

```bash
dvbt-rx --rate 10000000 < multiplex.cf32
```

- `--rate HZ`: the input sample rate. Required. The library resamples to the rate each channel
  width needs, so any rate from about 0.8 to 2.2 times that works (10 MS/s covers all four widths).
- `--bandwidth auto|HZ`: the channel width. `auto` (the default) tries 8, 7, 6 and 5 MHz in turn,
  half a second of signal each, until one acquires.
- `--every SECONDS`: the time between lines while locked (default 2).
- `--no-filter`: skip the channel filter. By default the input is low-passed to the width being
  tried (cut at half the channel width, 80 dB down inside the gap to the next channel), so an input
  wider than the multiplex may carry its neighbours: unfiltered, a neighbour 20 dB or more above
  the wanted multiplex takes the receiver's timing. The filter costs about 7 % of a core at 10 MS/s.
- `--verbose`: the library's acquisition log on stderr.

A line is printed when lock is gained or lost, when the signalled parameters change, every
`--every` seconds of signal while locked and every 5 seconds while searching:

```json
{"locked":true,"lock":"locked","bandwidth_hz":7000000,"mode":"8K","guard":"1/16","constellation":"64-QAM","code_rate":"3/4","hierarchy":null,"lp_code_rate":null,"cell_id":null,"offset_hz":690,"clock_ppm":3.8,"snr_db":12.3,"ripple_db":16.6,"resyncs":0}
```

`offset_hz` is the multiplex's carrier offset from the centre of the input and `clock_ppm` the
sample clock error. `snr_db` is the mean error of the equalised pilot carriers against their known
values, unweighted across the multiplex: carriers in a fade count as much as the rest, so a strong
multiplex with deep multipath reads low, and below 0 dB while still locked. It is not the level
over the noise floor. `ripple_db` is the spread of the channel's frequency response across the
multiplex (flat is a few dB; multipath or a same-frequency network's second site makes it tens).
`resyncs` counts the timing jumps the receiver recovered from without losing lock since it last
acquired; each is usually samples lost on the way in. `cell_id` is `null` when
the multiplex does not signal one, and `hierarchy` and `lp_code_rate` are `null` for a
non-hierarchical multiplex. While nothing is locked every field but `locked` and `lock` is `null`.

## Licence

GPL-3.0-or-later, as the library it links ([LICENSE](LICENSE)).
