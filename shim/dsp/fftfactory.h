// SPDX-License-Identifier: GPL-3.0-or-later
// Part of dvbt-rx: a stand-in for the header of this name that
// SDRangel's modemdvb library includes, so the library builds without Qt.
#pragma once
#include "dsp/fftengine.h"
class FFTFactory {
public:
    unsigned int getEngine(unsigned int, bool, FFTEngine **, const QString& = QString()) { return 0; }
    void releaseEngine(unsigned int, bool, unsigned int) {}
};
