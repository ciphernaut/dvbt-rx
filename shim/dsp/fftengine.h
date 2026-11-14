// SPDX-License-Identifier: GPL-3.0-or-later
// Part of dvbt-rx: a stand-in for the header of this name that
// SDRangel's modemdvb library includes, so the library builds without Qt.
#pragma once
#include <complex>
#include <QString>
typedef std::complex<float> Complex;
class FFTEngine {
public:
    virtual ~FFTEngine() {}
    virtual void configure(int n, bool inverse) = 0;
    virtual void transform() = 0;
    virtual Complex* in() = 0;
    virtual Complex* out() = 0;
    static FFTEngine* create(const QString& wisdom, const QString& preferred = QString());
};
