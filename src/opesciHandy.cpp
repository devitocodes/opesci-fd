/*  Copyright (C) 2015 Imperial College London and others.
 *
 *  Please see the AUTHORS file in the main source directory for a
 *  full list of copyright holders.
 *
 *  Gerard Gorman
 *  Department of Earth Science and Engineering
 *  Imperial College London
 *
 *  g.gorman@imperial.ac.uk
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above
 *  copyright notice, this list of conditions and the following
 *  disclaimer in the documentation and/or other materials provided
 *  with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 *  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 *  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 *  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
 *  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 *  SUCH DAMAGE.
 */

#include "opesciHandy.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

// Abort execution.
void opesci_abort(std::string msg)
{
    std::cerr<<msg<<std::endl;
    exit(-1);
}

// Compute the field of Lame constants from subsurface model - i.e. the p-wave velocity, s-wave velocity and density fields.
// See http://scienceworld.wolfram.com/physics/LameConstants.html
void opesci_calculate_lame_costants(const std::vector<float> &vp, const std::vector<float> &vs, const std::vector<float> &rho,
                                    std::vector<float> &mu, std::vector<float> &lam)
{
    size_t size=vp.size();
    assert(size==vs.size());
    assert(size==rho.size());
    mu.resize(size);
    lam.resize(size);

    #pragma omp parallel for
    for(size_t i=0; i<size; i++) {
        mu[i] = rho[i]*vs[i]*vs[i];
        lam[i] = rho[i]*(vp[i]*vp[i]-2.0*vs[i]*vs[i]);
    }
}

float opesci_calculate_dt(const std::vector<float> &vp, float h)
{
    float maxv = 0;
    size_t size = vp.size();
#if defined _OPENMP && _OPENMP >= 200711
    #pragma omp parallel for reduction(max:maxv)
    for(size_t i=0; i<size; ++i) {
        if(vp[i]>maxv) {
            maxv = vp[i];
        }
    }
#else
    for(size_t i=0; i<size; ++i) {
        if(vp[i]>maxv) {
            maxv = vp[i];
        }
    }
#endif

    return (6.0/7.0)*(1./sqrt(3.0))*(h/maxv);
}

/*
 * Discrete Fourier transform
 * by Project Nayuki, 2014. Public domain.
 * http://www.nayuki.io/page/how-to-implement-the-discrete-fourier-transform
 *
 * Computes the discrete Fourier transform (DFT) of the given vector.
 * All the array arguments must have the same length.
 */

void opesci_dft(const float inreal[], const float inimag[], float outreal[], float outimag[], int n)
{
    #pragma omp parallel for
    for(int k=0; k<n; ++k) { /* For each output element */
        float sumreal = 0;
        float sumimag = 0;
        for(int t=0; t<n; ++t) { /* For each input element */
            float angle = 2 * M_PI * t * k / n;
            sumreal +=  inreal[t] * cos(angle) + inimag[t] * sin(angle);
            sumimag += -inreal[t] * sin(angle) + inimag[t] * cos(angle);
        }
        outreal[k] = sumreal;
        outimag[k] = sumimag;
    }
}

/*
 * If dt < sdt, the source is resampled by perfoming a discrete
 * fourier transfor, filling the resulting spectrum with as many
 * zeroes as needed and inverse fourier transforming to get the
 * desired sampling without changing the frequency content.
 *
 * If sdt < dt, the source is resampled as above but removing samples
 * of the spectrum until the desired sample rate is desired.
 */
void opesci_resample_timeseries(const std::vector<float> &src, float dt, double sdt, std::vector<float> &resampled_src)
{
    // Break out early if necessary.
    if(fabs(dt-sdt)<std::numeric_limits<float>::epsilon()*(dt+sdt)) {
        resampled_src = src;
        return;
    }

    // Calculate new source size.
    int snt = src.size();
    int snt2 = (int)round(snt*sdt/dt);

    std::vector<float> imsrc(snt, 0.0);

    std::vector<float> ft_resrc(snt);
    std::vector<float> ft_imsrc(snt);

    opesci_dft(src.data(), imsrc.data(), ft_resrc.data(), ft_imsrc.data(), snt);

    // Normalize
    float nrm = 1./sqrt(snt);
    for(auto &i : ft_resrc)
        i*=nrm;
    for(auto &i : ft_imsrc)
        i*=nrm;

    std::vector<float> ft_resrc2(snt2, 0.0);
    std::vector<float> ft_imsrc2(snt2, 0.0);

    if(dt<sdt) {
        // Add zeroes in the middle of the real and imaginary part of the
        // transform. This is the same as adding zeroes at the end of the
        // spectrum.
        int midpoint = snt/2;
        for(int i=0; i<midpoint; i++) {
            ft_resrc2[i] = ft_resrc[i];
            ft_imsrc2[i] = ft_imsrc[i];
        }
        int offset = (snt2-snt);
        for(int i=midpoint; i<snt; i++) {
            ft_resrc2[offset+i] = ft_resrc[i];
            ft_imsrc2[offset+i] = ft_imsrc[i];
        }
    } else if(dt>sdt) {
        // Substract the middle of the real and imaginary part of the transform.
        int midpoint = snt2/2;
        for(int i=0; i<midpoint; i++) {
            ft_resrc2[i] = ft_resrc[i];
            ft_imsrc2[i] = ft_imsrc[i];
        }
        int offset = (snt-snt2);
        for(int i=midpoint; i<snt2; i++) {
            ft_resrc2[i] = ft_resrc[offset+i];
            ft_imsrc2[i] = ft_imsrc[offset+i];
        }
    }

    // Inverse fourier transform.
    for(auto &i : ft_imsrc2)
        i*=-1;

    resampled_src.resize(snt2);
    std::vector<float> imsrc2(snt2);

    opesci_dft(ft_resrc2.data(), ft_imsrc2.data(), resampled_src.data(), imsrc2.data(), snt2);

    // Normalise
    float nrm2 = 1.0/sqrt((float)snt2);
    for(auto &i : resampled_src)
        i*=nrm2;
}
