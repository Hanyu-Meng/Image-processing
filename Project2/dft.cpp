/*****************************************************************************/
// File: dft.cpp
// Author: David Taubman
// Last Revised: 28 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#define _USE_MATH_DEFINES  // Makes `M_PI' available as a constant.

#include <stdlib.h>
#include <math.h>
#include "dft.h"
#include <complex>

/*****************************************************************************/
/*                            my_direct_dft::init                            */
/*****************************************************************************/

void my_direct_dft::init(int N, bool is_forward)
{
    cleanup(); // Delete any pre-existing buffers.
    this->N = N;
    real_buf = new float[N];
    imag_buf = new float[N];
    real_trig = new double[N];
    imag_trig = new double[N];
    for (int n = 0; n < N; n++)
    {
        real_trig[n] = cos(n * 2.0 * M_PI / N);
        imag_trig[n] = sin(n * 2.0 * M_PI / N);
        if (is_forward)
            imag_trig[n] = -imag_trig[n];
    }
}

/*****************************************************************************/
/*                   my_direct_dft::perform_transform                        */
/*****************************************************************************/

void my_direct_dft::perform_transform(float* real, float* imag, int stride)
{
    // First copy the input values to the real and imaginary temporary buffers.
    int n, k;
    float* rp, * ip;
    for (rp = real, ip = imag, n = 0; n < N; n++, rp += stride, ip += stride)
    {
        real_buf[n] = *rp;  imag_buf[n] = *ip;
    }

    // Now compute each output coefficient in turn
    for (rp = real, ip = imag, k = 0; k < N; k++, rp += stride, ip += stride)
    {
        int index = 0; // This holds n*k mod N; it indexes the trig tables
        double real_sum = 0.0, imag_sum = 0.0;
        for (n = 0; n < N; n++, index += k)
        {
            if (index >= N)
                index -= N;
            real_sum += real_buf[n] * real_trig[index]
                - imag_buf[n] * imag_trig[index];
            imag_sum += real_buf[n] * imag_trig[index]
                + imag_buf[n] * real_trig[index];
        }
        *rp = (float)real_sum;
        *ip = (float)imag_sum;
    }
}

void my_direct_fft::init(int N, bool is_forward)
{
    cleanup(); // Delete any pre-existing buffers.
    this->N = N; // block size N
    real_trig = new double[N];
    imag_trig = new double[N];
    for (int n = 0; n < N; n++)
    {
        real_trig[n] = cos(n * 2.0 * M_PI / N);
        imag_trig[n] = sin(n * 2.0 * M_PI / N);
        if (is_forward)
            imag_trig[n] = -imag_trig[n];
    }
}

void my_direct_fft::FFT(float* x_real, float* x_imag, int N, int mode, int stride, float* res_real, float* res_imag)
{
    // _complex* X = new _complex[N];
     // First copy the input values to the real and imaginary temporary buffers.
    float* xo_real;
    float* xo_imag;
    float* xe_real;
    float* xe_imag;
    // original signal: first half: even; second half: odd;
    xe_real = res_real;
    xe_imag = res_imag;
    xo_real = res_real + (N / 2) * stride;
    xo_imag = res_imag + (N / 2) * stride;

    //xo = xe + (N / 2) * stride;

    if (N == 1) {
        return;
    }

    for (int n = 0; n < (N / 2); n ++) {
        xe_real[n * stride] = x_real[2 * n * stride];
        xe_imag[n * stride] = x_imag[2 * n * stride];
        xo_real[n * stride] = x_real[(2 * n + 1) * stride];
        xo_imag[n * stride] = x_imag[(2 * n + 1) * stride];
    }

    //_complex* Xe = new _complex[N/2 * N/2];
    //_complex* Xo = new _complex[N / 2 * N / 2];
    this->FFT(xo_real, xo_imag, N / 2, mode, stride, x_real, x_imag);
    this->FFT(xe_real, xe_imag, N / 2, mode, stride, x_real, x_imag);
    this->init(N, mode);
    //int k, i;
    for (int i = 0, k = 0; k < (N / 2) * stride; k += stride, i++) {
        x_real[k] = xo_real[k] * real_trig[i] - imag_trig[i] * xo_imag[k] + xe_real[k];
        x_imag[k] = xo_imag[k] * real_trig[i] + xo_real[k] * imag_trig[i] + xe_imag[k];
        //  x_imag[k] = Xe[k].y + tmp.y;
        x_real[k + N / 2 * stride] = xe_real[k] - xo_real[k] * real_trig[i] + imag_trig[i] * xo_imag[k];
        x_imag[k + N / 2 * stride] = xe_imag[k] - xo_imag[k] * real_trig[i] - xo_real[k] * imag_trig[i];
        //   res[k + N / 2 * stride].x = Xe[k].x - tmp.x;
        //   res[k + N / 2 * stride].y = Xe[k].y - tmp.y;
    }
    //  return X;
}