/*This header file includes the fuctions name and packages will be used in filters.cpp*/

#include <stdio.h>
#include <math.h>
#define Pi 3.141592653589793F
#include <malloc.h>
#include "image_comps.h"

// apply the 2D convolution directly
void directly_filter(my_image_comp* in, my_image_comp* out, int filter_length);

// 1D horizontal convolution
void horizontal_filter(my_image_comp* in, my_image_comp* out, int filter_length);

// 1D vertical convolution
void vertical_filter(my_image_comp* in, my_image_comp* out, int filter_length);

// Create a 1D sinc with hanning window
void sinc_windowing(int FILTER_EXTENT, int FILTER_DIM, float* mirror_psf);

// Create a 2D sinc kernel with hanning window
void direct_sinc_windowing(int FILTER_EXTENT, int FILTER_DIM, float* filter_buf);