#include <stdio.h>
#include <math.h>
#define Pi 3.141592653589793F
#include <malloc.h>
#include "image_comps.h"

// apply the 2D convolution directly
void directly_filter(my_image_comp* in, my_image_comp* out, int filter_length);

// 1D horizontal convolution
void vertical_filter(my_image_comp* in, my_image_comp* out, int filter_length);
void horizontal_filter(my_image_comp* in, my_image_comp* out, int filter_length);

// 1D vertical convolution
//void horizontal_filter(my_image_comp* in, my_image_comp* out, int filter_length);
