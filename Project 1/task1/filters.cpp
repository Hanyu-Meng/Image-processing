#include "filters.h"


void sinc_windowing(int FILTER_EXTENT, int FILTER_DIM, float* mirror_psf) {

	float tau = FILTER_EXTENT; // Half length of the hanning window
	// 'mirror_psf' points to the central tap in the filter

	// Build windowed sinc filter separably
	for (int m = -FILTER_EXTENT; m <= FILTER_EXTENT; m++)
	{
		// Special case: the 1-D filter only has one tap
		if (m == 0)
			mirror_psf[m] = (1.0F / 3.0F);
		else
			mirror_psf[m] = (1.0F / 3.0F) * (sinf((1.0F / 3.0F) * Pi * m) / ((1.0F / 3.0F) * Pi * m)) * (0.5F) * (1 + cosf(Pi * m / tau));
	}


	/*normalisation*/

	float gain0 = 0.0F;
	for (int i = -FILTER_EXTENT; i < FILTER_EXTENT; i++)
	{
		gain0 += mirror_psf[i];
	}

	for (int i = -FILTER_EXTENT; i < FILTER_EXTENT; i++)
	{
		mirror_psf[i] = mirror_psf[i]/gain0;
	}

}

void direct_sinc_windowing(int FILTER_EXTENT, int FILTER_DIM, float* filter_buf) {
	// Create a buffer of the desired filter
	float* filter_handle = new float[FILTER_DIM];
	float* filter_tmp = filter_handle + FILTER_EXTENT; // points to the center
	float tau = FILTER_EXTENT;

	for (int m = -FILTER_EXTENT; m <= FILTER_EXTENT; m++) {
		if (m == 0)
			filter_tmp[m] = 1.0f / 3.0f;
		else
			filter_tmp[m] = (1.0F / 3.0F) * (sinf((1.0F / 3.0F) * Pi * m) / ((1.0F / 3.0F) * Pi * m)) * (0.5F) * (1 + cosf(Pi * m / tau));
	}

	/*dc gain normalisation*/
	float gain = 0.0f;
	for (int i = 0; i < FILTER_DIM; i++) {
		gain += filter_handle[i];
	}

	for (int i = 0; i < FILTER_DIM; i++) {
		filter_handle[i] /= gain;
	}

	// create a 2D sinc kernel
	for (int r = 0; r < FILTER_DIM; r++) {
		for (int c = 0; c < FILTER_DIM; c++) {
			filter_buf[r * FILTER_DIM + c] = filter_handle[c] * filter_handle[r];
		}
	}

}

void directly_filter (my_image_comp *in, my_image_comp *out, int filter_length) {
	int FILTER_EXTENT = filter_length;
	int FILTER_DIM = (2 * FILTER_EXTENT + 1);
	int FILTER_TAPS = FILTER_DIM * FILTER_DIM;
	float* filter_buf = new float[FILTER_TAPS];  // No shift sinc kernel
	// Assign the sinc kernel to the filter_buf
	direct_sinc_windowing(FILTER_EXTENT, FILTER_DIM, filter_buf);
	float* mirror_psf = filter_buf + FILTER_EXTENT * FILTER_DIM + FILTER_EXTENT;

	// r1->input row index; r2->output row index;
	// c1->input column index; c2->output column index;

	for (int r1 = 0, r2 = 0; (r1 < in->height) && (r2 < out->height); r1 += 3, r2 += 1) {
		for (int c1 = 0, c2 = 0; (c1 < in->width) && (c2 < out->width); c1 += 3, c2 += 1) {
			float* ip = in->buf + r1 * in->stride + c1;
			float* op = out->buf + r2 * out->stride + c2;
			float sum = 0.0F; 
			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) 
				for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) 
					sum += ip[y * in->stride + x] * mirror_psf[y * FILTER_DIM + x]; // dot product  
			*op = sum;
		}
	}
}

void vertical_filter(my_image_comp *in, my_image_comp *out, int filter_length) {
	int FILTER_EXTENT = filter_length;
	int FILTER_DIM = (2 * FILTER_EXTENT + 1);
	float* filter_buf = new float[FILTER_DIM];  // No shift sinc kernel
	float* mirror_psf = filter_buf + FILTER_EXTENT;
	// 'mirror_psf' points to the central tap in the filter
	sinc_windowing(FILTER_EXTENT, FILTER_DIM, mirror_psf); //windowing

	for (int r1 = 0, r2 = 0; (r1 < in->height) && (r2 < out->height); r1 += 3, r2 += 1) {
		for (int c1 = 0, c2 = 0; c2 < out->width; c1 += 3, c2 += 1)
		{
			float* ip = in->buf + r1 * in->stride + c2;
			float* op = out->buf + r2 * out->stride + c2;
			float sum = 0.0F;
			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
				sum += ip[y * in->stride] * mirror_psf[y]; // dot product
			}
			*op = sum;
		}
	}
}

void horizontal_filter(my_image_comp* in, my_image_comp* out, int filter_length) {
	int FILTER_EXTENT = filter_length;
	int FILTER_DIM = (2 * FILTER_EXTENT + 1); 
	float* filter_buf = new float[FILTER_DIM]; 
	float* mirror_psf = filter_buf + FILTER_EXTENT;
	sinc_windowing(FILTER_EXTENT, FILTER_DIM, mirror_psf);
	for (int r1 = 0, r2 = 0; r2 < out->height; r1 += 3, r2 += 1) {
		for (int c1 = 0, c2 = 0; (c1 < in->width) && (c2 < out->width); c1 += 3, c2 += 1) {
			float* ip = in->buf + r2 * in->stride + c1;
			float* op = out->buf + r2 * out->stride + c2;
			float sum = 0.0F;
			for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++)
				sum += ip[x] * mirror_psf[x]; // dot product
			*op = sum;
		}
	}

}