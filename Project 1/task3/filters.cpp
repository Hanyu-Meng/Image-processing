#include "filters.h"



void sinc_windowing(int FILTER_EXTENT, int FILTER_DIM, float* mirror_psf, float shift) {

	int tau = FILTER_EXTENT; // Half length of the hanning window
	// Build windowed sinc filter separably
	for (int m = -FILTER_EXTENT; m <= FILTER_EXTENT; m++)
	{
		// Special case: the 1-D filter only has one tap
		if (m == shift)
			mirror_psf[m] = 1.0F * (0.5F) * (1 + cosf(Pi * (1.0F*m - shift) / tau));
		else
			mirror_psf[m] = (sinf(Pi * (1.0F*m-shift)) / (Pi * (1.0F*m-shift))) * (0.5F) * (1 + cosf(Pi * (1.0F*m) / tau));
	}

	/*normalisation*/

	float gain0 = 0.0F;
	for (int i = -FILTER_EXTENT; i < FILTER_EXTENT; i++)
	{
		gain0 += mirror_psf[i];
	}

	for (int i = -FILTER_EXTENT; i < FILTER_EXTENT; i++)
	{
		mirror_psf[i] = mirror_psf[i] / gain0;
	}
}

void direct_sinc_windowing(int FILTER_EXTENT, int FILTER_DIM, float* filter_buf, float shift) {

	float* filter_handle = new float[FILTER_DIM];
	float* filter_tmp = filter_handle + FILTER_EXTENT; // points to the center
	int tau = FILTER_EXTENT;

	for (int m = -FILTER_EXTENT; m <= FILTER_EXTENT; m++) {
		if (m == 0)
			filter_tmp[m] = 1;
		else
			filter_tmp[m] = (sinf(Pi * (m - shift)) / (Pi * (m - shift))) * (0.5F) * (1 + cosf(Pi * (m - shift) / tau));
	}

	/*normalisation*/
	float gain = 0.0F;
	for (int i = 0; i < FILTER_DIM; i++) {
		gain += filter_handle[i];
	}

	for (int c = 0; c < FILTER_DIM; c++) {
		filter_handle[c] /= gain;
	}

	for (int r = 0; r < FILTER_DIM; r++) {
		for (int c = 0; c < FILTER_DIM; c++) {
			filter_buf[r * FILTER_DIM + c] = filter_handle[c] * filter_handle[r];
		}
	}

}

void directly_filter(my_image_comp* in, my_image_comp* out, int filter_length) {

	int FILTER_EXTENT = filter_length;
	int FILTER_DIM = (2 * FILTER_EXTENT + 1);
	int FILTER_TAPS = FILTER_DIM * FILTER_DIM;
	float* filter_buf1 = new float[FILTER_TAPS];
	float* filter_buf2 = new float[FILTER_TAPS];
	float* filter_buf3 = new float[FILTER_TAPS];
	//float* mirror_psf = filter_buf + FILTER_EXTENT*FILTER_DIM + FILTER_EXTENT;
	direct_sinc_windowing(FILTER_EXTENT, FILTER_DIM, filter_buf1, 0.0F);
	direct_sinc_windowing(FILTER_EXTENT, FILTER_DIM, filter_buf2, 1.0F/3.0F);
	direct_sinc_windowing(FILTER_EXTENT, FILTER_DIM, filter_buf2, 2.0F/3.0F);
	float* mirror_psf1 = filter_buf1 + FILTER_EXTENT * FILTER_DIM + FILTER_EXTENT;
	float* mirror_psf2 = filter_buf2 + FILTER_EXTENT * FILTER_DIM + FILTER_EXTENT;
	float* mirror_psf3 = filter_buf3 + FILTER_EXTENT * FILTER_DIM + FILTER_EXTENT;

	// r1->input,row; r2->output,row;
	// c1->input,column; c2->output,column;

	for (int r1 = 0, r2 = 0; (r1 < in->height) && (r2 < out->height); r1 += 1, r2 += 3) {
		for (int c1 = 0, c2 = 0; (c1 < in->width) && (c2 < out->width); c1 += 1, c2 += 3) {
			float* ip = in->buf + r1 * in->stride + c1;
			float* op = out->buf + r2 * out->stride + c2;
			float sum = 0.0F;
			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
				for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) {
					sum += ip[y * in->stride + x] * mirror_psf1[y * FILTER_DIM + x]; // dot product
				}
			}

			op[0] = sum;
			sum = 0.0F;
			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
				for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) {
					sum += ip[y * in->stride + x] * mirror_psf2[y * FILTER_DIM + x]; // dot product
				}
			}
			op[1] = sum;
			sum = 0.0F;
			
			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
				for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) {
					sum += ip[y * in->stride + x] * mirror_psf3[y * FILTER_DIM + x]; // dot product
				}
			}
			op[2] = sum;
			
		}
	}
}

void vertical_filter(my_image_comp* in, my_image_comp* out, int filter_length){
	int FILTER_EXTENT = filter_length;
	int FILTER_DIM = (2 * FILTER_EXTENT + 1);
	float* filter_buf1 = new float[FILTER_DIM];  
	float* filter_buf2 = new float[FILTER_DIM];
	float* filter_buf3 = new float[FILTER_DIM];
	float* mirror_psf1 = filter_buf1 + FILTER_EXTENT;
	float* mirror_psf2 = filter_buf2 + FILTER_EXTENT;
	float* mirror_psf3 = filter_buf3 + FILTER_EXTENT;
	sinc_windowing(FILTER_EXTENT, FILTER_DIM, mirror_psf1, 0);
	sinc_windowing(FILTER_EXTENT, FILTER_DIM, mirror_psf2, 1.0F / 3.0F);
	sinc_windowing(FILTER_EXTENT, FILTER_DIM, mirror_psf3, 2.0F / 3.0F);
	// r1->input row, r2->output, row
	// c1->input col, c2->output col
	for (int r1 = 0, r2 = 0; (r1 < in->height); r1 += 1, r2 += 3) {
		for (int c2 = 0, c1 = 0; (c2 < out->width); c2 ++, c1++)
		{
			float* ip = in->buf + r1 * in->stride + c1;
			float* op = out->buf + r2 * out->stride + c2;
			float sum1 = 0.0F;
			float sum2 = 0.0F;
			float sum3 = 0.0F;

			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
				sum1 += ip[y * in->stride] * mirror_psf1[y];
				op[0] = sum1;
				sum2 += ip[y * in->stride] * mirror_psf2[y];
				op[out->stride] = sum1;
				sum3 += ip[y * in->stride] * mirror_psf3[y];
				op[out->stride * 2] = sum3;

			}
		}
	}

}

void horizontal_filter(my_image_comp* in, my_image_comp* out, int filter_length) {
	int FILTER_EXTENT = filter_length;
	int FILTER_DIM = (2 * FILTER_EXTENT + 1); // No shift sinc kernel
	float* filter_buf1 = new float[FILTER_DIM];  // No shift sinc kernel
	float* filter_buf2 = new float[FILTER_DIM];
	float* filter_buf3 = new float[FILTER_DIM];
	float* mirror_psf1 = filter_buf1 + FILTER_EXTENT;
	float* mirror_psf2 = filter_buf2 + FILTER_EXTENT;
	float* mirror_psf3 = filter_buf3 + FILTER_EXTENT;
	sinc_windowing(FILTER_EXTENT, FILTER_DIM, mirror_psf1, 0);
	sinc_windowing(FILTER_EXTENT, FILTER_DIM, mirror_psf2, 1.0F/3.0F);
	sinc_windowing(FILTER_EXTENT, FILTER_DIM, mirror_psf3, 2.0F/3.0F);
	for (int r1 = 0, r2 = 0; (r2 < in->height); r1 += 1, r2 += 1) {
		for (int c1 = 0, c2 = 0; (c1 < in->width); c1 += 1, c2 += 3)
		{
			float* ip = in->buf + r2 * in->stride + c1;
			float* op = out->buf + r2 * out->stride + c2;
			float sum1 = 0.0F;
			float sum2 = 0.0F;
			float sum3 = 0.0F;

			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
				sum1 += ip[y] * mirror_psf1[y];
				op[0] = sum1;
				sum2 += ip[y] * mirror_psf2[y];
				op[1] = sum1;
				sum3 += ip[y] * mirror_psf3[y];
				op[2] = sum3;

			}

		}
	}
}