/*****************************************************************************/
// File: dft_main.cpp
// Author: David Taubman
// Last Revised: 28 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include <math.h>
#include "io_bmp.h"
#include "image_comps.h"
#include "dft.h"


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */
/* ========================================================================= */
/*                              Boundary Extension                           */
/* ========================================================================= */
void my_image_comp::perform_boundary_extension()
{
    int r, c;

    // First extend upwards
    float* first_line = buf;
    for (r = 1; r <= border; r++)
        for (c = 0; c < width; c++)
            first_line[-r * stride + c] = first_line[c];

    // Now extend downwards
    float* last_line = buf + (height - 1) * stride;
    for (r = 1; r <= border; r++)
        for (c = 0; c < width; c++)
            last_line[r * stride + c] = last_line[c];

    // Now extend all rows to the left and to the right
    float* left_edge = buf - border * stride;
    float* right_edge = left_edge + width - 1;
    for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
        for (c = 1; c <= border; c++)
        {
            left_edge[-c] = left_edge[0];
            right_edge[c] = right_edge[0];
        }
}

/*****************************************************************************/
/*                                correlation_filter                         */
/*****************************************************************************/
void phaseCOR_dft_domain(my_image_comp* input_comps1, my_image_comp* padding_pattern, int w1, int h1) {
    // initialized
    int in_height = input_comps1->height, in_width = input_comps1->width;
    int in_stride = input_comps1->stride;

    int p_height = padding_pattern->height;
    int p_width = padding_pattern->width;
    int p_stride = padding_pattern->stride;

    // Perform dft
    // DFT for input image
    int r, c;
    float* dft_real_in = new float[in_height * in_width];
    float* dft_imag_in = new float[in_height * in_width];
    // First copy all samples to the `dft_real' buffer
    for (r = 0; r < in_height; r++)
        for (c = 0; c < in_width; c++)
        {
            dft_real_in[r * in_width + c] = input_comps1->buf[r * in_stride + c];
            dft_imag_in[r * in_width + c] = 0.0F;
        }

    // Next, perform the 2D DFT seperately
    my_direct_dft hor_dft_in, vert_dft_in;
    hor_dft_in.init(in_width, true);
    vert_dft_in.init(in_height, true);

    for (r = 0; r < in_height; r++)
        hor_dft_in.perform_transform(dft_real_in + r * in_width, dft_imag_in + r * in_width, 1);
    for (c = 0; c < in_width; c++)
        vert_dft_in.perform_transform(dft_real_in + c, dft_imag_in + c, in_width);

    // DFT for pattern
    float* dft_real_p = new float[p_height * p_width];
    float* dft_imag_p = new float[p_height * p_width];
    // First copy all samples to the `dft_real' buffer
    for (r = 0; r < p_height; r++)
        for (c = 0; c < p_width; c++)
        {
            dft_real_p[r * p_width + c] = padding_pattern->buf[r * p_stride + c];
            dft_imag_p[r * p_width + c] = 0.0F;
        }

    // Next, perform the 2D DFT seperately
    my_direct_dft hor_dft_p, vert_dft_p;
    hor_dft_p.init(p_stride, true);
    vert_dft_p.init(p_stride, true);

    for (r = 0; r < p_height; r++)
        hor_dft_in.perform_transform(dft_real_p + r * p_width, dft_imag_p + r * p_width, 1);
    for (c = 0; c < p_width; c++)
        vert_dft_in.perform_transform(dft_real_p + c, dft_imag_p + c, p_width);

    // perform multiplication
    float* y_real = new float[in_height * in_width];
    float* y_imag = new float[in_height * in_width];
    float* y_mag = new float[in_height * in_width];
    for (int i = 0; i < in_height * in_width; i++) {
        y_real[i] = (dft_real_in[i] * dft_real_p[i]) + (dft_imag_in[i] * dft_imag_p[i]);
        y_imag[i] = dft_real_p[i] * dft_imag_in[i] - dft_real_in[i] * dft_imag_p[i];
        y_mag[i] = sqrt(pow(dft_real_in[i], 2) + pow(dft_imag_in[i], 2)) * sqrt(pow(dft_real_p[i], 2) + pow(dft_imag_p[i], 2));
        y_real[i] = y_real[i] / y_mag[i];
        y_imag[i] = y_imag[i] / y_mag[i];
    }

    // printf("imaginary part:%f, real part:%f\n", y_real[1], y_imag[1]);
    // Next, perform the 2D DFT seperately
    my_direct_dft hor_idft, vert_idft;
    hor_idft.init(in_width, false);
    vert_idft.init(in_height, false);
    for (r = 0; r < in_height; r++)
        hor_idft.perform_transform(y_real + r * in_width, y_imag + r * in_width, 1);
    for (c = 0; c < in_width; c++)
        vert_idft.perform_transform(y_real + c, y_imag + c, in_width);

    // find local maxima
    float max = 0;
    int plc = 0;
    for (int i = 0; i < in_height * in_width; i++) {
        float mag = (float)sqrt(pow(y_real[i], 2) + pow(y_imag[i], 2));
        //printf("mag:%f\n", mag);
        if (mag > max) {
            max = mag;
            plc = i;
        }
    }

    int x = (plc % in_height + w1) % in_height;
    int y = (plc / in_height + h1) % in_height;

    printf("location in the array: %d\n", plc);
    printf("padding: %d, %d\n", w1, h1);
    printf("The pattern is located at (column, row): (%d, %d)\r\n", x, y);
    printf("Maximum correlation is : %f\r\n", max);
}

/*****************************************************************************/
/*                                zero_padding                               */
/*****************************************************************************/
void zero_padding(my_image_comp* in, my_image_comp* out, int w1, int w2, int h1, int h2) {
    int r, c;
    int height = out->height;
    int width = out->width;

    for (r = 0; r < height; r++) {
        for (c = 0; c < width; c++) {
            float* ip = in->buf + (r - h1) * in->stride + (c - w1);
            float* op = out->buf + r * out->stride + c;
            if (r < h1 || r >= h1 + in->height || c < w1 || c >= w1 + in->width) {
                *op = 0.0F;
            }
            else {
                *op = *ip;
            }

        }
    }

}

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/
int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s <in bmp file> <in bmp pattern>\n", argv[0]);
        return -1;
    }

    int err_code = 0;
    try {
        // Read the input image
        bmp_in in;
        if ((err_code = bmp_in__open(&in, argv[1])) != 0)
            throw err_code;

        int width = in.cols;
        int height = in.rows;
        int n, num_comps = in.num_components;
        my_image_comp* input_comps1 = new my_image_comp[num_comps];
        for (n = 0; n < num_comps; n++)
            input_comps1[n].init(height, width, 0); // Leave a border of 0

        int r; // Declare row index
        io_byte* line = new io_byte[width * num_comps];
        for (r = height - 1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
            if ((err_code = bmp_in__get_line(&in, line)) != 0)
                throw err_code;
            for (n = 0; n < num_comps; n++)
            {
                io_byte* src = line + n; // Points to first sample of component n
                float* dst = input_comps1[n].buf + r * input_comps1[n].stride;
                for (int c = 0; c < width; c++, src += num_comps)
                    dst[c] = (float)*src; // The cast to type "float" is not
                          // strictly required here, since bytes can always be
                          // converted to floats without any loss of information.
            }
        }
        bmp_in__close(&in);

        int err_code = 0;
        // Read the input image
        if ((err_code = bmp_in__open(&in, argv[2])) != 0)
            throw err_code;

        int pattern_width = in.cols;
        int pattern_height = in.rows;
        int padding_amount = (height - pattern_height) / 2;
        n, num_comps = in.num_components;
        my_image_comp* input_comps2 = new my_image_comp[num_comps];
        for (n = 0; n < num_comps; n++)
            input_comps2[n].init(pattern_height, pattern_width, padding_amount); // Leave a border of 0

        for (r = pattern_height - 1; r >= 0; r--)
        {
            // "r" holds the true row index we are reading, since the image is
            // stored upside down in the BMP file.
            if ((err_code = bmp_in__get_line(&in, line)) != 0)
                throw err_code;
            for (n = 0; n < num_comps; n++)
            {
                io_byte* src = line + n; // Points to first sample of component n
                float* dst = input_comps2[n].buf + r * input_comps2[n].stride;
                for (int c = 0; c < pattern_width; c++, src += num_comps)
                    dst[c] = (float)*src; // The cast to type "float" is not
                            // strictly required here, since bytes can always be
                            // converted to floats without any loss of information.
            }
        }

        bmp_in__close(&in);

        // Allocate storage for the output image
        int w1, w2, h1, h2;
        h1 = (input_comps1->height - input_comps2->height) / 2;
        h2 = input_comps1->height - h1 - input_comps2->height;

        w1 = (input_comps1->width - input_comps2->width) / 2;
        w2 = input_comps1->width - w1 - input_comps2->width;

        my_image_comp* padding_pattern = new my_image_comp[num_comps];
        /************* zero padding ***************/
        for (int n = 0; n < num_comps; n++) {
            padding_pattern[n].init(height, width, 0);
            zero_padding(input_comps2, padding_pattern, w1, w2, h1, h2);
        }

        /*************** Find phase correlation in DFT domain ***************/
        for (int n = 0; n < num_comps; n++) {
            phaseCOR_dft_domain(input_comps1, padding_pattern, w1, h1);
        }

    }
    catch (int exc) {
        if (exc == IO_ERR_NO_FILE)
            fprintf(stderr, "Cannot open supplied input or output file.\n");
        else if (exc == IO_ERR_FILE_HEADER)
            fprintf(stderr, "Error encountered while parsing BMP file header.\n");
        else if (exc == IO_ERR_UNSUPPORTED)
            fprintf(stderr, "Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
        else if (exc == IO_ERR_FILE_TRUNC)
            fprintf(stderr, "Input or output file truncated unexpectedly.\n");
        else if (exc == IO_ERR_FILE_NOT_OPEN)
            fprintf(stderr, "Trying to access a file which is not open!(?)\n");
        return -1;
    }
    return 0;
}
