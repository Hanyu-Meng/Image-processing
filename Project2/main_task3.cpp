/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007

/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"
#include <math.h>
/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/* ========================================================================= */
/*                              Global Functions                             */
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
/*                Symmetry and non-negative filter                           */
/*****************************************************************************/
void Pre_filter(my_image_comp* in, my_image_comp* out)
{
#define FILTER_EXTENT 2
#define FILTER_DIM (2*FILTER_EXTENT+1)
#define FILTER_TAPS (FILTER_DIM*FILTER_DIM)

    // Create the filter kernel as a local array on the stack, which can accept
    // row and column  indices in the range -FILTER_EXTENT to +FILTER_EXTENT.
    // `mirror_psf' points to the central tap in the filter
    // The HPF designed by the expression of unsharp masking substract the guassian LPF
    float h1[FILTER_TAPS] = {-0.0635F, -0.0750F, -0.0793F, -0.0750F, -0.0635F,
        -0.0750F, -0.0886F, -0.0937F, -0.0886F, -0.0750F,
        -0.0793F, -0.0937F, 2.9009F, -0.0937F, -0.0793F,
        -0.0750F, -0.0886F, -0.0937F, -0.0886F, -0.0750F,
        -0.0635F, -0.0750F, -0.0793F, -0.0750F, -0.0635F};

    float* mirror_hpf = h1 + (FILTER_DIM * FILTER_EXTENT) + FILTER_EXTENT;

    // `mirror_psf' points to the central tap in the filter
    // Check for consistent dimensions
    assert(in->border >= FILTER_EXTENT);
    assert((out->height <= in->height) && (out->width <= in->width));

    int r, c;
    // Perform the convolution
    for (r = 0; r < out->height; r++)
        for (c = 0; c < out->width; c++)
        {
            float* ip = in->buf + r * in->stride + c;
            float* op = out->buf + r * out->stride + c;
            float sum = 0.0F;
            for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++)
                for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++)
                    sum += ip[y * in->stride + x] * mirror_hpf[y * FILTER_DIM + x];
            *op = sum;
        }
}
/*****************************************************************************/
/*                               zero_intensity                              */
/*****************************************************************************/

void zero_intensity(my_image_comp* out) {
    float sum = 0.0F;
    float average = 0.0F;
    int r, c;
    for (r = 0; r < out->height; r++) {
        for (c = 0; c < out->width; c++) {
            float* op = out->buf + r * out->stride + c;
            sum += op[0];
        }
    }
    int size_out = out->height * out->width;
    average = sum / size_out;

    for (r = 0; r < out->height; r++) {
        for (c = 0; c < out->width; c++) {
            float* op = out->buf + r * out->stride + c;
            op[0] = op[0] - average;
            //printf("Intensity:%f\r\n", *op);
        }
    }

}

/*****************************************************************************/
/*                                correlation_filter                              */
/*****************************************************************************/
void COR_pattern_matching(my_image_comp* in1, my_image_comp* in2) {
    // in1 -> big picture; in2 -> pattern;
    int N1 = in1->height, N2 = in1->width;
    int P1 = in2->height, P2 = in2->width;
    int k1, k2;
    int p1, p2;
    int x = 0, y = 0;
    float max = 0;
    float multiple;
    printf("Pattern size P1: %d, P2:%d\n", P1, P2);
    printf("Input size N1: %d, N2:%d\n", N1, N2);
    for (p1 = 0; p1 <= (N1 - P1); p1 += 1) {
        for (p2 = 0; p2 <= (N2 - P2); p2 += 1) {
            float sum = 0.0F;
            for (k1 = 0; k1 < P1; k1++) {
                for (k2 = 0; k2 < P2; k2++) {
                    float* ip1 = in1->buf + (p1 + k1) * in1->stride + (p2 + k2);
                    float* ip2 = in2->buf + k1 * in2->stride + k2;
                    multiple = ip1[0] * ip2[0];
                    sum = sum + multiple;
                }
            }
            if (sum > max) {
                max = sum;
                x = p2;
                y = p1;
            }
        }
    }
    printf("The pattern is located at (%d, %d)\r\n", x, y);
    printf("Maximum correlation is : %f\r\n", max);
}



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
        n, num_comps = in.num_components;
        my_image_comp* input_comps2 = new my_image_comp[num_comps];
        for (n = 0; n < num_comps; n++)
            input_comps2[n].init(pattern_height, pattern_width, 2); // Leave a border of 0

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

        /***********Process the image**************/
        // Allocate storage for the filtered output
        my_image_comp* output_comps = new my_image_comp[num_comps];
        my_image_comp* filtered_pattern = new my_image_comp[num_comps];
        for (int n = 0; n < num_comps; n++) {
            //output_comps[n].init(height, width, 0); // Don't need a border for output
            filtered_pattern[n].init(pattern_height, pattern_width, 0);
            input_comps2[n].perform_boundary_extension();
        }
        // Process the image, all in floating point (easy)
        for (int n = 0; n < num_comps; n++) {
            Pre_filter(input_comps2 + n, filtered_pattern + n);
            zero_intensity(filtered_pattern + n);
        }

        for (int n = 0; n < num_comps; n++) {
            COR_pattern_matching(input_comps1, filtered_pattern + n);
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