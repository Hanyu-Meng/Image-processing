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
#include <time.h>

/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/
void my_image_comp::perform_symmetry_extension() {
    int r, c;

    // First extend upwards by border
    float* first_line = buf;
    for (r = 1; r <= border; r++)
        for (c = 0; c < width; c++)
            first_line[-r * stride + c] = first_line[r * stride + c];

    // Now extend downwards by border
    float* last_line = buf + (height - 1) * stride;
    for (r = 1; r <= border; r++)
        for (c = 0; c < width; c++)
            last_line[r * stride + c] = last_line[-r * stride + c];

    // Now extend all rows to the left and to the right
    float* left_edge = buf - border * stride;
    float* right_edge = left_edge + width - 1;
    for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
        for (c = 1; c <= border; c++)
        {
            left_edge[-c] = left_edge[c];
            right_edge[c] = right_edge[-c];
        }
}

/*****************************************************************************/
/*                              Erosion                                      */
/*****************************************************************************/
void erosion(my_image_comp* in, my_image_comp* out, int N_a, int* a_off) {
    int r, c;
    for (r = 0; r < out->height; r++) {
        for (c = 0; c < out->width; c++) {
            float* pn = in->buf + r * in->stride + c;
            float* op = out->buf + r * out->stride + c;
            int val = 255;
            for (int i = 0; i < N_a; i++) {
                val &= (int)pn[a_off[i]];
            }
            *op = val;
            //printf("%f\n", *op);
        }
    }
}
/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/


int main(int argc, char* argv[])
{
    if (argc < 5 || argc % 2 == 0)
    {
        fprintf(stderr, "Usage: %s <in bilevel bmp file> <out bmp file> <vertical displacement> <horizontal displacement>\n", argv[0]);
        return -1;
    }

    int err_code = 0;
    /***************Find largest border**************/
    // border = max{a1, a2}
    int max_border = 0;
    for (int i = 3; i < argc; i++) {
        if (abs(atoi(argv[i])) > max_border) {
            max_border = abs(atoi(argv[i]));
        }
    }
    printf("Border: %d\n", max_border);

    try {
        // Read the input image
        bmp_in in;
        if ((err_code = bmp_in__open(&in, argv[1])) != 0)
            throw err_code;

        int width = in.cols, height = in.rows;
        int n, num_comps = in.num_components;
        my_image_comp* input_comps = new my_image_comp[num_comps];
        for (n = 0; n < num_comps; n++)
            input_comps[n].init(height, width, max_border); // Do not need to extent

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
                float* dst = input_comps[n].buf + r * input_comps[n].stride;
                for (int c = 0; c < width; c++, src += num_comps)
                    dst[c] = (float)*src; // The cast to type "float" is not
                          // strictly required here, since bytes can always be
                          // converted to floats without any loss of information.
            }
        }
        bmp_in__close(&in);

        int N_a = (argc - 3) / 2;
        printf("N_a:%d\n", N_a);
        int* a_off = new int[N_a];
        // vector: a_off contains actual displacement of the image
        for (int i = 0; i < N_a; i++) {
            printf("(a1[%d], a2[%d]): (%d, %d)\n", i, i, atoi(argv[3 + i * 2]), atoi(argv[4 + i * 2]));
            a_off[i] = atoi(argv[3 + i * 2]) * input_comps->stride + atoi(argv[4 + i * 2]);
            printf("a_off[%d]: %d\n", i, a_off[i]);
        }

        // Allocate storage for the filtered output
        my_image_comp* output_comps = new my_image_comp[num_comps];

        // initialize the output that has the same size as input
        for (n = 0; n < num_comps; n++)
            output_comps[n].init(height, width, 0);

        //boundary extent input image
        for (n = 0; n < num_comps; n++)
            input_comps[n].perform_symmetry_extension(); //boundary extent input image

        /*********************Erosion*********************/
        for (n = 0; n < num_comps; n++)
            erosion(input_comps + n, output_comps + n, N_a, a_off);



       // clock_t start_time = clock();


        //   clock_t end_time = clock();
        //   float computation_time = ((float)(end_time - start_time)) / CLOCKS_PER_SEC;
        //   printf("The computation time for filtering:%lf seconds.\n", computation_time);

        
        // Write the image back out again
        bmp_out out;
        if ((err_code = bmp_out__open(&out, argv[2], width, height, num_comps)) != 0)
            throw err_code;
        io_byte* output_line = new io_byte[width * num_comps];  // shrinked line size
        for (r = height - 1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
            for (n = 0; n < num_comps; n++)
            {
                io_byte* dst = output_line + n; // Points to first sample of component n
                float* src = output_comps[n].buf + r * output_comps[n].stride;
                for (int c = 0; c < width; c++, dst += num_comps) {
                    if (src[c] > 255.5) {
                        src[c] = 255;
                    }
                    else if (src[c] < 0.5) {
                        src[c] = 0;
                    }
                    else {
                        src[c] = int(src[c] + 0.5F);
                    }

                    *dst = (io_byte)src[c];
                }
            }
            bmp_out__put_line(&out, output_line);
        }
        bmp_out__close(&out);
        delete[] line;
        delete[] input_comps;
        delete[] output_line;
        delete[] output_comps;
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
