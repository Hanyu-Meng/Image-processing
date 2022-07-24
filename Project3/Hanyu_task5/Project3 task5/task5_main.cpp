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
#include <vector>

/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_symmetry_extension                */
/*****************************************************************************/
void my_image_comp::perform_symmetry_extension() {
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
            left_edge[-c] = left_edge[c];
            right_edge[c] = right_edge[-c];
        }
}


/*****************************************************************************/
/*                  my_image_comp::perform_zero_extension                    */
/*****************************************************************************/
void my_image_comp::perform_zero_extension() {
    int r, c;

    // First extend upwards
    float* first_line = buf;
    for (r = 1; r <= border; r++)
        for (c = 0; c < width; c++)
            first_line[-r * stride + c] = 0;

    // Now extend downwards
    float* last_line = buf + (height - 1) * stride;
    for (r = 1; r <= border; r++)
        for (c = 0; c < width; c++)
            last_line[r * stride + c] = 0;

    // Now extend all rows to the left and to the right
    float* left_edge = buf - border * stride;
    float* right_edge = left_edge + width - 1;
    for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
        for (c = 1; c <= border; c++)
        {
            left_edge[-c] = 0;
            right_edge[c] = 0;
        }
}

/*****************************************************************************/
/*                              set output to background                     */
/*****************************************************************************/
void output_init(my_image_comp* g) {
    int r, c;
    for (r = 0; r < g->height; r++) {
        for (c = 0; c < g->width; g++) {
            float* op = g->buf + r * g->stride + c;
            *op = 0.0F;
        }
    }
}

/*****************************************************************************/
/*                              Opening                                      */
/*****************************************************************************/
void opening(my_image_comp* in, my_image_comp* out, int N_a, int* a_off) {
    int r, c;
    int b = in->border;
    for (r = b/2; r < in->height + b + b/2; r++) {
        for (c = b/2; c < in->width + b + b/2; c++) {
            float* pn = in->handle + r * in->stride + c;
            float* op = out->handle + r * out->stride + c;
            int val = 255;
            for (int i = 0; i < N_a; i++) {
                val &= (int)pn[a_off[i]];
            }
            if (val != 0) {
                for (int i = 0; i < N_a; i++) {
                    op[a_off[i]] = 255;
                }   
            }
        }
    }
}
/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/
int main(int argc, char* argv[])
{
    if (argc < 5)
    {
        fprintf(stderr, "Usage: %s <in bilevel bmp file> <out bmp file> <radius r/a1> <c1/a2> <c2/a1'> <a2'...>\n", argv[0]);
        return -1;
    }

    int N_a;
    int* a1;
    int* a2;
    if (argc == 6) {
        printf("circle mode\n");
        /***************Get parameters from command line**************/
        int radius = atoi(argv[3]);
        int c1 = atoi(argv[4]);
        int c2 = atoi(argv[5]);
        float R = radius + 0.5;
        printf("center:(%d, %d)\n", c1, c2);
        printf("Radius: %d\n", radius);
        /***********Creat the array*************/
        int size = (2 * radius + 1) * (2 * radius + 1);
        a1 = new int[size];
        a2 = new int[size];
        int x, y;
        int n = 0;
        for (y = c2 - radius; y <= radius + c2; y++) {
            for (x = c1 - radius; x <= radius + c1; x++) {
                if ((x - c1) * (x - c1) + (y - c2) * (y - c2) < R * R) {
                    a1[n] = y;
                    a2[n] = x;
                    printf("a1[%d],a2[%d] : (%d, %d)\n", n, n, a1[n], a2[n]);
                    n++;
                }
            }
        }
        N_a = n;
    }
    else {
        printf("vector mode\n");
        N_a = (argc - 3) / 2;
        a1 = new int[N_a];
        a2 = new int[N_a];
        for (int i = 0; i < N_a; i++) {
            a1[i] = atoi(argv[3 + i * 2]);
            a2[i] = atoi(argv[4 + i * 2]);
            printf("a1[%d]: %d\n", i, *a1);
            printf("a2[%d]: %d\n", i, *a2);
        }
    }

    printf("N_a:%d\n", N_a);


    /**********Find Max**********/
    int max = abs(a1[0]);

    for (int i = 0; i < N_a; i++) {
        if (abs(a1[i]) > max) {
            max = abs(a1[i]);
        }

        if (a2[i] > max) {
            max = abs(a2[i]);
        }
    }

    int border = 2 * max;
    printf("Border: %d\n", border);

    int err_code = 0;
    try {
        // Read the input image
        bmp_in in;
        if ((err_code = bmp_in__open(&in, argv[1])) != 0)
            throw err_code;

        int width = in.cols, height = in.rows;
        int n, num_comps = in.num_components;
        my_image_comp* input_comps = new my_image_comp[num_comps];
        for (n = 0; n < num_comps; n++)
            input_comps[n].init(height, width, border); // Do not need to extent

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

        /**************Construct displacement vector --> a_off***************/
        int S = input_comps->stride;
        int* a_off = new int[N_a];
        for (int i = 0; i < N_a; i++) {
            a_off[i] = a1[i] * S + a2[i];
            printf("a_off_mirror[%d]: %d\n", i, a_off[i]);
        }
        printf("stride:%d\n", S);

        /**************Allocate storage for the filtered output**************/
        my_image_comp* output_comps = new my_image_comp[num_comps];

        /************set the output image as zero***************/
        for (n = 0; n < num_comps; n++) {
            output_comps[n].init(height, width, border);
            output_init(output_comps + n);
        }

        for (n = 0; n < num_comps; n++)
            input_comps[n].perform_symmetry_extension();
            //input_comps[n].perform_zero_extension();

        /*********************Opening*********************/
        for (n = 0; n < num_comps; n++)
            opening(input_comps + n, output_comps + n, N_a, a_off);


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
