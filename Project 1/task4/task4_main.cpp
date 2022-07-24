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

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::perform_boundary_extension()
{
  int r, c;

  // First extend upwards
  float *first_line = buf;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      first_line[-r*stride+c] = first_line[c];

  // Now extend downwards
  float *last_line = buf+(height-1)*stride;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      last_line[r*stride+c] = last_line[c];

  // Now extend all rows to the left and to the right
  float *left_edge = buf-border*stride;
  float *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++) {
        left_edge[-c] = left_edge[c];
        right_edge[c] = right_edge[-c];
    }
}

void difference_filter(my_image_comp* in1, my_image_comp* in2, my_image_comp* out) {
    float sum1 = 0.0F;
    float sum2 = 0.0F;
    float delta = 0.0F;
    int height = 0, width = 0;
    if (in1->height > in2->height) {
        height = in2->height;
    }
    else {
        height = in1->height;
    }

    if (in1->width > in2->width) {
        width = in2->width;
    }
    else {
        width = in1->width;
    }

    int r, c;
    for (r = 0; r < height; r++) {
        for (c = 0; c < width; c++) {
            float* ip1 = in1->buf + r * in1->stride + c;
            float* ip2 = in2->buf + r * in2->stride + c;
            float* op = out->buf + r * out->stride + c;
            delta = ip1[0] - ip2[0];
            *op = 0.5F * delta; //difference between 2 images
            sum1 = sum1 + delta;
            sum2 = sum2 + delta * delta;
        }
    }

    // calculate parameters
    float ME = sum1 / (height*width);
    float MSE = sum2 / (height*width);
    float PSNR = 10.0F * log10(255.0F * 255.0F / MSE);

    // print the parameters
    printf("Mean Error (ME) is : %f\r\n", ME);
    printf("MSE is : %f\r\n", MSE);
    printf("PSNR is : %f\r\n", PSNR);

}

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        fprintf(stderr, "Usage: %s <in bmp file1> <in bmp file1> <out bmp file> \n", argv[0]);
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
            input_comps1[n].init(height, width, 4); // Leave a border of 4

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

        width = in.cols;
        height = in.rows;
        n, num_comps = in.num_components;
        my_image_comp* input_comps2 = new my_image_comp[num_comps];
        for (n = 0; n < num_comps; n++)
            input_comps2[n].init(height, width, 4); // Leave a border of 4

        for (r = height - 1; r >= 0; r--)
        { 
            // "r" holds the true row index we are reading, since the image is
            // stored upside down in the BMP file.
            if ((err_code = bmp_in__get_line(&in, line)) != 0)
                throw err_code;
            for (n = 0; n < num_comps; n++)
            {
                io_byte* src = line + n; // Points to first sample of component n
                float* dst = input_comps2[n].buf + r * input_comps2[n].stride;
                for (int c = 0; c < width; c++, src += num_comps)
                    dst[c] = (float)*src; // The cast to type "float" is not
                            // strictly required here, since bytes can always be
                            // converted to floats without any loss of information.
            }
        }
        bmp_in__close(&in);


        /***********Process the image**************/
        // Allocate storage for the filtered output
        my_image_comp* output_comps = new my_image_comp[num_comps];

        for (int n = 0; n < num_comps; n++)
            output_comps[n].init(height, width, 0); // Don't need a border for output

            // Process the image, all in floating point (easy)
        for (int n = 0; n < num_comps; n++) {
            difference_filter(input_comps1 + n, input_comps2 + n, output_comps + n);
        }



        //write back
            // Write the image back out again
        bmp_out out;
        if ((err_code = bmp_out__open(&out, argv[3], width, height, num_comps)) != 0)
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
        delete[] input_comps1;
        delete[] input_comps2;
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