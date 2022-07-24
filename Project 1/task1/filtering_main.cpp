/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
//#include "image_comps.h"
#include "filters.h"
#include <math.h>
#include <time.h>

/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::perform_boundary_extension() {
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
    for (c=1; c <= border; c++)
    {
        left_edge[-c] = left_edge[0];
        right_edge[c] = right_edge[0];
    }
}


/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/


int main(int argc, char *argv[])
{
  if (argc != 5)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file> <filter length> <mode: 0 for direct, 1 for separate>\n",argv[0]);
      return -1;
    }

  int err_code=0;

  int filter_length = atoi(argv[3]);
  int mode = atoi(argv[4]);

  try {
      // Read the input image
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

      int width = in.cols, height = in.rows;
      int n, num_comps = in.num_components;
      my_image_comp *input_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        input_comps[n].init(height,width,filter_length); // Leave a border of 15
      
      int r; // Declare row index
      io_byte *line = new io_byte[width*num_comps];
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
          if ((err_code = bmp_in__get_line(&in,line)) != 0)
            throw err_code;
          for (n=0; n < num_comps; n++)
            {
              io_byte *src = line+n; // Points to first sample of component n
              float *dst = input_comps[n].buf + r * input_comps[n].stride;
              for (int c=0; c < width; c++, src+=num_comps)
                dst[c] = (float) *src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
            }
        }
      bmp_in__close(&in);

      /*show modes*/
      if (mode == 0) {
          printf("Direct Resampling Mode\n");
      }
      else if (mode == 1) {
          printf("Seperate Resampling Mode\n");
      }
      else {
          printf("Please chose implementation mode!\n");
          return 0;
      }


      // Change the output image size
      int output_width, output_height;
      output_width = (int)floor(width/3.0F);
      output_height = (int)floor(height/3.0F);

      // Allocate storage for the filtered output
      my_image_comp *output_comps = new my_image_comp[num_comps];

      // check the input and output size
      printf("Input image height: % d\n", height);
      printf("Input image width: %d\n", width);
      printf("Output image height: %d\n", output_height);
      printf("Output image width: %d\n", output_width);

      // initialize the output
      for (n = 0; n < num_comps; n++)
          output_comps[n].init(output_height, output_width, filter_length);

      //boundary extent input image
      for (n = 0; n < num_comps; n++)
          input_comps[n].perform_boundary_extension(); //boundary extent input image

      clock_t start_time = clock();
      if (mode == 0) {
          printf("Direct Fitering starts!\n");
          for (n = 0; n < num_comps; n++) {
              directly_filter(input_comps + n, output_comps + n, filter_length);
          }
      }
      else if (mode == 1) {
          printf("Separable Filtering starts!\n"); // perform the separable filtering, first horizontally then ve
          // allocate intermedia image
          my_image_comp* inter_comps = new my_image_comp[num_comps]; 
          // initialize the intermedia image
          for (n = 0; n < num_comps; n++)
            inter_comps[n].init(height, output_width, filter_length);

          for (n = 0; n < num_comps; n++) {
            horizontal_filter(input_comps + n, inter_comps + n, filter_length);
            // vertical_filter(input_comps + n, inter_comps + n, filter_length);
            inter_comps[n].perform_boundary_extension();
          }

          for (n = 0; n < num_comps; n++) {
           // horizontal_filter(inter_comps + n, output_comps + n, filter_length);
            vertical_filter(inter_comps + n, output_comps + n, filter_length);
          }
          delete[] inter_comps;
      }
      else {
          printf("Error: Can not find the mode!\n");
      }

      clock_t end_time = clock();
      float computation_time = ((float)(end_time - start_time)) / CLOCKS_PER_SEC;
      printf("The computation time for filtering:%lf seconds.\n", computation_time);

     // Write the image back out again
      bmp_out out;
      if ((err_code = bmp_out__open(&out, argv[2], output_width, output_height, num_comps)) != 0)
          throw err_code;
      io_byte* output_line = new io_byte[output_width * num_comps];  // shrinked line size
      for (r = output_height - 1; r >= 0; r--)
      { // "r" holds the true row index we are writing, since the image is
        // written upside down in BMP files.
          for (n = 0; n < num_comps; n++)
          {
              io_byte* dst = output_line + n; // Points to first sample of component n
              float* src = output_comps[n].buf + r * output_comps[n].stride;
              for (int c = 0; c < output_width; c++, dst += num_comps) {
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
