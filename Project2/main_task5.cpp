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
#include <time.h>


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
/*                Symmetry and non-negative filter                           */
/*****************************************************************************/
void Pre_filter(my_image_comp* in, my_image_comp* out)
{
#define FILTER_EXTENT 2
#define FILTER_DIM (2*FILTER_EXTENT+1)
#define FILTER_TAPS (FILTER_DIM*FILTER_DIM)

    // Create the filter kernel as a local array on the stack, which can accept
    // row and column indices in the range -FILTER_EXTENT to +FILTER_EXTENT.
    // `mirror_psf' points to the central tap in the filter
    // The HPF designed by the expression of unsharp masking substract the guassian LPF

    // Symmetry and non-negative
    float h1[FILTER_TAPS] = { -0.0571F, -0.0726F, -0.0786F, -0.0726F, -0.0571,
        -0.0726F, -0.0922F, -0.0999F, -0.0922F, -0.0726F,
        -0.0786F, -0.00999F, 2.8918F, -0.0786F, -0.00999F,
        -0.0726F, -0.0922F, -0.0999F, -0.0922F, -0.0726F,
        -0.0571F, -0.0726F, -0.0786F, -0.0726F, -0.0571F };

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
    int size_out = (out->height) * (out->width);
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
/*                                correlation_filter                         */
/*****************************************************************************/

void COR_fft_domain(my_image_comp* input_comps1, my_image_comp* padding_comps, int w1, int h1) {
    // initialized
    int in_height = input_comps1->height, in_width = input_comps1->width;
    int in_stride = input_comps1->stride;

    int p_height = padding_comps->height;
    int p_width = padding_comps->width;
    int p_stride = padding_comps->stride;

    // Perform dft
    // DFT for input image

    /***************************FFT for pattern************************/
    int r, c;
    _complex* p = new _complex[p_height * p_width];
    my_direct_fft hor_fft_p;
    hor_fft_p.init(p_width, true);
    _complex* line = new _complex[in_width];
    for (r = 0; r < p_height; r++) {
        for (c = 0; c < p_width; c++) {
            line[c].x = padding_comps->buf[r * p_stride + c];
            line[c].y = 0.0F;
        }
        line = hor_fft_p.FFT(line, in_height, 1);
        for (c = 0; c < p_width; c++) {
            p[r * in_height + c].x = line[c].x;
            p[r * in_height + c].y = line[c].y;
        }

    }

    my_direct_fft ver_fft_p;
    ver_fft_p.init(p_width, true);
    for (c = 0; c < p_width; c++) {
        for (r = 0; r < p_height; r++) {
            line[r].x = p[r * p_width + c].x;
            line[r].y = p[r * p_width + c].y;
        }
        line = ver_fft_p.FFT(line, in_width, 1);
        for (r = 0; r < p_height; r++) {
            p[r * in_width + c].x = line[r].x;
            p[r * in_width + c].y = line[r].y;
        }

    }

    /*******************************FFT for input image********************************/

    _complex* in = new _complex[in_height * in_width];
    my_direct_fft hor_fft_in, vert_fft_in;
    hor_fft_in.init(in_width, true);
    vert_fft_in.init(in_height, true);

    for (r = 0; r < p_height; r++) {
        for (c = 0; c < p_width; c++) {
            line[c].x = input_comps1->buf[r * p_stride + c];
            line[c].y = 0.0F;
        }
        line = hor_fft_p.FFT(line, in_height,1);
        for (c = 0; c < p_width; c++) {
            in[r * in_height + c].x = line[c].x;
            in[r * in_height + c].y = line[c].y;
        }

    }
    for (c = 0; c < in_width; c++) {
        for (r = 0; r < in_height; r++) {
            line[r].x = in[r * p_width + c].x;
            line[r].y = in[r * p_width + c].y;
        }
        line = ver_fft_p.FFT(line, in_width,1);
        for (r = 0; r < in_height; r++) {
            in[r * in_width + c].x = line[r].x;
            in[r * in_width + c].y = line[r].y;
        }

    }
    /*******************Take multiplication******************/
    _complex* res = new _complex[in_height * in_width];
    for (int i = 0; i < in_height * in_width; i++) {
        res[i].x = (in[i].x * p[i].x) + (in[i].y * p[i].y);
        res[i].y = p[i].x * in[i].y - in[i].x * p[i].y;
    }


    /********************inverse FFT*****************/
    // separable
    my_direct_fft hor_ifft, vert_ifft;
    hor_ifft.init(in_width, false);
    vert_ifft.init(in_height, false);

    for (r = 0; r < p_height; r++) {
        for (c = 0; c < p_width; c++) {
            line[c].x = res[r * p_width + c].x;
            line[c].y = res[r * p_width + c].y;
        }
        line = hor_ifft.FFT(line, in_height,0);
        for (c = 0; c < p_width; c++) {
            res[r * in_height + c].x = line[c].x;
            res[r * in_height + c].y = line[c].y;
        }

    }

    for (c = 0; c < in_width; c++) {
        for (r = 0; r < in_height; r++) {
            line[r].x = res[r * p_width + c].x;
            line[r].y = res[r * p_width + c].y;
        }

        line = vert_ifft.FFT(line, in_width,0);
        for (r = 0; r < p_height; r++) {
            res[r * p_width + c].x = line[r].x;
            res[r * p_width + c].y = line[r].y;
        }

    }

    /**********************Find Local Maxima**************************/
    float max = 0.0F;
    int plc = 0;
    int column = 0, row = 0;

    for (int i = 0; i < in_height * in_width; i++) {
        float mag = (float)sqrt(pow(res[i].x, 2) + pow(res[i].y, 2));
        if (mag > max) {
            max = mag;
            plc = i;
        }
    }

    // map to corresponding position and wrap around
    column = (plc % in_height + w1) % in_height;
    row = (plc / in_width + h1) % in_height;

    /************Print Results*******************/
    printf("padding: %d, %d\n", w1, h1);
    printf("height:%d\n", in_height);
    printf("location in the array: %d\n", plc);
    printf("The pattern is located at (column, row): (%d, %d)\n", column, row);
    printf("Maximum correlation is : %f\n", max);
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

        // Allocate storage for the output image


        /************* Pre-filtered the pattern ***************/
        my_image_comp* filtered_pattern = new my_image_comp[num_comps];
        my_image_comp* padding_pattern = new my_image_comp[num_comps];


        for (int n = 0; n < num_comps; n++) {
            filtered_pattern[n].init(pattern_height, pattern_width, 0);
            padding_pattern[n].init(height, width, 0);
            input_comps2[n].perform_boundary_extension();
        }

        for (int n = 0; n < num_comps; n++) {
            Pre_filter(input_comps2 + n, filtered_pattern + n);
            zero_intensity(filtered_pattern + n);
        }

        /************* zero padding ***************/

        // Calculate parameters for padding
        int w1, w2, h1, h2;
        h1 = (input_comps1->height - input_comps2->height) / 2;
        h2 = input_comps1->height - h1 - input_comps2->height;

        w1 = (input_comps1->width - input_comps2->width) / 2;
        w2 = input_comps1->width - w1 - input_comps2->width;

        // zero padding
        for (int n = 0; n < num_comps; n++) {
            zero_padding(filtered_pattern + n, padding_pattern + n, w1, w2, h1, h2);
        }


        /*************** Find correlation in DFT domain ***************/
        clock_t start_time = clock();
        COR_fft_domain(input_comps1, padding_pattern, w1, h1);
        clock_t end_time = clock();
        float computation_time = ((float)(end_time - start_time)) / CLOCKS_PER_SEC;
        printf("The computation time for FFT pattern matching:%lf seconds.\n", computation_time);


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
