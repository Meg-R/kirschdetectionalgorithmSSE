#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h> ///AVX
#include "bmp.h"
#include "measure.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#ifdef _WIN32
#include "Windows.h"
#include "malloc.h"
#define posix_memalign(address, alignment, size) *(address)=_aligned_malloc((size), (alignment))
#define sleep(s) Sleep(1000*s)
#else
#include <cstring>
#include "unistd.h"
#define memcpy_s(d, n, s, c) memcpy(d, s, c)
#endif
/* just used for time measurements */
#define REP 1
#define MIN(X, Y) (((X)<(Y))? X:Y)


#define DEBUG
using namespace std;

void kirsch_operator_sse(unsigned char *data_out,
	unsigned char *data_in, unsigned height,
	unsigned width)
{
    //your code 
	 /*
       Fill in your code here.
     */

  unsigned int size, i, j, lay;
  unsigned num_of_ints_per_row = width;
  while ((num_of_ints_per_row << 2) % 16)
    num_of_ints_per_row++;
  size = height * width;
  unsigned short *color_rows[3];

  /*
      We divide the 3 colours into 3 separate arrays. 
      We malloc an amount of memory based on number of ints 
      per row x height and set it to 64 bit aligned.
  */
	
	//int storeval =0;
  /*posix_memalign((void **)&color_rows[0], 64,
    sizeof(unsigned short) * num_of_ints_per_row * height);
  posix_memalign((void **)&color_rows[1], 64,
    sizeof(unsigned short) * num_of_ints_per_row * height);
  posix_memalign((void **)&color_rows[2], 64,
    sizeof(unsigned short) * num_of_ints_per_row * height);*/

#ifdef _WIN32
	void* ret = posix_memalign((void **)&color_rows[0], 64,
    sizeof(unsigned short) * num_of_ints_per_row * height);
	if (ret == nullptr) 	throw std::runtime_error("posix memalign");

	ret = posix_memalign((void **)&color_rows[1], 64,
    sizeof(unsigned short) * num_of_ints_per_row * height);
	if (ret == nullptr) 	throw std::runtime_error("posix memalign");

	ret = posix_memalign((void **)&color_rows[2], 64,
    sizeof(unsigned short) * num_of_ints_per_row * height);
	if (ret == nullptr) 	throw std::runtime_error("posix memalign");

#else
	int ret = posix_memalign((void **)&color_rows[0], 64,
    sizeof(unsigned short) * num_of_ints_per_row * height);
	if (ret) 	throw std::runtime_error("posix memalign");

	ret = posix_memalign((void **)&color_rows[1], 64,
    sizeof(unsigned short) * num_of_ints_per_row * height);
	if (ret) 	throw std::runtime_error("posix memalign");

	ret = posix_memalign((void **)&color_rows[2], 64,
    sizeof(unsigned short) * num_of_ints_per_row * height);
	if (ret) 	throw std::runtime_error("posix memalign");
#endif

  /*
      We make the the 3 arrays points to the data that is 
      supplied by the bitmap that is being passed in.
  */
  unsigned short *temp_rows;
  unsigned char *my_data_temp_in = data_in;
  for (lay = 0; lay < 3; lay++) {
    temp_rows = color_rows[lay];
    for (i = 0; i < height; i++) {
      for (j = 0; j < width; j++) {
        temp_rows[i * width + j] = my_data_temp_in[i * width + j];
      }
    }
    my_data_temp_in += size;
  }

  /*
      A struct definition which will be used later on to
      assign values.
  */
  typedef union {
    __m128i value;
    unsigned char vector[16];
  } vec_or_val;

  // Loop the 3 arrays
  for (lay = 0; lay < 3; lay++) {
    unsigned short *my_data_in_rmp = color_rows[lay];
    for (i = 1; i < height - 1; ++i) 
	{
      /* 
          Retrieve the following pointers, pointing to the 
          3 consecutive rows od sections upper, middle, lower
      */
      __m128i *upper_row_ptr =
        (__m128i *) & my_data_in_rmp[(i - 1) * num_of_ints_per_row];
      __m128i *middle_row_ptr =
        (__m128i *) & my_data_in_rmp[i * num_of_ints_per_row];
      __m128i *lower_row_ptr =
        (__m128i *) & my_data_in_rmp[(i + 1) * num_of_ints_per_row];

      // Load values
      __m128i upper_row, lower_row, middle_row;
      upper_row  = _mm_load_si128(upper_row_ptr);
      middle_row = _mm_load_si128(middle_row_ptr);
      lower_row  = _mm_load_si128(lower_row_ptr);

      // Retrieve value for next 8 elements

      upper_row_ptr++;
      middle_row_ptr++;
      lower_row_ptr++;

      // Load next 8 values
      __m128i next_upper_row, next_lower_row, next_middle_row;
      next_upper_row  = _mm_load_si128(upper_row_ptr);
      next_middle_row = _mm_load_si128(middle_row_ptr);
      next_lower_row  = _mm_load_si128(lower_row_ptr);

      // Defination of variables which will be used in the loop
      __m128i thrice_upper_row, thrice_middle_row, thrice_lower_row;
      __m128i hex_upper_row,    hex_middle_row,    hex_lower_row;
      __m128i oct_upper_row,    oct_middle_row,    oct_lower_row;

      __m128i next_thrice_upper_row, next_thrice_middle_row, next_thrice_lower_row;
      __m128i next_hex_upper_row,    next_hex_middle_row,    next_hex_lower_row;
      __m128i next_oct_upper_row,    next_oct_middle_row,    next_oct_lower_row;

      __m128i inter_value_2, inter_value_3;
      __m128i inter_value_6, inter_value_8, inter_value_9;

      __m128i upper_shift_once, upper_shift_twice, middle_shift_twice, lower_shift_once, lower_shift_twice;

      /*
          Iterate the array, we are calculating with 16 values at a time
          and every loop we get the values of 2 sums which will later be packed
          and store into the output file.
      */
      for (j = 1; j < width - 17; j = j + 16) {
        __m128i total_sum1, total_sum2;
        for (int k = 0; k < 2; k++) {

          // Calculate 3x the value for all 3 rows
          thrice_upper_row  = _mm_add_epi16(upper_row, _mm_add_epi16(upper_row, upper_row));
          thrice_middle_row = _mm_add_epi16(middle_row, _mm_add_epi16(middle_row, middle_row));
          thrice_lower_row  = _mm_add_epi16(lower_row, _mm_add_epi16(lower_row, lower_row));

          // Calculate 5x the value for all 5 rows
          hex_upper_row  = _mm_add_epi16(_mm_add_epi16(upper_row, upper_row), thrice_upper_row);
          hex_middle_row = _mm_add_epi16(_mm_add_epi16(middle_row, middle_row), thrice_middle_row);
          hex_lower_row  = _mm_add_epi16(_mm_add_epi16(lower_row, lower_row), thrice_lower_row);

          // Calculate 8x the value for all 8 rows
          oct_upper_row  = _mm_add_epi16(hex_upper_row, thrice_upper_row);
          oct_middle_row = _mm_add_epi16(hex_middle_row, thrice_middle_row);
          oct_lower_row  = _mm_add_epi16(hex_lower_row, thrice_lower_row);

          // Calculate 3x the value for all 3 rows of the next 8 values
          next_thrice_upper_row  = _mm_add_epi16(next_upper_row, _mm_add_epi16(next_upper_row, next_upper_row));
          next_thrice_middle_row = _mm_add_epi16(next_middle_row, _mm_add_epi16(next_middle_row, next_middle_row));
          next_thrice_lower_row  = _mm_add_epi16(next_lower_row, _mm_add_epi16(next_lower_row, next_lower_row));

          // Calculate 5x the value for all 5 rows of the next 8 values
          next_hex_upper_row  = _mm_add_epi16(_mm_add_epi16(next_upper_row, next_upper_row), next_thrice_upper_row);
          next_hex_middle_row = _mm_add_epi16(_mm_add_epi16(next_middle_row, next_middle_row), next_thrice_middle_row);
          next_hex_lower_row  = _mm_add_epi16(_mm_add_epi16(next_lower_row, next_lower_row), next_thrice_lower_row);

          // Calculate 8x the value for all 8 rows
          next_oct_upper_row  = _mm_add_epi16(next_hex_upper_row, next_thrice_upper_row);
          next_oct_middle_row = _mm_add_epi16(next_hex_middle_row, next_thrice_middle_row);
          next_oct_lower_row  = _mm_add_epi16(next_hex_lower_row, next_thrice_lower_row);

          /*
             8 versions of intermediate values that corresponds to the offsets between each matrix which
             will be used as addition or subtraction
             Intermediate values of 1, 4, 7 need not be calcualated as they refer to 
             oct_upper_row, oct_middle_row and oct_lower_row respectively which were calculated earlier
          */

          inter_value_2 = _mm_or_si128(_mm_slli_si128(next_oct_upper_row, 14), _mm_srli_si128(oct_upper_row, 2));
          inter_value_3 = _mm_or_si128(_mm_slli_si128(next_oct_upper_row, 12), _mm_srli_si128(oct_upper_row, 4));
          inter_value_6 = _mm_or_si128(_mm_slli_si128(next_oct_middle_row, 12), _mm_srli_si128(oct_middle_row, 4));
          inter_value_8 = _mm_or_si128(_mm_slli_si128(next_oct_lower_row, 14), _mm_srli_si128(oct_lower_row, 2));
          inter_value_9 = _mm_or_si128(_mm_slli_si128(next_oct_lower_row, 12), _mm_srli_si128(oct_lower_row, 4));

          // rotation 1
          /*
            { 5,  5,  5},
            {-3,  0, -3},
            {-3, -3, -3}
          */

          // oct_upper_row,                                                                                                     // 5u0
          upper_shift_once    = _mm_or_si128(_mm_slli_si128(next_hex_upper_row, 14), _mm_srli_si128(hex_upper_row, 2));         // 5u1
          upper_shift_twice   = _mm_or_si128(_mm_slli_si128(next_hex_upper_row, 12), _mm_srli_si128(hex_upper_row, 4));         // 5u2
          // oct_middle_row,                                                                                                    // 3m0
          middle_shift_twice  = _mm_or_si128(_mm_slli_si128(next_thrice_middle_row, 12), _mm_srli_si128(thrice_middle_row, 4)); // 3m2
          // oct_lower_row,                                                                                                     // 3l0
          lower_shift_once    = _mm_or_si128(_mm_slli_si128(next_thrice_lower_row, 14), _mm_srli_si128(thrice_lower_row, 2));   // 3l1
          lower_shift_twice   = _mm_or_si128(_mm_slli_si128(next_thrice_lower_row, 12), _mm_srli_si128(thrice_lower_row, 4));   // 3l2

          // We zero the starting value of the maximum sum
          __m128i max_sum = _mm_setzero_si128();

          /* 
              Perform the following operations to get the sum of the first rotation matrix
              It follows the following formula: 5u0 + 5u1 + 5u2 - 3m0 - 3m2 - 3l0 - 3l1 - 3l2
              After sum is calulated, we retrieve the max values between max_sum and rotation 
              matrix's sum and do this for the remaining 7 rotation matrices.
          */

          __m128i sum_rot_1 = _mm_add_epi16(hex_upper_row,  _mm_add_epi16(upper_shift_once, upper_shift_twice));
          sum_rot_1 = _mm_sub_epi16(_mm_sub_epi16(sum_rot_1, thrice_middle_row), middle_shift_twice);
          sum_rot_1 = _mm_sub_epi16(_mm_sub_epi16(sum_rot_1, thrice_lower_row), lower_shift_once);
          sum_rot_1 = _mm_sub_epi16(sum_rot_1, lower_shift_twice);
          max_sum = _mm_max_epi16(max_sum, sum_rot_1); // Get max

          /*
              From here on, the subsequent sum of the remaining matrices are calclated through
              addition and subtraction from the previous sum by the offsets
          */

          // rotation 2

          __m128i sum_rot_2 = _mm_add_epi16(oct_middle_row, _mm_sub_epi16(sum_rot_1, inter_value_3));
          max_sum = _mm_max_epi16(max_sum, sum_rot_2);

          // rotation 3

          __m128i sum_rot_3 = _mm_add_epi16(oct_lower_row, _mm_sub_epi16(sum_rot_2, inter_value_2));
          max_sum = _mm_max_epi16(max_sum, sum_rot_3); 

          // rotation 4

          __m128i sum_rot_4 = _mm_add_epi16(inter_value_8, _mm_sub_epi16(sum_rot_3, oct_upper_row));
          max_sum = _mm_max_epi16(max_sum, sum_rot_4); 

          // rotation 5

          __m128i sum_rot_5 = _mm_add_epi16(inter_value_9, _mm_sub_epi16(sum_rot_4, oct_middle_row));
          max_sum = _mm_max_epi16(max_sum, sum_rot_5); 

          // rotation 6

          __m128i sum_rot_6 = _mm_add_epi16(inter_value_6, _mm_sub_epi16(sum_rot_5, oct_lower_row));
          max_sum = _mm_max_epi16(max_sum, sum_rot_6); 

          // rotation 7

          __m128i sum_rot_7 = _mm_add_epi16(inter_value_3, _mm_sub_epi16(sum_rot_6, inter_value_8));
          max_sum = _mm_max_epi16(max_sum, sum_rot_7);

          // rotation 8

          __m128i sum_rot_8 = _mm_add_epi16(inter_value_2, _mm_sub_epi16(sum_rot_7, inter_value_9));
          max_sum = _mm_max_epi16(max_sum, sum_rot_8);
          
          max_sum = _mm_srai_epi16(max_sum, 3);                   // divide each by 8
          max_sum = _mm_max_epi16(max_sum, _mm_setzero_si128());  // max(sum, 0)
          max_sum = _mm_min_epi16(max_sum, _mm_set1_epi16(255));  // min(sum, 255)

          // Store in sum1 and sum2 respectively
          if (k == 0)
            total_sum1 = max_sum;
          else
            total_sum2 = max_sum;

          // Set the variables to prepare for the next loop
          upper_row  = next_upper_row;
          middle_row = next_middle_row;
          lower_row  = next_lower_row;

          upper_row_ptr++;
          middle_row_ptr++;
          lower_row_ptr++;

          next_upper_row = _mm_load_si128(upper_row_ptr);
          next_middle_row = _mm_load_si128(middle_row_ptr);
          next_lower_row = _mm_load_si128(lower_row_ptr);

        }

        /*
        Packing the two rows of values we have computed and write into the
        data out buffer.
        */

        __m128i total_sum =
          _mm_packus_epi16(total_sum1, total_sum2);
        vec_or_val real_sum;
          _mm_store_si128(&real_sum.value, total_sum);

        memcpy_s(&data_out[lay * size + i * width + j], 16, real_sum.vector, 16);

      }


		/* Cleanup for the remaining pixels that did not get processed earlier */
		//This is the clean up part! 
		//The image size might not be multiple of 16, so we need to deal with the remaining pixels manually.
		//Notice that this code here is pretty much copy pasted from the basis algorithm

      /* Kirsch matrices for convolution */
      static int kirsch[8][3][3] = 
      {
        {{  5,  5,  5 }, { -3,  0, -3 }, { -3, -3, -3 }}, /*rotation 1 */
        {{  5,  5, -3 }, {  5,  0, -3 }, { -3, -3, -3 }}, /*rotation 2 */
        {{  5, -3, -3 }, {  5,  0, -3 }, {  5, -3, -3 }}, /*rotation 3 */
        {{ -3, -3, -3 }, {  5,  0, -3 }, {  5,  5, -3 }}, /*rotation 4 */
        {{ -3, -3, -3 }, { -3,  0, -3 }, {  5,  5,  5 }}, /*rotation 5 */
        {{ -3, -3, -3 }, { -3,  0,  5 }, { -3,  5,  5 }}, /*rotation 6 */
        {{ -3, -3,  5 }, { -3,  0,  5 }, { -3, -3,  5 }}, /*rotation 7 */
        {{ -3,  5,  5 }, { -3,  0,  5 }, { -3, -3, -3 }}  /*rotation 8 */
      
      };


      for (; j < width - 1; j++) {
        int max_sum;
        max_sum = 0;

        /*
        Perform convolutions for
        all 8 masks in succession. Compare and find the one
        that has the highest value. The one with the
        highest value is stored into the final bitmap.
        */
        for (unsigned m = 0; m < 8; ++m) {
          int sum;
          sum = 0;
          /* Convolution part */
          for (int k = -1; k < 2; k++)
            for (int l = -1; l < 2; l++) {
              sum =
                sum + kirsch[m][k + 1][l + 1] *
                (int)data_in[lay * size +
                (i + k) * width + (j + l)];
            }
          if (sum > max_sum)
            max_sum = sum;
        }
        data_out[lay * size + i * width + j] = min(max(max_sum / 8, 0), 255);
      }
    }
  }
}

//original code in c 
void
kirsch_operator_basic(unsigned char *data_out,
                      unsigned char *data_in,
                      unsigned height, unsigned width)
{
    /* Kirsch matrices for convolution */
    int kirsch[8][3][3] = {
        {
         {5, 5, 5},
         {-3, 0, -3},           /*rotation 1 */
         {-3, -3, -3}
         },
        {
         {5, 5, -3},
         {5, 0, -3},            /*rotation 2 */
         {-3, -3, -3}
         },
        {
         {5, -3, -3},
         {5, 0, -3},            /*rotation 3 */
         {5, -3, -3}
         },
        {
         {-3, -3, -3},
         {5, 0, -3},            /*rotation 4 */
         {5, 5, -3}
         },
        {
         {-3, -3, -3},
         {-3, 0, -3},           /*rotation 5 */
         {5, 5, 5}
         },
        {
         {-3, -3, -3},
         {-3, 0, 5},            /*rotation 6 */
         {-3, 5, 5}
         },
        {
         {-3, -3, 5},
         {-3, 0, 5},            /*rotation 7 */
         {-3, -3, 5}
         },
        {
         {-3, 5, 5},
         {-3, 0, 5},            /*rotation 8 */
         {-3, -3, -3}
         }
    };

    unsigned int size, y, x, lay;

    size = height * width;

    for (lay = 0; lay < 3; lay++) {
        for (y = 1; y < height - 1; ++y) {
            for (x = 1; x < width - 1; x++) {
                int max_sum;
                max_sum = 0;

                /*
                   Perform convolutions for 
                   all 8 masks in succession. Compare and find the one
                   that has the highest value. The one with the
                   highest value is stored into the final bitmap.
                 */
                for (unsigned m = 0; m < 8; ++m) {
                    int sum;
                    sum = 0;
                    /* Convolution part */
                    for (int k = -1; k < 2; k++)
                        for (int l = -1; l < 2; l++) {
                            sum =
                                sum + kirsch[m][k + 1][l + 1] *
                                (int) data_in[lay * size +
                                              (y + k) * width + (x + l)];
                        }
                    if (sum > max_sum)
                        max_sum = sum;
                }
				max_sum = max_sum/8 > 255 ? 255: max_sum/8;
				data_out[lay * size + y * width + x] = max_sum;
            }
        }
    }
}

