/* sobel.c */
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "mypgm.h"
#include <omp.h>

#define num_threads 4 


void sobel_filtering( )
     /* Spatial filtering of image data */
     /* Sobel filter (horizontal differentiation */
     /* Input: image1[y][x] ---- Outout: image2[y][x] */
{
  /* Definition of Sobel filter in horizontal direction */
  int weight[3][3] = {{ -1,  0,  1 },
		      { -2,  0,  2 },
		      { -1,  0,  1 }};

  omp_set_num_threads(num_threads);


  double pixel_value;
  double min, max;
  int x, y, i, j;  /* Loop variable */
  
  /* Maximum values calculation after filtering*/
  printf("Now, filtering of input image is performed\n\n");
  min = DBL_MAX;
  max = -DBL_MAX;

  printf("The ysize and xsize is %d and %d", y_size1,x_size1 );
  pixel_value = 0.0;
  
  #pragma omp parallel for collapse(2) private(x,y,i,j,pixel_value)\
                                       shared(min,max) \
                                       
  for (y = 1; y < y_size1 - 1; y++) {
    for (x = 1; x < x_size1 - 1; x++) {
        pixel_value = 0.0;
        for (j = -1; j <= 1; j++) {
	    for (i = -1; i <= 1; i++) {
	      pixel_value += weight[j + 1][i + 1] * image1[y + j][x + i];
	    }
      }
            if (pixel_value < min) min = pixel_value;
            if (pixel_value > max) max = pixel_value;
    }
  }


  if ((int)(max - min) == 0) {
    printf("Nothing exists!!!\n\n");
    exit(1);
  }

  /* Initialization of image2[y][x] */
  x_size2 = x_size1;
  y_size2 = y_size1;

  #pragma omp parallel for collapse(2) private(x,y)\
                                       shared(image2)
  for (y = 0; y < y_size2; y++) {
    for (x = 0; x < x_size2; x++) {
      image2[y][x] = 0;
    }
  }

   pixel_value = 0;
  /* Generation of image2 after linear transformtion */
  #pragma omp parallel for collapse(2) private(x,y,i,j,pixel_value)\
                                       shared(image2)
  for (y = 1; y < y_size1 - 1; y++) {
    for (x = 1; x < x_size1 - 1; x++) {
      pixel_value = 0.0;
      for (j = -1; j <= 1; j++) {
	    for (i = -1; i <= 1; i++) {
	      pixel_value += weight[j + 1][i + 1] * image1[y + j][x + i];
	    }
      }
      pixel_value = MAX_BRIGHTNESS * (pixel_value - min) / (max - min);
      image2[y][x] = (unsigned char)pixel_value;
    }
  }
}

main( )
{    
  load_image_data( );   /* Input of image1 */ 
  double start_time = omp_get_wtime();
  sobel_filtering( );   /* Sobel filter is applied to image1 */
  printf("The total time for filering is %f", omp_get_wtime() - start_time);
  save_image_data( );   /* Output of image2 */ 
  return 0;
}
