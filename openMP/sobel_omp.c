/* sobel.c */
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "mypgm.h"
#include <omp.h>

#define num_threads 12

void sobel_filtering( )
     /* Spatial filtering of image data */
     /* Sobel filter (horizontal differentiation */
     /* Input: image1[y][x] ---- Outout: image2[y][x] */
{
  /* Definition of Sobel filter in horizontal direction */
  int weight[3][3] = {{ -1,  0,  1 },
		      { -2,  0,  2 },
		      { -1,  0,  1 }};
  double min[num_threads], max[num_threads];
  double min_value, max_value;
  /* Maximum values calculation after filtering*/
  printf("Now, filtering of input image is performed\n\n");
  min_value = DBL_MAX;
  max_value = -DBL_MAX;
  
  int v;
  for(v = 0; v < num_threads;v++)
  {
	min[v] = DBL_MAX;
	max[v] = -DBL_MAX;
  }
  
  #pragma omp parallel
  {
 
  double pixel_value;
  
  int x, y, i, j;  /* Loop variable */
  int thr_id = omp_get_num_threads();
  for (y = (1+thr_id); y < y_size1 - 1; y= y+num_threads) {
    for (x = (1+thr_id); x < x_size1 - 1; x=x+num_threads) {
      pixel_value = 0.0;
      for (j = -1; j <= 1; j++) {
	    for (i = -1; i <= 1; i++) {
	      pixel_value += weight[j + 1][i + 1] * image1[y + j][x + i];
	    }
      }
      if (pixel_value < min[thr_id]) min[thr_id] = pixel_value;
      if (pixel_value > max[thr_id]) max[thr_id] = pixel_value;
    }
  }
  }
  
  int m,x,y,i,j;
  double pixel_value;

  for(m = 0;m<num_threads;m++)
  {
      if (min[m] < min_value) min_value = min[m];
      if (max[m] > max_value) max_value = max[m];
  }
  
  printf("The Value of Minimum Shared  %f",min_value);
  if ((int)(max - min) == 0) {
    printf("Nothing exists!!!\n\n");
    exit(1);
  }

  /* Initialization of image2[y][x] */
  x_size2 = x_size1;
  y_size2 = y_size1;
  for (y = 0; y < y_size2; y++) {
    for (x = 0; x < x_size2; x++) {
      image2[y][x] = 0;
    }
  }
  /* Generation of image2 after linear transformtion */
  for (y = 1; y < y_size1 - 1; y++) {
    for (x = 1; x < x_size1 - 1; x++) {
      pixel_value = 0.0;
      for (j = -1; j <= 1; j++) {
	    for (i = -1; i <= 1; i++) {
	      pixel_value += weight[j + 1][i + 1] * image1[y + j][x + i];
	    }
      }
      pixel_value = MAX_BRIGHTNESS * (pixel_value - min_value) / (max_value - min_value);
      image2[y][x] = (unsigned char)pixel_value;
    }
  }
}

main( )
{
  load_image_data( );   /* Input of image1 */ 
  sobel_filtering( );   /* Sobel filter is applied to image1 */
  save_image_data( );   /* Output of image2 */
  return 0;
}
