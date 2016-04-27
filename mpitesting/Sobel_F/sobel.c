/* sobel.c */
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "mypgm.h"
#include "mpi.h"

#define MASTER 0


unsigned char** imageLocal;
unsigned char** imageOpLocal;
int global_x_size, global_y_size;

//Create a local array
int malloc2dchar(unsigned char ***array, int n, int m) {
    
    int i = 0;
    unsigned char *p = (unsigned char *)malloc(n*m*sizeof(unsigned char));
    if (!p) return -1;

    (*array) = (unsigned char **)malloc(n*sizeof(unsigned char*));
    if (!(*array)) {
    free(p);
    return -1;
                                  }
    for (i=0; i<n; i++)
    (*array)[i] = &(p[i*m]);
    return 0;
}

int free2dchar(unsigned char ***array) {
        /* free the memory - the first element of the array is at the start */
        free(&((*array)[0][0]));

            /* free the pointers into the memory */
            free(*array);

                return 0;
}

void sobel_filtering(int rank)
     /* Spatial filtering of image data */
     /* Sobel filter (horizontal differentiation */
     /* Input: image1[y][x] ---- Outout: image2[y][x] */
{
  /* Definition of Sobel filter in horizontal direction */
  int weight[3][3] = {{ -1,  0,  1 },
		      { -2,  0,  2 },
		      { -1,  0,  1 }};
  double pixel_value;
  double min, max;
  int x, y, i, j;  /* Loop variable */
  /* Maximum values calculation after filtering*/
  printf("Now, filtering of input image is performed\n\n");
  min = DBL_MAX;
  max = -DBL_MAX;
  for (y = 1; y < ((global_y_size+2)/2 - 1); y++) {
    for (x = 1; x < ((global_x_size+2)/2 - 1); x++) {
      pixel_value = 0.0;
      for (j = -1; j <= 1; j++) {
	    for (i = -1; i <= 1; i++) {
	      pixel_value += weight[j + 1][i + 1] * imageLocal[y + j][x + i];
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

  for (y = 0; y < (global_y_size + 2)/2; y++) {
    for (x = 0; x < (global_y_size +2)/2; x++) {
      imageOpLocal[y][x] = 0;
    }
  }

  /* Generation of image2 after linear transformtion */
  for (y = 1; y < ((global_y_size+2)/2 - 1); y++) {
    for (x = 1; x < ((global_x_size+2)/2 - 1); x++) {
      pixel_value = 0.0;
      for (j = -1; j <= 1; j++) {
	    for (i = -1; i <= 1; i++) {
	      pixel_value += weight[j + 1][i + 1] * imageLocal[y + j][x + i];
	    }
      }
      pixel_value = MAX_BRIGHTNESS * (pixel_value - min) / (max - min);
      imageOpLocal[y][x] = (unsigned char)pixel_value;
    }
  }

}

int main(int argc, char** argv)
{
  double start_time, time;
  int i, j, disp;
  int rank, size;
  unsigned char* ipImageptr = NULL;
  unsigned char* opImageptr = NULL;
  int totPixCount;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
  ipImageptr = NULL;
  opImageptr = NULL;

  if(rank == MASTER)
  {
    load_image_file("image25.pgm");      
    global_x_size = x_size1;
    global_y_size = y_size1;
    ipImageptr = &(image1[0][0]);
    opImageptr = &(image2[0][0]);
  }
  

  //Broadcast the image size to all the nodes

  MPI_Bcast(&global_x_size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&global_y_size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);

  totPixCount   = (global_x_size+2)*(global_y_size+2); 
  
  int ipImageDim[2] = {(global_y_size+2), (global_x_size+2)};
  int scatterImageDim[2] = {(global_y_size+2)/2, (global_x_size+2)/2};
  int scatterStart[2]   = {0,0};
  int sendcount[4];
  int displ[4];

  int opImageDim[2] = {(global_y_size), (global_x_size)};
  int gatherImageDim[2] = {(global_y_size)/2, (global_x_size)/2};
  int gatherStart[2]   = {1,1};
  int recvcount[4];
  int recvdispl[4];
  
  //Scatter dimensions for input image from Master process

  MPI_Datatype scatter_type, scatter_subarraytype;
  MPI_Type_create_subarray(2, ipImageDim, scatterImageDim, scatterStart, MPI_ORDER_C, MPI_UNSIGNED_CHAR, &scatter_type);
  MPI_Type_create_resized(scatter_type, 0, (global_y_size+2)/2, &scatter_subarraytype);
  MPI_Type_commit(&scatter_subarraytype);

  //Gather dimesions for output image from all the porcess

  MPI_Datatype gather_type, gather_subarraytype;
  MPI_Type_create_subarray(2, ipImageDim, scatterImageDim, scatterStart, MPI_ORDER_C, MPI_UNSIGNED_CHAR, &gather_type);
  MPI_Type_create_resized(gather_type, 0, (global_y_size)/2, &gather_subarraytype);
  MPI_Type_commit(&gather_subarraytype);

    
  malloc2dchar(&imageLocal, (global_x_size+2)/2, (global_y_size+2)/2); 
  malloc2dchar(&imageOpLocal, (global_x_size+2)/2, (global_y_size+2)/2); 

  if(rank == MASTER)
  {
    for(i = 0; i < 4; i++) {
        sendcount[i] = 1;
        recvcount[i] = 1;
    }
      
    for(i = 0; i < 2; i++) {
        disp = 0;
        for(j = 0; j < 2; j++) {
            displ[(i*2)+j] = disp;    
            disp += 1;
        }
        disp += (((global_x_size+2)/2) - 1)*2;
    } 
  
    
    for(i = 0; i < 2; i++) {
        disp = 0;
        for(j = 0; j < 2; j++) {
            recvdispl[(i*2)+j] = disp;    
            disp += 1;
        }
        disp += ((global_x_size/2) - 1)*2;
    }
}


  MPI_Scatterv(ipImageptr, sendcount, displ, scatter_subarraytype, &(imageLocal[0][0]), totPixCount/4, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);  
  
  
  MPI_Barrier(MPI_COMM_WORLD);
  sobel_filtering(rank);
  
  
  printf("Process %d  \n\n",rank);
  
  MPI_Gatherv(&(imageOpLocal[1][1]), (global_x_size*global_y_size)/4, MPI_UNSIGNED_CHAR, &(image2[0][0]), recvcount, recvdispl, gather_subarraytype ,0, MPI_COMM_WORLD);  
  //MPI_Gatherv(&(imageOpLocal[0][0]), totPixCount/4, MPI_UNSIGNED_CHAR, &(image2[0][0]), sendcount, displ, scatter_subarraytype, 0, MPI_COMM_WORLD);  
 
  free2dchar(&imageLocal);

  MPI_Type_free(&scatter_subarraytype);
  MPI_Type_free(&gather_subarraytype);
  
  if(rank == MASTER){   
    x_size2 = global_x_size+2;
    y_size2 = global_y_size+2;
    save_image_file("image25_par.pgm");   /* Output of image2 */
  }

  MPI_Finalize();

}
