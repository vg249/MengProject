#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

int malloc2dchar(int ***array, int n, int m) {
int i;
        /* allocate the n*m contiguous items */
        int *p = (int *)malloc(n*m*sizeof(int));
            if (!p) return -1;

                /* allocate the row pointers into the memory */
                (*array) = (int **)malloc(n*sizeof(int*));
                    if (!(*array)) {
                               free(p);
                                      return -1;
                                          }

                        /* set up the pointers into the contiguous memory */
                        for (i=0; i<n; i++)
                                   (*array)[i] = &(p[i*m]);

                            return 0;
}

int malloc2dint(int ***array, int n, int m) {
int i;
        /* allocate the n*m contiguous items */
        int *p = (int *)malloc(n*m*sizeof(int));
            if (!p) return -1;

                /* allocate the row pointers into the memory */
                (*array) = (int **)malloc(n*sizeof(int*));
                    if (!(*array)) {
                               free(p);
                                      return -1;
                                          }

                        /* set up the pointers into the contiguous memory */
                        for (i=0; i<n; i++)
                                   (*array)[i] = &(p[i*m]);

                            return 0;
}

int free2dint(int ***array) {
        /* free the memory - the first element of the array is at the start */
        free(&((*array)[0][0]));

            /* free the pointers into the memory */
            free(*array);

                return 0;
}
int free2dchar(char ***array) {
        /* free the memory - the first element of the array is at the start */
        free(&((*array)[0][0]));

            /* free the pointers into the memory */
            free(*array);

                return 0;
}

int main(int argc, char **argv) {
        int **global, **local;
        int **local1;
        int i ,j;
            const int gridsize=8; // size of grid
                const int procgridsize=4;  // size of process grid
                    int rank, size;        // rank of current process and no. of processes

                        MPI_Init(&argc, &argv);
                            MPI_Comm_size(MPI_COMM_WORLD, &size);
                                MPI_Comm_rank(MPI_COMM_WORLD, &rank);


                                    if (size != procgridsize) {
                                                fprintf(stderr,"%s: Only works with np=%d for now\n", argv[0], procgridsize);
                                                        MPI_Abort(MPI_COMM_WORLD,1);
                                                            }
int d = 0;

if (rank == 0) {
            /* fill in the array, and prit */
malloc2dint(&global, gridsize, gridsize);
for (i=0; i<gridsize; i++) {
for (j=0; j<gridsize; j++){
    global[i][j] = d;
    d++;                        
}
}


printf("Global array is:\n");
        for (i=0; i<gridsize; i++) {
        for (j=0; j<gridsize; j++){
                    printf("%d\t",global[i][j]);
        }
        printf("\n");
    }
}
                            /* create the local array which we'll process */
                            malloc2dint(&local, (gridsize/procgridsize)+1, gridsize);
                            malloc2dint(&local1, (gridsize/procgridsize)+1, gridsize);

/* create a datatype to describe the subarrays of the global array */

        int sizes[2]    = {gridsize, gridsize};         /* global size */
        int subsizes[2] = {(gridsize/procgridsize)+1, gridsize};     /* local size */
        int starts[2]   = {0,0};                        /* where this one starts */
        MPI_Datatype type, subarrtype;
        MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &type);
        MPI_Type_create_resized(type, 0, (gridsize)*sizeof(int), &subarrtype);
        MPI_Type_commit(&subarrtype);

        int sizes1[2]    = {8, 8};         /* global size */
        int subsizes1[2] = {4, 4};     /* local size */
        int starts1[2]   = {1,1};                        /* where this one starts */
        MPI_Datatype type1, subarrtype1;
        MPI_Type_create_subarray(2, sizes1, subsizes1, starts1, MPI_ORDER_C, MPI_INT, &type1);
        MPI_Type_create_resized(type1, 0, 2*sizeof(int), &subarrtype1);
        MPI_Type_commit(&subarrtype1);
  
int *globalptr=NULL;
if (rank == 0) globalptr = &(global[0][0]);

    /* scatter the array to all processors */
    int sendcounts[procgridsize*procgridsize];
        int displs[procgridsize*procgridsize];
         int k;
        if (rank == 0) {
            int dest = 0;
            for ( k=0; k<size; k++) {

            sendcounts[k] = 1;
            displs[k] = dest; 
            dest += ((gridsize/procgridsize) - 1);
        }
    }

MPI_Scatterv(globalptr, sendcounts, displs, subarrtype, &(local[0][0]), ((gridsize/procgridsize)+1)*(gridsize), MPI_INT, 0, MPI_COMM_WORLD);

                    /* now all processors print their local data: */
                    int p = 0;
                    for (p=0; p<size; p++) {
                                if (rank == p) {
            printf("Local array is %d:\n",p);
            for (i=0; i<(gridsize/procgridsize)+1; i++) {
                                printf("|");
            for (j=0; j<(gridsize); j++) {
            printf("%d\t",local[i][j]);
                                                                                                                                                                        }                                                                                                                                                                                                printf("|\n");
                                                                                                                                                                                                            }
MPI_Barrier(MPI_COMM_WORLD);
                                }
                    }
                                                                                                        /* now each processor has its local array, and can process it */
        for (i=0; i<gridsize/procgridsize; i++) {
                    for (j=0; j<gridsize/procgridsize; j++) {
                                    local1[i][j] = local[i][j];
                                            }
                        }

            /* it all goes back to process 0 */
            MPI_Gatherv(&(local1[1][1]), 4,  MPI_CHAR,
                                     globalptr, sendcounts, displs, subarrtype,
                                                      0, MPI_COMM_WORLD);

                /* don't need the local data anymore */
                free2dint(&local);

                    /* or the MPI data type */
                    MPI_Type_free(&subarrtype);

                        if (rank == 0) {
                                    printf("Processed grid:\n");
                                            for (i=0; i<gridsize; i++) {
                                                            for (j=0; j<gridsize; j++) {
                                                                                printf("%d\t",global[i][j]);
                                                                                            }
                                                                        printf("\n");
                                                                                }

                                                    free2dint(&global);
                                                                                                                                                        }


                                                                                                                            MPI_Finalize();

                                                                                                                                return 0;
}
