#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stddef.h>

#include "RPS_MPI.h"
#include <time.h>

void initialize();
void exchange_borders(bool inverted);
void iterate_CA(int i);
void gather_petri();
void create_types();
void free_stuff();


int rank;
int size;

// I denote mpi process specific values with hungarian notation, adding a p

// The dimensions of the processor grid. Same for every process
int p_x_dims;
int p_y_dims;

// The location of a process in the process grid. Unique for every process
int p_my_x_dim;
int p_my_y_dim;

int p_north, p_south, p_east, p_west;

// The dimensions for the process local petri
int p_local_petri_x_dim;
int p_local_petri_y_dim;

//local buffer size
int p_local_petri_size;

MPI_Comm cart_comm;

// some datatypes, useful for sending data with somewhat less primitive semantics
MPI_Datatype border_row_t;  // TODO: Implement this
MPI_Datatype border_col_t;  // TODO: Implement this
MPI_Datatype local_petri_t; // Already implemented
MPI_Datatype mpi_cell_t;    // Already implemented
MPI_Datatype inner_petri;
MPI_Datatype final_petri_part;

// Each process is responsible for one part of the petri dish.
// Since we can't update the petri-dish in place each process actually
// gets two petri-dishes which they update in a lockstep fashion.
// dish A is updated by writing to dish B, then next step dish B updates dish A.
// (or you can just swap them inbetween iterations)
cell* local_petri_A;
cell* local_petri_B;

cell* final_petri;

int main(int argc, char** argv){

  srand(1234);

  clock_t before = clock();

  // Ask MPI what size (number of processors) and rank (which process we are)
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ////////////////////////////////
  // Create cartesian communicator
  int dims[2];
  dims[0] = p_y_dims;
  dims[1] = p_x_dims;

  int periods[2]; // we set these to 0 because we are not interested in wrap-around
  periods[0] = 0;
  periods[1] = 0;

  int coords[2];
  coords[0] = p_my_y_dim;
  coords[1] = p_my_x_dim;

  MPI_Dims_create(size, 2, dims);
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
  MPI_Cart_coords(cart_comm, rank, 2, coords);

  MPI_Cart_shift(cart_comm, 0, 1, &p_north, &p_south);
  MPI_Cart_shift(cart_comm, 1, 1, &p_west, &p_east);

  p_y_dims = dims[0];
  p_x_dims = dims[1];

  p_my_y_dim = coords[0];
  p_my_x_dim = coords[1];

  ////////////////////////////////
  ////////////////////////////////

  initialize();

  create_types();

  for(int i = 0; i < ITERATIONS; i++){
    iterate_CA(i);
  }

  gather_petri();

  MPI_Finalize();

  if(rank==0){
    clock_t after = clock();
    par_make_bmp(final_petri);
    printf("Time taken: %ldms\n", (after - before)* 1000 / CLOCKS_PER_SEC);
  }

  // You should probably make sure to free your memory here
  // We will dock points for memory leaks, don't let your hard work go to waste!
  // free_stuff()
  free_stuff();

  exit(0);
}

void free_stuff(){
  free(final_petri);
  free(local_petri_B);
  free(local_petri_A);
}


void create_types(){

  ////////////////////////////////
  ////////////////////////////////
  // cell type
  const int    nitems=2;
  int          blocklengths[2] = {1,1};
  MPI_Datatype types[2] = {MPI_INT, MPI_INT};
  MPI_Aint     offsets[2];

  offsets[0] = offsetof(cell, color);
  offsets[1] = offsetof(cell, strength);

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cell_t);
  MPI_Type_commit(&mpi_cell_t);
  ////////////////////////////////
  ////////////////////////////////



  ////////////////////////////////
  ////////////////////////////////
  // A message for a local petri-dish
  MPI_Type_contiguous(p_local_petri_x_dim * p_local_petri_y_dim,
                      mpi_cell_t,
                      &local_petri_t);
  MPI_Type_commit(&local_petri_t);
  ////////////////////////////////
  ///////////////////////////////

  //TODO: Create MPI types for border exchange
  //border_row_t, border_col_t
  MPI_Type_vector(p_local_petri_x_dim, 1, 1, mpi_cell_t, &border_row_t);
  MPI_Type_commit(&border_row_t);

  MPI_Type_vector(p_local_petri_y_dim, 1, p_local_petri_x_dim + 2, mpi_cell_t, &border_col_t);
  MPI_Type_commit(&border_col_t);


  MPI_Type_vector(p_local_petri_y_dim, p_local_petri_x_dim, p_local_petri_x_dim+2, mpi_cell_t, &inner_petri);
  MPI_Type_commit(&inner_petri);

  MPI_Type_vector(p_local_petri_y_dim, p_local_petri_x_dim, IMG_X, mpi_cell_t, &final_petri_part);
  MPI_Type_commit(&final_petri_part);

}


void initialize(){
  //TODO: assign the following to something more useful than 0

  //Calculate the dimensions of the local petri

  p_local_petri_x_dim = IMG_X/p_x_dims;
  p_local_petri_y_dim = IMG_Y/p_y_dims;

  //Calculate total size of local petri with the borders

  p_local_petri_size = (p_local_petri_x_dim+2)*(p_local_petri_y_dim+2);

  // TODO: When allocating these buffers, keep in mind that you might need to allocate a little more
  // than just your piece of the petri.

  //Allocate the full local petri with the borders also

  local_petri_A = malloc(p_local_petri_size*sizeof(cell));
  local_petri_B = malloc(p_local_petri_size*sizeof(cell));

  //Set seed values

  for(int ii = 0; ii < 100; ii++){

    int currX = (rand() % p_local_petri_x_dim)+1;
    int currY = (rand() % p_local_petri_y_dim)+1;

    int rx = currY*(p_local_petri_x_dim+2) + currX;

    int rt = rand() % 4;

    local_petri_A[rx].color = rt;
    local_petri_A[rx].strength = 1;

  }
}

void exchange_borders(bool inverted){

  //Exchange borders with the correct neighbors

  if(inverted){
    if(p_north > -1){

      MPI_Send(&local_petri_B[p_local_petri_x_dim + 3], 1, border_row_t, p_north, 0, MPI_COMM_WORLD);
      MPI_Recv(&local_petri_A[1], 1, border_row_t, p_north, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(p_south > -1){

      MPI_Send(&local_petri_B[p_local_petri_size - 3 - p_local_petri_x_dim*2], 1, border_row_t, p_south, 0, MPI_COMM_WORLD);
      MPI_Recv(&local_petri_A[p_local_petri_size - 1 - p_local_petri_x_dim], 1, border_row_t, p_south, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(p_west > -1){

      MPI_Send(&local_petri_B[p_local_petri_x_dim + 3], 1, border_col_t, p_west, 0, MPI_COMM_WORLD);
      MPI_Recv(&local_petri_A[p_local_petri_x_dim + 2], 1, border_col_t, p_west, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(p_east > - 1){

      MPI_Send(&local_petri_B[2*p_local_petri_x_dim + 2], 1, border_col_t, p_east, 0, MPI_COMM_WORLD);
      MPI_Recv(&local_petri_A[2*p_local_petri_x_dim + 3], 1, border_col_t, p_east, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
  else{
    if(p_north > -1){

      MPI_Send(&local_petri_A[p_local_petri_x_dim + 3], 1, border_row_t, p_north, 0, MPI_COMM_WORLD);
      MPI_Recv(&local_petri_B[1], 1, border_row_t, p_north, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(p_south > -1){

      MPI_Send(&local_petri_A[p_local_petri_size - 3 - p_local_petri_x_dim*2], 1, border_row_t, p_south, 0, MPI_COMM_WORLD);
      MPI_Recv(&local_petri_B[p_local_petri_size - 1 - p_local_petri_x_dim], 1, border_row_t, p_south, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(p_west > -1){

      MPI_Send(&local_petri_A[p_local_petri_x_dim + 3], 1, border_col_t, p_west, 0, MPI_COMM_WORLD);
      MPI_Recv(&local_petri_B[p_local_petri_x_dim + 2], 1, border_col_t, p_west, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(p_east > - 1){

      MPI_Send(&local_petri_A[2*p_local_petri_x_dim + 2], 1, border_col_t, p_east, 0, MPI_COMM_WORLD);
      MPI_Recv(&local_petri_B[2*p_local_petri_x_dim + 3], 1, border_col_t, p_east, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
}

//Iterate the CA, calculate next image, if iteration is a square number, send from A to B, else send from B to A

void iterate_CA(int i){
  if(i%2 == 0){
    par_iterate_image(local_petri_A,local_petri_B, p_local_petri_size, p_local_petri_x_dim + 2, p_local_petri_y_dim + 2);
    exchange_borders(false);
  }
  else{
    par_iterate_image(local_petri_B,local_petri_A, p_local_petri_size, p_local_petri_x_dim + 2, p_local_petri_y_dim + 2);
    exchange_borders(true);
  }
}

//Gather petris

void gather_petri(){

  if(rank == 0) {

    int j = 0;
    final_petri = malloc(IMG_X*IMG_Y*sizeof(cell));

    for(int i = 0; i < p_local_petri_size; i++){
      if(!is_border(i, p_local_petri_x_dim + 2, p_local_petri_y_dim + 2)){

        //Fancy calculation to find right position of local_petri in main petri (for 0)

        final_petri[j%p_local_petri_x_dim + IMG_X*(j/p_local_petri_x_dim)] = local_petri_B[i];

        j++;
      }
    }

    for(int i = 1; i < size; i++) {
      MPI_Recv(&final_petri[(p_local_petri_x_dim * p_local_petri_y_dim * p_x_dims * (i/p_x_dims)) + (p_local_petri_x_dim * (i%p_x_dims))], 1, final_petri_part, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  }
  else {
    MPI_Send(&local_petri_B[p_local_petri_x_dim + 3], 1, inner_petri, 0, 0, MPI_COMM_WORLD);
  }
}
