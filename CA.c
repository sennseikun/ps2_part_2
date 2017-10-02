#include "CA.h"

cell pick_neighbor(int x, int y, cell** image);
int par_pick_neighbor(int rows, int cols, int pos);
int MPI_pick_neighbor(int position, int width, int height);

cell** alloc_img(cell* buffer, int index){
cell** image = malloc(IMG_X*sizeof(cell*));

  for(int ii = 0; ii < IMG_X; ii++){
    image[ii] = &buffer[(IMG_X*IMG_Y*index) + IMG_X*ii];
  }

  return image;
}

void free_img(cell** image){
  for (int ii = 0; ii < IMG_X; ii++){
    free(image[ii]);
  }
  free(image);
}


bool beats(cell me, cell other){
  return
    (((me.color == SCISSOR) && (other.color == PAPER)) ||
     ((me.color == PAPER) && (other.color == ROCK))    ||
     ((me.color == ROCK) && (other.color == SCISSOR))  ||
     (me.color == other.color));
}

cell next_cell(int x, int y, cell** image){

  cell neighbor_cell = pick_neighbor(x, y, image);
  cell my_cell = image[x][y];
  if(neighbor_cell.color == WHITE){
    return my_cell;
  }
  cell next_cell = my_cell;

  if(my_cell.color == WHITE){
    next_cell.strength = 1;
    next_cell.color = neighbor_cell.color;
    return next_cell;
  }
  else {
    if(beats(my_cell, neighbor_cell)){
      next_cell.strength++;
    }
    else{
      next_cell.strength--;
    }
  }

  if(next_cell.strength == 0){
    next_cell.color = neighbor_cell.color;
    next_cell.strength = 1;
  }

  if(next_cell.strength > 4){
    next_cell.strength = 4;
  }

  return next_cell;
}

//Parallel version of the code picking the next cell
//I split it up to allow having the serial code working as well

cell par_next_cell(int rows, int cols, int pos, cell* image){

  int newPos = par_pick_neighbor(rows, cols, pos);
  cell neighbor_cell = image[newPos];
  cell my_cell = image[pos];

  if(neighbor_cell.color == WHITE){
    return my_cell;
  }
  cell next_cell = my_cell;

  if(my_cell.color == WHITE){
    next_cell.strength = 1;
    next_cell.color = neighbor_cell.color;
    return next_cell;
  }
  else {
    if(beats(my_cell, neighbor_cell)){
      next_cell.strength++;
    }
    else{
      next_cell.strength--;
    }
  }

  if(next_cell.strength == 0){
    next_cell.color = neighbor_cell.color;
    next_cell.strength = 1;
  }

  if(next_cell.strength > 4){
    next_cell.strength = 4;
  }

  return next_cell;
}

int par_pick_neighbor(int rows, int cols ,int pos){

  int y = pos/rows;
  int x = pos%rows;

  int chosen = rand() % 8;

  if(chosen == 4){ chosen++; } // a cell cant choose itself
  int c_x = chosen % 3;
  int c_y = chosen / 3;

  int new_x = c_x + x - 1;
  int new_y = c_y + y - 1;

  int new_pos = rows*new_y + new_x;

  //printf("%d\n", rows);

  return new_pos;
}

//Checks if position is a border point

bool is_border(int pos, int rows, int cols){
  if(pos/rows == 0){
    return true;
  }
  else if(pos%rows == 0){
    return true;
  }
  else if(pos/rows == rows - 1){
    return true;
  }
  else if(pos%rows == rows  - 1){
    return true;
  }

  return false;
}


cell pick_neighbor(int x, int y, cell** image){
  int chosen = rand() % 8;

  if(chosen == 4){ chosen++; } // a cell cant choose itself
  int c_x = chosen % 3;
  int c_y = chosen / 3;

  return image[x + c_x - 1][y + c_y - 1];
}

//Function that iterates over image, taking one dimensional arrays as input

void par_iterate_image(cell* old_image, cell* next_image, int imgsize, int rows, int cols){

  for(int pos = rows + 1; pos < rows*cols - rows - 2; pos++){

      if (is_border(pos, rows, cols)) {
	       pos += 2;
      }
      next_image[pos] = par_next_cell(rows, cols, pos, old_image);
  }
  }

void iterate_image(cell** old_image, cell** next_image){

  for(int xx = 1; xx < IMG_X - 2; xx++){
    for(int yy = 1; yy < IMG_Y - 2; yy++){
      // printf("%d %d\n", xx, yy);
      next_image[xx][yy] = next_cell(xx, yy, old_image);
    }
  }
}
