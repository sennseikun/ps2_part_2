#ifndef CA_H
#define CA_H

#include "RPS.h"


void free_img(cell** image);
bool neighborhood_contains(int x, int y, cell** image, int color);
int next_color(int x, int y, cell** image);
bool is_border(int pos, int rows, int cols);
void iterate_image(cell** old_image, cell** next_image);
void par_iterate_image(cell* old_image, cell* next_image, int imgsize, int rows, int cols);
cell** alloc_img();

#endif
