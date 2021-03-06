#ifndef TRANSFER_TO_GRID_H
#define TRANSFER_TO_GRID_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <particles.h>
#include <grid.h>

// Transfer velocity from particle to grid
void transfer_to_grid(Grid& grid, Particle& particles);

#endif //TRANSFER_TO_GRID_H