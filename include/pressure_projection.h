#ifndef PRESSURE_PROJ_H
#define PRESSURE_PROJ_H

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <grid.h>
#include <particles.h>

// Performs pressure projection update
void pressure_projection(Grid& grid);

#endif // !PRESSURE_PROJ_H