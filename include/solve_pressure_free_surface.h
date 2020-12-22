#ifndef SOLVE_PRESSURE_FREE_H
#define SOLVE_PRESSURE_FREE_H

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <grid.h>
#include <particles.h>

// Solve for pressure Ap=d with free surface, builds A and d
void solve_pressure_free_surface(Grid& grid);

#endif // !SOLVE_PRESSURE_FREE_H
