#ifndef SOLVE_PRESSURE_H
#define SOLVE_PRESSURE_H

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <grid.h>
#include <particles.h>

void solve_pressure(Grid& grid);

#endif // !SOLVE_PRESSURE_H
