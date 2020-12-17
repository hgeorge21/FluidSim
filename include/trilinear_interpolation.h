#ifndef TRILINERP_H
#define TRILINERP_H

#include <Eigen/Core>

void trilinear_interpolation(const Eigen::RowVector3d& d, Eigen::VectorXd &w);

#endif // !TRILINERP_H
