#include <Eigen/Core>
#include <Eigen/Sparse>

#include <grid.h>
#include <particles.h>

void finite_difference(
    const int& nx,
    const int& ny,
    const int& nz,
    const int& dim,
    const double& h,
    const Eigen::VectorXd& grid_values,
    Eigen::SparseMatrix<double>& D);