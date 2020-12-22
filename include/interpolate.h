#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <Eigen/Core>
#include <Eigen/Sparse>

// trilinear interpolation based on distance from corner of a cell
void trilinear_interpolation(
    const Eigen::RowVector3d& d, 
    Eigen::VectorXd& w);

// computes distance from corner of the cell
void bary(
    const double& x, 
    const double& h, 
    int& i, 
    double& d);

// computes distance from corner of the cell with clamping
void bary_center(
    const int& n, 
    const double& x, 
    const double& h, 
    int& i, 
    double& d);

// Compute the trilinear interpolation matrix between
// staggered grid and particles based on position
void interpolate(
    const int& nx,
    const int& ny,
    const int& nz,
    const int& dim,
    const Eigen::RowVector3d& h,
    const Eigen::RowVector3d& corner,
    const Eigen::MatrixXd& q,
    Eigen::VectorXd& sum,
    Eigen::SparseMatrix<double>& W);

#endif //INTERPOLATE_H