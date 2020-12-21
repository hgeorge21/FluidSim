#include <Eigen/Core>
#include <Eigen/Sparse>

void trilinear_interpolation(const Eigen::RowVector3d& d, Eigen::VectorXd& w);

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