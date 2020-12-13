#include <Eigen/Core>
#include <Eigen/Sparse>

void interpolate(
    const int& nx,
    const int& ny,
    const int& nz,
    const int& dim,
    const double& h,
    const Eigen::RowVector3d& corner,
    const Eigen::MatrixXd& q,
    Eigen::VectorXd& sum,
    Eigen::SparseMatrix<double>& W);