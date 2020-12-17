#include <Eigen/Core>
#include <Eigen/Sparse>

//Input:
// nx, ny, nz: numbers of cells per row in these directions in the grid
// dim: the grid(x only / y only / z only) to calculate
// h: width of the cell in the given dimension
// corner: the location of the left_lower_corner
//Output:
// sum: 
// W: weight matrix for all cells in that dimension grid
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