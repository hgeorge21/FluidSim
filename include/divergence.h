#include <Eigen/Core>
#include <Eigen/Sparse>

void divergence_op(
	const int& nx,
	const int& ny,
	const int& nz,
	const int& dim,
	const Eigen::Vector3d& h,
	const Eigen::VectorXi& markers,
	Eigen::SparseMatrix<double>& D);