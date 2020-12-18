#include <divergence.h>
#include <grid.h>
#include <vector>

void divergence_op(
	const int& nx,
	const int& ny,
	const int& nz,
	const int& dim,
	const Eigen::Vector3d& h,
	const Eigen::VectorXi& markers,
	Eigen::SparseMatrix<double>& D) {

	int dim0 = (dim == 0) ? nx + 1 : nx;
	int dim1 = (dim == 1) ? ny + 1 : ny;
	int dim2 = (dim == 2) ? nz + 1 : nz;

	D.resize(nx * ny * nz, dim0 * dim1 * dim2);

	typedef Eigen::Triplet<double> T;
	std::vector<T> trip;

	const auto& get_idx = [&](int xi, int yi, int zi) { return xi * ny * nz + yi * nz + zi; };
	const auto& get_idx2 = [&](int xi, int yi, int zi) { return xi * dim1 * dim2 + yi * dim2 + zi; };

	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			for (int k = 1; k < nz - 1; k++) {
				int row = get_idx(i, j, k);
				if (markers(row) == FLUIDCELL) {
					if (dim < 3) {
						trip.emplace_back(T(row, get_idx2(i + 1, j, k), 1. / h(dim)));
						trip.emplace_back(T(row, get_idx2(i, j, k), -1. / h(dim)));
					}
					/*else {
						trip.emplace_back(T(row, get_idx(i - 1, j, k), -1. / h(0)));
						trip.emplace_back(T(row, get_idx(i + 1, j, k),  1. / h(0)));
						trip.emplace_back(T(row, get_idx(i, j - 1, k), -1. / h(1)));
						trip.emplace_back(T(row, get_idx(i, j + 1, k),  1. / h(1)));
						trip.emplace_back(T(row, get_idx(i, j, k - 1), -1. / h(2)));
						trip.emplace_back(T(row, get_idx(i, j, k + 1), -1. / h(2)));
						trip.emplace_back(T(row, get_idx(i - 1, j, k), -1. / h(0)));
					}*/
				}
			}
		}
	}

	D.setFromTriplets(trip.begin(), trip.end());
}