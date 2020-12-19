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
	D.setZero();

	typedef Eigen::Triplet<double> T;
	std::vector<T> trip1;

	// for cell
	const auto& get_idx = [&](int xi, int yi, int zi) { return xi * ny * nz + yi * nz + zi; };
	// for different grids
	const auto& get_idx2 = [&](int xi, int yi, int zi) { return xi * dim1 * dim2 + yi * dim2 + zi; };

	for (int i = 1; i < nx - 2; i++) {
		for (int j = 1; j < ny - 2; j++) {
			for (int k = 1; k < nz - 2; k++) {
				int index = get_idx(i, j, k);
				if (markers(index) == FLUIDCELL) {
					if (dim < 3) {
						trip1.emplace_back(T(index, get_idx2(i + 1, j, k), 1. / h(dim)));
						trip1.emplace_back(T(index, get_idx2(i, j, k), -1. / h(dim)));
					}
				}
			}
		}
	}

	D.setFromTriplets(trip1.begin(), trip1.end());
}

void gradient_op(
	const int& nx,
	const int& ny,
	const int& nz,
	const int& dim,
	const Eigen::Vector3d& h,
	const Eigen::VectorXi& markers,
	Eigen::SparseMatrix<double>& G) {

	int ll0 = (dim == 0) ? 2 : 1;
	int ll1 = (dim == 1) ? 2 : 1;
	int ll2 = (dim == 2) ? 2 : 1;
	int ul0 = (dim == 0) ? nx - 1 : nx;
	int ul1 = (dim == 1) ? ny - 1 : ny;
	int ul2 = (dim == 2) ? nz - 1 : nz;
	int dim0 = (dim == 0) ? nx + 1 : nx;
	int dim1 = (dim == 1) ? ny + 1 : ny;
	int dim2 = (dim == 2) ? nz + 1 : nz;

	G.resize(dim0 * dim1 * dim2, nx * ny * nz);
	typedef Eigen::Triplet<double> T;
	std::vector<T> trip1;

	// for cell
	const auto& get_idx = [&](int xi, int yi, int zi) { return xi * ny * nz + yi * nz + zi; };
	// for different grids
	const auto& get_idx2 = [&](int xi, int yi, int zi) { return xi * dim1 * dim2 + yi * dim2 + zi; };

	for (int i = ll0; i < ul0; i++) {
		for (int j = ll1; j < ul1; j++) {
			for (int k = ll2; k < ul2; k++) {
				int index = get_idx(i, j, k); // index of grid edge
				int index_prev;
				if (dim == 0)
					index_prev = get_idx(i - 1, j, k);
				else if (dim == 1)
					index_prev = get_idx(i, j - 1, k);
				else
					index_prev = get_idx(i, j, k - 1);

				int index_grid = get_idx2(i, j, k);
				// at least one of the two cells is fluid; none is solid
				if ((markers(index) == FLUIDCELL || markers(index_prev) == FLUIDCELL) && 
					(markers(index) != SOLIDCELL && markers(index_prev) != SOLIDCELL)) {
					trip1.emplace_back(T(index_grid, index, 1. / h(dim)));
					trip1.emplace_back(T(index_grid, index_prev, -1. / h(dim)));

				}
			}
		}
	}

	G.setFromTriplets(trip1.begin(), trip1.end());
}
