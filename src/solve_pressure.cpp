#include <solve_pressure.h>

void solve_pressure(Grid& grid) {
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
	Eigen::Vector3d h = grid.h;
	Eigen::Vector3d inv_h = h.unaryExpr([](double x) { return 1. / pow(x, 2); });

	const auto& get_idx = [&](const int& xi, const int& yi, const int& zi) {
		return xi * ny * nz + yi * nz + zi;
	};
	const auto& get_idx2 = [&](const int& xi, const int& yi, const int& zi, const int& dim) {
		int dim0 = (dim == 0) ? nx + 1 : nx;
		int dim1 = (dim == 1) ? ny + 1 : ny;
		int dim2 = (dim == 2) ? nz + 1 : nz;
		return xi * dim1 * dim2 + yi * dim2 + zi;
	};

	Eigen::VectorXd div = Eigen::VectorXd::Zero(grid.n_grids);
	Eigen::SparseMatrix<double> A;
	A.resize(grid.n_grids, grid.n_grids);

	typedef Eigen::Triplet<double> T;
	std::vector<T> trip;

	int idx;
	int i, j, k;

	const auto& set_entry = [&](const int &xi, const int &yi, const int &zi, const int &dir) {
		int index2 = get_idx(xi, yi, zi);
		if (grid.markers[index2] != SOLIDCELL)
			trip.push_back(T(idx, index2, inv_h(dir)));
		else
			trip.push_back(T(idx, idx, inv_h(dir)));
	};

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				// build RHS
				idx = get_idx(i, j, k);
				if (grid.markers[idx] == FLUIDCELL) {
					div[idx] =
						(grid.Vx[get_idx2(i + 1, j, k, 0)] - grid.Vx[get_idx2(i, j, k, 0)]) / h(0) +
						(grid.Vy[get_idx2(i, j + 1, k, 1)] - grid.Vy[get_idx2(i, j, k, 1)]) / h(1) +
						(grid.Vz[get_idx2(i, j, k + 1, 2)] - grid.Vz[get_idx2(i, j, k, 2)]) / h(2);
				}

				if (grid.markers[idx] == FLUIDCELL) {
					trip.push_back(T(idx, idx, -2. * inv_h.sum()));
					set_entry(i - 1, j, k, 0);
					set_entry(i + 1, j, k, 0);
					set_entry(i, j - 1, k, 1);
					set_entry(i, j + 1, k, 1);
					set_entry(i, j, k - 1, 2);
					set_entry(i, j, k + 1, 2);
				}
			}
		}
	}

	A.setFromTriplets(trip.begin(), trip.end());
	div = grid.density / grid.dt * div;

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
	cg.setTolerance(1e-5);
	cg.compute(A.transpose() * A);
	grid.pressure = cg.solve(A.transpose() * div);
}