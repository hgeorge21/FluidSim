#include <solve_pressure_free_surface.h>

void solve_pressure_free_surface(Grid& grid) {
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
	int n_fluid_cells = grid.num_fluid_cells;
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

	Eigen::VectorXd div = Eigen::VectorXd::Zero(n_fluid_cells);
	Eigen::SparseMatrix<double> A;
	A.resize(n_fluid_cells, n_fluid_cells);
	typedef Eigen::Triplet<double> T;
	std::vector<T> trip;

	int idx;
	int i, j, k;
	int fcell_num = 0;
	int marker;
	double coeff;

	const auto& set_entry = [&](const int &xi, const int &yi, const int &zi, const int &dir) {
		int index2 = get_idx(xi, yi, zi);
		if (grid.markers[index2] != SOLIDCELL) { // make sure it's not solid cell
			trip.push_back(T(fcell_num, fcell_num, -1. * inv_h(dir)));

			// add the ghost pressure if adjacent cell is air cell
			if (grid.fluid_map(index2) > -1)
				trip.push_back(T(fcell_num, grid.fluid_map(index2), 1. * inv_h(dir)));
			else {
				double alpha = grid.phi(index2) / (grid.phi(index2) - grid.phi(idx));
				coeff = (alpha / (1 - alpha)) * -1. * inv_h(dir);
				trip.push_back(T(fcell_num, fcell_num, coeff));
			}
		}
	};

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				idx = get_idx(i, j, k);
				if (grid.fluid_map(idx) > -1) {
					div[fcell_num] =
						(grid.Vx[get_idx2(i + 1, j, k, 0)] - grid.Vx[get_idx2(i, j, k, 0)]) / h(0) +
						(grid.Vy[get_idx2(i, j + 1, k, 1)] - grid.Vy[get_idx2(i, j, k, 1)]) / h(1) +
						(grid.Vz[get_idx2(i, j, k + 1, 2)] - grid.Vz[get_idx2(i, j, k, 2)]) / h(2);

					// all 6 neighbors
					set_entry(i - 1, j, k, 0);
					set_entry(i + 1, j, k, 0);
					set_entry(i, j - 1, k, 1);
					set_entry(i, j + 1, k, 1);
					set_entry(i, j, k - 1, 2);
					set_entry(i, j, k + 1, 2);
					fcell_num += 1;
				}
			}
		}
	}


	A.setFromTriplets(trip.begin(), trip.end());
	div = grid.density / grid.dt * div;

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
	cg.setTolerance(1e-5);
	cg.compute(A.transpose() * A);
	Eigen::VectorXd temp_p = cg.solve(A.transpose() * div);

	// set the pressure for fluid cells
	grid.pressure.setZero();
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				int map_index = grid.fluid_map(get_idx(i, j, k));
				if(map_index > -1)
					grid.pressure(get_idx(i, j, k)) = temp_p(map_index);
			}
		}
	}
}