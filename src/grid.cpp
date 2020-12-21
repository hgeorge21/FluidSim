#include <grid.h>
#include <chrono>
#include <random>
#include <iostream>

void Grid::init() {
	// calculate number of cells in each dimension
	Eigen::Vector3d ns = (right_upper_corner - left_lower_corner).cwiseQuotient(h);
	nx = ceil(ns(0));
	ny = ceil(ns(1));
	nz = ceil(ns(2));

	// initialization
	n_grids = nx * ny * nz;
	num_fluid_cells = 0;
	pressure = Eigen::VectorXd::Zero(n_grids);
	phi = Eigen::VectorXd::Zero(n_grids);
	fluid_map = Eigen::VectorXi::Zero(n_grids);
	markers = Eigen::VectorXi::Constant(n_grids, AIRCELL);
}

int Grid::get_idx(const int& xi, const int& yi, const int& zi) {
	return xi * ny * nz + yi * nz + zi;
}


void Grid::apply_boundary_condition() {
	int i, j, k;
	int dim0, dim1, dim2;
	const auto& get_idx2 = [&](const int& xi, const int& yi, const int& zi) {
		return xi * dim1 * dim2 + yi * dim2 + zi;
	};

	// boundary cells and filters out velocity
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			int k1 = 0; int k2 = nz - 1;
			markers(get_idx(i, j, k1)) = SOLIDCELL;
			markers(get_idx(i, j, k2)) = SOLIDCELL;

			dim0 = nx; dim1 = ny; dim2 = nz + 1;
			Vz(get_idx2(i, j, 0)) = 0.;
			Vz(get_idx2(i, j, 1)) = 0.;
			Vz(get_idx2(i, j, nz - 1)) = 0.;
			Vz(get_idx2(i, j, nz)) = 0.;
		}
		for (k = 0; k < nz; k++) {
			int j1 = 0; int j2 = ny - 1;
			markers(get_idx(i, j1, k)) = SOLIDCELL;
			markers(get_idx(i, j2, k)) = SOLIDCELL;

			dim0 = nx; dim1 = ny + 1; dim2 = nz;
			Vy(get_idx2(i, 0, k)) = 0.;
			Vy(get_idx2(i, 1, k)) = 0.;
			Vy(get_idx2(i, ny - 1, k)) = 0.;
			Vy(get_idx2(i, ny, k)) = 0.;
		}
	}
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			int i1 = 0; int i2 = nx - 1;
			markers(get_idx(i1, j, k)) = SOLIDCELL;
			markers(get_idx(i2, j, k)) = SOLIDCELL;

			dim0 = nx + 1; dim1 = ny; dim2 = nz;
			Vx(get_idx2(0, j, k)) = 0.;
			Vx(get_idx2(1, j, k)) = 0.;
			Vx(get_idx2(nx - 1, j, k)) = 0.;
			Vx(get_idx2(nx, j, k)) = 0.;
		}
	}
}


void Grid::save_grids() {
	Vx_ = Vx;
	Vy_ = Vy;
	Vz_ = Vz;
}


void Grid::create_free_boundary() {
	phi.setZero();
	fluid_map.setConstant(-1);

	// get smoothing kernel coefficients
	double c[4];
	c[0] = 1.0;
	c[1] = exp(-1.);
	c[2] = exp(-2.);
	c[3] = exp(-3.);

	// performing smoothing on  +-1 marker cells with 3x3x3 Guassian cube
	int index;
	int i, j, k, di, dj, dk;
	int counter = 0;
	double weight, result;
	for (i = 1; i < nx - 1; i++)
		for (j = 1; j < ny - 1; j++)
			for (k = 1; k < nz - 1; k++) {
				result = 0.;
				weight = 0.;

				index = get_idx(i, j, k);
				for (di = -1; di <= 1; di++)
					for (dj = -1; dj <= 1; dj++)
						for (dk = -1; dk <= 1; dk++) {
							int total = abs(di) + abs(dj) + abs(dk);
							if (i + di > 0 && i + di < nx - 1 && j + dj > 0 && j + dj < ny - 1 && k + dk > 0 && k + dk < nz - 1) {
								result += c[total] * markers(get_idx(i + di, j + dj, k + dk));
								weight += c[total];
							}
						}
				phi(index) = result / weight;
				if (result < 0) {
					fluid_map(index) = counter;
					counter++;
				}
			}
	num_fluid_cells = counter;
}