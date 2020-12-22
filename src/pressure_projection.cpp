#include <pressure_projection.h>

void pressure_projection(Grid& grid) {
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
	Eigen::Vector3d h = grid.h;

	const auto& get_idx = [&](const int& xi, const int& yi, const int& zi) {
		return xi * ny * nz + yi * nz + zi;
	};
	const auto& get_idx2 = [&](const int& xi, const int& yi, const int& zi, const int& dim) {
		int dim0 = (dim == 0) ? nx + 1 : nx;
		int dim1 = (dim == 1) ? ny + 1 : ny;
		int dim2 = (dim == 2) ? nz + 1 : nz;
		return xi * dim1 * dim2 + yi * dim2 + zi;
	};
	
	int i, j, k;
	int index, index_prev;
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				index = get_idx(i, j, k);
				if (i > 1 && i < nx - 1) {
					index_prev = get_idx(i - 1, j, k);
					grid.Vx(get_idx2(i, j, k, 0)) -= (grid.dt / grid.density) * (grid.pressure(index) - grid.pressure(index_prev)) / h(0);
				}
				if (j > 1 && j < ny - 1) {
					index_prev = get_idx(i, j - 1, k);
					grid.Vy(get_idx2(i, j, k, 1)) -= (grid.dt / grid.density) * (grid.pressure(index) - grid.pressure(index_prev)) / h(1);
				}
				if (k > 1 && k < nz - 1) {
					index_prev = get_idx(i, j, k - 1);
					grid.Vz(get_idx2(i, j, k, 2)) -= (grid.dt / grid.density) * (grid.pressure(index) - grid.pressure(index_prev)) / h(2);
				}
			}
		}
	}
}