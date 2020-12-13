#include <finite_difference.h>

void finite_difference(
    const int &nx,
    const int &ny,
    const int &nz,
    const int &dim, 
    const double &h, 
    const Eigen::VectorXd &grid_values, 
    Eigen::SparseMatrix<double> &D) {
    
    int dim0 = (dim == 0) ? nx : nx + 1;
    int dim1 = (dim == 1) ? ny : ny + 1;
    int dim2 = (dim == 2) ? nz : nz + 1;

    typedef Eigen::Triplet<double> T;
    std::vector<T> triplet;

    for (int i = 0; i < dim0; i++) {
        for (int j = 0; j < dim1; j++) {
            for (int k = 0; k < dim2; k++) {
                int idx = i * dim1 * dim2 + j * dim2 + dim2;
                int row = i * ny * nz + j * nz + k;

                triplet.emplace_back(T(row, idx, -1./h));
                if (dim == 0)
                    triplet.emplace_back(T(row, idx + ny * nz, 1./h));
                else if (dim == 1)
                    triplet.emplace_back(T(row, idx + nz, 1./h));
                else
                    triplet.emplace_back(T(row, idx + 1, 1.));
            }
        }
    }

    D.resize(nx * ny * nz, grid_values.rows());
    D.setFromTriplets(triplet.begin(), triplet.end());
}