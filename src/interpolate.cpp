#include <interpolate.h>
#include <omp.h>
#include <iostream>

#include <trilinear_interpolation.h>

void interpolate(
    const int &nx,
    const int &ny,
    const int &nz,
    const int &dim,
    const double &h,
    const Eigen::RowVector3d& corner,
    const Eigen::MatrixXd& q,
    Eigen::VectorXd& sum,
    Eigen::SparseMatrix<double>& W) {

    int tmp = (dim == 0) ? nx : (dim == 1) ? ny : nz;
    int dim0 = (dim == 0) ? nx + 1 : nx;
    int dim1 = (dim == 1) ? ny + 1 : ny;
    int dim2 = (dim == 2) ? nz + 1 : nz;

    int n_particles = q.rows();
    int i, j, k;
    double xd, yd, zd;
    Eigen::VectorXd w;

    W.resize(n_particles, dim0*dim1*dim2);
    W.reserve(8 * n_particles);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> trip;

    for (int l = 0; l < n_particles; l++) {
        Eigen::RowVector3d x = q.row(l) - corner;
        Eigen::RowVector3d d = x / h;
        Eigen::RowVector3d ids = d.unaryExpr([](double x){ return floor(std::max(x, 0.)); });
        ids(dim) = std::min(int(ids(dim)), tmp - 2);

        d = d - ids;
        d(dim) = std::min(d(dim), 1.0);
        
        i = int(ids(0)); j = int(ids(1)); k = int(ids(2));
        int ind = i * ny * nz + j * nz + k;

        trilinear_interpolation(d, w);

        trip.emplace_back(T(l, ind,                    w[0]));
        trip.emplace_back(T(l, ind + ny*nz,            w[1]));
        trip.emplace_back(T(l, ind + nz,               w[2]));
        trip.emplace_back(T(l, ind + ny*nz + nz,       w[3]));
        trip.emplace_back(T(l, ind + 1,                w[4]));
        trip.emplace_back(T(l, ind + nz * ny + 1,      w[5]));
        trip.emplace_back(T(l, ind + nz + 1,           w[6]));
        trip.emplace_back(T(l, ind + nz + nz * ny + 1, w[7]));
    }

    W.setFromTriplets(trip.begin(), trip.end());
  
    sum = Eigen::VectorXd::Ones(n_particles);
    sum = W.transpose() * sum;
    sum = sum.unaryExpr([](double x) { return (x == 0.) ? 1.0 : x; });
}