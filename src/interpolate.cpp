#include <interpolate.h>

#include <iostream>

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
    int dim0 = (dim == 0) ? nx : nx + 1;
    int dim1 = (dim == 1) ? ny : ny + 1;
    int dim2 = (dim == 2) ? nz : nz + 1;

    int n_particles = q.rows();
    int i, j, k;
    double xd, yd, zd;
    double c_000, c_001, c_010, c_011, c_100, c_101, c_110, c_111;

    sum = Eigen::VectorXd::Zero(dim0*dim1*dim2);
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
        
        xd = d(0); yd = d(1); zd = d(2);
        i = int(ids(0)); j = int(ids(1)); k = int(ids(2));

        c_000 = (1 - xd) * (1 - yd) * (1 - zd);
        c_001 = xd * (1 - yd) * (1 - zd);
        c_010 = (1 - xd) * yd * (1 - zd);
        c_011 = xd * yd * (1 - zd);
        c_100 = (1 - xd) * (1 - yd) * zd;
        c_101 = xd * (1 - yd) * zd;
        c_110 = (1 - xd) * yd * zd;
        c_111 = xd * yd * zd;

        int ind = i * ny * nz + j * nz + k;

        trip.emplace_back(T(l, ind, c_000));
        trip.emplace_back(T(l, ind + ny*nz, c_001));
        trip.emplace_back(T(l, ind + nz, c_010));
        trip.emplace_back(T(l, ind + ny*nz + nz, c_011));
        trip.emplace_back(T(l, ind + 1, c_100));
        trip.emplace_back(T(l, ind + nz * ny + 1, c_101));
        trip.emplace_back(T(l, ind + nz + 1, c_110));
        trip.emplace_back(T(l, ind + nz + nz * ny + 1, c_111));

        sum(ind) += c_000;
        sum(ind + ny * nz) += c_001;
        sum(ind + nz + 1) += c_010;
        sum(ind + ny * nz + nz) += c_011;
        sum(ind + 1) += c_100;
        sum(ind + nz * ny + 1) += c_101;
        sum(ind + nz + 1) += c_110;
        sum(ind + nz + nz * ny + 1) += c_111;
    }
    sum = sum.unaryExpr([](double x) { return (x == 0.) ? 1.0 : x; });
    W.setFromTriplets(trip.begin(), trip.end());
}