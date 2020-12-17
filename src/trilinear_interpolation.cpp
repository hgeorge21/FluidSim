#include <trilinear_interpolation.h>

void trilinear_interpolation(const Eigen::RowVector3d& d, Eigen::VectorXd& w) {
    w = Eigen::VectorXd(8);

    double xd, yd, zd;
    xd = d(0); yd = d(1); zd = d(2);

    w(0) = (1 - xd) * (1 - yd) * (1 - zd);
    w(1) = xd * (1 - yd) * (1 - zd);
    w(2) = (1 - xd) * yd * (1 - zd);
    w(3) = xd * yd * (1 - zd);
    w(4) = (1 - xd) * (1 - yd) * zd;
    w(5) = xd * (1 - yd) * zd;
    w(6) = (1 - xd) * yd * zd;
    w(7) = xd * yd * zd;
}