#ifndef RECT_H
#define RECT_H

#include <Eigen/Core>

void create_rectangle(Eigen::MatrixXd& V, Eigen::MatrixXi &F) {
	V = Eigen::MatrixXd::Zero(14, 3);
	V.row(0) << -0.5, -0.3, -0.5;
	V.row(1) << 0.5, -0.3, -0.5;
	V.row(2) << 0, 0, -0.5;
	V.row(3) << -0.5, 0.3, -0.5;
	V.row(4) << 0.5, 0.3, -0.5;

	V.row(5) << -0.5, -0.3, 0.5;
	V.row(6) << 0.5, -0.3, 0.5;
	V.row(7) << 0, 0, 0.5;
	V.row(8) << -0.5, 0.3, 0.5;
	V.row(9) << 0.5, 0.3, 0.5;

	V.row(10) << 0, 0.3, 0;
	V.row(11) << 0.5, 0, 0;
	V.row(12) << 0, -0.3, 0;
	V.row(13) << -0.5, 0, 0;

	F = Eigen::MatrixXi(24, 3);
	F.row(0) << 2, 1, 0;
	F.row(1) << 3, 2, 0;
	F.row(2) << 4, 1, 2;
	F.row(3) << 3, 4, 2;

	F.row(4) << 3, 0, 13;
	F.row(5) << 3, 13, 8;
	F.row(6) << 0, 5, 13;
	F.row(7) << 13, 5, 8;

	F.row(8) << 5, 7, 8;
	F.row(9) << 7, 9, 8;
	F.row(10) << 5, 6, 7;
	F.row(11) << 6, 9, 7;

	F.row(12) << 4, 11, 1;
	F.row(13) << 4, 9, 11;
	F.row(14) << 6, 11, 9;
	F.row(15) << 6, 1, 11;

	F.row(16) << 3, 4, 10;
	F.row(17) << 3, 8, 10;
	F.row(18) << 4, 10, 9;
	F.row(19) << 8, 9, 10;

	F.row(20) << 1, 12, 0;
	F.row(21) << 1, 6, 12;
	F.row(22) << 5, 12, 6;
	F.row(23) << 5, 0, 12;
}

#endif // !RECT_H
