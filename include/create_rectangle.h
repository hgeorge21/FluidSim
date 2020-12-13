#ifndef RECT_H
#define RECT_H

#include <Eigen/Core>

void create_rectangle(
	const Eigen::Vector3d &left_lower_corner,
	const Eigen::Vector3d &right_upper_corner,
	Eigen::MatrixXd &V, 
	Eigen::MatrixXi &F) {
	
	double llx = left_lower_corner(0);
	double lly = left_lower_corner(1);
	double llz = left_lower_corner(2);

	double ulx = right_upper_corner(0);
	double uly = right_upper_corner(1);
	double ulz = right_upper_corner(2);

	double mdx = 0.5 * (llx + ulx);
	double mdy = 0.5 * (lly + uly);
	double mdz = 0.5 * (llz + ulz);

	V = Eigen::MatrixXd::Zero(14, 3);
	V.row(0) << llx, lly, llz;
	V.row(1) << ulx, lly, llz;
	V.row(2) << mdx, mdy, llz;
	V.row(3) << llx, uly, llz;
	V.row(4) << ulx, uly, llz;

	V.row(5) << llx, lly, ulz;
	V.row(6) << ulx, lly, ulz;
	V.row(7) << mdx, mdy, ulz;
	V.row(8) << llx, uly, ulz;
	V.row(9) << ulx, uly, ulz;

	V.row(10) << mdx, uly, mdz;
	V.row(11) << ulx, mdy, mdz;
	V.row(12) << mdx, lly, mdz;
	V.row(13) << llx, mdy, mdz;

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
