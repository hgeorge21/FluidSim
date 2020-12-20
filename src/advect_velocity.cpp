#include <advect_velocity.h>
#include <interpolate.h>

void advect_velocity(Grid& grid, Particle& particles, double dt) {
	particles.q = particles.q + dt * particles.v;
	
	// clamping to boundary	
	Eigen::Vector3d v1 = grid.left_lower_corner + grid.h;
	Eigen::Vector3d v2 = grid.right_upper_corner - grid.h;
	particles.q.col(0) = particles.q.col(0).unaryExpr([&](double x) { return std::max(v1(0), std::min(v2(0),x)); });
	particles.q.col(1) = particles.q.col(1).unaryExpr([&](double x) { return std::max(v1(1), std::min(v2(1),x)); });
	particles.q.col(2) = particles.q.col(2).unaryExpr([&](double x) { return std::max(v1(2), std::min(v2(2),x)); });
}