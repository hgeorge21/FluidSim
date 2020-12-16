#include <advect_velocity.h>
#include <interpolate.h>

void advect_velocity(Grid& grid, Particle& particles, double dt) {
	if (particles.method == advection_method::SIMPLE) {
		particles.q = particles.q + dt * particles.v;
	}
	else {
		// 2nd Order Runge Kutta method
	}
}