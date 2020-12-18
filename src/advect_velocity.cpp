#include <advect_velocity.h>
#include <interpolate.h>

void advect_velocity(Grid& grid, Particle& particles, double dt) {
	particles.q = particles.q + dt * particles.v;
}