#ifndef PARTICLE_H
#define PARTICLE_H

#include <Eigen/Core>
#include <grid.h>

class Particle
{
public:
	Eigen::MatrixXd q;
	Eigen::MatrixXd v;

	Particle();
	~Particle();

	void generate_particles(Grid &grid);

private:
	
};




#endif // !PARTICLE_H
