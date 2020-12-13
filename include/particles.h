#ifndef PARTICLE_H
#define PARTICLE_H

#include <Eigen/Core>

class Particle
{
public:
	Eigen::MatrixXd q;
	Eigen::MatrixXd v;

	Particle();
	~Particle();

private:
	
};




#endif // !PARTICLE_H
