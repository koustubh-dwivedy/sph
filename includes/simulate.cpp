#include "sph.h"

void environment::simulate(){
	if(ready == 1){
		/********
		Initializing SPH System
		********/
		//// Creating n particles and setting their initial positions, velocities and mass
		float**** particles = new float***[int(b*particle_factor)];
		for(int i=0; i<b*particle_factor; i++){
			particles[i] = new float**[int(c*particle_factor)];
			for(int j=0; j<c*particle_factor; j++){
				particles[i][j] = new float*[int(w*particle_factor)];
				for(int k=0; k<w*particle_factor; k++){
					particles[i][j][k] = new float[7];	
					//This contains velocity and position
					//[x, y, z, Vx, Vy, Vz]
					//As mass is same, therefore it is stored in a single float value
					for(int l=0; l<6; l++){
						particles[i][j][k][0] = 1;//specifying that this is a particle
						particles[i][j][k][1] = i/float(particle_factor);
						particles[i][j][k][2] = j/float(particle_factor);
						particles[i][j][k][3] = k/float(particle_factor);
						particles[i][j][k][4] = 0;
						particles[i][j][k][5] = 0;
						particles[i][j][k][6] = 0;
					}

				}
			}
		}
		////

		//// Initialize the smoothing kernels using (5.14) to compute the compact support radius
		// x has been hard-coded as 13
		float compact_support_radius = pow(3*b*c*w*13/(4*M_PI*number_of_particles),1/3.);



















	}
}