#include <stdio.h>
#include <CL/cl.h>
#include <math.h>
#include <iostream>
#include "include.h"
#include "hash.h"

using namespace std;


int main(int argv, char ** argc){

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
				particles[i][j][k] = new float[6];	
				//This contains velocity and position
				//[x, y, z, Vx, Vy, Vz]
				//As mass is same, therefore it is stored in a single float value
				for(int l=0; l<6; l++){
					particles[i][j][k][0] = i/float(particle_factor);
					particles[i][j][k][1] = j/float(particle_factor);
					particles[i][j][k][2] = k/float(particle_factor);
					particles[i][j][k][3] = 0;
					particles[i][j][k][4] = 0;
					particles[i][j][k][5] = 0;
				}

			}
		}
	}
	////

	//// Initialize the smoothing kernels using (5.14) to compute the compact support radius
	// x has been hard-coded as 13
	float compact_support_radius = pow(3*b*c*w*13/(4*M_PI*number_of_particles),1/3.);



















	/****freeing up memory*******/
	for(int i=0; i<b*particle_factor; i++){
		for(int j=0; j<c*particle_factor; j++){
			delete(particles[i][j]);
		}
		delete(particles[i]);
	}
	delete(particles);
	/****freeing up memory*******/
		
}


