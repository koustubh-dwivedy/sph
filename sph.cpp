//// 	WE ARE SIMULATING A DAM BREAK
	/*
	http://www.ijser.org/paper/Simulation-of-Dam-Break-Using-Modified-Incompressible-Smoothed/Image_034.jpg
	
		****w corresponds to depth****
	\						\
	\						\
	\ <-b->					\
	\\\\\\\\\				\
	\		\ 				\  |
	\		\ |				\  d
	\		\ c				\  |
	\		\ |				\
	\		\ 				\
	\		\				\
	\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
			<--a-->

			!!!!!This does not follow vector right hand rule!!!!!
			y     z
			|   /
			|  /
			| /
			--->x

	*/
////

#include <stdio.h>
#include <CL/cl.h>
#include <math.h>
#include <iostream>
#include "includes/include.h"

using namespace std;


int main(int argv, char ** argc){

	/********
	Initializing SPH System
	********/
	//// Environment Parameters
	const float gravity = 9.81;
	const float time_step = 0.01;
	const float temperature =  283.15;
	const int atm_pressure = 101325;
	const float density = 1000;
	const float a  = 10;
	const float b = 4;
	const float d = 10;
	const float c = 6;
	const float w = 10;
	const float particle_factor = 10;
	const int number_of_particles = b*c*w*particle_factor*particle_factor*particle_factor; // 40*60*100 ~ b*c*w
	const float mass = density*b*c*w/float(number_of_particles);
	////

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

		
}


