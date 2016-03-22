//// 	WE ARE SIMULATING A DAM BREAK
	/*
	http://www.ijser.org/paper/Simulation-of-Dam-Break-Using-Modified-Incompressible-Smoothed/Image_034.jpg
	
		****w corresponds to depth****
	\										\
	\										\
	\ <-b->									\
	\\\\\\\\\								\
	\		\ 								\  |
	\		\ |								\  d
	\		\ c								\  |
	\		\ |								\
	\		\ 								\
	\		\								\
	\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
			<--a-->

		  z|   /y
		   |  /
		   | /
		   |/_____x

	*/
////

#include "includes/sph.h"

int main(int argc, char **argv){
	int part[3] = {20, 50, 30};//THIS IS NOT NUMBER OF PARTICLES, BUT THE DIMENSION OF particles ARRAY. (WHICH TURNS OUT TO BE EQUAL TO NUMBER OF PARTICLES)
	float**** particles;
	//Define your particles here
	/********
	Initializing SPH System
	********/
	//// Creating n particles and setting their initial positions, velocities and mass
	particles = new float***[(part[0])];
	for(int i=0; i<(part[0]); i++){
		particles[i] = new float**[(part[1])];
		for(int j=0; j<(part[1]); j++){
			particles[i][j] = new float*[(part[2])];
			for(int k=0; k<(part[2]); k++){
				particles[i][j][k] = new float[7];	
				//This contains velocity and position and whether particle is present or not
				//[isPresent, x, y, z, Vx, Vy, Vz]
				for(int l=0; l<7; l++){
					particles[i][j][k][0] = 1;//specifying that this is a particle
					particles[i][j][k][1] = 4*i/part[0];
					particles[i][j][k][2] = 10*j/part[1];
					particles[i][j][k][3] = 6*k/part[2];
					particles[i][j][k][4] = 0;
					particles[i][j][k][5] = 0;
					particles[i][j][k][6] = 0;
				}

			}
		}
	}
	////

	environment environ1;
	environ1.gravity(9.81);
	environ1.time_step(0.01);
	environ1.temperature(283.15);
	environ1.atmosphericPressure(101325);
	environ1.density(1000);
	environ1.numParticles(20*50*30);//AT PRESENT, IT IS LIMITED BY MAXIMUM VALUE OF int
	environ1.fluidVolume(4*10*6);//THIS IS THE ACTUAL OCCUPIED VOLUME OF FLUID (x*y*z)
	environ1.environmentLength(10);
	environ1.environmentHeight(10);
	environ1.environmentWidth(10);
	environ1.particlesDimen(part);//THIS  IS THE ARRAY CONTAINING THE DIMENSIONS OF particles ARRAY
	environ1.addParticles(particles);//THIS IS THE ARRAY CONTAINING [0 OR 1](PARTICLE OR NOT) [X, Y, Z](INITIAL POSITIONS) AND [Vx, Vy, Vz](INITIAL VELOCITIES)
	environ1.envInit();
	environ1.simulate();
	environ1.environmentFree();
	delete(part);//this is needed to be done manually	
}