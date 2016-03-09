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

using namespace std;


float default_kernel(float position[3], float particle[6], float radius);
void default_kernel_gradient(float position[3], float particle[6], float radius, float* receiver);
float default_kernel_laplacian(float position[3], float particle[6], float radius);
float pressure_kernel(float position[3], float particle[6], float radius);
void pressure_kernel_gradient(float position[3], float particle[6], float radius, float* receiver);
float pressure_kernel_laplacian(float position[3], float particle[6], float radius);
float viscosity_kernel(float position[3], float particle[6], float radius);
void viscosity_kernel_gradient(float position[3], float particle[6], float radius, float* receiver);
float viscosity_kernel_laplacian(float position[3], float particle[6], float radius);




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


// radius -> compact support radius
// position -> position where to calculate kernel value
// particle -> particle for which kernel is being calculated
float default_kernel(float position[3], float particle[6], float radius){
	float separation = pow(pow(position[0] - particle[0], 2) + pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2), 1/2.);
	if (separation <= radius){
		float value = 315*(pow(pow(radius, 2) - pow(separation, 2), 3))/(64*M_PI*pow(radius, 9));
		return value;
	}
	else{
		return 0;
	}
}
//receiver should be first initialized by the following-> float receiver[3]
void default_kernel_gradient(float position[3], float particle[6], float radius, float* receiver){
	float separation = pow(pow(position[0] - particle[0], 2) + pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2), 1/2.);
	if (separation <= radius){
		float value = -945*(pow(pow(radius, 2) - pow(separation, 2), 3))/(32*M_PI*pow(radius, 9));
		receiver[0] = value*(position[0] - particle[0]);
		receiver[1] = value*(position[1] - particle[1]);
		receiver[2] = value*(position[2] - particle[2]);
			
	}
}

float default_kernel_laplacian(float position[3], float particle[6], float radius){
	float separation = pow(pow(position[0] - particle[0], 2) + pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2), 1/2.);
	if (separation <= radius){
		float value = -945*(pow(radius, 2) - pow(separation, 2))*(3*pow(radius, 2) - 7*pow(separation, 2))/(32*M_PI*pow(radius, 9));
		return value;
	}
	else{
		return 0;
	}	
}



float pressure_kernel(float position[3], float particle[6], float radius){
	float separation = pow(pow(position[0] - particle[0], 2) + pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2), 1/2.);
	if (separation <= radius){
		float value = 15*(pow(radius - separation, 3))/(M_PI*pow(radius, 6));
		return value;
	}
	else{
		return 0;
	}
}


/////////////THIS FUNCTION MAY CAUSE PROBLEMS DUE TO LIMITING VALUES////////////
void pressure_kernel_gradient(float position[3], float particle[6], float radius, float* receiver){
	float separation = pow(pow(position[0] - particle[0], 2) + pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2), 1/2.);
	if (separation <= radius){
		float value = -45*(pow(radius - separation, 2))/(M_PI*pow(radius, 6));
		receiver[0] = value*(position[0] - particle[0])/separation;
		receiver[1] = value*(position[1] - particle[1])/separation;
		receiver[2] = value*(position[2] - particle[2])/separation;
			
	}
}
////////////////////////////////////////////////////////////////////////////////
/////////////THIS FUNCTION MAY CAUSE PROBLEMS DUE TO LIMITING VALUES////////////
float pressure_kernel_laplacian(float position[3], float particle[6], float radius){
	float separation = pow(pow(position[0] - particle[0], 2) + pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2), 1/2.);
	if (separation <= radius){
		float value = -90*(radius - separation)*(radius - 2*separation)/(M_PI*separation*pow(radius, 6));
		return value;
	}
	else{
		return 0;
	}	
}
////////////////////////////////////////////////////////////////////////////////



/////////////THIS FUNCTION MAY CAUSE PROBLEMS DUE TO LIMITING VALUES////////////
float viscosity_kernel(float position[3], float particle[6], float radius){
	float separation = pow(pow(position[0] - particle[0], 2) + pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2), 1/2.);
	if (separation <= radius){
		float value1 = 15/(2*M_PI*pow(radius, 3));
		float value2 = (-pow(separation, 3)/(2*pow(radius, 3))) + (pow(separation, 2)/pow(radius, 2)) + (radius/(2*separation)) -1;
		return value1*value2;
	}
	else{
		return 0;
	}
}
////////////////////////////////////////////////////////////////////////////////
/////////////THIS FUNCTION MAY CAUSE PROBLEMS DUE TO LIMITING VALUES////////////
void viscosity_kernel_gradient(float position[3], float particle[6], float radius, float* receiver){
	float separation = pow(pow(position[0] - particle[0], 2) + pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2), 1/2.);
	if (separation <= radius){
		float value1 = -15/(2*M_PI*pow(radius, 3));
		float value2 = (-3*separation/(2*pow(radius, 3))) + 2/pow(radius, 2) - (radius/(2*pow(separation, 3)));
		receiver[0] = value1*value2*(position[0] - particle[0]);
		receiver[1] = value1*value2*(position[1] - particle[1]);
		receiver[2] = value1*value2*(position[2] - particle[2]);
			
	}
}
////////////////////////////////////////////////////////////////////////////////

float viscosity_kernel_laplacian(float position[3], float particle[6], float radius){
	float separation = pow(pow(position[0] - particle[0], 2) + pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2), 1/2.);
	if (separation <= radius){
		float value = 45*(radius - separation)/(M_PI*separation*pow(radius, 6));
		return value;
	}
	else{
		return 0;
	}	
}

