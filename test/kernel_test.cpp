#include "../includes/kernels.cpp"
#include <iostream>
using namespace std;

int main(void){
	int part[3] = {4, 10, 6};//THIS IS NOT NUMBER OF PARTICLES, BUT THE DIMENSION OF particles ARRAY. (WHICH TURNS OUT TO BE EQUAL TO NUMBER OF PARTICLES)
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
					particles[i][j][k][1] = 4.0*i/part[0];
					particles[i][j][k][2] = 10.0*j/part[1];
					particles[i][j][k][3] = 6.0*k/part[2];
					particles[i][j][k][4] = 0;
					particles[i][j][k][5] = 0;
					particles[i][j][k][6] = 0;
				}

			}
		}
	}
	////
	float radius = 10.0;
	float receiver[3];
	/*
		If the second parameter in default_kernel (rather all the kernels) is particles[1][1][1] then the values outputted are -inf, inf, nan and -nan
	*/
	cout << default_kernel(particles[1][1][1], particles[1][1][1], radius) << endl;
	default_kernel_gradient(particles[1][1][1], particles[1][1][1], radius, receiver);
	for(int i=0; i<3; i++){
		cout << receiver[i] << "		";
	}
	cout << endl;
	cout << default_kernel_laplacian(particles[1][1][1], particles[1][1][1], radius) << endl;
	cout << pressure_kernel(particles[1][1][1], particles[1][1][1], radius) << endl;
	pressure_kernel_gradient(particles[1][1][1], particles[1][1][1], radius, receiver);
	for(int i=0; i<3; i++){
		cout << receiver[i] << "		";
	}
	cout << endl;
	cout << pressure_kernel_laplacian(particles[1][1][1], particles[1][1][1], radius) << endl;
	cout << viscosity_kernel(particles[1][1][1], particles[1][1][1], radius) << endl;
	viscosity_kernel_gradient(particles[1][1][1], particles[1][1][1], radius, receiver);
	for(int i=0; i<3; i++){
		cout << receiver[i] << "		";
	}
	cout << endl;
	cout << viscosity_kernel_laplacian(particles[1][1][1], particles[1][1][1], radius) << endl;
}