#include "../includes/hash.h"
#include "../includes/prime.cpp"
#include <iostream>
using namespace std;

int main(void){

	int part[3] = {4, 10, 6};
	float**** particles;
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
	
	hashTable table;

	//checked that program is releasing/freeing up memory properly
	float radius = 1;//radius above 0.3 is OK. BELOW THAT IS GIVING SEGFAULT (irrespective of search argument of 0/1 in particleQuery)
	//the above bug is not giving pains in intel's compiler
	int hash_size = NextPrime(part[0]*part[1]*part[2]);
	cout << hash_size << endl;
	table.insertParticles(particles, hash_size, part, radius);

	//cout << "OK" << endl;
	
	int actualNearestNeighbours;
	float* nearestNeighbourAdresses[40];

	for(int i=0; i<5; i++){
		table.particleQuery(particles[3][i][5][1], particles[3][i][5][2], particles[3][i][5][3], nearestNeighbourAdresses, 40, &actualNearestNeighbours, radius, 0);//the last argument of 0/1 is working
		cout << particles[3][i][5] << "		" << particles[3][i][5][1] << " " << particles[3][i][5][2] << " " << particles[3][i][5][3] << endl;
		cout << actualNearestNeighbours << endl;
		for(int j=0; j<actualNearestNeighbours; j++){
			cout << "		" << nearestNeighbourAdresses[j] << " " << nearestNeighbourAdresses[j][1] << " " << nearestNeighbourAdresses[j][2] << " " << nearestNeighbourAdresses[j][3] << " " << endl;
		}
		cout << endl << endl;
	}

	
	table.free();


	//checked that program is releasing/freeing up memory properly
	radius = 1;//radius above 0.3 is OK. BELOW THAT IS GIVING SEGFAULT (irrespective of search argument of 0/1 in particleQuery)
	//the above bug is not giving pains in intel's compiler
	hash_size = NextPrime(part[0]*part[1]*part[2]);
	cout << hash_size << endl;
	table.insertParticles(particles, hash_size, part, radius);

	for(int i=0; i<5; i++){
		table.particleQuery(particles[3][i][5][1], particles[3][i][5][2], particles[3][i][5][3], nearestNeighbourAdresses, 40, &actualNearestNeighbours, radius, 0);//the last argument of 0/1 is working
		cout << particles[3][i][5] << "		" << particles[3][i][5][1] << " " << particles[3][i][5][2] << " " << particles[3][i][5][3] << endl;
		cout << actualNearestNeighbours << endl;
		for(int j=0; j<actualNearestNeighbours; j++){
			cout << "		" << nearestNeighbourAdresses[j] << " " << nearestNeighbourAdresses[j][1] << " " << nearestNeighbourAdresses[j][2] << " " << nearestNeighbourAdresses[j][3] << " " << endl;
		}
		cout << endl << endl;
	}
}