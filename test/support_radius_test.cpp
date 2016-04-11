#include <iostream>
#include <math.h>
using namespace std;

int main(void){
	float compactsupportrad;
	float volume_of_fluid = 240;
	float number_of_particles = 240;
	compactsupportrad = pow(3*volume_of_fluid*40/((4*M_PI*number_of_particles)),1/3.0);
	cout << M_PI << endl;
	cout << compactsupportrad << endl;
}