#include <math.h>
#define M_PI 3.141592653589
// radius -> compact support radius
// position -> position where to calculate kernel value
// particle -> particle for which kernel is being calculated
float default_kernel(float position[12], float particle[12], float radius){
	float separation = pow(pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2) + pow(position[3] - particle[3], 2), 1/2.0);
	if (separation <= radius){
		float value = 315*(pow(pow(radius, 2) - pow(separation, 2), 3))/(64.0*M_PI*pow(radius, 9));
		return value;
	}
	else{
		return 0;
	}
}
//receiver should be first initialized by the following-> float receiver[3]
void default_kernel_gradient(float position[12], float particle[12], float radius, float* receiver){
	float separation = pow(pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2) + pow(position[3] - particle[3], 2), 1/2.0);
	if (separation <= radius){
		float value = -945*(pow(pow(radius, 2) - pow(separation, 2), 3))/(32.0*M_PI*pow(radius, 9));
		receiver[0] = value*(position[1] - particle[1]);
		receiver[1] = value*(position[2] - particle[2]);
		receiver[2] = value*(position[3] - particle[3]);

	}
}

float default_kernel_laplacian(float position[12], float particle[12], float radius){
	float separation = pow(pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2) + pow(position[3] - particle[3], 2), 1/2.0);
	if (separation <= radius){
		float value = -945*(pow(radius, 2) - pow(separation, 2))*(3*pow(radius, 2) - 7*pow(separation, 2))/(32.0*M_PI*pow(radius, 9));
		return value;
	}
	else{
		return 0;
	}
}



float pressure_kernel(float position[12], float particle[12], float radius){
	float separation = pow(pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2) + pow(position[3] - particle[3], 2), 1/2.0);
	if (separation <= radius){
		float value = 15.0*(pow(radius - separation, 3))/(M_PI*pow(radius, 6));
		return value;
	}
	else{
		return 0;
	}
}


/////////////THIS FUNCTION MAY CAUSE PROBLEMS DUE TO LIMITING VALUES////////////
void pressure_kernel_gradient(float position[12], float particle[12], float radius, float* receiver){
	float separation = pow(pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2) + pow(position[3] - particle[3], 2), 1/2.0);
	if (separation <= radius){
		float value = -45.0*(pow(radius - separation, 2))/(M_PI*pow(radius, 6));
		receiver[0] = value*(position[1] - particle[1])/separation;
		receiver[1] = value*(position[2] - particle[2])/separation;
		receiver[2] = value*(position[3] - particle[3])/separation;

	}
}
////////////////////////////////////////////////////////////////////////////////
/////////////THIS FUNCTION MAY CAUSE PROBLEMS DUE TO LIMITING VALUES////////////
float pressure_kernel_laplacian(float position[12], float particle[12], float radius){
	float separation = pow(pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2) + pow(position[3] - particle[3], 2), 1/2.0);
	if (separation <= radius){
		float value = -90.0*(radius - separation)*(radius - 2*separation)/(M_PI*separation*pow(radius, 6));
		return value;
	}
	else{
		return 0;
	}
}
////////////////////////////////////////////////////////////////////////////////



/////////////THIS FUNCTION MAY CAUSE PROBLEMS DUE TO LIMITING VALUES////////////
float viscosity_kernel(float position[12], float particle[12], float radius){
	float separation = pow(pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2) + pow(position[3] - particle[3], 2), 1/2.0);
	if (separation <= radius){
		float value1 = 15/(2.0*M_PI*pow(radius, 3));
		float value2 = (-pow(separation, 3)/(2*pow(radius, 3))) + (pow(separation, 2)/pow(radius, 2)) + (radius/(2*separation)) -1;
		return value1*value2;
	}
	else{
		return 0;
	}
}
////////////////////////////////////////////////////////////////////////////////
/////////////THIS FUNCTION MAY CAUSE PROBLEMS DUE TO LIMITING VALUES////////////
void viscosity_kernel_gradient(float position[12], float particle[12], float radius, float* receiver){
	float separation = pow(pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2) + pow(position[3] - particle[3], 2), 1/2.0);
	if (separation <= radius){
		float value1 = -15/(2.0*M_PI*pow(radius, 3));
		float value2 = (-3.0*separation/(2.0*pow(radius, 3))) + 2.0/pow(radius, 2) - (radius/(2.0*pow(separation, 3)));
		receiver[0] = value1*value2*(position[1] - particle[1]);
		receiver[1] = value1*value2*(position[2] - particle[2]);
		receiver[2] = value1*value2*(position[3] - particle[3]);

	}
}
////////////////////////////////////////////////////////////////////////////////

float viscosity_kernel_laplacian(float position[12], float particle[12], float radius){
	float separation = pow(pow(position[1] - particle[1], 2) + pow(position[2] - particle[2], 2) + pow(position[3] - particle[3], 2), 1/2.0);
	if (separation <= radius){
		float value = 45.0*(radius - separation)/(M_PI*separation*pow(radius, 6));
		return value;
	}
	else{
		return 0;
	}
}
