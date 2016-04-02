#include <cstddef>
#include <math.h>

#define P1 73856093
#define P2 19349663
#define P3 83492791

struct list{
public:
	float* value;
	list* nextNode;
};

struct hashTable
{
public:
	long int size;
	list** table;
	bool* boolArray; // This informs if any list has ever been inserted into hashTable.table[value]

	void createHashTable(long int numberOfElements){
		this->table = new list*[numberOfElements];
		this->size = numberOfElements;
		this->boolArray = new bool[numberOfElements];
		for(int i=0; i<numberOfElements; i++){
			boolArray[i] = 0;
		}
	}

	long int hashValue(float hashKeyX, float hashKeyY, float hashKeyZ, float gridCellSize){//gridCellSize is same as compact support radius
		int x_mod, y_mod, z_mod;
		x_mod = 1 + (long int)(hashKeyX/((float)gridCellSize));//http://matthias-mueller-fischer.ch/publications/tetraederCollision.pdf
		y_mod = 1 + (long int)(hashKeyY/((float)gridCellSize));//The x, y, z queries have been divided by gridsize and rounded off to next integer
		z_mod = 1 + (long int)(hashKeyZ/((float)gridCellSize));
		return (((long int)((x_mod*P1)^(y_mod*P2)^(z_mod*P3)))%size);
	}

	void insertElement(float* address, long int value){


		//is there a requirement of "this->" ??


		list* new_item;
		new_item = new list;
		new_item->value = address;
		new_item->nextNode = NULL;

		if(this->boolArray[value] == 0){
			this->table[value] = new_item;
			this->boolArray[value] = 1;
		}
		else if(this->boolArray[value] == 1){
			list* temp = this->table[value];
			while(1){
				if(temp->nextNode == NULL){
					temp->nextNode = new_item;
					break;
				}
				else{
					temp = temp->nextNode;
				}
			}
		}
	}
	/*//Will be implemented when required
	void lookup(){

	}
	*/

	/*
	For particle search->
	 _________+_________+_________+
	|		  |  		|		  |
	|		  |			|		  |
	|		  |			|		  |
	|_________+_________o_________+
	|		  |			|		  |
	|		  |			|		  |
	|		  |  *      |		  |
	|_________+_________+_________+
	|		  |			|		  |
	|		  |			|		  |
	|		  |         |		  |
	|_________|_________|_________|

	* -> location of particle
	o -> location of x_mod, y_mod (1+hashKeyX/gridCellSize, 1+hashKeyY/gridCellSize)
	+ -> locations we want to search

	*/
	void insertParticles(float**** particles, int tablesize, int* number_of_particles_array, float compactSupportRadius){
		createHashTable(tablesize);
		float x, y, z;
		long int hashvalue;
		for(int i=0; i<number_of_particles_array[0]; i++){
			for(int j=0; j<number_of_particles_array[1]; j++){
				for(int k=0; k<number_of_particles_array[2]; k++){
					if(particles[i][j][k][0] == 1){//if this is a particle
						x = particles[i][j][k][1];
						y = particles[i][j][k][2];
						z = particles[i][j][k][3];
						hashvalue = hashValue(x, y, z, compactSupportRadius);
						insertElement(particles[i][j][k], hashvalue);
					}
				}
			}
		}
	}
	void particleQuery(float x, float y, float z, float** nearestNeighbourList, int estimatedNumNearestNeighbours, int* numNearestNeighbours, float compactSupportRadius){
		int BB_min_x, BB_min_y, BB_min_z;
		int BB_max_x, BB_max_y, BB_max_z;//BB stands for bounding box. see page number 43 in the paper
		BB_min_x = (int)(x/((float)compactSupportRadius));
		BB_max_x = 2 + (int)(x/((float)compactSupportRadius));
		BB_min_y = (int)(y/((float)compactSupportRadius));
		BB_max_y = 2 + (int)(y/((float)compactSupportRadius));
		BB_min_z = (int)(z/((float)compactSupportRadius));
		BB_max_z = 2 + (int)(z/((float)compactSupportRadius));

		long int hvalue;
		float* location;
		list* tempNode;
		float tempDist;

		int tempNumNearestNeighbours = 0;

		for(int i=BB_min_x; i<=BB_max_x; i++){
			for(int j=BB_min_y; j<=BB_max_y; j++){
				for(int k=BB_min_z; k<=BB_max_z; k++){
					hvalue = (((long int)((i*P1)^(j*P2)^(k*P3)))%size);
					if(boolArray[hvalue] == 1){
						tempNode = table[hvalue];
						while(1){
							if(tempNode->nextNode == NULL){								
								location = tempNode->value;
								tempDist = sqrt((x - location[1])*(x - location[1]) + (y - location[2])*(y - location[2]) + (z - location[3])*(z - location[3]));
								if((tempDist <= compactSupportRadius) && (tempNumNearestNeighbours < estimatedNumNearestNeighbours)){
									nearestNeighbourList[tempNumNearestNeighbours] = location;
									tempNumNearestNeighbours += 1;
								}
								break;
							}
							else{
								location = tempNode->value;
								tempDist = sqrt((x - location[1])*(x - location[1]) + (y - location[2])*(y - location[2]) + (z - location[3])*(z - location[3]));
								if((tempDist <= compactSupportRadius) && (tempNumNearestNeighbours < estimatedNumNearestNeighbours)){
									nearestNeighbourList[tempNumNearestNeighbours] = location;
									tempNumNearestNeighbours += 1;
								}
								tempNode = tempNode->nextNode;
							}
						}
					}
				}
			}
		}
		*numNearestNeighbours = tempNumNearestNeighbours;

	}
	void free(void){
		//free up table
		// NOTE-> instead of delete, you can use free as well
		list** tempPointer1, **tempPointer2;
		for(int i=0; i<size; i++){
			if(boolArray[i] == 1){
				tempPointer1 = &(table[i]);
				if((*tempPointer1)->nextNode != NULL){
					tempPointer2 = &((*tempPointer1)->nextNode);
				}
				while(1){
					if((*tempPointer1)->nextNode == NULL){
						delete (*tempPointer1);
						break;
					}
					else{
						delete (*tempPointer1);
						tempPointer1 = tempPointer2;
						if((*tempPointer1)->nextNode != NULL){
							tempPointer2 = &((*tempPointer1)->nextNode);
						}
					}
				}

			}
		}
		delete [] table;
		//free up bool array
		delete [] boolArray;
	}
};







// For testing hash table implementation
/* 
#include <iostream>
using namespace std;

int main(){
	hashTable a;
	float x, y, z;
	x = 1.1;
	y = 1.2;
	z = 1.3;
	a.createHashTable(3);
	a.insertElement(&x, 1);
	a.insertElement(&y, 1);
	a.insertElement(&z, 1);
	
	cout << a.size << endl;
	cout << a.boolArray[0] << endl << a.boolArray[1] << endl << a.boolArray[2] << endl;
	cout << a.table[1]->value << endl << a.table[1]->nextNode->value << endl << a.table[1]->nextNode->nextNode->value << endl;
	cout << *(a.table[1]->value) << endl << *(a.table[1]->nextNode->value) << endl << *(a.table[1]->nextNode->nextNode->value) << endl;
}
*/