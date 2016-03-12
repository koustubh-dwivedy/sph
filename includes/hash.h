#include <cstddef>

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
	int size;
	list** table;
	bool* boolArray; // This informs if any list has ever been inserted into hashTable.table[value]

	void createHashTable(int numberOfElements){
		this->table = new list*[numberOfElements];
		this->size = numberOfElements;
		this->boolArray = new bool[numberOfElements];
	}

	unsigned int hashValue(int hashKeyX, int hashKeyY, int hashKeyZ, int tableSize){
		return (((hashKeyX*P1)^(hashKeyY*P2)^(hashKeyZ*P3))%tableSize);
	}

	void insertElement(float* address, int value){
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
	void free(){
		
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