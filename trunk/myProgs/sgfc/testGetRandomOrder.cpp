#include <iostream>
#include <fstream>
#include "RandomNumberGenerator.h"

using namespace std;
using namespace MSL;

int main () {
	RandomNumberGenerator RNG;
	cout << "========getRandomOrder(start, end)===============" << endl;
	vector <unsigned int> wector = RNG.getRandomOrder(6,12);
	for (std::vector<unsigned int>::iterator i = wector.begin(); i != wector.end(); ++i) {
		cout << *i << std::endl;
	}
*/

/*	
	cout << "========getRandomOrder(size)=====================" << endl;
	vector <unsigned int> wectorWector = getRandomOrder(10);
	for (std::vector<unsigned int>::iterator i = wectorWector.begin(); i != wectorWector.end(); ++i) {
		cout << *i << std::endl;
	}
*/
	return 0;
}
