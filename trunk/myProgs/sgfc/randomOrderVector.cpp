#include "randomOrderVector.h"
#include "RandomNumberGenerator.h"
#include <iostream>


using namespace std;
using namespace MSL;

/*
This program will generate a vector of integers in a random order.
This will be useful for algorithms like "greedy", which loop through positions of a protein randomly
*/


std::vector <unsigned int> getRandomOrder (unsigned int _size) {
	unsigned int _start = 1;
	return getRandomOrder (_start, _size); 
}


std::vector<unsigned int> getRandomOrder (unsigned int _start, unsigned int _end) {
	/*Debug comments to terminal
	cout << "==============================" << endl;
	cout << "Ordered vector" << endl;
	cout << "==============================" << endl;
	*/

	RandomNumberGenerator RNG;
	std::vector <unsigned int> ordered;
	std::vector <unsigned int> random;
	for (unsigned int i = _start; i <= _end; i++) {
		ordered.push_back (i);

	//	cout << ordered.back() << endl;
	}

	/*Debug comments to terminal

	cout << "==============================" << endl;
	cout << "Random vector" << endl;
	cout << "==============================" << endl;
	*/

	while (!ordered.empty()) {
		double randNum = RNG.getRandomDouble() * ordered.size();
		int randInt = randNum;
		random.push_back (ordered[randInt]);
		ordered.erase(ordered.begin() + randInt);
	//	cout << random.back() << endl;
	}

	return random;

}

