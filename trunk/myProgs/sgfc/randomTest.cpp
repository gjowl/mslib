#include "RandomNumberGenerator.h"

using namespace std;
using namespace MSL;

int main() {
	RandomNumberGenerator RNG;
	RNG.setTimeBasedSeed();
	for (uint i = 0; i < 100; i++) {
		vector<uint> random = RNG.getRandomOrder(10);
		for (uint j = 0; j < random.size(); j++) {
			cout << random[j] << " ";
		}
		cout << endl;
		cout << "================================="  << endl;
	}
	return 0;
}
