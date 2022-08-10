#ifndef RANDOMORDERVECTOR_H
#define RANDOMORDERVECTOR_H

#include <vector>

namespace MSL {
	std::vector <unsigned int> getRandomOrder (unsigned int _size);
	std::vector <unsigned int> getRandomOrder (unsigned int _start, unsigned int _end);
}

#endif
