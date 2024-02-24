/**
 * @Author: Gilbert Loiseau
 * @Date:   2022/02/12
 * @Email:  gjowl04@gmail.com
 * @Filename: design.cpp
 * @Last modified by:   Gilbert Loiseau
 * @Last modified time: 2022/02/22
 */

#include <sstream>
#include <iterator>
#include <unistd.h>
#include "designFunctions.h"
#include "functions.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;

map<string,string> getChainSeqMap(System &_sys) {
	map<string,string> sequences;
	for (uint i=0; i<_sys.getChainSize(); i++){
		Chain &chain = _sys.getChain();
		string chainId = chain.getChainId();
		string seq = chain.toString();
		sequences[chainId] = seq;
	}
	return sequences;
}

void setActiveSequence(System &_sys, string _sequence, string _chainId){
	// Set the active sequence for the system
	Chain &chain = _sys.getChain(_chainId);
	for (uint i=0; i<chain.chainSize(); i++){
		Position &pos = _sys.getPosition(i);
		string posId = posA.getPositionId();
		string aa = MslTools::getThreeLetterCode(_sequence.substr(i, 1));
		_sys.setActiveIdentity(posId, aa);
	}
}