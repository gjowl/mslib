#ifndef HELIXDIMERCLUSTER_H
#define HELIXDIMERCLUSTER_H

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "DeadEndElimination.h"
#include "MonteCarloManager.h"
#include "SelfConsistentMeanField.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "AtomSelection.h"
#include "AtomContainer.h"
#include "FormatConverter.h"
#include "CRDReader.h"
#include "CRDWriter.h"
#include "SysEnv.h"
#include "HelixDimer.h"


class HelixDimerCluster {
	public:
	HelixDimerCluster(HelixDimer * _structure);

	void addHelixDimer(HelixDimer * _structure);

	AtomPointerVector& getAtomPointers();

	void setDetails(string _origSeq, string _seq, string _uName, string _uAccession,string _outputDir,int _resStart, int _resEnd, int _tmStart, int _tmEnd);

	void convertToPdbNames();
	void printHelixDimerClusterPdbs(int id);
	void printHelixDimerClusterCrds(int id, bool allStructures, bool writeAxis);

	void printDetails(int id);

	void makePse(int id);

	vector<HelixDimer*>& getMembers() {return members;}

	string getOrigSequence() {return origSeq;}
	string getModeledSequence() {return seq;}

	string getUniprotName() {return uniprotName;}
	string getUniprotAccession(){return uniprotAccession;}

	int getResStart() {return resStart;}
	int getResEnd() {return resEnd;}

	int getTmStart() {return tmStart;}
	int getTmEnd() {return tmEnd;}

	void printTermEnergies(int id);

	private:
	vector<HelixDimer*> members; // should be added in sorted order
	string origSeq;
	string seq;
	string uniprotName;
	string uniprotAccession;
	string outputDir;
	int resStart;
	int resEnd;
	int tmStart;
	int tmEnd;

};

# endif