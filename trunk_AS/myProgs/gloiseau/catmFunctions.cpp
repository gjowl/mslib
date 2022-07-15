#include <sstream>
#include <iterator>
#include <unistd.h>
#include "catmFunctions.h"
#include "functions.h"

using namespace std;
using namespace MSL;

static SysEnv ENV;

void printOptions(catmOptions & _op, ofstream & _fout) {
	//TODO: when reorganizing I had to get rid of the below: can I put that info into a file that all of my code for this program inherits from?
	//_fout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;

	_fout << "backboneCrd " << _op.backboneCrd << endl;
	_fout << "pdbOutputDir " << _op.pdbOutputDir << endl;

	_fout << "fullSequence " << _op.fullSequence << endl;
	_fout << "tmStart " << _op.tmStart << endl;
	_fout << "tmEnd " << _op.tmEnd << endl;

	_fout << "helixGeoFile " << _op.helixGeoFile << endl;
	_fout << "rulesFile " << _op.rulesFile << endl;

	_fout << "topFile " << _op.topFile << endl;
	_fout << "parFile " << _op.parFile << endl;
	_fout << "solvFile " << _op.solvFile << endl;
	_fout << "hBondFile " << _op.hBondFile << endl;
	_fout << "rotLibFile " << _op.rotLibFile << endl;
	_fout << "monoRotLibFile " << _op.monoRotLibFile << endl;

	_fout << "MCCycles " << _op.MCCycles << endl;
	_fout << "MCMaxRejects " << _op.MCMaxRejects << endl;
	_fout << "MCStartTemp " << _op.MCStartTemp << endl;
	_fout << "MCEndTemp " << _op.MCEndTemp << endl;
	_fout << "MCCurve " << _op.MCCurve << endl;

	_fout << "deltaZ " << _op.deltaZ << endl;
	_fout << "deltaAx " << _op.deltaAx << endl;
	_fout << "deltaCross " << _op.deltaCross << endl;
	_fout << "deltaX " << _op.deltaX << endl;

	_fout << "verbose " << _op.verbose << endl;
	_fout << "greedyCycles " << _op.greedyCycles << endl;
	_fout << "seed " << _op.seed << endl;

	_fout << "numberOfStructuresToMCRepack " << _op.numberOfStructuresToMCRepack << endl;
	_fout << "energyCutOff " << _op.energyCutOff << endl;

	_fout << "uniprotName " << _op.uniprotName << endl;
	_fout << "uniprotAccession " << _op.uniprotAccession << endl;

	_fout << "monoE_vdw " << _op.monoE_vdw << endl;
	_fout << "monoE_solv " << _op.monoE_solv << endl;
	_fout << "monoE_solvRef" << _op.monoE_solvRef << endl;
	_fout << "monoE_hbond" << _op.monoE_hbond << endl;

	_fout << "rmsdCutoff " << _op.rmsdCutoff << endl;
	_fout << "clusterSolutions " << _op.clusterSolutions << endl;
	_fout << "printAllCrds " << _op.printAllCrds << endl;
	_fout << "printAxes " << _op.printAxes << endl;
	_fout << "printTermEnergies " << _op.printTermEnergies << endl;
	_fout << "deleteTerminalBonds " << _op.deleteTerminalBonds << endl;
	_fout << "deleteTerminalInteractions: ";
	for (uint i=0; i<_op.deleteTerminalInteractions.size(); i++){
		_fout << _op.deleteTerminalInteractions[i] << " ";
	}
	_fout << endl;

	_fout << "fullSequenceStart " << _op.fullSequenceStart << endl;

	_fout << "startResNum " << _op.startResNum << endl;
	_fout << "endResNum " << _op.endResNum << endl;

	_fout << "threadStart " << _op.threadStart << endl;
	_fout << "threadEnd " << _op.threadEnd << endl;

	_fout << "weight_elec " << _op.weight_elec << endl;
	_fout << "weight_vdw " << _op.weight_vdw << endl;
	_fout << "weight_hbond " << _op.weight_hbond << endl;
	_fout << "weight_solv " << _op.weight_solv << endl;

	if(_op.configfile != "") {
		_fout << "configfile " << _op.configfile << endl;
	}

	_fout << endl;

}

unsigned int CAOfClosestApproach(Chain & _chainA, Chain & _chainB) {
	double minCAtoCAdist = MslTools::doubleMax;
	string closestCA_A = "";
	string closestCA_B = "";

	unsigned int closestAtomIndex = 0;
	for (unsigned int i=0; i < _chainA.positionSize(); i++) {
		Position& posA = _chainA.getPosition(i);
		Position& posB = _chainB.getPosition(i);
		if (posA.atomExists("CA") && posB.atomExists("CA")) {
			Atom& posACA = posA.getLastFoundAtom();
			Atom& posBCA = posB.getLastFoundAtom();
			double dist = posACA.distance(posBCA);
			if (dist < minCAtoCAdist) {
				minCAtoCAdist = dist;
				closestAtomIndex = i;
			}
		}
	}

	return _chainA(closestAtomIndex).getResidueNumber();
}

map<string, unsigned int> interfaceResidueCheck(AtomPointerVector & _chainA, AtomPointerVector & _chainB) {
	map<string, unsigned int> atomsWithIn4AOfPosition;
	for (uint i=0; i < _chainA.size(); i++) {
		if (_chainA[i]->getName() != "CA" && _chainA[i]->getName() != "C" && _chainA[i]->getName() != "N" && _chainA[i]->getName() != "HN" && _chainA[i]->getName() != "O") {
			for (uint j=0; j < _chainB.size(); j++) {
				if (_chainB[j]->getName() != "CA" && _chainB[j]->getName() != "C" && _chainB[j]->getName() != "N" && _chainB[j]->getName() != "HN" && _chainB[j]->getName() != "O") {
					if (_chainA[i]->distance(*_chainB[j]) < 4.0) {
						if (atomsWithIn4AOfPosition.find(_chainA[i]->getIdentityId()) != atomsWithIn4AOfPosition.end()) {
							atomsWithIn4AOfPosition[_chainA[i]->getIdentityId()]++;
						}
						else {
							atomsWithIn4AOfPosition[_chainA[i]->getIdentityId()] = 1;
						}
					}
				}
			}
		}
	}

	//for(map<string, unsigned int>::iterator it=atomsWithIn4AOfPosition.begin(); it != atomsWithIn4AOfPosition.end(); it++) {
	//	_fout << "pos: " << it->first << " count: " << it->second << endl;
	//}

	return atomsWithIn4AOfPosition;
}

void reThreadResidues(vector<Position*> & positions, int offset) {
	for(int i = 0; i < positions.size(); i++) {
		int resNum = positions[i]->getResidueNumber();
		positions[i]->setResidueNumber(resNum + offset);
	}
}

bool hydrogenBondCheck(System & _sys, catmOptions &_opt, vector<string> _parsedGeoInformation, double & _xShiftStart) {

	vector<HbondInfo*> hbondList;
	for(unsigned k = 5; k < _parsedGeoInformation.size(); k+=5 ) {
		HbondInfo * t = new HbondInfo(_parsedGeoInformation[k],_parsedGeoInformation[k+1],_parsedGeoInformation[k+2],_parsedGeoInformation[k+3],_parsedGeoInformation[k+4]);
		if(t->isValid(_sys)) {
			hbondList.push_back(t);
		} else {
			delete t;
		}
	}


	/*
	for(int i = 0; i < hbondList.size(); i++) {
		hbondList[i]->print(fout);
	}
	*/
	if(hbondList.size() < _opt.hbondCheckNumber) {
		for(int i = 0; i < hbondList.size(); i++) {
			delete hbondList[i];
		}
		return false;
	} else {
		_xShiftStart = hbondList[hbondList.size() - _opt.hbondCheckNumber]->distBreak - 0.1;
		for(int i = 0; i < hbondList.size(); i++) {
			delete hbondList[i];
		}
		return true;
	}

}

void renumberResidues(System& _sys, int _startResNum) {
	// just set the index + startResNum of each position as its residue number
	vector<Position*>& positions = _sys.getPositions();
	for(int i = 0; i < positions.size(); i++) {
		int resNum = positions[i]->getIndexInChain();
		positions[i]->setResidueNumber(resNum + _startResNum);
	}
}

bool rulesCheck(System & _sys, string _geoIndex, map<int, string> _rulesMap) {

	int idx = MslTools::toInt(_geoIndex);
	if (_rulesMap.find(idx) == _rulesMap.end()) {
		return true;
	}

	vector<string> parsedRules = MslTools::tokenize(_rulesMap[idx], ",");
	// read over all the rules for a given model
	for (uint m = 1; m < parsedRules.size(); m++) {
		int delimiter  = parsedRules[m].find_first_of(":");
		string position = parsedRules[m].substr(0,1)+","+parsedRules[m].substr(1,delimiter-1); // makes string A,35
		string rule = parsedRules[m].substr(delimiter+1,(parsedRules[m].length()-delimiter-1));
		if(!_sys.positionExists(position)) {
			continue;
		}
		//fout << "position " << position << " must be " << rule << " , " << " is actually " << _sys.getPosition(position).getResidueName() << endl;

		// Check if there is a ! "not"
		size_t findNot = rule.find("!");

		// if there is a "not"
		if (findNot == 0) { // ! will only appear in position 0
			rule = rule.substr(2);
			rule = MslTools::trim(rule, "]"); // residues that cannot exist in given position
			// look at position given in rule
			string resName = MslTools::getOneLetterCode(_sys.getPosition(position).getResidueName());
			for (uint n=0; n < rule.length(); n++) {
				if (rule.substr(n, 1) == resName) {
					return false;
				}
			}

		}
		else { // no "not" found
			rule = rule.substr(1);
			rule = MslTools::trim(rule, "]"); // residues that must exist in given position
			// look at position given in rule
			string resName = MslTools::getOneLetterCode(_sys.getPosition(position).getResidueName());
			bool found = false;
			for (uint n=0; n < rule.length(); n++) {
				if (rule.substr(n, 1) == resName) {
					found = true;
					break;
				}
			}
			if(!found) {
				return false;
			}
		}
	}

	return true;
}

void repackSideChains(SelfPairManager & _spm, int _greedyCycles, vector<vector<vector<vector<bool> > > > _savedEnergyFlagTable) {

	_spm.setOnTheFly(1);
	//_spm.recalculateNonSavedEnergies(_savedEnergyFlagTable);
	_spm.calculateEnergies(); // CHANGE BACK!!!
	_spm.runGreedyOptimizer(_greedyCycles);
}

void addStructureToTopMCRepackList(vector<vector<double> > & _topMCRepackVector, vector<double> & _currentStuctureInfo, ofstream & _fout, unsigned int _numberOfStructuresToMCRepack) {
	//cout << "size of top 5 vector: " << _topMCRepackVector.size() << ", energy to be added: " << _currentStuctureInfo[0] << endl;
	if(_topMCRepackVector.size() < _numberOfStructuresToMCRepack) {
		_topMCRepackVector.push_back(_currentStuctureInfo);
		_fout << "Structure added to Monte Carlo Repack List" << endl;
	}
	else{
		double currentWorstEnergy = -1.0*MslTools::doubleMax;
		uint positionOfWorstEnergy = _topMCRepackVector.size() + 1;
		for(uint k = 0; k<_topMCRepackVector.size(); k++) {
			//cout << "pos: " << k << "\tcurrent energies: " << _topMCRepackVector[k][0];
			if (_topMCRepackVector[k][0] > _currentStuctureInfo[0] && _topMCRepackVector[k][0] > currentWorstEnergy) {
				positionOfWorstEnergy = k;
				currentWorstEnergy = _topMCRepackVector[k][0];
				//cout << "*";
			}
			//cout << endl;
		}
		if (positionOfWorstEnergy <= _topMCRepackVector.size()) {
			_topMCRepackVector.erase(_topMCRepackVector.begin()+positionOfWorstEnergy);
			_topMCRepackVector.push_back(_currentStuctureInfo);
			_fout << "Structure added to Monte Carlo Repack List" << endl;
		}
	}

}

void deleteTerminalBondInteractions(System &_sys, catmOptions &_opt) {
	// look at all hbond interactions and remove those
	// remove any interaction that has a donor or acceptor from residues 1 2 3 and n-2 n-1 n on each chain that are not part of the TM

	EnergySet* pESet = _sys.getEnergySet();
	int chainSize = _sys.chainSize(); // number of positions in chain A = chain B
	// compute the extensions at the beginning and the end
	int frontExt = _opt.tmStart - _opt.startResNum;
	int endExt = _opt.endResNum - _opt.tmEnd;
	AtomPointerVector atoms;
	for(int i = 0; i < chainSize; i++) {
		Chain& thisChain = _sys.getChain(i);
		vector<Position*>& positions = thisChain.getPositions();
		for(int i = 0; i < 3; i++) {
			if(frontExt > i) {
				atoms += positions[i]->getAtomPointers();
				//cout << "Removing Hbonds from " << positions[i]->getPositionId()  << endl;
			}
			if(endExt > i) {
				atoms += positions[positions.size() - 1 - i]->getAtomPointers();
				//cout << "Removing Hbonds from " << positions[positions.size() - 1 - i]->getPositionId()  << endl;
			}
		}
	}
	for (uint i=0; i<_opt.deleteTerminalInteractions.size(); i++){
		pESet->deleteInteractionsWithAtoms(atoms,_opt.deleteTerminalInteractions[i]);
	}
}

void readRulesFile(string _fileName, map<int, string> & _rulesFileMap) {
	ifstream file;
	file.open(_fileName.c_str());
	if(!file.is_open()) {
		cerr << "Unable to open " << _fileName << endl;
		exit(0);
	}

	string tmpRulesLine;

	while(file) {
		getline(file, tmpRulesLine);
		if (tmpRulesLine.length() > 1) {
			vector<string> token = MslTools::tokenizeAndTrim(tmpRulesLine,",");
			_rulesFileMap[MslTools::toInt(token[0])] = tmpRulesLine;
		}
	}
	file.close();
}

vector<string> getInterHelicalHbonds(EnergySet* & _ESet) {
	unsigned int numHbonds = 0;
	// Why are we doing this?
	//_ESet->setAllTermsInactive();
	//_ESet->setTermActive("SCWRL4_HBOND", true);
	vector<Interaction*> hbondInteractions = (*(_ESet->getEnergyTerms()))["SCWRL4_HBOND"];

	vector<string> hbonds;
	for(int i = 0; i < hbondInteractions.size(); i++) {
		vector<Atom*> atoms =  hbondInteractions[i]->getAtomPointers();
		if(atoms[0]->getChainId() == atoms[2]->getChainId()) {
			continue;
		}
		double e = hbondInteractions[i]->getEnergy();
		if(e < 0) {
			//_fout << atoms[0]->getAtomId() << " " << atoms[2]->getAtomId() << endl;
			hbonds.push_back(atoms[0]->getAtomOfIdentityId() + ":" + atoms[2]->getAtomOfIdentityId() + "=" + MslTools::doubleToString(atoms[0]->distance(*atoms[2])) + "," + MslTools::doubleToString(e));
			numHbonds++;
		}
	}
	// Why are we doing this?
	//_ESet->setAllTermsActive();
	//_ESet->setTermActive("CHARMM_ELEC", false);
	return hbonds;
}

void readGeometryFile(string _filename, vector<string>& _fileVec) {
	ifstream file;
	file.open(_filename.c_str());
	if(!file.is_open()) {
		cerr << "Unable to open " << _filename << endl;
		exit(0);
	}

	string parameterList;

	while(file) {
		getline(file, parameterList);
		_fileVec.push_back(parameterList);
	}
	file.close();
}


void clusterSolutions(vector<HelixDimerCluster*>& clusters,vector<HelixDimer*>& _structures, double _rmsdCutoff, string _origSeq, string _builtSeq, catmOptions &_opt) {
	for(int i  = 0; i < _structures.size(); i++) {
		bool diff = true;
		AtomPointerVector& a1 = _structures[i]->getAtomPointers();
		AtomSelection sel1(a1);
		AtomPointerVector& ca1 = sel1.select("name CA");
		//cout << i << endl;
		for(int j = 0; j < clusters.size(); j++) {
			AtomPointerVector& a2 = clusters[j]->getAtomPointers();
			AtomSelection sel2(a2);
			AtomPointerVector& ca2 = sel2.select("name CA");
			//cout << _structures[i]->getName() << " " << _structures[j]->getName() << " ";
			double r = ca1.rmsd(ca2);
			if(r < _rmsdCutoff) {
				clusters[j]->addHelixDimer(_structures[i]);
				diff = false;
				break;
			}
		}
		if(diff) {
			clusters.push_back(new HelixDimerCluster(_structures[i]));
			clusters.back()->setDetails(_origSeq, _builtSeq, _opt.uniprotName,_opt.uniprotAccession,_opt.pdbOutputDir,_opt.startResNum,_opt.endResNum,_opt.tmStart,_opt.tmEnd);
		}
	}

}

vector<vector<vector<vector<bool> > > > createSavedEnergyFlagTable(System & _sys) {
	vector<vector<vector<vector<bool> > > > savedEnergyFlagTable; // true if the energy is computed already and available
	vector<unsigned int> variablePos = _sys.getVariablePositions();

	for (int i=0; i<variablePos.size(); i++) {
		savedEnergyFlagTable.push_back(vector<vector<vector<bool> > >());
		//cout << i;

		for (int j=0; j < _sys.getPosition(variablePos[i]).getTotalNumberOfRotamers(); j++) {
			savedEnergyFlagTable[i].push_back(vector<vector<bool> >());
			//cout << " " << i << "/" << j << "-";

			for (int k=0; k < i; k++) {
				if (i==k) {
					continue;
				}
				savedEnergyFlagTable[i][j].push_back(vector<bool>());

				for (int l=0; l < _sys.getPosition(variablePos[k]).getTotalNumberOfRotamers(); l++) {
			//		cout << " " << k << "/" << l;
					if (_sys.getPosition(variablePos[i]).getChainId() == _sys.getPosition(variablePos[k]).getChainId()) {
						savedEnergyFlagTable[i][j][k].push_back(true);
			//			cout << "=T";
					}
					else {
						savedEnergyFlagTable[i][j][k].push_back(false);
			//			cout << "=F";
					}
				}
			}
		}
		//cout << endl;
	}

	return savedEnergyFlagTable;
}

double computeMonomerEnergy(System & _sys, Transforms & _trans, catmOptions & _opt, System & _helicalAxis, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _fout, int _greedyCycles, int _MCCycles, int _MCMaxRejects) {

	time_t startTimeMono, endTimeMono;
	double diffTimeMono;
	time(&startTimeMono);

	Chain & inputChain = _sys.getChain(0);
	AtomPointerVector &axisA = _helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = _helicalAxis.getChain("B").getAtomPointers();

	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile, _opt.solvFile);
	CSBMono.setBuildTerm("CHARMM_ELEC", true);
	CSBMono.setBuildTerm("CHARMM_ANGL", false);
	CSBMono.setBuildTerm("CHARMM_BOND", false);
	CSBMono.setBuildTerm("CHARMM_DIHE", false);
	CSBMono.setBuildTerm("CHARMM_IMPR", false);
	CSBMono.setBuildTerm("CHARMM_U-BR", false);

	CSBMono.setSolvent("MEMBRANE");
	CSBMono.setIMM1Params(15, 10);
	CSBMono.buildSystemFromPDB(inputChain.getAtomPointers());

	SystemRotamerLoader monoRot(monoSys, _opt.monoRotLibFile);
	monoRot.defineRotamerSamplingLevels();

	// Add hydrogen bond term
	HydrogenBondBuilder monohb(monoSys, _opt.hBondFile);
	monohb.buildInteractions(30);

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* monoEset = monoSys.getEnergySet();
	monoEset->setAllTermsActive();
	monoEset->setTermActive("CHARMM_ELEC", true);
	monoEset->setTermActive("CHARMM_ANGL", false);
	monoEset->setTermActive("CHARMM_BOND", false);
	monoEset->setTermActive("CHARMM_DIHE", false);
	monoEset->setTermActive("CHARMM_IMPR", false);
	monoEset->setTermActive("CHARMM_U-BR", false);

	deleteTerminalBondInteractions(monoSys, _opt);

	/******************************************************************************
	 *              === LOAD ROTAMERS FOR MONOMER & SET-UP SPM ===
	 ******************************************************************************/
	SelfPairManager monoSpm;
	monoSpm.seed(_opt.seed);
	monoSpm.setVerbose(_opt.verbose);

	for (uint k=0; k < monoSys.positionSize(); k++) {
		Position &pos = monoSys.getPosition(k);

		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!monoRot.loadRotamers(&pos, pos.getResidueName(), "SL97.00")) {
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}

	monoEset->setWeight("CHARMM_ELEC", _opt.weight_elec);
	monoEset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	monoEset->setWeight("SCWRL4_HBOND", _opt.weight_hbond);
	monoEset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
	monoEset->setWeight("CHARMM_IMM1", _opt.weight_solv);
	_fout << "Monomer - ELEC weight: " << monoEset->getWeight("CHARMM_ELEC") << " VDW weight: " << monoEset->getWeight("CHARMM_VDW") << " HB weight: " << monoEset->getWeight("SCWRL4_HBOND") << " IMM1REF weight: " << monoEset->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << monoEset->getWeight("CHARMM_IMM1") << endl;

	monoSpm.setSystem(&monoSys);
	monoSpm.updateWeights();

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	AtomPointerVector &chainA = monoSys.getAtomPointers();

	/******************************************************************************
	 *                     === SHIFT HELICES INTO MEMBRANE ===
	 ******************************************************************************/
	CartesianPoint moveAxisBOneAngstrom;
	moveAxisBOneAngstrom.setCoor(1.0, 0.0, 0.0);
	_trans.translate(axisB, moveAxisBOneAngstrom);

	monoSys.calcEnergy();

	// move center of mass to origin
	//moveZCenterOfCAMassToOrigin(chainA, _helicalAxis.getAtomPointers(), _trans);
	AtomSelection sel(chainA);
	AtomPointerVector & caApV = sel.select("name CA");
	double centerHelix = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		centerHelix += (caApV[i]->getCoor()).getZ();
	}
	centerHelix = -1.0 * centerHelix/double(caApV.size());

	CartesianPoint interDistVect;
	interDistVect.setCoor(0.0, 0.0, centerHelix);
	_trans.translate(chainA, interDistVect);


	// Initial Z Shift move -5A down
	CartesianPoint zUnitVector;
	zUnitVector.setCoor(0.0, 0.0, 1.0);

	CartesianPoint move5Down = zUnitVector * -5.0;
	_trans.translate(chainA, move5Down);
	double bestZ = -5.0;

	monoSys.calcEnergy();

	// Repack side chains
	monoSpm.setOnTheFly(1);
	monoSpm.calculateEnergies();
        monoSpm.runGreedyOptimizer(_greedyCycles);

	double currentEnergy = monoSpm.getMinBound()[0];
	double bestEnergy = currentEnergy;
	monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
	monoSys.saveAltCoor("savedBestMonomer");
	_helicalAxis.saveAltCoor("BestMonomerAxis");
	_fout << "current Z: -5 Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

	// Test -5 to +5A shifts in Membrane
	for(int i=0; i<=10; i++) {

		_trans.translate(chainA, zUnitVector);

		double currentZ = -5.0 + ((i+1)*1.0);
		monoSpm.calculateEnergies();
		monoSpm.runGreedyOptimizer(_greedyCycles);
		currentEnergy = monoSpm.getMinBound()[0];
		_fout << "current Z: " << currentZ << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix

		if(currentEnergy < bestEnergy) {
			bestEnergy = currentEnergy;
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			bestZ = -5.0 + ((i+1)*1.0);
		}
	}

	// Test at different tilts and rotations
	monoSys.applySavedCoor("savedBestMonomer");
	_helicalAxis.applySavedCoor("BestMonomerAxis");

	monoSys.saveAltCoor("bestZ");
	_helicalAxis.saveAltCoor("bestZ");

	double bestTilt = 0.0;
	double bestRotation = 0.0;
	double monoTilt = 0.0;
	double monoAxialRotation = 0.0;
	for(int i=1; i<=3; i++) { // test at 3 tilts: 15, 30 and 45 degrees
		//==================================
		//====== Membrane Tilt ======
		//==================================
		monoSys.applySavedCoor("bestZ");
		_helicalAxis.applySavedCoor("bestZ");

		monoTilt = i * 15;
		_trans.rotate(chainA, monoTilt, axisA(0).getCoor(), axisB(0).getCoor());
		_trans.rotate(axisA, monoTilt, axisA(0).getCoor(), axisB(0).getCoor());
		for(int j=0; j<=3; j++) { // test at 4 rotations 0, 90, 180 and 270 degrees
			//==================================
			//====== Axial Rot ======
			//==================================
			monoAxialRotation = j * 90.0;

			monoSpm.calculateEnergies();
			monoSpm.runGreedyOptimizer(_greedyCycles);
			currentEnergy = monoSpm.getMinBound()[0];
			_fout << "current tilt: " << monoTilt << " current rotation: " << monoAxialRotation << " Energy: " << currentEnergy*2.0 << endl; // must double the energy, as only computed energy for 1 helix
			//monoSys.writePdb("mono_" + MslTools::doubleToString(monoTilt) + "_" + MslTools::doubleToString(monoAxialRotation) + ".pdb");

			if(currentEnergy < bestEnergy) {
				bestEnergy = currentEnergy;
				bestTilt = monoTilt;
				bestRotation = monoAxialRotation;
				monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
				monoSys.saveAltCoor("savedBestMonomer");
				_helicalAxis.saveAltCoor("BestMonomerAxis");
			}

			_trans.rotate(chainA, 90.0, axisA(0).getCoor(), axisA(1).getCoor());

		}
	}

	MonteCarloManager MCMngr(1000.0, 0.5, _MCCycles, MonteCarloManager::EXPONENTIAL, _MCMaxRejects);
	MCMngr.setEner(bestEnergy);

	double zShift = bestZ;
	double crossingAngle = bestTilt;
	double axialRotation = bestRotation;
	unsigned int counter = 0;

	while(!MCMngr.getComplete()) {

		monoSys.applySavedCoor("savedBestMonomer");
		_helicalAxis.applySavedCoor("BestMonomerAxis");

		int moveToPreform = _RNG.getRandomInt(2);

		double deltaZShift = 0.0;
		double deltaTilt = 0.0;
		double deltaAxialRotation = 0.0;

		//======================================
		//====== Z Shift ======
		//======================================
		if (moveToPreform == 0) {
			deltaZShift = getStandardNormal(_RNG) * 1.0;
			CartesianPoint translateA = axisA(1).getCoor() - axisA(0).getCoor(); // vector minus helical center
			translateA = translateA.getUnit() * deltaZShift; // unit vector of helical _axis times the amount to shift by
			_trans.translate(chainA, translateA);
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "Zshift: " << deltaZShift << endl;

		} else if (moveToPreform == 1) {
		//==================================
		//====== Axial Rot ======
		//==================================
			deltaAxialRotation = getStandardNormal(_RNG) * 20.0;
			_trans.rotate(chainA, deltaAxialRotation, axisA(0).getCoor(), axisA(1).getCoor());
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "axial: " << deltaAxialRotation << endl;

		} else if (moveToPreform == 2) {
		//==================================
		//====== Membrane Tilt ======
		//==================================
			deltaTilt = getStandardNormal(_RNG) * 10;
			_trans.rotate(chainA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			_trans.rotate(axisA, deltaTilt, axisA(0).getCoor(), axisB(0).getCoor());
			//_fout << setiosflags(ios::fixed) << setprecision(3)<< "tilt: " << deltaTilt << endl;
		}

		// Run Optimization
		// Run repack every N steps
		if (counter % 10 == 0) {
			//_fout << "repack." << endl;
			monoSpm.calculateEnergies();
			monoSpm.runGreedyOptimizer(_greedyCycles);

			currentEnergy = monoSpm.getMinBound()[0];
		} else {
			currentEnergy = monoSys.calcEnergy();
			//_fout << monoEset->getSummary() << endl;
		}

		if (!MCMngr.accept(currentEnergy)) {
			//_fout << "state rejected   energy: " << currentEnergy << endl;
		}
		else {
			monoSys.setActiveRotamers(monoSpm.getMinStates()[0]);
			monoSys.saveAltCoor("savedBestMonomer");
			_helicalAxis.saveAltCoor("BestMonomerAxis");
			bestEnergy = currentEnergy;

			crossingAngle = crossingAngle + deltaTilt;
			axialRotation = axialRotation + deltaAxialRotation;
			zShift = zShift +  deltaZShift;

			_fout << setiosflags(ios::fixed) << setprecision(3) << "MCAccept   axial Tilt: " << crossingAngle << " zShift: " << zShift << " axialRot: " << axialRotation << " energy: " << currentEnergy*2 << endl;
		}

		counter++;
	}

	time(&endTimeMono);
	diffTimeMono = difftime (endTimeMono, startTimeMono);
	_fout << endl << "Total Monomer Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTimeMono << " seconds" << endl;


	/******************************************************************************
	 *               === PRINT OUT MONOMER / STORE ENERGIES ===
	 ******************************************************************************/
	monoSys.applySavedCoor("savedBestMonomer");
	_helicalAxis.applySavedCoor("BestMonomerAxis");
	_fout << monoEset->getSummary();
	_fout << endl;

	// print the monomer
	string monoOutCrdFile  = _opt.pdbOutputDir + "/" + _opt.uniprotAccession + "_monomer.crd";
	CRDWriter monoCrd;
	monoCrd.open(monoOutCrdFile);
	if(!monoCrd.write(monoSys.getAtomPointers())) {
		cerr << "Unable to write " << monoOutCrdFile << endl;
		exit(0);
	}

	string monoOutPdbFile  = _opt.pdbOutputDir + "/" + _opt.uniprotAccession + "_monomer.pdb";
	PDBWriter monoPdb;
	monoPdb.setConvertFormat("CHARMM22","PDB2");
	monoPdb.open(monoOutPdbFile);
	if(!monoPdb.write(monoSys.getAtomPointers())) {
		cerr << "Unable to write " << monoOutPdbFile << endl;
		exit(0);
	}

	// Store monomer energy by term
	if(_opt.printTermEnergies) {
		monoSys.calcEnergy();
		_monomerEnergyByTerm = getEnergyByTermDoubled(monoSys.getEnergySet()); // must double the energy, as only computed energy for 1 helix

		ofstream meOut;
		string meOutName  = _opt.pdbOutputDir + "/" + _opt.uniprotAccession + "_monomer.energy";
		meOut.open(meOutName.c_str());
		if(!meOut.is_open()) {
			cerr << "Unable to open " << meOutName << endl;
			exit(0);
		}
		for(map<string,double>::iterator it = _monomerEnergyByTerm.begin(); it != _monomerEnergyByTerm.end(); it++) {
			meOut << it->first << " " << it->second << endl;
		}
		meOut.close();

	}

	double finalEnergy = 2.0 * monoSpm.getMinBound()[0]; // double the energy for 2 helices
	return finalEnergy;

}

void setupOutputDirectory(catmOptions &_opt, string &_logFile){
	string cmd = "mkdir -p " + _opt.pdbOutputDir;
	_logFile = _opt.pdbOutputDir + "/" + _opt.uniprotAccession + ".log";
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
}


/****************************************
 *
 *  ======= CONFIG FILE OPTIONS =======
 *
 ****************************************/

