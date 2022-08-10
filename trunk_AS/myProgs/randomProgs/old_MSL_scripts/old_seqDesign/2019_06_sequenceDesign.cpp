#include <iostream>
#include <fstream>
#include <sstream>

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
#include "ResidueSelection.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;
string programName = "GASrightSequenceDesign";
string programDescription = "This program designs GASright sequences onto the interfacial residues of a given GASright geometry";
string programAuthor = "Gilbert Loiseau";
string programVersion = "1";
string programDate = "22 May 2019";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

/******************************************************************************************************************************************************************************/


class HelixDimer : public AtomContainer {
	public:
	HelixDimer(string _name, double _energy, int _thread);

	void print();
	void setHelixDimerDetails(double _x, double _z, double _a, double _c, string _interface, string _prolineMask, vector<string> _hbonds, int _centroidIDNum);

	void setDeltaEnergyByTerm(map<string,double> _deltaEbyTerm) {deltaEByTerm = _deltaEbyTerm;}

	map<string,double>& getDeltaEnergyByTerm(){return deltaEByTerm;}
	void addAxisAtom(Atom* _atom){axisAtoms.push_back(_atom);}
	AtomPointerVector& getAxes(){return axisAtoms;}

	bool operator< (HelixDimer& _b) const {
		if(energy < _b.energy) {
			return true;
		} else {
			return false;
		}
	}

	void printHelixDimerDetails(ofstream & _fout);
	int getNumHbonds() {return hbonds.size();}
	int getThread() {return thread;}
	int getCentId() {return centroidIDNum;}
	
	string getName() {return name;}
	double getEnergy() {return energy;}
	double getXShift() {return xshift;}
	double getZShift() {return zshift;}
	double getAxialRotation() {return axialrot;}
	double getRelZShift() {return relZshift;}
	double getRelAxialRotation() {return relAxialrot;}
	double getCrossingAngle() {return crossang;}
	string getInterface() {return interface;}
	string getProlineMask() {return prolineMask;}

	// Format: A,11,CA:B,12,O=2.6;
	vector<string>& getHbonds() {return hbonds;}
	

	private:

	string name;
	double energy;
	double xshift;
	double zshift;
	double relZshift;
	double axialrot;
	double relAxialrot;
	double crossang;
	string interface;
	int thread;
	string prolineMask;
	int centroidIDNum;
	vector<string> hbonds;
	map<string,double> deltaEByTerm;
	AtomPointerVector axisAtoms; 
};

struct compareHelixDimers {
	bool operator () (HelixDimer* lhs, HelixDimer* rhs) {return *lhs < *rhs;}
};


HelixDimer::HelixDimer(string _name,double _energy, int _thread) : AtomContainer() {
	name = _name;
	energy = _energy;
	thread = _thread;
}

void HelixDimer::print() {
	cout << name << " " << energy << " " << hbonds.size() << " "  << thread << endl;
}

void HelixDimer::setHelixDimerDetails(double _x, double _z, double _a, double _c, string _interface, string _prolineMask, vector<string> _hbonds, int _centroidIDNum) {
	xshift = _x;
	zshift = _z;
	axialrot = _a;
	crossang = _c;
	interface = _interface;
	centroidIDNum = _centroidIDNum;
	
	prolineMask = _prolineMask;
	hbonds = _hbonds;
	relAxialrot = 10.0/9.0 * axialrot + 200.0/27.0 * zshift;
	relZshift = 10.0/9.0 * zshift + axialrot/60.0;

	// set the 35 - thread - {0,1,4,5} as 2 - to mark the trapezoid
	int trapCorner = 35-thread;
	if(relAxialrot > 0) {
		trapCorner -= 1;
	} else if (relAxialrot < -100) {
		trapCorner += 1;
	}
	if(relZshift > 6) {
		trapCorner -= 4;
	} else if (relZshift < 0) {
		trapCorner += 4;
	}
	//cout << "TRAPCORNER " << trapCorner << endl;
	//cout << "interface " << interface << endl;
	if(trapCorner >= 0) {
		interface.replace(trapCorner,1,"2");
	}
	if(trapCorner >= 1) {
		interface.replace(trapCorner-1,1,"2");
	}
	if(trapCorner >= 4) {
		interface.replace(trapCorner-4,1,"2");
	}
	if(trapCorner >= 5) {
		interface.replace(trapCorner-5,1,"2");
	}

}

void HelixDimer::printHelixDimerDetails(ofstream & _fout) {
	_fout << "interfacial " << interface << endl;
	_fout << "prolineMask " << prolineMask << endl;
	_fout << "name " << name << endl;
	_fout << "thread " << thread << endl;
	_fout << "Xshift " << xshift << endl;
	_fout << "Zshift " << zshift << endl;
	_fout << "relZshift " << relZshift << endl;
	_fout << "axialRot " << axialrot << endl;
	_fout << "relAxialRot " << relAxialrot << endl;
	_fout << "crossAng " << crossang << endl;
	_fout << "energy " << energy << endl;
	string hbondsList = "";
	for(int h = 0; h < hbonds.size(); h++) {
		if(h == 0) {
			hbondsList += hbonds[h];
		} else {
			hbondsList += ";" + hbonds[h];
		}
	}
	_fout << "hBonds " << hbondsList << endl;
}


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


HelixDimerCluster::HelixDimerCluster(HelixDimer * _structure) {
	members.push_back(_structure);
}

void HelixDimerCluster::addHelixDimer(HelixDimer * _structure) {
	if(members.size() < 1) {
		members.push_back(_structure);
		return;
	}
	if(_structure->getEnergy() < 0) {
		members.push_back(_structure);
	}
}
AtomPointerVector& HelixDimerCluster::getAtomPointers() {
	if(members.size() > 0) {
		return members[0]->getAtomPointers();
	} else {
		cerr << "ERROR 13543: No members in cluster" << endl;
		exit(0);
	}
}
void HelixDimerCluster::setDetails(string _origSeq, string _seq, string _uName, string _uAccession, string _outputDir, int _resStart, int _resEnd, int _tmStart, int _tmEnd) {
	uniprotName = _uName;
	uniprotAccession = _uAccession;
	seq = _seq;
	origSeq = _origSeq;
	outputDir = _outputDir;
	resStart = _resStart;
	resEnd = _resEnd;
	tmStart = _tmStart;
	tmEnd = _tmEnd;
}

void HelixDimerCluster::printTermEnergies(int id) {
	char name[1000];
	sprintf(name,"%s/%s_%02d.energy",outputDir.c_str(),uniprotAccession.c_str(),id);
	ofstream eOut;
	eOut.open(name);
	map<string,double>& deltaEByTerm = members[0]->getDeltaEnergyByTerm();
	for(map<string,double>::iterator it = deltaEByTerm.begin(); it != deltaEByTerm.end(); it++) {
		eOut << it->first << " " << it->second << endl;
	}
	eOut.close();
}

void HelixDimerCluster::printHelixDimerClusterCrds(int id, bool allStructures,bool writeAxis) {
	//members[0]->print();
	CRDWriter writer;
	char name[1000];
	sprintf(name,"%s/%s_%02d.crd",outputDir.c_str(),uniprotAccession.c_str(),id);
	writer.open(string(name));
	writer.addRemark("Energy " + MslTools::doubleToString(members[0]->getEnergy()) );
	writer.addRemark("Name - " + members[0]->getName() + " HydrogenBonds " + MslTools::intToString(members[0]->getNumHbonds()) );
	AtomPointerVector& ats = members[0]->getAtomPointers();
	writer.write(ats);
	writer.close();

	if(writeAxis) {
		sprintf(name,"%s/%s_%02d_axis.crd",outputDir.c_str(),uniprotAccession.c_str(),id);
		writer.open(string(name));
		writer.write(members[0]->getAxes());
		writer.close();
	}

	if(allStructures) {
		for(int i = 0; i < members.size(); i++) {
			sprintf(name,"%s/%s_%02d_%03d.crd",outputDir.c_str(),uniprotAccession.c_str(),id,i);
			writer.open(string(name));
			writer.clearRemarks();
			writer.addRemark("Energy " + MslTools::doubleToString(members[i]->getEnergy()) );
			writer.addRemark("Name - " + members[i]->getName() + " HydrogenBonds " + MslTools::intToString(members[i]->getNumHbonds()) );
			AtomPointerVector& ats = members[i]->getAtomPointers();
			writer.write(ats);
			writer.close();
		}
	}
}

void HelixDimerCluster::convertToPdbNames() {
	FormatConverter fc;
	fc.setNamespaces("CHARMM22","PDB2");

	for(int i  = 0; i < members.size(); i++) {
		AtomPointerVector& ats = members[i]->getAtomPointers();
		fc.convert(ats);
	}
}

void HelixDimerCluster::printHelixDimerClusterPdbs(int id) {
	//members[0]->print();
	PDBWriter writer;
	writer.setConvertFormat("CHARMM22","PDB2");
	char name[1000];
	sprintf(name,"%s/%s_%02d.pdb",outputDir.c_str(),uniprotAccession.c_str(),id);
	writer.open(string(name));
	char helixA[1000];
	char helixB[1000];

	
	// Format is here http://www.wwpdb.org/documentation/format33/sect5.html
	//COLUMNS        DATA  TYPE     FIELD         DEFINITION
	//-----------------------------------------------------------------------------------
	// 1 -  6        Record name    "HELIX "
	// 8 - 10        Integer        serNum        Serial number of the helix. This starts
	//                                            at 1  and increases incrementally.
	//12 - 14        LString(3)     helixID       Helix  identifier. In addition to a serial
	//                                            number, each helix is given an 
	//                                            alphanumeric character helix identifier.
	//16 - 18        Residue name   initResName   Name of the initial residue.
	//20             Character      initChainID   Chain identifier for the chain containing
	//                                            this  helix.
	//22 - 25        Integer        initSeqNum    Sequence number of the initial residue.
	//26             AChar          initICode     Insertion code of the initial residue.
	//28 - 30        Residue  name  endResName    Name of the terminal residue of the helix.
	//32             Character      endChainID    Chain identifier for the chain containing
	//                                            this  helix.
	//34 - 37        Integer        endSeqNum     Sequence number of the terminal residue.
	//38             AChar          endICode      Insertion code of the terminal residue.
	//39 - 40        Integer        helixClass    Helix class (see below).
	//41 - 70        String         comment       Comment about this helix.
	//72 - 76        Integer        length        Length of this helix.
	//HELIX    1  HA ALA A   19  ILE A   41  1                                  23
	//HELIX    2  HA ALA B   19  ILE B   41  1                                  23


	
	FormatConverter fc;
	fc.setNamespaces("CHARMM22","PDB2");

	for(int i = 0; i < members.size(); i++) {
		AtomPointerVector& ats = members[i]->getAtomPointers();

		if(i == 0) {
			// BEWARE: hacky code - will work only for Homo
			Atom& startAtom = *ats[0];
			Atom& endAtom = *ats[ats.size() -1];

			int helixSize = resEnd - resStart + 1;

			sprintf(helixA,"HELIX    1  HA %3s A %4d  %3s A %4d  1                              %5d\n",fc.getResidueName(startAtom.getResidueName()).c_str(),startAtom.getResidueNumber(),fc.getResidueName(endAtom.getResidueName()).c_str(),endAtom.getResidueNumber(),helixSize);
			sprintf(helixB,"HELIX    2  HB %3s B %4d  %3s B %4d  1                              %5d\n",fc.getResidueName(startAtom.getResidueName()).c_str(),startAtom.getResidueNumber(),fc.getResidueName(endAtom.getResidueName()).c_str(),endAtom.getResidueNumber(),helixSize);

			string line = string(helixA) + string(helixB);
			writer.writeln(line); 

		}
		writer.clearRemarks();
		writer.addRemark("Energy " + MslTools::doubleToString(members[i]->getEnergy()) );
		writer.addRemark("Name - " + members[i]->getName() + " HydrogenBonds " + MslTools::intToString(members[i]->getNumHbonds()) );
		writer.addRemark("MSL " + mslVersion + ", CATM " + programVersion);

		writer.writeREMARKS();
		writer.write(ats,true,false,true);
	}
	writer.close();
}

void HelixDimerCluster::printDetails(int id) {
	ofstream fout;
	char filename[1000];
	sprintf(filename,"%s/%s_%02d.txt",outputDir.c_str(),uniprotAccession.c_str(),id);
	fout.open(filename);
	fout << "originalSeq " << origSeq << endl;
	fout << "modelledSeq " << seq << endl;
	fout << "interfacial " << members[0]->getInterface() << endl;
	fout << "prolineMask " << members[0]->getProlineMask() << endl;
	fout << "resStart " << resStart << endl;
	fout << "resEnd " << resEnd << endl;
	fout << "tmStart " << tmStart << endl;
	fout << "tmEnd " << tmEnd << endl;
	fout << "numModels " << members.size() << endl;
	fout << "thread " << members[0]->getThread() << endl;
	fout << "XShift " << members[0]->getXShift() << endl;
	fout << "ZShift " << members[0]->getZShift() << endl;
	fout << "relZshift " << members[0]->getRelZShift() << endl;
	fout << "axialRot " << members[0]->getAxialRotation() << endl;
	fout << "relAxialRot " << members[0]->getRelAxialRotation() << endl;
	fout << "crossAng " << members[0]->getCrossingAngle() << endl;
	fout << "energy " << members[0]->getEnergy() << endl;
	fout << "hbonds ";
	vector<string> & hbonds = members[0]->getHbonds();
	for(int i = 0; i < hbonds.size(); i++) {
		fout << hbonds[i] << ";";
	}
	fout << endl;
	// print details of cluster members
	fout << "MODEL num thread centID   XShift   ZShift relZShift    axRot relAxRot crossAng numhb     energy " << endl;
	for(int i = 0; i < members.size(); i++) {
		char line[1000];
		sprintf(line,"MODEL %3d %6d %6d %8.3f %8.3f %9.3f %8.3f %8.3f %8.3f %5d %f",i,members[i]->getThread(), members[i]->getCentId(), members[i]->getXShift(), members[i]->getZShift(), members[i]->getRelZShift(), members[i]->getAxialRotation(),members[i]->getRelAxialRotation(),members[i]->getCrossingAngle(),members[i]->getNumHbonds(),members[i]->getEnergy());
		//fout << "MODEL " << i << " " << members[i]->getThread() << " " << members[i]->getXShift() << " " << members[i]->getZShift() << " " << members[i]->getAxialRotation() << " " << members[i]->getCrossingAngle() << " " << members[i]->getEnergy() << " " << members[i]->getNumHbonds() << endl; 
		fout << line << endl;
	}
	fout.close();
}

void HelixDimerCluster::makePse(int id) {

	char filename[1000];
	///sprintf(filename,"%s/%s_%02d",outputDir.c_str(),uniprotAccession.c_str(),id);
	sprintf(filename,"%s_%02d",uniprotAccession.c_str(),id);


	ofstream fout;

	char scriptfilename[1000];
	sprintf(scriptfilename,"%s/%s.inp",outputDir.c_str(),filename);
	
	char pdbfilename[1000];
	sprintf(pdbfilename,"%s.pdb",filename);

	char psefilename[1000];
	sprintf(psefilename,"%s.pse",filename);

	// PyMOL script
	/*
	bg_color white
	Load MODEL_000.pdb
	Rotate X, 90
	Color green, chain B and element C
	Color cyan, chain B and element C
	Show cartoon
	Show stick, resi 28+31+32+35+36+39+40
	Distance hb000, ///A/35/HA1, ///B/32/O
	Color black, hb000
	Distance hb001, ///B/35/HA1, ///A/32/O
	Color black, hb001
	Distance hb002, ///A/36/HA, ///B/35/O
	Color black, hb002
	Distance hb003, ///B/36/HA, ///A/35/O
	Color black, hb003
	Distance hb004, ///A/39/HG1, ///B/36/O
	Color black, hb004
	Distance hb005, ///B/39/HG1, ///A/36/O
	Color black, hb005
	Zoom all, 3
	Set label_color, black
	Save MODEL_000.pse
	Hide labels
	Ray 400,600
	Png MODEL_000_view1.png
	Rotate Y, 90
	Ray 500,600
	Png MODEL_000_view2.png
*/

	fout.open(scriptfilename);
	fout << "bg_color white" << endl;
	fout << "load " << pdbfilename << endl;
	fout << "rotate X, 90, all, 0" << endl;
	fout << "color green, chain B and element C" << endl;
	fout << "color cyan, chain B and element C" << endl;
	fout << "show cartoon" << endl;

	// find the interface residue numbers from the seq, resStart and interface mask
	string interface = members[0]->getInterface();

	vector<int> interfaceRes;
	for(int i = 0; i < interface.length();i++) {
		if(interface[i] == '1' || interface[i] == '2') {
			interfaceRes.push_back(resStart + i);
		}
	}
	//fout << "show stick, resi " << interfaceRes[0];

	// form the strings
	/*
	select interfaceA, chain A and resi 6+9+10+13+14+16+18+20+21
	select interfaceB, chain B and resi 6+9+10+13+14+16+18+20+21
	show stick, interfaceA interfaceB
	*/

	stringstream ss;
	ss <<  interfaceRes[0];
	for(int i = 1; i < interfaceRes.size(); i++) {
		ss << "+" << interfaceRes[i];
	}

	fout << "select interfaceA, chain A and resi " << ss.str() << endl;
	fout << "select interfaceB, chain B and resi " << ss.str() << endl;
	fout << "show stick, interfaceA interfaceB" << endl;

	// A,12,ALA,HA1:B,34,LEU,O=2.6;A,12,ALA,HA1:B,34,LEU,O=3.5
	vector<string> hbondTokens  = members[0]->getHbonds();
	
	FormatConverter fc;
	fc.setNamespaces("CHARMM22","PDB2");
	for(int i = 0; i < hbondTokens.size(); i++) {
		vector<string> bondData = MslTools::tokenize(hbondTokens[i],":=");
		string chain1,resNum1,resName1,atom1;
		vector<string> donorData = MslTools::tokenize(bondData[0],",");
		chain1 = donorData[0];
		resNum1 = donorData[1];
		resName1 = donorData[2];
		if(resName1 == "HSD" || resName1 == "HSE" || resName1 == "HSP") {
			resName1 = "HIS";
		}
		atom1 = donorData[3];
		atom1 = fc.getAtomName(atom1,resName1);
		string chain2,resNum2,resName2,atom2;
		vector<string> acceptorData = MslTools::tokenize(bondData[1],",");
		chain2 = acceptorData[0];
		resNum2 = acceptorData[1];
		resName2 = acceptorData[2];
		if(resName2 == "HSD" || resName2 == "HSE" || resName2 == "HSP") {
			resName2 = "HIS";
		}
		atom2 = acceptorData[3];
		atom2 = fc.getAtomName(atom2,resName2);
		char distLine[500];
		sprintf(distLine,"distance hb%03d, ///%s/%s/%s, ///%s/%s/%s\ncolor black, hb%03d",i,chain1.c_str(),resNum1.c_str(),atom1.c_str(),chain2.c_str(),resNum2.c_str(),atom2.c_str(),i);
		fout << distLine << endl;
	}

	fout << "zoom all, 3" << endl;
	fout << "set label_color, black" << endl;
	fout << "save " << psefilename << endl;
	fout << "hide labels" << endl;
	fout << "set ray_opaque_background, false" << endl;
	fout << "ray 400,600" << endl;
	fout << "png " <<  filename << "_view1.png" << endl;
	fout << "rotate Y, 90, all, 0" << endl;
	fout << "ray 500,600" << endl;
	fout << "png " << filename << "_view2.png" << endl;
//TODO : once pymol is fixed uncomment this part - update the correct path before invoking the .inp
/*
	char command [1000];
	sprintf(command,"pymol -cqd @%s",scriptfilename);
	system(command);
	*/
}

/******************************************************************************************************************************************************************************/



struct Options{
	string sequence;
	string backboneAA;
	int backboneLength;

	// optional
	int tmStart;
	int tmEnd;

	// the actual AAs being modeled
	// TODO: think I have to change this to a vector for specifically the interfacial residues 
	int startResNum;
	int endResNum;
	int sequenceStart;

	int threadStart;
	int threadEnd;
	bool threadBool;
	
	bool deleteTerminalHbonds;
	
	vector<string> AAs;
	string SL; //number of rotamers

	// transformation
	double xShift;
	double zShift;
	double crossingAngle;
	double axialRotation;
	bool transform;
	int thread;
	int bbThread;
	
	// datafile options
	string originalSeq;
	string modelledSeq;
	int resStart;
	int resEnd;
	int seqLength;
	double XShift;
	double ZShift;
	double axialRot;
	double crossAng;
	double energy;

	// input files
	string helixGeoFile;
	string backboneCrd;	
	string pdbOutputDir;
	string topFile;
	string parFile;
	string solvFile;
	string hBondFile;
	string rotLibFile;
	string monoRotLibFile;
	string infile;
	string rulesFile;

	// side-chain repack variable
	int mcCycles;
	int mcMaxRejects;
	double mcStartTemp;
	double mcEndTemp;
	int mcCurve;

	double deltaZ;
	double deltaAx;
	double deltaCross;
	double deltaX;

	bool verbose;
	int greedyCycles;
	int seed;

	int numberOfStructuresToMCRepack;
	double energyCutOff;
	
	// weights
	double weight_vdw;
	double weight_hbond;
	double weight_solv;
	
	// input monomerEnergy
	bool inputMonomerE;
	int monoE_vdw;
	int monoE_hbond;
	int monoE_solv;
	int monoE_solvRef;

	// clustering options (optional)
	double rmsdCutoff;
	bool clusterSolutions;
	bool printAllCrds;
	bool printAxes;
	bool printTermEnergies;

	int start;
	int end;
	double ener;
	vector<int> ivalues;

	// alternate identities
	vector<string> Ids;

	// protein information (optional)
	string uniprotName;
	string uniprotAccession;

	/***** MANAGEMENT VARIABLES ******/
	string pwd; // the present working directory obtained with a getenv
	string host; // the host name obtained with a getenv
	bool version; // ask for the program version
	bool help; // ask for the program help

	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages

	vector<string> allowed; //list of allowed options
	vector<string> required; //list of required options

	vector<string> disallowed;  // disallowed options that were given
	vector<string> missing; // required options that were not given
	vector<string> ambiguous; // required options that were not given
	vector<string> defaultArgs; // the default arguments can be specified in command line without "--option"
	vector<vector<string> > equivalent; // this links short options to long ones (for example -x can be given for --extended)

	string OPerrors; //the errors from the option parser
	string rerunConf; // data for a configuration file that would rerun the job as the current run

	string configfile;
	string datafile;
};

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}

void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
	
	/* Faster code 
	for (uint i=0; i < _apvA.size(); i++) {
			_apvB[i]->copyAllCoor(*_apvA[i]);
			vector<CartesianPoint*>& bCoors = _apvB[i]->getAllCoor();

			for(uint j = 0; j < bCoors.size(); j++) {
				bCoors[j]->setX(0 - bCoors[j]->getX());
				bCoors[j]->setY(0 - bCoors[j]->getY());
			}
					
		}
	*/

	// Set coordinates of chain A to chain B
	for (uint i=0; i < _apvA.size(); i++) {
		_apvB[i]->copyAllCoor(*_apvA[i]);
	}

	// Rotation matrix for 180 degrees
	// flips the sign on the x and y coordinates
	Matrix m(3,3,0.0);
	m[0][0] = -1.0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = -1.0;
	m[1][2] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = 1.0;

	// Rotate chain B around Z axis
	Transforms trans; 
	trans.rotate(_apvB, m);
	
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

string convertToPolymerSequence(string _seq, int _startResNum) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		if(resName == "HIS") {
			ps = ps + " HSE";
		} else {
			ps = ps + " " + resName;
		}
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

string convertToPolymerSequence(string _seq, int _startResNum, vector<string> _alternateIds, vector<int> _varPos) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	int counter = 0;
	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
		counter++;
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		if (_varPos[counter] == 1){
			ps = ps + " [";
			if(resName == "HIS") {
				ps = ps + " HSE";
			} else {
				ps = ps + " " + resName;
			}
			for (uint i=0; i<_alternateIds.size(); i++){
				if(_alternateIds[i] == "HIS") {
					ps = ps + " HSE";
				} else {
					ps = ps + " " + _alternateIds[i];
				}
			}
			ps = ps + "] ";
		} else {
			if(resName == "HIS") {
				ps = ps + " HSE";
			} else {
				ps = ps + " " + resName;
			}
		}
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans) {
	
	//====== Z Shift (Crossing Point) ======
	CartesianPoint zShiftCP(0.0, 0.0, _zShift);
	_trans.translate(_chainA, zShiftCP);

	//===== Axial Rotation ======
	_trans.rotate(_chainA, _axialRotation, _ori, _zAxis);

	//====== Local Crossing Angle ======
	_trans.rotate(_chainA, (_crossingAngle/2.0), _ori, _xAxis);
	_trans.rotate(_axisA, (_crossingAngle/2.0), _ori, _xAxis);

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((-1.0*_xShift/2.0), 0.0, 0.0);
	_trans.translate(_chainA, interDistVect);
	_trans.translate(_axisA, interDistVect);

	c2Symmetry(_chainA, _chainB);
	c2Symmetry(_axisA, _axisB);
}

void repackSideChains(SelfPairManager & _spm, int _greedyCycles, vector<vector<vector<vector<bool> > > > _savedEnergyFlagTable) {

	_spm.setOnTheFly(1);
	//_spm.recalculateNonSavedEnergies(_savedEnergyFlagTable);
	_spm.calculateEnergies(); // CHANGE BACK!!!
	_spm.runGreedyOptimizer(_greedyCycles);
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

void clusterSolutions(vector<HelixDimerCluster*>& clusters,vector<HelixDimer*>& _structures, double _rmsdCutoff, string _origSeq, string _builtSeq, Options& opt) {
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
			clusters.back()->setDetails(_origSeq, _builtSeq, opt.uniprotName,opt.uniprotAccession,opt.pdbOutputDir,opt.startResNum,opt.endResNum,opt.tmStart,opt.tmEnd);
		}
	}

}

struct HbondInfo{
	string chain;
	string resA;
	string resB;
	string atomName;
	double distBreak;

	HbondInfo(string _c, string _resA, string _resB, string _aName, string _dist) {
		chain = _c;
		resA = _resA;
		resB = _resB;
		atomName = _aName;
		distBreak = MslTools::toDouble(_dist);
	}

	bool isValid(System & _sys) {
		string posA = chain + "," + resA;
		string posB = "";
		if(chain == "A") {
			posB = "B," + resB;
		} else {
			posB = "A," + resB;
		}
		if(!_sys.positionExists(posA) || !_sys.positionExists(posB) ) {
			return false;
		}

		if(atomName == "HA1") {
			if(!_sys.identityExists(posA + ",GLY" )) {
				return false;
			}
		}
		return true;
	}

	void print(ofstream & pout) {
		pout << chain << " " << resA << " " << resB << " " << atomName << " " << distBreak << endl;
	}
};

bool hydrogenBondCheck(System & _sys, vector<string> _parsedGeoInformation, double & _xShiftStart) {
	
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
	if(hbondList.size() < 4) {
		for(int i = 0; i < hbondList.size(); i++) {
			delete hbondList[i];
		}
		return false;
	} else {
		_xShiftStart = hbondList[hbondList.size() - 4]->distBreak - 0.1;
		for(int i = 0; i < hbondList.size(); i++) {
			delete hbondList[i];
		}
		return true;
	}

}//*/

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

void deleteTerminalHydrogenBondInteractions(System &_sys, Options& _opt) {
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
	pESet->deleteInteractionsWithAtoms(atoms,"SCWRL4_HBOND");

}

void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV, AtomPointerVector& _axis, Transforms & _trans) {

	AtomSelection sel(_apV);
	AtomPointerVector & caApV = sel.select("name CA");
	double zShift = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		zShift += (caApV[i]->getCoor()).getZ();
	}
	zShift = -1.0 * zShift/double(caApV.size());
	//fout << x << " " << y << " " << pt << " " << caApV.size() << endl;

	// old code
	//for(int i = 0; i < _apV.size(); i++) {
	//	CartesianPoint& pt = _apV[i]->getCoor();
	//	pt.setZ(pt.getZ() +  zShift);
	//}

	CartesianPoint interDistVect;
	interDistVect.setCoor(0.0, 0.0, zShift);
	_trans.translate(_apV, interDistVect);
	_trans.translate(_axis, interDistVect);
	

}

map<string,double> getEnergyByTerm(EnergySet* _eSet) {
	// get all terms
	map<string,double> eByTerm;
	map<string,vector<Interaction*> > * allTerms = _eSet->getEnergyTerms();
	for(map<string,vector<Interaction*> >::iterator it = allTerms->begin(); it != allTerms->end(); it++) {
		if(_eSet->isTermActive(it->first)) {
			eByTerm[it->first] =  _eSet->getTermEnergy(it->first);
		}
	}
	return eByTerm;
}

map<string,double> getEnergyByTermDoubled(EnergySet* _eSet) {
	// get all terms
	map<string,double> eByTerm;
	map<string,vector<Interaction*> > * allTerms = _eSet->getEnergyTerms();
	for(map<string,vector<Interaction*> >::iterator it = allTerms->begin(); it != allTerms->end(); it++) {
		if(_eSet->isTermActive(it->first)) {
			eByTerm[it->first] =  2.0* _eSet->getTermEnergy(it->first);
		}
	}
	return eByTerm;
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

void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType) {

	 if (moveType == 0) {
		// Z Shift
		CartesianPoint translateA = _axisA(1).getCoor() - _axisA(0).getCoor(); // vector minus helical center 
		translateA = translateA.getUnit() * _deltaMove; // unit vector of helical _axis times the amount to shift by

		_trans.translate(_chainA, translateA);

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else if (moveType == 1) {
		// Axial Rotation
		_trans.rotate(_chainA, (_deltaMove), _axisA(0).getCoor(), _axisA(1).getCoor());

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else 	if (moveType == 2) {
		// Crossing Angle 
		_trans.rotate(_chainA, (_deltaMove * 0.5), _axisA(0).getCoor(), _axisB(0).getCoor());
		_trans.rotate(_axisA, (_deltaMove * 0.5), _axisA(0).getCoor(), _axisB(0).getCoor());

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else if (moveType == 3) {
		// XShift
		// Helix A interhelical distance
		CartesianPoint translateA = _axisB(0).getCoor() - _axisA(0).getCoor(); // vector minus helical center 
		translateA = translateA.getUnit() * _deltaMove * -0.5; // unit vector of helical axis times the amount to shift by

		_trans.translate(_chainA, translateA);
		_trans.translate(_axisA, translateA);

		// Helix B interhelical distance
		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else {
		cerr << "Unknown moveType " << moveType << " in backboneMovement. Should be 0-3 " << endl;
	}
}

double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
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

double computeMonomerEnergy(System & _sys, Transforms & _trans, Options& _opt, System & _helicalAxis, RandomNumberGenerator & _RNG, map<string,double> & _monomerEnergyByTerm, ofstream & _fout, int _greedyCycles, int _MCCycles, int _MCMaxRejects) {

	time_t startTimeMono, endTimeMono;
	double diffTimeMono;
	time(&startTimeMono);	
	
	Chain & inputChain = _sys.getChain(0);
	AtomPointerVector &axisA = _helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = _helicalAxis.getChain("B").getAtomPointers();

	// Declare new system
	System monoSys;
	CharmmSystemBuilder CSBMono(monoSys, _opt.topFile, _opt.parFile, _opt.solvFile);
	CSBMono.setBuildTerm("CHARMM_ELEC", false);
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
	monoEset->setTermActive("CHARMM_ELEC", false);
	monoEset->setTermActive("CHARMM_ANGL", false);
	monoEset->setTermActive("CHARMM_BOND", false);
	monoEset->setTermActive("CHARMM_DIHE", false);
	monoEset->setTermActive("CHARMM_IMPR", false);
	monoEset->setTermActive("CHARMM_U-BR", false);

	deleteTerminalHydrogenBondInteractions(monoSys, _opt);
	
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


	monoEset->setWeight("CHARMM_VDW", _opt.weight_vdw);
	monoEset->setWeight("SCWRL4_HBOND", _opt.weight_hbond);
	monoEset->setWeight("CHARMM_IMM1REF", _opt.weight_solv);
	monoEset->setWeight("CHARMM_IMM1", _opt.weight_solv);
	_fout << "Monomer - VDW weight: " << monoEset->getWeight("CHARMM_VDW") << " HB weight: " << monoEset->getWeight("SCWRL4_HBOND") << " IMM1REF weight: " << monoEset->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << monoEset->getWeight("CHARMM_IMM1") << endl;

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


Options parseOptions(int _argc, char * _argv[], Options defaults);

//skeleton code for adding given sequence into different portions of a polyLeu chain
void generateSequence(int _startResNum, int _sequenceLength, int _backboneLength, string _backboneAA, string &_sequence, bool _backbone){
	string seq = "";
	for(uint i=0; i<_sequenceLength; i++){
		//somehow multiply the amount of _backboneAA, then add in the _sequence, then add in leftover _backboneAA
		//might be able to switch this to a while loop like below?
		seq += _backboneAA;	
	}
	//seq += _sequence; use at some point if I want to do something like threading through the sequence
	if (_backbone == true){
		while (seq.length() < _backboneLength){
			seq += _backboneAA;
			//cout << seq << endl;
		}
	}
	_sequence = seq;
}

void printOptions(Options & _op, ofstream & _fout) {
	_fout << "Datafile: " << _op.datafile << endl << endl;

	_fout << "Options from Datafile" << endl;
	_fout << "Original Sequence: " << _op.originalSeq << endl;
	_fout << "Modelled Sequence: " << _op.modelledSeq << endl;
	_fout << "Residue Start: " << _op.resStart << " Residue End: " << _op.resEnd << endl;
	_fout << "Sequence Length: " << _op.seqLength << endl;
	_fout << setiosflags(ios::fixed) << setprecision(3) << "XShift: " << _op.XShift << " ZShift: " << _op.ZShift << " Axial Rotation: " << _op.axialRot << " Crossing Angle: " << _op.crossAng << " Energy: " << _op.energy << endl;
	
	_fout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;

	_fout << "Warning Messages: " << _op.warningMessages << endl << endl;

	_fout << "Other Parameters" << endl;
	_fout << "backboneCrd " << _op.backboneCrd << endl;
	//_fout << "logFile " << _op.logFile << endl;
	_fout << "pdbOutputDir " << _op.pdbOutputDir << endl;

	//_fout << "fullSequence " << _op.fullSequence << endl;
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

	_fout << "MCCycles " << _op.mcCycles << endl;
	_fout << "MCMaxRejects " << _op.mcMaxRejects << endl;
	_fout << "MCStartTemp " << _op.mcStartTemp << endl;
	_fout << "MCEndTemp " << _op.mcEndTemp << endl;
	_fout << "MCCurve " << _op.mcCurve << endl;

	//_fout << "deltaZ " << _op.deltaZ << endl;
	//_fout << "deltaAx " << _op.deltaAx << endl;
	//_fout << "deltaCross " << _op.deltaCross << endl;
	//_fout << "deltaX " << _op.deltaX << endl;

	_fout << "verbose " << _op.verbose << endl;
	_fout << "greedyCycles " << _op.greedyCycles << endl;
	_fout << "seed " << _op.seed << endl;

	_fout << "numberOfStructuresToMCRepack " << _op.numberOfStructuresToMCRepack << endl;
	_fout << "energyCutOff " << _op.energyCutOff << endl;

	//_fout << "uniprotName " << _op.uniprotName << endl;
	//_fout << "uniprotAccession " << _op.uniprotAccession << endl;

	_fout << "monoE_vdw " << _op.monoE_vdw << endl;
	_fout << "monoE_solv " << _op.monoE_solv << endl;
	_fout << "monoE_solvRef" << _op.monoE_solvRef << endl;
	_fout << "monoE_hbond" << _op.monoE_hbond << endl;

	_fout << "rmsdCutoff " << _op.rmsdCutoff << endl;
	_fout << "clusterSolutions " << _op.clusterSolutions << endl;
	_fout << "printAllCrds " << _op.printAllCrds << endl;
	_fout << "printAxes " << _op.printAxes << endl;
	_fout << "printTermEnergies " << _op.printTermEnergies << endl;
	_fout << "deleteTerminalHbonds " << _op.deleteTerminalHbonds << endl;

	_fout << "fullSequenceStart " << _op.sequenceStart << endl;

	_fout << "startResNum " << _op.startResNum << endl;
	_fout << "endResNum " << _op.endResNum << endl;

	_fout << "threadStart " << _op.threadStart << endl;
	_fout << "threadEnd " << _op.threadEnd << endl;

	_fout << "weight_vdw " << _op.weight_vdw << endl;
	_fout << "weight_hbond " << _op.weight_hbond << endl;
	_fout << "weight_solv " << _op.weight_solv << endl;

	if(_op.configfile != "") {
		_fout << "configfile " << _op.configfile << endl;
	}

	_fout << endl;

}

void loadRotamers(System &_sys, SystemRotamerLoader &_sysRot, string _SL){
	/******************************************************************************
	 *                     === LOAD ROTAMERS ===
	 ******************************************************************************/
	for (uint k=0; k<_sys.positionSize(); k++) {
		Position &pos = _sys.getPosition(k);
		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!_sysRot.loadRotamers(&pos, pos.getResidueName(), _SL)) {
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}	
}

void identifyInterface(System &_sys, vector<int> &_pos, vector<double> &_dists){
	int count = 0;
	for (uint k=0; k<_sys.positionSize()/2; k++) {
		//cout << _sys.positionSize() << endl;
		Position &pos1 = _sys.getPosition(k);
		Position &pos2 = _sys.getPosition(k+_sys.positionSize()/2);
		
		Atom &c1 = pos1.getAtom("CA");
		Atom &c2 = pos2.getAtom("CA");
		
		double dist;
		double most;
		double x = 0;
		dist = c1.distance(c2);
		//cout << "Dist " << k << ": " << dist << endl;
		if (_pos.size() < 12){
			_dists.push_back(dist);
			_pos.push_back(k);
		}
		else {
			count++;
			if (count <= 5){
				most = _dists[0];
				x = 0;
				//cout << "Current: " << endl;
				for (uint j=0; j<8; j++){
					//cout << _pos[j] << ": " << _dists[j] << endl;
					if (most < _dists[j]){
						most = _dists[j];
						x = j;
					}
				}
				if (dist <= most){
					//cout << "Before: " << x << ": " << _dists[x];
					_dists[x] = dist;
					_pos[x] = k;
					//cout << " ; After: " << _pos[x] << ": " << _dists[x] << endl;
					count = 0;
				}
			}
			else{
				cout << "Interface Identified!" << endl;
				k = _sys.positionSize()/2;
			}
		}
	}
}

void interface01(System &_sys, vector<int> &_pos, vector<int> &_varPos){
	for (uint k=0; k<_sys.positionSize(); k++){
		_varPos.push_back(0);
	}
	for (uint j=0; j<_pos.size(); j++){
		_varPos[_pos[j]] = 1;
		_varPos[_pos[j]+_sys.positionSize()/2] = 1;
	}
}

void addIdentities(CharmmSystemBuilder &_CSB, System &_sys, vector<string> &_ids, vector<int> &_varPos){
	for (uint k=0; k<_varPos.size(); k++){
		if (_varPos[k] == 1){
			Position &pos = _sys.getPosition(k);
			for (uint j=0; j<_ids.size(); j++){
				_CSB.addIdentity(pos, _ids[j]);
			}
		}
		else{
			continue;
		}
	}
	cout << "Ids added!" << endl;
}

void switchAA(CharmmSystemBuilder &_CSB, System &_sys, vector<int> &_varPos){
	for (uint k=0; k<_varPos.size(); k++){
		if (_varPos[k] == 1){
			Position &pos = _sys.getPosition(k);
			pos.setActiveIdentity(2);
			cout << "Id switched at position " << k << "!" << endl;
			//for (uint j=0; j<_ids.size(); j++){
			//	CSB.addIdentity(pos, _ids[j]);
			//}
		}
		else{
			continue;
		}
	}
}


//TODO: how can I set it so that the optimizer results in the same switches for each position the same rather than different?
//I can think of adding a boolean to the actual code, but is this the best option? Maybe there's already one somewhere in the code
//I think I found it in System.h: void setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions); So I need to transition the positions to this "A,19" "B,19" format!

vector<vector<string>> positionToString(System &_sys, vector<int> &_varPos){
	vector<vector<string>> stringPositions;
	for (uint k=0; k<_varPos.size()/2; k++){
		//cout << "string" << endl;
		if (_varPos[k] == 1){
			vector<string> tempPos;

			Position &posA = _sys.getPosition(k);
			Position &posB = _sys.getPosition(k+_varPos.size()/2);

			string A = posA.toString();
			string B = posB.toString();

			string delimiter = " ";
			
			size_t p = 0;
			p = A.find(delimiter);

			tempPos.push_back(A.substr(0, p));
			tempPos.push_back(B.substr(0, p));

			stringPositions.push_back(tempPos);
		}
		else{
			continue;
		}
	}
	return stringPositions;
}

void identifyInterface2(System &_sys, vector<int> &_pos, vector<double> &dists){
	/*int count = 0;
	for (uint k=0; k<39; k++){
		Position &pos1 = _sys.getChain("A").getPosition(k);
		Position &pos2 = _sys.getChain("B").getPosition(k);
		
		Atom &H1 = pos1.getAtom("HA");
		Atom &H2 = pos2.getAtom("HA");
		Atom &O1 = pos1.getAtom("O");
		Atom &O2 = pos2.getAtom("O");

		//dist;
		//TODO: Think of what I can do with this code for identifying interface by CalphaH and carboxyl O
		//a bit harder to do with this than the other because they're not symmetrical hydrogen bonds and not all of the atoms have them?
		//This might be better as a complementary interfacial finder to the previous: maybe as a check for hydrogen bond distances after changing those positions to glycine?
		
	}*/
}
/******************************************
 *  
 *  =======  BEGIN MAIN =======
 *
 ******************************************/
int main(int argc, char *argv[]){

	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];
	
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	
	strftime(buffer,sizeof(buffer),"%m_%d_%Y",timeinfo);
	string date(buffer);
	
	cout << date << endl;

	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;

	Options opt = parseOptions(argc, argv, defaults);
	if (opt.errorFlag) {
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << opt.errorMessages << endl;
		cerr << endl;
		cerr << opt.OPerrors << endl;

		usage();
		exit(1);
	}

	//const char *dirname = opt.datafile + "_8_20_18";
	//_tmkdir(dirname);
	/******************************************************************************
	 *                     === READ IN GEOMETRY FILE ===
	 ******************************************************************************/
	vector<string> fileVec;
	readGeometryFile(opt.helixGeoFile, fileVec);

	// Import sequence, determine length
	//string originalTMSeq1 = opt.originalSeq.substr(opt.startResNum-opt.sequenceStart,opt.endResNum - opt.startResNum + 1);
	unsigned int sequenceLength1 = opt.originalSeq.length(); 
	//sequenceLength = 69;
	// cant do sequences of size less than 4
	cout << "Seq: " << opt.originalSeq << endl;
	cout << "SeqLength:" << sequenceLength1 << endl;
	if(sequenceLength1 < 4) {
		cerr << "Sequence " << opt.originalSeq << "is too small (should be >= 4 AA long)" << endl; 
		exit(0);
	}

	//generateSequence(opt.startResNum, sequenceLength, opt.backboneLength, opt.backboneAA, opt.sequence);
	sequenceLength1 = opt.originalSeq.length();
	
	string modelledTMSeq1 =""; // replace P with A
	string prolineMask1 = "";
	for(int i = 0; i < sequenceLength1; i++) {
		char AA = opt.originalSeq[i];
		if(AA == 'P') {
			modelledTMSeq1 += "A";
			prolineMask1 += "1";
		} else {
			modelledTMSeq1 += AA;
			prolineMask1 += "0";
		}
	}

	/**********************************************************************************
	*
	*    printProteinOutFile
	*
	**********************************************************************************/
	ofstream pout;
	opt.pdbOutputDir = opt.pdbOutputDir + "/" + date + "/" + opt.datafile;
	string dir = opt.pdbOutputDir;
	string cmd = "mkdir -p " + dir;
	if (system(cmd.c_str())){
		cout << "Unable to make directory" << endl;
		exit(0);
	}
	string poutName = dir + "/" + opt.datafile + ".out";
	
	pout.open(poutName.c_str());
	
	printOptions(opt, pout);
	pout << "Original sequence: " << opt.originalSeq << endl;
	pout << "Original length: " << sequenceLength1 << endl;

	string sequence1 = convertToPolymerSequence(modelledTMSeq1,1); // so that the 4th residue will be the middle one (35th) on the GLY 69 backbone
	PolymerSequence PS1(sequence1); 

	// Create system with sequence - a string with the following format
	// A:{startingResNum} ALA ILE ...\n
	// B:{startingResNum} ALA ILE ...
		
	/******************************************************************************
	 *                     === DECLARE SYSTEM ===
	 ******************************************************************************/
	System sys1;
	CharmmSystemBuilder CSB1(sys1,opt.topFile,opt.parFile);
	CSB1.setBuildTerm("CHARMM_ELEC", false);
	CSB1.setBuildTerm("CHARMM_ANGL", false);
	CSB1.setBuildTerm("CHARMM_BOND", false);
	CSB1.setBuildTerm("CHARMM_DIHE", false);
	CSB1.setBuildTerm("CHARMM_IMPR", false);
	CSB1.setBuildTerm("CHARMM_U-BR", false);

	CSB1.setSolvent("MEMBRANE");
	CSB1.setIMM1Params(15, 10);

	//CSB.buildSystemFromPDB(opt.infile);
	//I think I used the above so that I didn't have to use all of the CATM stuff...whoops...well if it works then cool, but otherwise, I might just want to go back to that
	PDBWriter writer;
	PDBWriter writer2;
	PDBWriter writer3;
	PDBWriter writer4;
	PDBWriter writer5;
	PDBWriter writer6;
	sys1.saveAltCoor("start");

	if(!CSB1.buildSystem(PS1)) {
		cerr << "Unable to build system from " << sequence1 << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}
	
	writer.open(dir + "/helixAlign.pdb");
	writer2.open(dir + "/helixAlign2.pdb");
	writer3.open(dir + "/IDswitch.pdb");
	writer4.open(dir + "/MCOFinal.pdb");
	writer5.open(dir + "/SCMF.pdb");
	writer6.open(dir + "/DEE.pdb");
	Chain & chain1A = sys1.getChain("A");
	Chain & chain1B = sys1.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChain1A = chain1A.getAtomPointers();
	AtomPointerVector & apvChain1B = chain1B.getAtomPointers();

	//vector<Position*>& positions = sys.getPositions();

	vector<Position*>& positionsChain1A = sys1.getChain("A").getPositions();
	vector<Position*>& positionsChain1B = sys1.getChain("B").getPositions();


	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(opt.backboneCrd); 
	if(!cRead.read()) {
		cout << "Unable to read " << opt.backboneCrd << endl;
		exit(0);
	}
	cRead.close();

	AtomPointerVector& glyAPV = cRead.getAtomPointers();//*/

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);

	// Objects used for transformations
	Transforms trans1; 
	trans1.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans1.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// Random Number Generator
	RandomNumberGenerator RNG1;
	RNG1.setSeed(opt.seed);

	// Read Rotamer Library File
	SystemRotamerLoader sysRot1(sys1, opt.rotLibFile);
	sysRot1.defineRotamerSamplingLevels();

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sys1.assignCoordinates(glyAPV,false);
	sys1.buildAtoms();
	
	// Add hydrogen bond term
	HydrogenBondBuilder hb1(sys1, opt.hBondFile);
	hb1.buildInteractions(30);//maybe the +charges?

	sys1.saveAltCoor("copiedBackbone");
	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* Eset1 = sys1.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset1->setAllTermsActive();
	Eset1->setTermActive("CHARMM_ELEC", false);
	Eset1->setTermActive("CHARMM_ANGL", false);
	Eset1->setTermActive("CHARMM_BOND", false);
	Eset1->setTermActive("CHARMM_DIHE", false);
	Eset1->setTermActive("CHARMM_IMPR", false);
	Eset1->setTermActive("CHARMM_U-BR", false);

	/******************************************************************************
	*             Delete hydrogen bonds in the 3 (+3) terminal residues
	* 	      Residue numbers 1,2,3 and last,last-1,last-2
	*******************************************************************************/
	if(opt.deleteTerminalHbonds) {
		pout << "\nDeleting hbonds from the N- and C-terminal extension residues\n" << endl;
		deleteTerminalHydrogenBondInteractions(sys1,opt);
	}

	/******************************************************************************
	 *                     === HELICAL AXIS SET UP ===
	 ******************************************************************************/

	string axis = "\
ATOM      1  O   DUM A   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      2  Z   DUM A   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
ATOM      3  O   DUM B   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      4  Z   DUM B   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
END";
	

	PDBReader readAxis;
	if(!readAxis.read(axis)) {
		cerr << "Unable to read axis" << endl;
		exit(0);
	}
	//cout << "4" << endl;

	System helicalAxis;
	//helicalAxis.readPdb(opt.helicalAxisPdbFile);
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	helicalAxis.saveCoor("originState");

	// Declare SelfPairManager and Set Seed
	SelfPairManager spm1;
	spm1.seed(RNG1.getSeed()); 
	
	// Set weights
	Eset1->setWeight("CHARMM_VDW", opt.weight_vdw);
	Eset1->setWeight("SCWRL4_HBOND", opt.weight_hbond);
	Eset1->setWeight("CHARMM_IMM1REF", opt.weight_solv);
	Eset1->setWeight("CHARMM_IMM1", opt.weight_solv);
	pout << "VDW weight: " << Eset1->getWeight("CHARMM_VDW") << " HB weight: " << Eset1->getWeight("SCWRL4_HBOND") << " IMM1REF weight: " << Eset1->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << Eset1->getWeight("CHARMM_IMM1") << endl;
	pout << endl;

	// TEMPORARY SET VARIABLE
	vector<string> variablePositionIdList1;
	for (uint k=0; k < sys1.positionSize(); k++) {
		Position &pos = sys1.getPosition(k);

		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			variablePositionIdList1.push_back(pos.getPositionId());
		}
	}
	sys1.setVariablePositions(variablePositionIdList1);

	spm1.setSystem(&sys1);
	spm1.updateWeights();
	//spm1.setVerbose(true);
	//spm.getMinStates()[0];
	sys1.saveAltCoor("beforeInterfaceLoop");
	//cout << "5" << endl;
	//adds in little positive charges?
	
	//This portion of the code takes all of the parameters from the datafile and sets the system to it, allowing it to be used for system and CSB later (why didn't I just transform those?)
	if (opt.transform == true){
					
		chain1A.renumberChain(opt.thread);
		chain1B.renumberChain(opt.thread);
		
		/******************************************************************************
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		sys1.wipeAllCoordinates();
		sys1.assignCoordinates(glyAPV,false);		
		sys1.buildAllAtoms();
		
		sys1.saveAltCoor("copiedBackbone2");
		//writer.write(sys.getAtomPointers(), true, false, true);
		loadRotamers(sys1, sysRot1, opt.SL);
		
		sys1.saveAltCoor("currentTMFragment");
		//needToCalcEnergiesForThread = true;
	
		transformation(apvChain1A, apvChain1B, axisA, axisB, ori, xAxis, zAxis, opt.ZShift, opt.axialRot, opt.crossAng, opt.XShift, trans1);
		//writer2.write(sys.getAtomPointers(), true, false, true);
		moveZCenterOfCAMassToOrigin(sys1.getAtomPointers(), helicalAxis.getAtomPointers(), trans1);
		sys1.buildAllAtoms();
		writer.write(sys1.getAtomPointers(), true, false, true);
		cout << "Finished Transforming!" << endl;
	}

	/******************************************************************************
	 *                         === EXTENDING HELIX ===
	 ******************************************************************************/
	cout << "Extending helix..." << endl;
	//string originalTMSeq = opt.sequence.substr(opt.startResNum-opt.sequenceStart,opt.endResNum - opt.startResNum + 1);
	unsigned int sequenceLength = opt.originalSeq.length(); 
	//sequenceLength = 69;
	// cant do sequences of size less than 4
	//if(sequenceLength < 4) {
	//	cerr << "Sequence " << originalTMSeq << "is too small (should be >= 4 AA long)" << endl; 
	//	exit(0);
	//}

	generateSequence(opt.startResNum, sequenceLength, opt.backboneLength, opt.backboneAA, opt.sequence, true);
	sequenceLength = opt.sequence.length();
	cout << "SEQlength: " << sequenceLength << endl;	
	string modelledTMSeq =""; // replace P with A
	string prolineMask = "";
	for(int i = 0; i < sequenceLength; i++) {
		char AA = opt.sequence[i];
		if(AA == 'P') {
			modelledTMSeq += "A";
			prolineMask += "1";
		} else {
			modelledTMSeq += AA;
			prolineMask += "0";
		}
	}

	cout << "Extended sequence: " << opt.sequence << endl;
	cout << "New length: " << sequenceLength << endl;

	/******************************************************************************
	 *                   === IDENTIFY INTERFACIAL POSITIONS ===
	 ******************************************************************************/
	vector<int> pos;
	vector<double> dist;
	vector<vector<vector<double>>> posDistVector;
	identifyInterface(sys1, pos, dist);
	cout << "Positions: " << endl;
	vector<int> organizedPos;
	vector<double>organizedDist;
	//Organizes the position and distance combo by from lowest position to largest
	for (uint i=0; i<pos.size(); i++){
		organizedPos.push_back(pos[i]);
		organizedDist.push_back(dist[i]);
		if (organizedPos.size() > 1){
			int count = organizedPos.size()-1;
			double tempDist = 0;
			int tempPos = 0;
			while (count > 0){
				if (organizedPos[count] < organizedPos[count-1]){
					tempPos = organizedPos[count];
					tempDist = organizedDist[count];
					organizedPos[count] = organizedPos[count-1];
					organizedDist[count] = organizedDist[count-1];
					organizedPos[count-1] = tempPos;
					organizedDist[count-1] = tempDist;
					count -= 1;
				}
				else{
					count -= 1;
				}
			}
		}
	}
	// Prints out organized positions and distance combos
	for (uint i=0; i < organizedPos.size(); i++){
		cout << organizedPos[i] << ": " << organizedDist[i] << endl;
	}

	//TODO: Check how accurate the above is with a distance checker for the specific interfacial residues of known sequences
	//TODO: figure out where to go with this position vector
	//how do I initialize these specific positions are amino acids to go through for the spm DEE?
	//TODO: Samson had an idea to set these interfacial positions as a specific vector with 0s and 1s where the 1s are interfacial for the positions
	//then I can add identities for each of these different positions by looping through this vector and adding identities for positions with 1s; might be a better long term way to code it
	vector<int> varPos;
	interface01(sys1, pos, varPos);
	//for (uint i=0; i<varPos.size(); i++){
	//	cout << varPos[i];
	//}
	//cout << endl;
	

	/******************************************************************************
	 *               === ADD IDENTITIES TO INTERFACIAL POSITIONS ===
	 ******************************************************************************/
	//addIdentities(CSB1, sys1, opt.Ids, varPos);
	//switchAA(CSB, sys, varPos);
	//sys1.buildAllAtoms();
	string seq = convertToPolymerSequence(modelledTMSeq, 1, opt.Ids, varPos);
	PolymerSequence PS(seq);
	
	//TODO: figure out where would be best to put these parts to add identities and have it still work
	if(!CSB1.buildSystem(PS)) {
		cerr << "Unable to build system from " << seq << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}


	// Reset the chains for the new PolymerSequence
	Chain & chain2A = sys1.getChain("A");
	Chain & chain2B = sys1.getChain("B");

	// Set up chain A and chain B atom pointer vectors for new PolymerSequence
	AtomPointerVector & apvChain2A = chain2A.getAtomPointers();
	AtomPointerVector & apvChain2B = chain2B.getAtomPointers();
	
	if (opt.transform == true){
					
		chain2A.renumberChain(opt.thread);
		chain2B.renumberChain(opt.thread);
		
		/******************************************************************************
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		sys1.wipeAllCoordinates();
		sys1.assignCoordinates(glyAPV,false);		
		sys1.buildAllAtoms();
		
		sys1.saveAltCoor("copiedBackbone2");
		//writer.write(sys.getAtomPointers(), true, false, true);
		loadRotamers(sys1, sysRot1, opt.SL);
		
		sys1.saveAltCoor("currentTMFragment");
		//needToCalcEnergiesForThread = true;
	
		transformation(apvChain2A, apvChain2B, axisA, axisB, ori, xAxis, zAxis, opt.ZShift, opt.axialRot, opt.crossAng, opt.XShift, trans1);
		//writer2.write(sys.getAtomPointers(), true, false, true);
		moveZCenterOfCAMassToOrigin(sys1.getAtomPointers(), helicalAxis.getAtomPointers(), trans1);
		sys1.buildAllAtoms();
		writer2.write(sys1.getAtomPointers(), true, false, true);
		cout << "Finished Transforming!" << endl;
	}
	/*vector<vector<string>> pos2String;
	pos2String = positionToString(sys1, varPos);

	for (uint i=0; i<pos2String.size(); i++){
		for (uint k=0; k<pos2String[i].size(); k++){
			cout << pos2String[i][k] << " ";
		}
		cout << endl;
	}*/

	//string sequence = convertToPolymerSequence(modelledTMSeq,1); // so that the 4th residue will be the middle one (35th) on the GLY 69 backbone
	//PolymerSequence PS(sequence); 

	// Create system with sequence - a string with the following format
	// A:{startingResNum} ALA ILE ...\n
	// B:{startingResNum} ALA ILE ...
		
	/******************************************************************************
	 *                === DECLARE SYSTEM FOR EXTENDED HELIX ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile);
	CSB.setBuildTerm("CHARMM_ELEC", false);
	CSB.setBuildTerm("CHARMM_ANGL", false);
	CSB.setBuildTerm("CHARMM_BOND", false);
	CSB.setBuildTerm("CHARMM_DIHE", false);
	CSB.setBuildTerm("CHARMM_IMPR", false);
	CSB.setBuildTerm("CHARMM_U-BR", false);

	CSB.setSolvent("MEMBRANE");
	CSB.setIMM1Params(15, 10);

	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from " << seq << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}
	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized
	
	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");
	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	vector<Position*>& positions = sys.getPositions();

	vector<Position*>& positionsChainA = sys.getChain("A").getPositions();
	vector<Position*>& positionsChainB = sys.getChain("B").getPositions();

	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sys.assignCoordinates(glyAPV,false);
	sys.buildAtoms();
	
	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, opt.hBondFile);
	hb.buildInteractions(30);//maybe the +charges?

	sys.saveAltCoor("copiedBackbone");
	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* Eset = sys.getEnergySet();
	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsActive();
	Eset->setTermActive("CHARMM_ELEC", false);
	Eset->setTermActive("CHARMM_ANGL", false);
	Eset->setTermActive("CHARMM_BOND", false);
	Eset->setTermActive("CHARMM_DIHE", false);
	Eset->setTermActive("CHARMM_IMPR", false);
	Eset->setTermActive("CHARMM_U-BR", false);

	// Declare SelfPairManager and Set Seed
	SelfPairManager spm;
	spm.seed(RNG1.getSeed()); 
	
	// Set weights
	Eset->setWeight("CHARMM_VDW", opt.weight_vdw);
	Eset->setWeight("SCWRL4_HBOND", opt.weight_hbond);
	Eset->setWeight("CHARMM_IMM1REF", opt.weight_solv);
	Eset->setWeight("CHARMM_IMM1", opt.weight_solv);
	pout << "VDW weight: " << Eset->getWeight("CHARMM_VDW") << " HB weight: " << Eset->getWeight("SCWRL4_HBOND") << " IMM1REF weight: " << Eset->getWeight("CHARMM_IMM1REF") << " IMM1 weight: " << Eset->getWeight("CHARMM_IMM1") << endl;
	pout << endl;

	// TEMPORARY SET VARIABLE
	// don't think I need this?
	/*vector<string> variablePositionIdList;
	for (uint k=0; k < sys.positionSize(); k++) {
		Position &pos = sys.getPosition(k);

		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			variablePositionIdList.push_back(pos.getPositionId());
		}
	}
	sys.setVariablePositions(variablePositionIdList);*/

	spm.setSystem(&sys);
	spm.updateWeights();
	spm.setVerbose(true);
	spm.getMinStates()[0];
	sys.saveAltCoor("beforeInterfaceLoop");
	//cout << "5" << endl;
	//adds in little positive charges?
	//trans1.smartRmsdAlignment(apvChain1A, apvChainA, 2);
	//TODO: Figure out why this compute monomer code crashes my computer; for now, just keep at since solution has to have lower energy that 0 anyways
	//works well not, but not sure why?
	/******************************************************************************
	 *              === COMPUTE MONOMER ENERGY ===
	 ******************************************************************************/
	map<string,double> monomerEnergyByTerm;
	double monomerEnergy = 0;

	if (opt.inputMonomerE == false) {
		monomerEnergy = opt.monoE_vdw + opt.monoE_hbond + opt.monoE_solv + opt.monoE_solvRef;
		monomerEnergyByTerm["CHARMM_VDW"] = opt.monoE_vdw;
		monomerEnergyByTerm["CHARMM_IMM1"] = opt.monoE_solv;
		monomerEnergyByTerm["CHARMM_IMM1REF"] = opt.monoE_solvRef;
		monomerEnergyByTerm["SCWRL4_HBOND"] = opt.monoE_hbond;
	}
	else {
		monomerEnergy = computeMonomerEnergy(sys, trans, opt, helicalAxis, RNG1, monomerEnergyByTerm, pout, opt.greedyCycles, opt.mcCycles, opt.mcMaxRejects);
	}

	cout << setiosflags(ios::fixed) << setprecision(4) << "Monomer Energy: " << monomerEnergy << endl;
	cout << setiosflags(ios::fixed) << setprecision(4) << "Current Energy: " << sys.calcEnergy() << endl;
	cout << endl;

	/******************************************************************************
	 *                  === READ IN RULES ===
	 ******************************************************************************/
	//map<int, string> rulesFileMap;
	//readRulesFile(opt.rulesFile, rulesFileMap);
	//TODO: Could be a good idea for me to make a new rules file? Or use this rules file to check my final sequence

	/******************************************************************************
	 *                  === NECESSARY VARIABLES ===
	 ******************************************************************************/
	// To store the structures in case we are clustering
	vector<HelixDimer*> structures;

	// We will need a format converter to convert to PDB names
	FormatConverter fc;

	map<string,double> energyByTerm;

	// Only calculate energies if necessary
	bool needToCalcEnergiesForThread = true;
	vector<vector<vector<vector<bool> > > > savedEnergyFlagTable;

	// Store info of top MC Repack structures
	vector<vector<double> > topMCRepackVector;

	// Global lowest energy found (if above monomer we won't save anyways)
	double globalLowestE = monomerEnergy;

	if (opt.transform == true){
		//if (sequenceLength % 2 == 0){
		chainA.renumberChain(opt.bbThread);
		chainB.renumberChain(opt.bbThread);
		//}	
		//else{
		//	chainA.renumberChain(opt.bbThread+1);
		//	chainB.renumberChain(opt.bbThread+1);
		//}
		/******************************************************************************
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		//cout << "Polymer Seq: \n" << PS << endl;
		sys.wipeAllCoordinates();
		sys.assignCoordinates(glyAPV,false);		
		sys.buildAllAtoms();

		//vector<int> vP;
		//vP.push_back(3);
	
		//addIdentities(CSB, sys, opt.Ids, vP);
		//cout << "Polymer Seq added: \n" << PS << endl;
		
		sys.saveAltCoor("copiedBackbone2");
		//writer.write(sys.getAtomPointers(), true, false, true);
		loadRotamers(sys, sysRot, opt.SL);
		
		sys.saveAltCoor("currentTMFragment");
		//needToCalcEnergiesForThread = true;
	
		//transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, opt.zShift, opt.axialRotation, opt.crossingAngle, opt.xShift, trans);
		trans.smartRmsdAlignment(apvChainA, apvChain2A, 2);
		trans.smartRmsdAlignment(apvChainB, apvChain2B, 2);
		cout << setiosflags(ios::fixed) << setprecision(4) << "Parameters: crossingAngle: " << opt.crossAng << " axialRotation: " << opt.axialRot << " zShift: " << opt.ZShift;
		cout << " xShift: " << opt.XShift << endl;
		//writer2.write(sys.getAtomPointers(), true, false, true);
		moveZCenterOfCAMassToOrigin(sys1.getAtomPointers(), helicalAxis.getAtomPointers(), trans);
		repackSideChains(spm, opt.greedyCycles, savedEnergyFlagTable);
		sys.buildAllAtoms();
		writer2.write(sys.getAtomPointers(), true, false, true);
		cout << "Finished Transforming!" << endl;
	}

	/******************************************************************************
	 *               === IDENTIFY INTERFACIAL POSITIONS -- OLD ===
	 ******************************************************************************/
	/*vector<int> pos;
	vector<double> dist;
	vector<vector<vector<double>>> posDistVector;
	identifyInterface(sys, pos, dist);
	cout << "Positions: " << endl;
	vector<int> organizedPos;
	//TODO: make below code print out in order for ease of understanding
	for (uint i=0; i < pos.size(); i++){
		cout << pos[i] << ": " << dist[i] << endl;
	}

	//TODO: Check how accurate the above is with a distance checker for the specific interfacial residues of known sequences
	//TODO: figure out where to go with this position vector
	//how do I initialize these specific positions are amino acids to go through for the spm DEE?
	//TODO: Samson had an idea to set these interfacial positions as a specific vector with 0s and 1s where the 1s are interfacial for the positions
	//then I can add identities for each of these different positions by looping through this vector and adding identities for positions with 1s; might be a better long term way to code it
	vector<int> varPos;
	interface01(sys, pos, varPos);
	//for (uint i=0; i<varPos.size(); i++){
	//	cout << varPos[i];
	//}
	//cout << endl;
	

	/******************************************************************************
	 *               === ADD IDENTITIES TO INTERFACIAL POSITIONS ===
	 ******************************************************************************/
	/*addIdentities(CSB, sys, opt.Ids, varPos);
	//switchAA(CSB, sys, varPos);
	sys.buildAllAtoms();
	//string seq = convertToPolymerSequence(modelledTMSeq, 1, opt.Ids, varPos);
	//PolymerSequence altPS(seq);
	
	//TODO: figure out where would be best to put these parts to add identities and have it still work
	//if(!CSB.buildSystem(altPS)) {
	//	cerr << "Unable to build system from " << sequence << endl;
	//	exit(0);
	//} else {
	//	//fout << "CharmmSystem built for sequence" << endl;
	//}*/

	vector<vector<string>> pos2String;
	pos2String = positionToString(sys, varPos);

	for (uint i=0; i<pos2String.size(); i++){
		for (uint k=0; k<pos2String[i].size(); k++){
			cout << pos2String[i][k] << " ";
		}
		cout << endl;
	}
	
	sys.setLinkedPositions(pos2String);
	//error 48923...different number of rotamers for identity at position 34, the one that is being difficult for me
	//loadRotamers(sys, sysRot, opt.SL);
	//TODO: figure out how to change the fixed energy before setting the system
	//the way to do that is to set only one rotamer for the rest of the positions
	
	// TEMPORARY SET VARIABLE
	/*vector<string> variablePositionIdList;
	for (uint k=0; k < sys.positionSize(); k++) {
		Position &pos = sys.getPosition(k);

		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			variablePositionIdList.push_back(pos.getPositionId());
		}
	}
	sys.setVariablePositions(variablePositionIdList);*/
	
	spm.setSystem(&sys);
	sys.buildAllAtoms();
	//spm.setVerbose(opt.verbose);
	spm.setRunDEE(true, true);
	repackSideChains(spm, opt.greedyCycles, savedEnergyFlagTable);
	spm.calculateEnergies();
	
	cout << "Current Energy: " << spm.getMinBound()[0];

	spm.setRunSCMF(false);
	spm.setRunUnbiasedMC(false);
	spm.runOptimizer();
	spm.runGreedyOptimizer(10);

	//TODO: Figure out how to set cutoffs to not only get the best state, but a variety of states
	//TODO: how can I set it so that the optimizer results in the same switches for each position the same rather than different?
	//I can think of adding a boolean to the actual code, but is this the best option? Maybe there's already one somewhere in the code
	//I think I found it in System.h: void setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions); So I need to transition the positions to this "A,19" "B,19" format!
	vector<vector<unsigned int>> rotamerStates = spm.getDEEAliveRotamers();
	vector<unsigned int> rotamerStates1;

	//TODO: I think I need to understand a bit more about the CATM code? For some reason I'm not moving the backbone the way I expect; get these xShift repacks to work

	/******************************************************************************
	 *                     === X SHIFT REPACKS ===
	 ******************************************************************************/
	/*double xShift = opt.XShift;
	double currentEnergy = 0;
	double bestEnergy = currentEnergy;
	double savedXShift = opt.XShift;
	double previousEnergy = monomerEnergy;
	double deltaXShift = 0.1;
	double xShiftEnd = -5;

	int counterx = 0;
	while (xShift >= xShiftEnd) {
	
		xShift += deltaXShift;

		// Move the helix
		backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, 3 );
	
		// Run Optimization
		repackSideChains(spm, opt.greedyCycles, savedEnergyFlagTable);
	
		vector<unsigned int> MCOFinal;
		MCOFinal = spm.getMinStates()[0];
		sys.setActiveRotamers(MCOFinal);
	
		currentEnergy = spm.getMinBound()[0];
	
		if (currentEnergy < bestEnergy) {
			bestEnergy = currentEnergy;
			savedXShift = xShift;
			sys.saveAltCoor("savedBestState");
			helicalAxis.saveAltCoor("BestAxis");
			sys.applySavedCoor("savedBestState");
			sys.buildAllAtoms();
			writer3.write(sys.getAtomPointers(), true, false, true);
		}
	
		cout << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;

		//TODO: Make this run until completion of movements or until it decreases below the global lowest at least once
		if (currentEnergy > globalLowestE){
			counterx++;
			previousEnergy = currentEnergy;
		}
		else{
			globalLowestE = currentEnergy;
			cout << "Energy finally below monomer energy! Current Energy: " << currentEnergy << endl;
			break;
		}

	}
	cout << "Best Energy at x shift: " << bestEnergy-monomerEnergy << " at " << savedXShift << endl;
	cout << "Counter: " << counterx << endl;*/
	//vector<unsigned int> MCOFinal;
	//for (uint i=0; i<spm.getMinStates().size(); i++){
	//	MCOFinal = spm.getMinStates()[i];
	//	sys.setActiveRotamers(MCOFinal);
	//	sys.buildAllAtoms();
	//	writer4.write(sys.getAtomPointers(), true, false, true);
	//}
	////spm.getSummary(rotamerStates);
	////spm.getStatePositionIdentityRotamerIndeces(rotamerStates1);
	////the above do not help with writing states!
	////cout << "Alive Rots: " << endl;
	////for (uint i=0; i<rotamerStates.size(); i++){
	////	cout << rotamerStates[i] << " ";
	////}
	//cout << endl;
	//cout << "Min Energy: " << spm.getMinBound()[0] << endl;
	//cout << "SCMF States: " << endl;
	//sys.setActiveRotamers(spm.getSCMFstate());
	//sys.buildAllAtoms();
	//writer5.write(sys.getAtomPointers(), true, false, true);
	//cout << "SPM Energy: " << spm.getStateEnergy(spm.getSCMFstate()) << endl;
	//cout << "MC States: " << endl;
	//for (uint i=0; i<spm.getBestUnbiasedMCState().size(); i++){
	//	cout << spm.getBestUnbiasedMCState()[i] << " ";
	//}
	//sys.setActiveRotamers(rotamerStates);
	sys.buildAllAtoms();
	writer6.write(sys.getAtomPointers(), true, false, true);
	cout << endl;
	cout << "System Energy: " << sys.calcEnergy() << endl;

	//TODO: Shouldn't system energy and spm energy be the same if in the proper state? Even when building into those states though, the energy is higher than the spm. Find our why this could be

	/******************************************************************************
	 *      === LOCAL BACKBONE MONTE CARLO REPACKS FOR TOP STRUCTURES ===
	 ******************************************************************************/
	/*xShift = opt.XShift;
	double crossingAngle = opt.crossAng;
	double axialRotation = opt.axialRot;
	double zShift = opt.ZShift;

	// Optimizatize Initial Starting Position
	repackSideChains(spm, opt.greedyCycles, savedEnergyFlagTable);

	sys.setActiveRotamers(spm.getMinStates()[0]);
	currentEnergy = spm.getMinBound()[0];
	sys.saveAltCoor("savedBestState");
	helicalAxis.saveAltCoor("BestAxis");

	/******************************************************************************
	 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
	 ******************************************************************************/
	/*bestEnergy = currentEnergy;
	time_t startTimeMC, endTimeMC;
	double diffTimeMC;
	time(&startTimeMC);	

	if (opt.mcCycles > 0) {
		//MonteCarloManager MCMngr(1000.0, 0.5, opt.MCCycles, MonteCarloManager::EXPONENTIAL, opt.MCMaxRejects);
		MonteCarloManager MCMngr(opt.mcStartTemp, opt.mcEndTemp, opt.mcCycles, opt.mcCurve, opt.mcMaxRejects);

		MCMngr.setEner(bestEnergy);
		
		while(!MCMngr.getComplete()) {

			sys.applySavedCoor("savedBestState");
			helicalAxis.applySavedCoor("BestAxis");

			int moveToPreform = RNG1.getRandomInt(3);

			double deltaXShift = 0.0;
			double deltaZShift = 0.0;
			double deltaCrossingAngle = 0.0;
			double deltaAxialRotation = 0.0; 

			//======================================
			//====== Z Shift (Crossing Point) ======
			//======================================
			if (moveToPreform == 0) {
				//deltaZShift = getStandardNormal(RNG1) * 0.1;
				deltaZShift = getStandardNormal(RNG1) * opt.deltaZ;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaZShift, moveToPreform);
			} else if (moveToPreform == 1) {
			//===========================
			//===== Axial Rotation ======
			//===========================
				//deltaAxialRotation = getStandardNormal(RNG1) * 1.0;
				deltaAxialRotation = getStandardNormal(RNG1) * opt.deltaAx;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaAxialRotation, moveToPreform);
			} else if (moveToPreform == 2) {
			//==================================
			//====== Local Crossing Angle ======
			//==================================
				//deltaCrossingAngle = getStandardNormal(RNG1) * 1.0;
				deltaCrossingAngle = getStandardNormal(RNG1) * opt.deltaCross;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaCrossingAngle, moveToPreform);
			} else if (moveToPreform == 3) {
			//==============================================
			//====== X shift (Interhelical Distance) =======
			//==============================================
				//deltaXShift = getStandardNormal(RNG1) * 0.1;
				deltaXShift = getStandardNormal(RNG1) * opt.deltaX;
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, moveToPreform);
			}

			// Run Optimization
			repackSideChains(spm, opt.greedyCycles, savedEnergyFlagTable);
	
			vector<unsigned int> MCOFinal = spm.getMinStates()[0];
			sys.setActiveRotamers(MCOFinal);
			currentEnergy = spm.getMinBound()[0];

			if (!MCMngr.accept(currentEnergy)) {
				//fout << "state rejected   energy: " << currentEnergy << endl;
			}
			else {
				// check number of hydrogen bonds and if it is less than 4 dont accept
				// TODO: hardcoded 4 
				vector<string> interHelicalHbonds =  getInterHelicalHbonds(Eset);
				if(interHelicalHbonds.size() < 4) {
					// dont accept if we lost hydrogen bonds
					pout << "state rejected, lost Hbonds" << endl;
				} else {
					bestEnergy = currentEnergy;
					sys.saveAltCoor("savedBestState");
					helicalAxis.saveAltCoor("BestAxis");

					xShift = xShift + deltaXShift;
					crossingAngle = crossingAngle + deltaCrossingAngle;
					axialRotation = axialRotation + deltaAxialRotation;
					zShift = zShift +  deltaZShift;

					pout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
				}
			}
		}
	}
	time(&endTimeMC);
	diffTimeMC = difftime (endTimeMC, startTimeMC);
	cout << "Monte Carlo repack complete. Time: " << diffTimeMC << endl << endl;

	sys.applySavedCoor("savedBestState");

	double finalEnergy = spm.getMinBound()[0]-monomerEnergy;
	if(opt.printTermEnergies) {
		energyByTerm = getEnergyByTerm(sys.getEnergySet());

		for(map<string,double>::iterator it = monomerEnergyByTerm.begin(); it != monomerEnergyByTerm.end(); it++) {
			if(energyByTerm.find(it->first) != energyByTerm.end()) {
				energyByTerm[it->first] -= monomerEnergyByTerm[it->first];
			} else {
				// impossible
				pout << "ERROR 1242 Missing Energy term " << it->first << " Exiting" << endl; 
				exit(0);
			}
		}
	}
	if (finalEnergy > 0) { // TODO Check number of Hbonds???
		// renumber residues - rethreading
		cout << "Energy " << finalEnergy << " above zero - not written" << endl;
		cout << endl;
		//reThreadResidues(positions,-1);

		//continue;
	}
	
	//TODO: Figure out what I need from this bottom portion: I think if I do it properly, I can use the whole of the CATM code but instead of threading, with alternate sequences and keeping say the 5 or 10 best ones?
	*//*/ Renumber structures (all should go 1 - N, or from starting point determined by user)
	//int currStartNum = positions[0]->getResidueNumber();
	chainA.renumberChain(opt.startResNum);
	chainB.renumberChain(opt.startResNum);
	//renumberResidues(sys,opt.startResNum);

	// get vector of Hbonds, check size
	vector<string> interHelicalHbonds =  getInterHelicalHbonds(Eset);

	// Store everything
	// create class structure to store data
	HelixDimer * st = new HelixDimer("",finalEnergy,opt.thread);
	
	st->addAtoms(sys.getAtomPointers());
	// save the structure interface residues, energy, hbond list, CA of closest approach mark it with 2 in the mask 
	map<string, unsigned int> interfaceMap = interfaceResidueCheck(apvChainA, apvChainB);

	//unsigned int closestCA = CAOfClosestApproach(chainA,chainB);
	string interfaceString = "";
	for (uint r=0; r < chainA.positionSize(); r++) {
		sequence = sequence + MslTools::getOneLetterCode(chainA.getPosition(r).getResidueName());

		if (interfaceMap.find(chainA.getIdentity(r).getIdentityId()) != interfaceMap.end()) {
			interfaceString = interfaceString + "1";
		}
		else {
			interfaceString = interfaceString + "0";
		}
	}

	//mark residue with CA of closest approach with 2
	//interfaceString.replace(closestCA - opt.startResNum,1,"2");

	st->setHelixDimerDetails(xShift,zShift,axialRotation,crossingAngle,interfaceString,prolineMask,interHelicalHbonds,0);
	if(opt.printTermEnergies) {
		st->setDeltaEnergyByTerm(energyByTerm);
	}

	if (opt.clusterSolutions) {
		// Save structure to vector of structures
		cout << "clustering structure with Energy: " << finalEnergy << endl;
		cout << endl;
		structures.push_back(st);
		
		AtomPointerVector& axisAtoms = helicalAxis.getAtomPointers();
		for(AtomPointerVector::iterator it =  axisAtoms.begin(); it != axisAtoms.end(); it++) {
			Atom* a = new Atom(* *it);
			st->addAxisAtom(a);
		}
		
	}
	else { 
		char tmp[1000];


		sprintf(tmp,"%s_A_%02d_B_%02d_x%05.3f_z%05.3f_a%06.3f_c%05.3f.crd", opt.uniprotAccession.c_str(), positions[0]->getResidueNumber(), positions[positions.size()/2]->getResidueNumber(), xShift, zShift, axialRotation, crossingAngle);
		CRDWriter crd;
		string crdName = opt.pdbOutputDir+"/"+string(tmp);
		crd.open(crdName);
		crd.write(sys.getAtomPointers());
		pout << "wrote file: " << string(tmp) << " Energy: " << finalEnergy << endl;
		crd.close();

		// TODO convert to PDB names before printing
		sprintf(tmp,"%s_A_%02d_B_%02d_x%05.3f_z%05.3f_a%06.3f_c%05.3f.pdb", opt.uniprotAccession.c_str(), positions[0]->getResidueNumber(), positions[positions.size()/2]->getResidueNumber(), xShift, zShift, axialRotation, crossingAngle);
		pout << "wrote file: " << string(tmp) << " Energy: " << finalEnergy << endl;
		sys.writePdb(opt.pdbOutputDir+"/"+string(tmp));

		sprintf(tmp,"%s_A_%02d_B_%02d_x%05.3f_z%05.3f_a%06.3f_c%05.3f.txt", opt.uniprotAccession.c_str(), positions[0]->getResidueNumber(), positions[positions.size()/2]->getResidueNumber(), xShift, zShift, axialRotation, crossingAngle);
		pout << "wrote text file: " << string(tmp) << endl;
		
		ofstream txt;
		txt.open(tmp);
		if(txt.is_open()) {
			txt << "resStart " << opt.startResNum << endl;
			txt << "resEnd " << opt.endResNum << endl;
			txt << "TMStart " << opt.startResNum << endl;
			txt << "TMEnd " << opt.endResNum << endl;
			txt << "originalSeq " << opt.originalSeq << endl;
			txt << "modelledSeq " << modelledTMSeq << endl;
			st->printHelixDimerDetails(txt);
			delete st;

		} else {
			cerr << "Unable to open " << tmp << endl;
		}
		

		cout << sys.getEnergySummary();
		cout << endl;
	}



	/******************************************************************************
	 *                   === CLUSTER CATM MODEL SOLUTIONS ===
	 ******************************************************************************/
	/*if (opt.clusterSolutions) {
		sort(structures.begin(),structures.end(),compareHelixDimers());
		// convert all structures to PDB names
		/*
		for(int i = 0; i < structures.size(); i++) {
			How do we handle changed hydrogen bond names?
			AtomPointerVector & a = structures[i]->getAtomPointers();
			for(int j = 0; j < a.size(); j++) {
				int resNum = a[j]->getResidueNumber();
				fc.setPdbFromCharmm(*a[j],"22",resNum == opt.startResNum,resNum == opt.endResNum);
			}
		}
		*/
		/*vector<HelixDimerCluster*> clusters;
		clusterSolutions(clusters,structures,opt.rmsdCutoff,opt.originalSeq,modelledTMSeq,opt);
		//cout << "Output " << clusters.size() << " HelixDimerClusters at rmsd " << _rmsdCutoff <<  endl;
		for(int i = 1; i <= clusters.size(); i++) {
			clusters[i-1]->printHelixDimerClusterCrds(i,opt.printAllCrds,opt.printAxes);
			clusters[i-1]->printDetails(i);
			if(opt.printTermEnergies) {
				clusters[i-1]->printTermEnergies(i);
			}
		}


		for(int i = 1; i <= clusters.size(); i++) {
			//clusters[i-1]->convertToPdbNames();
			clusters[i-1]->printHelixDimerClusterPdbs(i);
			clusters[i-1]->makePse(i);
		}

		// TODO delete HelixDimerCluster* form the vector clusters

	}*/

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	cout << endl << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;
	
	writer.close();
	writer2.close();
	writer3.close();
	writer4.close();
	writer5.close();
	writer6.close();
}

Options parseOptions(int _argc, char * _argv[], Options defaults){
	
	/******************************************
	 *  Pass the array of argument and the name of
	 *  a configuration file to the ArgumentParser
	 *  object.  Then ask for the value of the argument
	 *  and collect error and warnings.
	 *
	 *  This function returns a Options structure
	 *  defined at the head of this file 
	 ******************************************/
	
	Options opt;
	Options data;

	/******************************************
	 *  Set the allowed and required options:
	 *
	 *  Example of configuration file:
	 *
	 *  /exports/home/gloiseau/mslib/trunk_AS/myProgs/gloiseau/helixGenerator.config
	 *  
	 ******************************************/

	vector<string> required;
	vector<string> allowed;

	//opt.required.push_back("");
	//opt.allowed.push_back("");

	//opt.allowed.push_back("");
	
	opt.allowed.push_back("sequence");
	opt.allowed.push_back("backboneAA");
	opt.allowed.push_back("backboneLength");
	
	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("sequenceStart");
	opt.allowed.push_back("endResNum");
	opt.allowed.push_back("tmStart");
	opt.allowed.push_back("tmEnd");

	opt.allowed.push_back("threadStart");
	opt.allowed.push_back("threadEnd");
	opt.allowed.push_back("threadBool");

	opt.allowed.push_back("xShift");
	opt.allowed.push_back("zShift");
	opt.allowed.push_back("axialRotation");
	opt.allowed.push_back("crossingAngle");
	opt.allowed.push_back("transform");
	//opt.allowed.push_back("thread");
	//opt.allowed.push_back("bbThread");

	opt.allowed.push_back("mcCycles");
	opt.allowed.push_back("mcMaxRejects");
	opt.allowed.push_back("mcStartTemp");
	opt.allowed.push_back("mcEndTemp");
	opt.allowed.push_back("mcCurve");

	opt.allowed.push_back("deltaZ");
	opt.allowed.push_back("deltaAx");
	opt.allowed.push_back("deltaCross");
	opt.allowed.push_back("deltaX");

	opt.allowed.push_back("weight_vdw");
	opt.allowed.push_back("weight_hbond");
	opt.allowed.push_back("weight_solv");
	
	opt.allowed.push_back("AAs");
	opt.allowed.push_back("SL");
	
	opt.allowed.push_back("start");
	opt.allowed.push_back("end");
	
	opt.allowed.push_back("ener");
	
	opt.allowed.push_back("ivalues");

	opt.allowed.push_back("verbose");
	opt.allowed.push_back("greedyOptimizer");
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");
	opt.allowed.push_back("numberOfStructuresToMCRepack");
	opt.allowed.push_back("energyCutOff");
	
	opt.allowed.push_back("uniprotName");
	opt.allowed.push_back("uniprotAccession");

	opt.allowed.push_back("inputMonomerE");
	opt.allowed.push_back("monoE_vdw");
	opt.allowed.push_back("monoE_hbond");
	opt.allowed.push_back("monoE_solv");
	opt.allowed.push_back("monoE_solvRef");

	opt.allowed.push_back("rmsdCutoff");
	opt.allowed.push_back("clusterSolutions");
	opt.allowed.push_back("printAllCrds");
	opt.allowed.push_back("printAxes");
	opt.allowed.push_back("printTermEnergies");
	opt.allowed.push_back("deleteTerminalHbonds");

	opt.allowed.push_back("helixGeoFile");
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("rotLibFile");
	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("hbondFile");
	opt.allowed.push_back("pdbOutputDir");
	opt.allowed.push_back("infile");
	opt.allowed.push_back("rulesFile");
	opt.allowed.push_back("configfile");
	
	opt.allowed.push_back("datafile");
	opt.allowed.push_back("originalSeq");
	opt.allowed.push_back("modelledSeq");
	opt.allowed.push_back("resStart");
	opt.allowed.push_back("resEnd");
	opt.allowed.push_back("seqLength");
	opt.allowed.push_back("thread");
	opt.allowed.push_back("bbThread");
	opt.allowed.push_back("XShift");
	opt.allowed.push_back("ZShift");
	opt.allowed.push_back("axialRot");
	opt.allowed.push_back("crossAng");
	
	opt.allowed.push_back("Ids");

	OptionParser OP;
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions();

	if (OP.countOptions() == 0){
		usage();
		cerr << "No options given!" << endl;
		exit(0);
	}
	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
			exit(1);
		}
	}
	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return opt;
	}
	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.errorMessages = "";
	opt.warningMessages = "";

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.clusterSolutions = OP.getBool("clusterSolutions");
	if (OP.fail()) {
		opt.clusterSolutions = true;
		opt.warningMessages += "clusterSolutions not specified using true\n";
		opt.warningFlag = true;
	}
	opt.printAllCrds = OP.getBool("printAllCrds");
	if (OP.fail()) {
		opt.warningMessages += "printAllCrds not specified using false\n";
		opt.warningFlag = true;
		opt.printAllCrds = false;
	}
	opt.printAxes = OP.getBool("printAxes");
	if (OP.fail()) {
		opt.warningMessages += "printAxes not specified using false\n";
		opt.warningFlag = true;
		opt.printAxes = false;
	}
	opt.printTermEnergies = OP.getBool("printTermEnergies");
	if (OP.fail()) {
		opt.printTermEnergies = true;
		opt.warningMessages += "printTermEnergies not specified using true\n";
		opt.warningFlag = true;
	}
	opt.deleteTerminalHbonds = OP.getBool("deleteTerminalHbonds");
	if (OP.fail()) {
		opt.deleteTerminalHbonds = true;
		opt.warningMessages += "deleteTerminalHbonds not specified using true\n";
		opt.warningFlag = true;
	}
	opt.rmsdCutoff = OP.getDouble("rmsdCutoff");
	if (OP.fail()) {
		opt.rmsdCutoff = 2.0;
		opt.warningMessages += "rmsdCutoff not specified using 2.0\n";
		opt.warningFlag = true;
	}

	opt.sequence = OP.getString("sequence");
	if (OP.fail()) {
		opt.errorMessages += "sequence (1 letter aa) not specified\n";
		opt.errorFlag = true;
	}
	
	opt.tmStart = OP.getInt("tmStart");
	if(OP.fail()) {
		opt.warningMessages += "tmStart not specified using 1\n";
		opt.warningFlag = true;
		opt.tmStart = 1;
	}

	opt.tmEnd = OP.getInt("tmEnd");
	if(OP.fail()) {
		opt.tmEnd = opt.sequence.length();
		opt.warningMessages += "tmEnd not specified using " + MslTools::intToString(opt.tmEnd) + "\n";
		opt.warningFlag = true;
	}

	opt.rulesFile = OP.getString("rulesFile");
	if (OP.fail()) {
		opt.rulesFile = "/data01/sabs/tmRepacks/GLY_69_Homo_2/tmRules/rules_10kcals_vdw_only/tmRules.out";
		opt.warningMessages += "rulesFile not specified using " + opt.rulesFile + "\n";
		opt.warningFlag = true;
	}

	opt.sequenceStart = OP.getInt("sequenceStart");
	if (OP.fail()) {
		opt.warningMessages += "sequenceStart not specified using 1\n";
		opt.warningFlag = true;
		opt.sequenceStart = 1;
	}

	opt.startResNum = OP.getInt("startResNum");
	if (OP.fail()) {
		opt.warningMessages += "startResNum not specified using " + MslTools::intToString(opt.tmStart) + "\n";
		opt.warningFlag = true;
		opt.startResNum = opt.tmStart;
	}

	opt.endResNum = OP.getInt("endResNum");
	if (OP.fail()) {
		opt.warningMessages += "endResNum not specified using " + MslTools::intToString(opt.tmEnd) + "\n";
		opt.warningFlag = true;
		opt.endResNum = opt.tmEnd;
	}

	opt.threadStart = OP.getInt("threadStart");
	if (OP.fail()) {
		// excluded 2 residues in the beginning
		// opt.threadStart = 35 - ( (tmEnd-tmStart + 1) - 3 );  
		opt.threadStart = opt.startResNum + 37 - opt.endResNum;  
		opt.warningMessages += "threadStart not specified using " + MslTools::intToString(opt.threadStart) + "\n";
		opt.warningFlag = true;
	}
	opt.threadEnd = OP.getInt("threadEnd");
	if (OP.fail()) {
		opt.threadEnd = 33;
		opt.warningMessages += "threadEnd not specified using " + MslTools::intToString(opt.threadEnd) + "\n";
		opt.warningFlag = true;
	}
	opt.threadBool = OP.getBool("threadBool");
	if (OP.fail()) {
		opt.warningMessages += "threadBool not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.threadBool = false;
	}

	opt.xShift = OP.getDouble("xShift");
	if (OP.fail()) {
		opt.warningMessages += "xShift not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.xShift = 6.7;
	}
	opt.zShift = OP.getDouble("zShift");
	if (OP.fail()) {
		opt.warningMessages += "zShift not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.zShift = 0;
	}
	opt.axialRotation = OP.getDouble("axialRotation");
	if (OP.fail()) {
		opt.warningMessages += "axialRotation not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.axialRotation = 0;
	}
	opt.crossingAngle = OP.getDouble("crossingAngle");
	if (OP.fail()) {
		opt.warningMessages += "crossingAngle not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.crossingAngle = -40;
	}
	opt.transform = OP.getBool("transform");
	if (OP.fail()) {
		opt.warningMessages += "transform not specified, defaulting to false\n";
		opt.warningFlag = true;
		opt.transform = false;
	}
//	opt.thread = OP.getInt("thread");
//	if (OP.fail()) {
//		opt.warningMessages += "thread not specified, defaulting to 0\n";
//		opt.warningFlag = true;
//		opt.thread = 0;
//	}
//	opt.bbThread = opt.thread + opt.sequence.length() - opt.backboneLength + 1;
		
	opt.mcCycles = OP.getInt("mcCycles");
	if (OP.fail()) {
		opt.errorMessages += "Number of MC cycles not specified!\n";
		opt.errorFlag = true;
	}

	opt.mcMaxRejects = OP.getInt("mcMaxRejects");
	if (OP.fail()) {
		opt.mcMaxRejects = 10;
		opt.warningMessages += "Number of MC max rejects not specified, default to using 10\n";
		opt.warningFlag = true;
	}

	opt.mcStartTemp = OP.getDouble("MCStartTemp");
	if (OP.fail()) {
		opt.warningMessages += "MCStartTemp not specified using 1000.0\n";
		opt.warningFlag = true;
		opt.mcStartTemp = 1000.0;
	}
	opt.mcEndTemp = OP.getDouble("MCEndTemp");
	if (OP.fail()) {
		opt.warningMessages += "MCEndTemp not specified using 0.5\n";
		opt.warningFlag = true;
		opt.mcEndTemp = 0.5;
	}
	opt.mcCurve = OP.getInt("MCCurve");
	if (OP.fail()) {
		opt.warningMessages += "MCCurve not specified using EXPONENTIAL(2)\n";
		opt.warningFlag = true;
		opt.mcCurve = 2;
	}

	opt.deltaZ = OP.getDouble("deltaZ");
	if (OP.fail()) {
		opt.warningMessages += "deltaZ not specified using 0.1\n";
		opt.warningFlag = true;
		opt.deltaZ = 0.1;
	}
	opt.deltaAx = OP.getDouble("deltaAx");
	if (OP.fail()) {
		opt.warningMessages += "deltaAx not specified using 1.0\n";
		opt.warningFlag = true;
		opt.deltaAx = 1.0;
	}
	opt.deltaCross = OP.getDouble("deltaCross");
	if (OP.fail()) {
		opt.warningMessages += "deltaCross not specified using 1.0\n";
		opt.warningFlag = true;
		opt.deltaCross = 1.0;
	}
	opt.deltaX = OP.getDouble("deltaX");
	if (OP.fail()) {
		opt.warningMessages += "deltaX not specified using 0.1\n";
		opt.warningFlag = true;
		opt.deltaX = 0.1;
	}

	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages += "verbose not specified using false\n";
		opt.warningFlag = true;
		opt.verbose = false;
	}
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 1\n";
		opt.warningFlag = true;
		opt.greedyCycles = 1;
	}
	opt.numberOfStructuresToMCRepack = OP.getInt("numberOfStructuresToMCRepack");
	if (OP.fail()) {
		opt.warningMessages += "numberOfStructuresToMCRepack not specified using 20\n";
		opt.warningFlag = true;
		opt.numberOfStructuresToMCRepack = 20;
	}
	opt.energyCutOff = OP.getDouble("energyCutOff");
	if (OP.fail()) {
		opt.warningMessages += "energyCutOff not specified using 100.0\n";
		opt.warningFlag = true;
		opt.energyCutOff = 100.0;
	}
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.seed = 1;
		opt.warningMessages += "Seed not specified!\n";
		opt.warningFlag = true;
	}

	opt.weight_vdw = OP.getDouble("weight_vdw");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_vdw not specified, default 1.0\n";
		opt.weight_vdw = 1.0;
	}
	opt.weight_hbond = OP.getDouble("weight_hbond");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_hbond not specified, default 1.0\n";
		opt.weight_hbond = 1.0;
	}
	opt.weight_solv = OP.getDouble("weight_solv");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "weight_solv not specified, default 1.0\n";
		opt.weight_solv = 1.0;
	}

	opt.AAs = OP.getMultiString("AAs");
	if (OP.fail()) {
		opt.errorMessages += "Amino acids not specified!\n";
		opt.errorFlag = true;
	}
	if (opt.AAs.size() > 2){
		opt.errorMessages += "Too many amino acids!\n";
		opt.errorFlag = true;
	}
	
	opt.SL = OP.getString("SL");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "SL not specified, default to SL70\n";
		opt.SL = "SL70.00";
	}
	else {
		opt.SL = "SL" + opt.SL;
	}

	opt.backboneAA = OP.getString("backboneAA");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "backboneAA not specified, default to glycine\n";
		opt.backboneAA = "G";
	}
	opt.backboneLength = OP.getInt("backboneLength");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "backboneLength not specified, default to 35\n";
		opt.backboneLength = 35;
	}

	opt.start = OP.getInt("start");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "Start point not specified, default to 0\n";
		opt.start = 0;
	}
	opt.end = OP.getInt("end");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "End point not specified, default to 60\n";
		opt.end = 60;
	}
	
	opt.ener = OP.getInt("ener");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "Energy not specified, default to -15\n";
		opt.end = -15;
	}

	opt.ivalues = OP.getMultiInt("ivalues");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "i+x values not specified, default to 1, 3, and 4\n";
		opt.ivalues.push_back(1);
		opt.ivalues.push_back(3);
		opt.ivalues.push_back(4);
	}

	opt.helixGeoFile = OP.getString("helixGeoFile");
	if (OP.fail()) {
		opt.warningFlag = true;
		opt.warningMessages += "helixGeoFile not specified, default to /data03/CATM/files/CENTROIDS_0.5_cutoff-6_byHbond";
		opt.helixGeoFile = "/data03/CATM/files/CENTROIDS_0.5_cutoff-6_byHbond";
	}

	opt.inputMonomerE = OP.getBool("inputMonomerE");
	if (OP.fail()) {
		opt.warningMessages += "monomer energy will be calculated\n";
		opt.warningFlag = true;
		opt.inputMonomerE = true;
	}
	opt.monoE_vdw = OP.getDouble("monoE_vdw");
	if (OP.fail()) {
		opt.monoE_vdw = 1000000; //Default large, easy to spot error.
	}
	opt.monoE_hbond = OP.getDouble("monoE_hbond");
	if (OP.fail()) {
		opt.monoE_hbond = 1000000; //Default large, easy to spot error.
	}
	opt.monoE_solv = OP.getDouble("monoE_solv");
	if (OP.fail()) {
		opt.monoE_solv = 1000000; //Default large, easy to spot error.
	}
	opt.monoE_solvRef = OP.getDouble("monoE_solvRef");
	if (OP.fail()) {
		opt.monoE_solvRef= 1000000; //Default large, easy to spot error.
	}

	opt.topFile = OP.getString("topFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_TOP";
		if(SYSENV.isDefined(envVar)) {
			opt.topFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "topFile not specified using " + opt.topFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine topFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.parFile = OP.getString("parFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.parFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "parFile not specified using " + opt.parFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine parFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.rotLibFile = OP.getString("rotLibFile");
	if (OP.fail()) {
		string envVar = "MSL_ROTLIB";
		if(SYSENV.isDefined(envVar)) {
			opt.rotLibFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "rotLibFile not specified using " + opt.rotLibFile + ", defaulting to " + SYSENV.getEnv(envVar) + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine rotLibFile - " + envVar + " - not set\n";
			opt.errorFlag = true;
		}
	}
	opt.monoRotLibFile = OP.getString("monoRotLibFile");
	if (OP.fail()) {
		opt.warningMessages += "monoRotLibFile not specified using " + opt.rotLibFile + "\n";
		opt.warningFlag = true;
		opt.monoRotLibFile = opt.rotLibFile;
	}

	opt.backboneCrd = OP.getString("backboneCrd");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine backboneCrd";
		opt.errorFlag = true;
	}
	
	opt.hBondFile = OP.getString("hbondFile");
	if (OP.fail()) {
		string envVar = "MSL_HBOND_CA_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.hBondFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "hbondFile not specified using " + opt.hBondFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine hbondFile - MSL_HBOND_CA_PAR - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.infile = OP.getString("infile");
	if (OP.fail()) { 
		opt.warningMessages += "infile not specified, default to /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdbi\n";
		opt.warningFlag = true;
		opt.infile = "/data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.pdb";
	}
	
	opt.uniprotName = OP.getString("uniprotName");
	if (OP.fail()) {
		opt.uniprotName = "PROTEIN_UNK";
		opt.warningMessages += "uniprotName not specified using " + opt.uniprotName + "\n";
		opt.warningFlag = true;
	}

	opt.uniprotAccession = OP.getString("uniprotAccession");
	if (OP.fail()) {
		opt.uniprotAccession = "P00000";
		opt.warningMessages += "uniprotAccession not specified using " + opt.uniprotAccession + "\n";
		opt.warningFlag = true;
	}

	opt.pdbOutputDir = OP.getString("pdbOutputDir");
	if (OP.fail()) {
		opt.errorMessages += "Unable to determine pdbOutputDir";
		opt.errorFlag = true;
	}

	opt.Ids = OP.getStringVector("Ids");
	if (OP.fail()) {
		opt.errorMessages += "Unable to identify alternate AA identities, make sure they are space separated\n";
		opt.errorFlag = true;
	}


	OptionParser OPD;
	OPD.readArgv(_argc, _argv);
	OPD.setDefaultArguments(data.defaultArgs);
	OPD.setRequired(data.required);
	OPD.setAllowed(data.allowed);
	OPD.autoExtendOptions();
	
	/*******************************************
	 *   DATA OPTIONS
	 *******************************************/

	if (OPD.countOptions() == 0){
		usage();
		cerr << "No options given!" << endl;
		exit(0);
	}
	opt.datafile = OPD.getString("datafile");
	if (opt.datafile != "") {
		string dataPath = "/data07/sanderson/CATM/runs/CATM_v24_Hsapiens/" + opt.datafile + "/";
		OPD.readFile(dataPath + opt.datafile + ".txt");
		if (OPD.fail()) {
			opt.warningFlag = true;
			opt.warningMessages += "Cannot read data file" + dataPath + opt.datafile + ".txt\n";
			OPD.readFile(dataPath + opt.datafile + "_01.txt");
			if (OPD.fail()) {
				opt.errorFlag = true;
				opt.errorMessages += "Cannot read data file " + dataPath + opt.datafile + "_01.txt\n";
				exit(1);
			}
		}
	}
	
	/*****************************************
	 *  CHECK THE NUMBER OF DATA OPTIONS
	 *****************************************/
	if (!OPD.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return data;
	}

	/*****************************************
	 *  CHECK THE DATA OPTIONS
	 *****************************************/
	opt.originalSeq = OPD.getString("originalSeq");
	if (OPD.fail()) {
		opt.warningMessages += "originalSeq not specified, defaulting to GGGGGGGGGGGGGGGGGGGGG";
		opt.warningFlag = true;
		opt.originalSeq = "GGGGGGGGGGGGGGGGGGGGG";
	}
	opt.modelledSeq = OPD.getString("modelledSeq");
	if (OPD.fail()) {
		opt.warningMessages += "modelledSeq not specified, defaulting to GGGGGGGGGGGGGGGGGGGGG";
		opt.warningFlag = true;
		opt.modelledSeq = "GGGGGGGGGGGGGGGGGGGGG";
	}

	opt.resStart = OPD.getInt("resStart");
	if (OPD.fail()) {
		opt.warningMessages += "resStart not specified, ...";
		opt.warningFlag = true;
		opt.resStart = 0;
	}
	opt.resEnd = OPD.getInt("resEnd");
	if (OPD.fail()) {
		opt.warningMessages += "resEnd not specified, ...";
		opt.warningFlag = true;
		opt.resStart = 0;
	}
	opt.seqLength = opt.resEnd - opt.resStart + 1;


	opt.thread = OPD.getInt("thread");
	if (OPD.fail()) {
		opt.warningMessages += "thread not specified, defaulting to 0\n";
		opt.warningFlag = true;
		opt.thread = 0;
	}
	opt.bbThread = opt.thread;// + opt.seqLength - opt.backboneLength;//might not need this +1 anymore because it has the right number for the sequence?
		
	opt.XShift = OPD.getDouble("XShift");
	if (OPD.fail()) {
		opt.warningMessages += "XShift not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.xShift = 6.7;
	}
	opt.ZShift = OPD.getDouble("ZShift");
	if (OPD.fail()) {
		opt.warningMessages += "ZShift not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.zShift = 0;
	}
	opt.axialRot = OPD.getDouble("axialRot");
	if (OPD.fail()) {
		opt.warningMessages += "axialRot not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.axialRotation = 0;
	}
	opt.crossAng = OPD.getDouble("crossAng");
	if (OPD.fail()) {
		opt.warningMessages += "crossingAngle not specified, defaulting to ...\n";
		opt.warningFlag = true;
		opt.crossingAngle = -40;
	}

	opt.energy = OPD.getDouble("energy");
	if (OPD.fail()) {

	}



	opt.rerunConf = OP.getConfFile();

	return opt;
	return data;
}
