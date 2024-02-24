#include <iostream>
#include "System.h"
#include "Helanal.h"
#include "Line.h"
#include "Transforms.h"
#include "OptionParser.h"


using namespace std;
using namespace MSL;

time_t start_time, end_time;
double diffTime;

string programName = "interhelicalCoordinates";
string programDescription = "This program examines the interhelical giometry of pairs of interaction helices in structures";
string programAuthor = "Alessandro Senes";
string programVersion = "0.0.00.1";
string programDate = "14 January 2020";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;


struct ICOptions {
	// Required
	string filename;

	// optional
	string PDBcode; // the pdb code if different from the file name

	// helix filtering parameters
	double backboneBondMin;
	uint minHelicalSegmentSize;

	// parameters for helix classification, H-bonding distance
	double minHBondDistanceStrict;
	double minHBondDistanceLoose;

	// parameters for selecting membrane residues only by their Z coordinate
	bool imposeZlimits;
	double minZ;
	double maxZ;

	// pair filtering parameters
	// interacting helices need to have a min number of close contact CA (using minCAdistance) 
	double minCAdistance;
	uint minNumOfContacts;
	uint excludeEndCAs; // contacts between the N- and C- terminal CAs do not count

	// for the alingment to the standard atoms
	uint numberOfCAinAlignment;

	string logFile;
	string outputdir;
	string rerunConfFile;

	//bool verbose;

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
};

struct ResidueHelanal {
	Helanal h;
	bool hasHelanal = false;
	bool isHelanalLooselyHelical = false;
	bool isHelanalStrictlyHelical = false;
};

struct ResiduePhiPsi {
	vector<double> angle = vector<double>(2, 0.0);
	vector<bool> hasPhiPsi = vector<bool>(2, false);
	bool isPhiPsiStrictlyHelical = false;
	bool isPhiPsiLooselyHelical = false;
};

struct ResidueHbonding {
	vector<double> distance = vector<double>(2, 0.0);
	vector<bool> hasHbonding = vector<bool>(2, false);
	bool isHbondingStrictlyHelical = false;
	bool isHbondingLooselyHelical = false;
};

struct ResidueProperties {
	Atom * pCA = NULL;
	Position * pPosition = NULL;
	bool isBondedPosition = false;
	bool hasCA = false;
	bool isLooselyHelical = false;
	bool isStrictlyHelical = false;
	bool isInMembrane = true;
	int belongsToSegment = -1;
	ResidueHelanal rHelanal;
	ResiduePhiPsi rPhiPsi;
	ResidueHbonding rHbond;
	uint chainIndex = 0;
	uint posIndexInChain = 0;
};

struct StandardInteractingHelices {
	AtomPointerVector standardHelixCAs1;
	AtomPointerVector standardHelixCAs2;
	double risePerResidue;
	double twist;
	double fitRMSD1;
	double fitRMSD2;
	Line axis1;
	Line axis2;
	uint axisIndex1;
	uint axisIndex2;
	CartesianPoint PoCA1;
	CartesianPoint PoCA2;
	double axialDistance;
	double segmentDistance;
	double angle;
	double Z1;
	double omega1;
	double Zp1;
	double omegaP1;
	double Z2;
	double omega2;
	double Zp2;
	double omegaP2;
	bool PoCA_is_inner1 = false;
	bool PoCA_is_inner2 = false;
	uint interhelicalQuadrant1;
	uint interhelicalQuadrant2;
	vector<Helanal> helanalsHelix1;
	vector<Helanal> helanalsHelix2;
	vector<bool> hasHelanal1;
	vector<bool> hasHelanal2;
};

struct InteractingHelicalSegments {
	uint segment1Index;
	uint segment2Index;
	vector<ResidueProperties> * pSegment1 = NULL;
	vector<ResidueProperties> * pSegment2 = NULL;
	Line axisOfSegment1;
	Line axisOfSegment2;
	uint axisIndex1;
	uint axisIndex2;
	double axialDistance;
	double segmentDistance;
	double angle;
	CartesianPoint PoCA1;
	CartesianPoint PoCA2;
	uint totalContacts;
	StandardInteractingHelices standardHelices;
	bool hasStandardHelices = false;
	string sequence1;
	string sequence2;
	string helicalString1;
	string helicalString2;
};

// pre-declaration of functions
void extractPositions(System & _sys, double & _backboneBondMin, vector<vector<ResidueProperties> > & _chainResidues); 
void computeHelanals(vector<vector<ResidueProperties> > & _chainResidues);
void computePhiPsis(vector<vector<ResidueProperties> > & _chainResidues);
void computeHbonding(vector<vector<ResidueProperties> > & _chainResidues);
void classifyHelicalHelanals(vector<vector<ResidueProperties> > & _chainResidues, double _minStrictHelanalHeight, double _maxStrictHelanalHeight, double _minStrictHelanalTwist, double _maxStrictHelanalTwist, double _minStrictHelanalRadius, double _maxStrictHelanalRadius, double _minLooseHelanalHeight, double _maxLooseHelanalHeight, double _minLooseHelanalTwist, double _maxLooseHelanalTwist, double _minLooseHelanalRadius, double _maxLooseHelanalRadius);
void findHelicalSegments(vector<vector<ResidueProperties> > & _chainResidues, uint & _minSegmentSize, vector<vector<ResidueProperties> > & _helicalSegmentResidues);
void printProteinReport(string _PDBcode, vector<vector<ResidueProperties> > & _chainResidues, stringstream & _ss);
void printSegmentListReport(string _PDBcode, vector<vector<ResidueProperties> > & _helicalSegmentResidues, stringstream & _ss);
void findInteractingSegments(vector<vector<ResidueProperties> > & _helicalSegmentResidues, double _minCAdistance, uint _minNumOfContacts, uint _excludeEndCAs, vector<InteractingHelicalSegments> & _interactingSegments);
void calculateInterhelicalGeometry(vector<InteractingHelicalSegments> & _interactingSegments);
void printInterhelicalListReport(string _PDBcode, vector<InteractingHelicalSegments> & _interactingSegments, stringstream & _ss);
bool writeSegmentPDBs(vector<vector<ResidueProperties> > & _helicalSegmentResidues, string _outputdir);
bool writePairsPDBs(vector<InteractingHelicalSegments> & _interactingSegments, string _outputdir);
bool writeStandardPairsPDBs(vector<InteractingHelicalSegments> & _interactingSegments, string _outputdir);
void deleteStandardAtomFits(vector<InteractingHelicalSegments> & _interactingSegments);
void deleteAtoms(AtomPointerVector & _standardAtoms);
void quadrantCoordinates(CartesianPoint & _pointOnOtherHelixAxis, Line & _helAxis, Atom * _refCAs, double _risePerResidue, double _twist, vector<double> & _coordinates, int & _quadrant);
void generateStandardHelicalPairs(vector<InteractingHelicalSegments> & _interactingSegments, uint _alignNumber);
void generateStandardAtoms(AtomPointerVector & _originalCAs, AtomPointerVector & _standardAtoms, uint _alignStart, uint _alignNumber, double & _risePerResidue, double & _twist, double & _RMSD);
void computeStandardHelixHelanals(AtomPointerVector & _standardHelixCAs, vector<Helanal>  &_helanals, vector<bool> & _hasHelanal);
void calculateStandardInterhelicalGeometry(vector<InteractingHelicalSegments> & _interactingSegments);
void calculateZ_OmegaCoordinates(CartesianPoint & _pointOfMinDThisHelix, CartesianPoint & _pointOfMinDOtherHelix, AtomPointerVector & _standardHelix, Line & _helanalAxis, uint & _helanalIndex, double _risePerResidue, double _twist, vector<double> & _coordinates, int & _quadrant, bool & _PoCA_is_inner);
void printInterhelicalStandardListReport(string _PDBcode, vector<InteractingHelicalSegments> & _interactingSegments, stringstream & _ss);
void classifyHelicalPhiPsi(vector<vector<ResidueProperties> > & _chainResidues);
void classifyHelicalHbonding(vector<vector<ResidueProperties> > & _chainResidues, double _minHBondDistanceStrict, double _minHBondDistanceLoose);
void classifyByZcoordinate(vector<vector<ResidueProperties> > & _chainResidues, double _minZ, double _maxZ);
void generateInterhelicalString(vector<InteractingHelicalSegments> & _interactingSegments);
bool writeReport(string _filename, string _outputdir, stringstream & _ss);

ICOptions parseICOptions(int _argc, char * _argv[], ICOptions defaults);
void usage();
void version();
void help(ICOptions defaults);

int main(int argc, char* argv[]) {

	time(&start_time);	

	/******************************************************************************
	 *                          === SETTINGS THE DEFAULTS ===
	 *
	 *  Put here the defaults for some options that are
	 *  not always required
	 ******************************************************************************/
	ICOptions defaults;
	// minimal distance for considering two backbone atoms bonded
	defaults.backboneBondMin = 2.0;
	// minimal lenght of a helical segment
	defaults.minHelicalSegmentSize = 13;
	// parameters for H-bonding distance
	defaults.minHBondDistanceStrict = 3.2;
	defaults.minHBondDistanceLoose = 4.5;

	// interacting helices need to have a min number of close contact CA (using minCAdistance) 
	defaults.minCAdistance = 9;
	defaults.minNumOfContacts = 3;
	defaults.excludeEndCAs = 3; // contacts between the N- and C- terminal CAs do not count

	// for the alingment to the standard atoms
	defaults.numberOfCAinAlignment = 8;

	// for imposing membrane residues only by Z coordinate
	defaults.imposeZlimits = false;
	defaults.minZ = -32.0;
	defaults.maxZ = 32.0;

	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
  	ICOptions opt = parseICOptions(argc, argv, defaults);
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

	if(opt.warningFlag) {
		cout << opt.warningMessages;
	}
	MslTools::createDir(opt.outputdir);
	/************************************************************
	 *   Output a configuration file that would rerun the 
	 *   program with the same options of the current run (this is
	 *   useful when mixing options between the configuration file and
	 *   the command line
	 ************************************************************/
	opt.rerunConfFile = opt.outputdir + "/rerun_conf.txt";
	ofstream rerunConf_fs;
	rerunConf_fs.open(opt.rerunConfFile.c_str());
	if (rerunConf_fs.fail()) {
		cerr << "ERROR: Cannot write " << opt.rerunConfFile << endl;
		exit(1);
	}
	rerunConf_fs << opt.rerunConf;
	rerunConf_fs.close();

//	ofstream fout;
//	fout.open(opt.logFile.c_str());
//	
//	if(!fout.is_open()) {
//		cerr << "Unable to open " << opt.logFile << endl;
//		exit(0);
//	}
//	
//	printICOptions(opt, fout);


	// defined the range of helanal parameters to consider a quadruplet
	// of CAs as helical
	// ideal: 1.5, 100, 2.27
	double minStrictHelanalHeight = 1.25;
	double maxStrictHelanalHeight = 1.75;
	double minStrictHelanalTwist = 90.0;
	double maxStrictHelanalTwist = 110.0;
	double minStrictHelanalRadius = 2.12;
	double maxStrictHelanalRadius = 2.42;
	double minLooseHelanalHeight = 1.25;
	double maxLooseHelanalHeight = 1.75;
	double minLooseHelanalTwist = 90.0;
	double maxLooseHelanalTwist = 110.0;
	double minLooseHelanalRadius = 2.12;
	double maxLooseHelanalRadius = 2.42;

	//System molecule(readHelix.getAtomPointers());
	System molecule;
	molecule.readPdb(argv[1]);

	/******************************************************
	 *  Extract the position in a structure ResidueProperties
	 *  that will later include information on helicity (helanal,
	 *  phi/psi and backbone hydrogen bonding
	 ******************************************************/
	vector<vector<ResidueProperties> > chainResidues;
	extractPositions(molecule, opt.backboneBondMin, chainResidues); 
	
	/******************************************************
	 * Let's now compute the helanals for each chain.
	 * If a position does not have a CA, it gets added an
	 * empty helanal.
	 * Same for Phi-Psi and backbone H-bonding at i-4 and i+4
	 ******************************************************/
	computeHelanals(chainResidues);
	computePhiPsis(chainResidues);
	computeHbonding(chainResidues);


	/******************************************************
	 * Classify the positions as strictly or loosely helical
	 * based on Phi-Psi and helanal (at the moment no filter will
	 * be later applied for helanal)
	 ******************************************************/
	classifyHelicalPhiPsi(chainResidues);
	classifyHelicalHbonding(chainResidues, opt.minHBondDistanceStrict, opt.minHBondDistanceLoose);
	classifyHelicalHelanals(chainResidues, minStrictHelanalHeight, maxStrictHelanalHeight, minStrictHelanalTwist, maxStrictHelanalTwist, minStrictHelanalRadius, maxStrictHelanalRadius, minLooseHelanalHeight, maxLooseHelanalHeight, minLooseHelanalTwist, maxLooseHelanalTwist, minLooseHelanalRadius, maxLooseHelanalRadius);
	if (opt.imposeZlimits) {
		cout << "Imposing Z limits of " << opt.minZ << " to " << opt.maxZ << endl;
		classifyByZcoordinate(chainResidues, opt.minZ, opt.maxZ);
	}

	
	/**********************************************************
	 * Now we identify the helical segments, stored in the following
	 * vector of ResidueProperties
	 **********************************************************/
	vector<vector<ResidueProperties> > helicalSegmentResidues;
	findHelicalSegments(chainResidues, opt.minHelicalSegmentSize, helicalSegmentResidues);


	/**********************************************************
	 * print a report for the protein and a report for the
	 * helical segment that were found
	 **********************************************************/
	//  * * * *  ADD PDB CODE TO ALL REPORTS  * * * *
	stringstream ss;
	printProteinReport(opt.PDBcode, chainResidues, ss);
	cout << "=====================================" << endl;
	cout << "PROTEIN HELANAL REPORT" << endl;
	cout << ss.str();
	writeReport("proteinReport", opt.outputdir, ss);

	printSegmentListReport(opt.PDBcode, helicalSegmentResidues, ss);
	cout << "=====================================" << endl;
	cout << "HELICAL SEGMENTS REPORT" << endl;
	cout << ss.str();
	writeReport("segmentReport", opt.outputdir, ss);


	/**********************************************************
	 * Now we identify which helical segments contact each other
	 * by looking at CA-CA distances.
	 **********************************************************/
	vector<InteractingHelicalSegments> interactingSegments;
	findInteractingSegments(helicalSegmentResidues, opt.minCAdistance, opt.minNumOfContacts, opt.excludeEndCAs, interactingSegments);
	
	/**********************************************************
	 * Now we calculate the inter-helical geometry of the
	 * interacting segments using their helanals.
	 **********************************************************/
	calculateInterhelicalGeometry(interactingSegments);
	printInterhelicalListReport(opt.PDBcode, interactingSegments, ss);
	cout << "=====================================" << endl;
	cout << "INTERACTING HELICAL SEGMENTS REPORT" << endl;
	cout << ss.str();
	writeReport("pairReport", opt.outputdir, ss);

	/**********************************************************
	 * Now we fit the helices to standard helices.
	 * The RMSD alignment will be peformed with a subset of CA
	 * atoms around the segment of closest approac, vector<bool>()h
	 * Compute also the helanals of the standard helices.
	 **********************************************************/
	vector<vector<vector<Helanal> > >  standardHelixHelanals;
	vector<vector<vector<bool> > > standardHelixHasHelanal; // all quadruplet will have helanals but the last 3 atoms do not because the helix ends
	generateStandardHelicalPairs(interactingSegments, opt.numberOfCAinAlignment);

	/**********************************************************
	 * Now we measure the interhelical geometry of the standard helices.
	 * We will calcuate the interhelical coordinates based on the geometry
	 * of these standard helices
	 *
	 **********************************************************/
	calculateStandardInterhelicalGeometry(interactingSegments);
	generateInterhelicalString(interactingSegments);
	printInterhelicalStandardListReport(opt.PDBcode, interactingSegments, ss);
	cout << "=====================================" << endl;
	cout << "INTERACTING HELICAL SEGMENTS STANDARD GEOMETRY REPORT" << endl;
	cout << ss.str();
	writeReport("pairGeometryReport", opt.outputdir, ss);


	/**********************************************************
	 * Write out the pdbs
	 **********************************************************/
	writeSegmentPDBs(helicalSegmentResidues, opt.outputdir);
	writePairsPDBs(interactingSegments, opt.outputdir);
	writeStandardPairsPDBs(interactingSegments, opt.outputdir);

	// garbage collection for the standard atoms
	deleteStandardAtomFits(interactingSegments);

	return 0;
}

void extractPositions(System & _sys, double & _backboneBondMin, vector<vector<ResidueProperties> > & _chainResidues) {
	/****************************************************
	 * Return variables:
	 *  - _chainResidues structure
	 ****************************************************/
	_chainResidues.clear();	
	//_chainCAs.clear();
	//_chainPositions.clear();
	//_bondedPositions.clear();
	//_positionHasCA.clear();

	for (uint i=0; i<_sys.chainSize(); i++) {
		Chain & chain = _sys.getChain(i);
	//	cout << "Chain " << chain.getChainId() << " has " << chain.positionSize() << " positions" << endl;
		//_chainCAs.push_back(AtomPointerVector());
		//_chainPositions.push_back(vector<Position*>());
		//_bondedPositions.push_back(vector<bool>());
		//_positionHasCA.push_back(vector<bool>());

		if (chain.positionSize() == 0) {
			continue;
		}
		vector<Position*> chainPositions = chain.getPositions();
		_chainResidues.push_back(vector<ResidueProperties>(chainPositions.size(), ResidueProperties()));
		for (uint j=0; j<chainPositions.size(); j++) {
			//_chainPositions[i].push_back(chainPositions[j]); // we should remove the local chainPositions and use direcly _chainPositions
			_chainResidues[i][j].pPosition = chainPositions[j];
			_chainResidues[i][j].chainIndex = i;
			_chainResidues[i][j].posIndexInChain = j;

			// check if we have a CA
			if (chainPositions[j]->atomExists("CA")) {
				//cout << "UUU get CA" << endl;
				Atom & CA = chainPositions[j]->getLastFoundAtom();
				//_chainCAs[i].push_back(&CA);
				//_positionHasCA[i].push_back(true);
				_chainResidues[i][j].pCA = &CA;
				_chainResidues[i][j].hasCA = true;

				// check if the bonding of the backbone is OK for this position
				if (j<chainPositions.size()-1) {
					if (chainPositions[j]->atomExists("N") && chainPositions[j]->atomExists("C") && chainPositions[j+1]->atomExists("N")) {
						//cout << "UUU get N not end" << endl;
						Atom & N = chainPositions[j]->getAtom("N");
						//cout << "UUU get C not end" << endl;
						Atom & C = chainPositions[j]->getAtom("C");
						//cout << "UUU get Nplus not end" << endl;
						Atom & Nplus = chainPositions[j+1]->getAtom("N");
						if (N.distance(CA) < _backboneBondMin && CA.distance(C) < _backboneBondMin && C.distance(Nplus) < _backboneBondMin) {
							_chainResidues[i][j].isBondedPosition = true;
						//	cout << j << " Bonded " << _chainResidues[i][j].isBondedPosition << endl;
						} else {
							_chainResidues[i][j].isBondedPosition = false;
						//	cout << j << " Non bonded " << _chainResidues[i][j].isBondedPosition << endl;
							//cout << "Found break of the backbone at position " << _chainResidues[i][j].pPosition->getPositionId() << endl;
						}
					} else {
						_chainResidues[i][j].isBondedPosition = false;
					}
				} else {
					// check is different if this is the last position, do not check the N+1
					if (chainPositions[j]->atomExists("N") && chainPositions[j]->atomExists("CA") && chainPositions[j]->atomExists("C")) {
						//cout << "UUU get N end" << endl;
						Atom & N = chainPositions[j]->getAtom("N");
						//cout << "UUU get C end" << endl;
						Atom & C = chainPositions[j]->getAtom("C");
						if (N.distance(CA) < _backboneBondMin && CA.distance(C) < _backboneBondMin) {
							_chainResidues[i][j].isBondedPosition = true;
						//	cout << j << " Bonded " << _chainResidues[i][j].isBondedPosition << endl;
						} else {
							_chainResidues[i][j].isBondedPosition = false;
						//	cout << j << " Non bonded " << _chainResidues[i][j].isBondedPosition << endl;
							//cout << "Found break of the backbone at position " << _chainResidues[i][j].pPosition->getPositionId() << endl;
						}
					} else {
						_chainResidues[i][j].isBondedPosition = false;
					}
				}
			}
		}
	}
}

void computeHelanals(vector<vector<ResidueProperties> > & _chainResidues) {
	/****************************************************
	 ****************************************************/

	for (uint i=0; i<_chainResidues.size(); i++) {
		for (uint j=0; j<_chainResidues[i].size(); j++) {
			if (j>=_chainResidues[i].size()-3) {
				_chainResidues[i][j].rHelanal.hasHelanal = false;
				continue;
			}
			// check that the four positions are properly bonded and that the last one has a CA (proper bonding assumes a CA for positions 1, 2 and 3)
			if (_chainResidues[i][j].isBondedPosition && _chainResidues[i][j+1].isBondedPosition && _chainResidues[i][j+2].isBondedPosition && _chainResidues[i][j+3].hasCA) {
				//Helanal	h(_chainCAs[i][j]->getCoor(), _chainCAs[i][j+1]->getCoor(), _chainCAs[i][j+2]->getCoor(), _chainCAs[i][j+3]->getCoor());   
				//_helanals[i].push_back(h);
				//_hasHelanal[i].push_back(true);
				_chainResidues[i][j].rHelanal.h = Helanal(_chainResidues[i][j].pCA->getCoor(), _chainResidues[i][j+1].pCA->getCoor(), _chainResidues[i][j+2].pCA->getCoor(), _chainResidues[i][j+3].pCA->getCoor());
				_chainResidues[i][j].rHelanal.hasHelanal = true;
			} else {
				// since the positions are not bonded, create an empty Helanal and flag the position has invalid
				//_helanals[i].push_back(Helanal());
				//_hasHelanal[i].push_back(false);
				_chainResidues[i][j].rHelanal.hasHelanal = false;
			}	
		}
	}
}

void computePhiPsis(vector<vector<ResidueProperties> > & _chainResidues) {
	/****************************************************
	 *  This function adds the phi-psi value to each
	 *  position of the chains
	 ****************************************************/
	for (uint i=0; i<_chainResidues.size(); i++) {
		for (uint j=0; j<_chainResidues[i].size(); j++) {
			bool checkPhi = true;
			bool checkPsi = true;

			// no phi if at the beginning of the chain or after a break...
			if (j==0) {
				checkPhi = false;
			} else {
				if (!_chainResidues[i][j-1].isBondedPosition) {
					checkPhi = false;
				}
			}
			// ... and no psi if it is the end of the chain or before a break
			if (j==_chainResidues[i].size()-1 || !_chainResidues[i][j].isBondedPosition) {
				checkPsi = false;
			}
			if (!checkPhi and !checkPsi) {
				// nothing to do
				continue;
			}

			// get the backbone atoms
			bool foundPhiAtoms = true;
			bool foundPsiAtoms = true;
			Atom * pCminus = NULL;
			Atom * pN = NULL;
			Atom * pCA = NULL;
			Atom * pC = NULL;
			Atom * pNplus = NULL;
			if (checkPhi) {
				if (_chainResidues[i][j-1].pPosition->atomExists("C")) {
					pCminus = &_chainResidues[i][j-1].pPosition->getLastFoundAtom();
				} else {
					foundPhiAtoms = false;
				}
			} else {
				foundPhiAtoms = false;
			}
			if (_chainResidues[i][j].pPosition->atomExists("N")) {
				pN = &_chainResidues[i][j].pPosition->getLastFoundAtom();
			} else {
				foundPhiAtoms = false;
				foundPsiAtoms = false;
			}
			if (_chainResidues[i][j].pPosition->atomExists("CA")) {
				pCA = &_chainResidues[i][j].pPosition->getLastFoundAtom();
			} else {
				foundPhiAtoms = false;
				foundPsiAtoms = false;
			}
			if (_chainResidues[i][j].pPosition->atomExists("C")) {
				pC = &_chainResidues[i][j].pPosition->getLastFoundAtom();
			} else {
				foundPhiAtoms = false;
				foundPsiAtoms = false;
			}
			if (checkPsi) {
				if (_chainResidues[i][j+1].pPosition->atomExists("N")) {
					pNplus = &_chainResidues[i][j+1].pPosition->getLastFoundAtom();
				} else {
					foundPsiAtoms = false;
				}
			} else {
				foundPsiAtoms = false;
			}
			if (foundPhiAtoms) {
				_chainResidues[i][j].rPhiPsi.angle[0] = pCminus->dihedral(*pN, *pCA, *pC);
				_chainResidues[i][j].rPhiPsi.hasPhiPsi[0] = true;
			}
			if (foundPsiAtoms) {
				_chainResidues[i][j].rPhiPsi.angle[1] = pN->dihedral(*pCA, *pC, *pNplus);
				_chainResidues[i][j].rPhiPsi.hasPhiPsi[1] = true;
			}
		}
	}
}

void computeHbonding(vector<vector<ResidueProperties> > & _chainResidues) {
	/****************************************************
	 * This function computes the backbone hydrogen bonding
	 * at i, i-4 and i, i+4
	 ****************************************************/
	for (uint i=0; i<_chainResidues.size(); i++) {
		for (uint j=0; j<_chainResidues[i].size(); j++) {
			bool checkNbond = true;
			bool checkCbond = true;

			// no N-H H-bonding if at the beginning of the chain or after a break...
			if (j<4) {
				checkNbond = false;
			} else {
				if (!_chainResidues[i][j-4].isBondedPosition || !_chainResidues[i][j-3].isBondedPosition || !_chainResidues[i][j-2].isBondedPosition || !_chainResidues[i][j-1].isBondedPosition) {
					checkNbond = false;
				}
			}
			// ... and no C=O H-bondinh if it is the end of the chain or before a break
			if (j>_chainResidues[i].size()-5 || !_chainResidues[i][j].isBondedPosition || !_chainResidues[i][j+1].isBondedPosition || !_chainResidues[i][j+2].isBondedPosition || !_chainResidues[i][j+3].isBondedPosition) {
				checkCbond = false;
			}
			if (!checkNbond and !checkCbond) {
				// nothing to do
				continue;
			}

			// get the backbone atoms
			bool foundNbondAtoms = true;
			bool foundCbondAtoms = true;
			Atom * pOminus4 = NULL;
			Atom * pN = NULL;
			Atom * pO = NULL;
			Atom * pNplus4 = NULL;
			if (checkNbond) {
				if (_chainResidues[i][j-4].pPosition->atomExists("O")) {
					pOminus4 = &_chainResidues[i][j-4].pPosition->getLastFoundAtom();
				} else {
					foundNbondAtoms = false;
				}
				if (_chainResidues[i][j].pPosition->atomExists("N")) {
					pN = &_chainResidues[i][j].pPosition->getLastFoundAtom();
				} else {
					foundNbondAtoms = false;
				}
			} else {
				foundNbondAtoms = false;
			}
			if (checkCbond) {
				if (_chainResidues[i][j+4].pPosition->atomExists("N")) {
					pNplus4 = &_chainResidues[i][j+4].pPosition->getLastFoundAtom();
				} else {
					foundCbondAtoms = false;
				}
				if (_chainResidues[i][j].pPosition->atomExists("O")) {
					pO = &_chainResidues[i][j].pPosition->getLastFoundAtom();
				} else {
					foundCbondAtoms = false;
				}
			} else {
				foundCbondAtoms = false;
			}
			if (foundNbondAtoms) {
				_chainResidues[i][j].rHbond.distance[0] = pOminus4->distance(*pN);
				_chainResidues[i][j].rHbond.hasHbonding[0] = true;
			}
			if (foundCbondAtoms) {
				_chainResidues[i][j].rHbond.distance[1] = pNplus4->distance(*pO);
				_chainResidues[i][j].rHbond.hasHbonding[1] = true;
			}
		}
	}
}

void classifyHelicalHelanals(vector<vector<ResidueProperties> > & _chainResidues, double _minStrictHelanalHeight, double _maxStrictHelanalHeight, double _minStrictHelanalTwist, double _maxStrictHelanalTwist, double _minStrictHelanalRadius, double _maxStrictHelanalRadius, double _minLooseHelanalHeight, double _maxLooseHelanalHeight, double _minLooseHelanalTwist, double _maxLooseHelanalTwist, double _minLooseHelanalRadius, double _maxLooseHelanalRadius) {
	/****************************************************
	 * This function flags the helanals that are helical
	 * depending on the range of their height, twist and radius.
	 ****************************************************/
	for (uint i=0; i<_chainResidues.size(); i++) {
		for (uint j=0; j<_chainResidues[i].size(); j++) {
			if (_chainResidues[i][j].rHelanal.hasHelanal) {
				if (_chainResidues[i][j].rHelanal.h.getHeight() >= _minStrictHelanalHeight && _chainResidues[i][j].rHelanal.h.getHeight() <= _maxStrictHelanalHeight && _chainResidues[i][j].rHelanal.h.getTwist() >= _minStrictHelanalTwist && _chainResidues[i][j].rHelanal.h.getTwist() <= _maxStrictHelanalTwist && _chainResidues[i][j].rHelanal.h.getRadius() >= _minStrictHelanalRadius && _chainResidues[i][j].rHelanal.h.getRadius() <= _maxStrictHelanalRadius) {
					_chainResidues[i][j].rHelanal.isHelanalStrictlyHelical = true;
					_chainResidues[i][j].rHelanal.isHelanalLooselyHelical = true;
				} else {
					if (_chainResidues[i][j].rHelanal.h.getHeight() >= _minLooseHelanalHeight && _chainResidues[i][j].rHelanal.h.getHeight() <= _maxLooseHelanalHeight && _chainResidues[i][j].rHelanal.h.getTwist() >= _minLooseHelanalTwist && _chainResidues[i][j].rHelanal.h.getTwist() <= _maxLooseHelanalTwist && _chainResidues[i][j].rHelanal.h.getRadius() >= _minLooseHelanalRadius && _chainResidues[i][j].rHelanal.h.getRadius() <= _maxLooseHelanalRadius) {
						_chainResidues[i][j].rHelanal.isHelanalStrictlyHelical = false;
						_chainResidues[i][j].rHelanal.isHelanalLooselyHelical = true;
					} else {
						_chainResidues[i][j].rHelanal.isHelanalStrictlyHelical = false;
						_chainResidues[i][j].rHelanal.isHelanalLooselyHelical = false;
					}
				}
			} else {
				// if we do not have a helanal for this quadruplet, it is obviously false
				_chainResidues[i][j].rHelanal.isHelanalStrictlyHelical = false;
				_chainResidues[i][j].rHelanal.isHelanalLooselyHelical = false;
			}
		}
	}

}

/*
void findHelicalSegments(vector<vector<ResidueProperties> > & _chainResidues, uint & _minSegmentSize, vector<vector<ResidueProperties> > & _helicalSegmentResidues) {
	/ ****************************************************
	 * This function identify the start and end of each
	 * helical segment (a series of consecutive helical
	 * helanals) in each chain. To be valid, the segments
	 * need to have a minimal lenght.
	 * Return variables:
	 *  - _segmentList: a vector of vector of int
	 *                  dimension 1: each segment
	 *                  dimension 2: index of chain, index 
	 *                               of start position, index 
	 *                               of end position
	 *  - _helicalSegmentResidues: 
	 **************************************************** /
	
	//_segmentList.clear();
	_helicalSegmentResidues.clear();

	// remove or downgrade any position that is not helical according to phi-phi
	for (uint i=0; i<_chainResidues.size(); i++) {
		for (uint j=0; j<_chainResidues[i].size(); j++) {
			if (_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical) {
				_chainResidues[i][j].isLooselyHelical = true;
				_chainResidues[i][j].isStrictlyHelical = true;
			} else {
				if (_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical) {
					_chainResidues[i][j].isLooselyHelical = true;
				}
			}
			if (!_chainResidues[i][j].isLooselyHelical) {
				continue;
			}
			// remove or downgrade any position that is not helical according to H bonding
			if (!_chainResidues[i][j].rHbond.isHbondingLooselyHelical) {
				_chainResidues[i][j].isLooselyHelical = false;
				_chainResidues[i][j].isStrictlyHelical = false;
			} else {
				if (!_chainResidues[i][j].rHbond.isHbondingStrictlyHelical) {
					_chainResidues[i][j].isStrictlyHelical = false;
				}
			}
			// remove by Z limits
			if (!_chainResidues[i][j].isInMembrane) {
				_chainResidues[i][j].isLooselyHelical = false;
				_chainResidues[i][j].isStrictlyHelical = false;
			}
		}
	}

	// WE NEED A SECTION FOR FLAGGING NON HELICAL HELANALS
					
	// lets create the helical segments
	for (uint i=0; i<_chainResidues.size(); i++) {
	//	cout << "UUU adding segments for chain " << i << endl;
		bool helixFound = false;
		uint segStart = 0;
		uint segEnd = 0;
		uint segLength = 0;
		bool addHelix = false;
		for (uint j=0; j<_chainResidues[i].size(); j++) {
		//	cout << "   UUU scanning chain " << i << ", position" << j <<  endl;
			if (!helixFound) {
				if (_chainResidues[i][j].isLooselyHelical) {
					// the start of a helix
					helixFound = true;
					segStart = j;
					segEnd = j;
		//			cout << "UUU found helix start at " << segStart << endl;
				}
			} else {
				// we are within a helix
				if (_chainResidues[i][j].isLooselyHelical) {
					// increment the helix
					segEnd = j;
					segLength = segEnd - segStart + 1;
		//			cout << "UUU increasing helix lenght, start " << segStart << ", end " << segEnd << endl;
					if (segLength >= _minSegmentSize) {
						// if the helix is long enough, set to add it to the list
				//		cout << "UUUU set to add helix for start " << segStart << ", end " << segEnd << endl;
						addHelix = true;
					}
				} else {
					// we found the end of the helix
					helixFound = false;
				}
			}
			if ((!helixFound || j==_chainResidues.size()-1) && addHelix) {
		//		cout << "UUUU adding segment: chain " << i << ", start " << segStart << ", end " << segEnd << endl;
				// if we found the end of the helix or we are at the end of the chain
				// add the helix, if the addHelix flag is on
				//_segmentList.push_back(vector<uint>(3, 0));
				//_segmentList.back()[0] = i; // the index of the chain;
				//_segmentList.back()[1] = segStart; // the index of the start position;
				//_segmentList.back()[2] = segEnd; // the index of the end position;
				_helicalSegmentResidues.push_back(vector<ResidueProperties>());
				// let's record that these position belong to the last segment
				for (uint k=segStart; k <= segEnd; k++) {
					//_chainResidues[i][k].belongsToSegment = _segmentList.size()-1;
					_chainResidues[i][k].belongsToSegment = _helicalSegmentResidues.size()-1;
					_helicalSegmentResidues.back().push_back(_chainResidues[i][k]);
					if (k>segEnd-3) {
						// mark the last 3 helanal of the helix as non existing
						_helicalSegmentResidues.back().back().rHelanal.hasHelanal = false;
					}
				}
				j = segEnd; // move forward the position past the last one we included (because of the +3 addition)
				segStart = 0;
				segEnd = 0;
				addHelix = false;
				helixFound = false;
			}
		}
	}
}
*/

void findHelicalSegments(vector<vector<ResidueProperties> > & _chainResidues, uint & _minSegmentSize, vector<vector<ResidueProperties> > & _helicalSegmentResidues) {
	// new version that only considers the i,i-4 H-bond

	/****************************************************
	 * This function identify the start and end of each
	 * helical segment (a series of consecutive helical
	 * helanals) in each chain. To be valid, the segments
	 * need to have a minimal lenght.
	 * Return variables:
	 *  - _segmentList: a vector of vector of int
	 *                  dimension 1: each segment
	 *                  dimension 2: index of chain, index 
	 *                               of start position, index 
	 *                               of end position
	 *  - _helicalSegmentResidues: 
	 ****************************************************/
	
	//_segmentList.clear();
	_helicalSegmentResidues.clear();

//	// remove or downgrade any position that is not helical according to phi-phi
//	for (uint i=0; i<_chainResidues.size(); i++) {
//		for (uint j=0; j<_chainResidues[i].size(); j++) {
//			if (_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical) {
//				_chainResidues[i][j].isLooselyHelical = true;
//				_chainResidues[i][j].isStrictlyHelical = true;
//			} else {
//				if (_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical) {
//					_chainResidues[i][j].isLooselyHelical = true;
//				}
//			}
//			if (!_chainResidues[i][j].isLooselyHelical) {
//				continue;
//			}
//			// remove or downgrade any position that is not helical according to H bonding
//			if (!_chainResidues[i][j].rHbond.isHbondingLooselyHelical) {
//				_chainResidues[i][j].isLooselyHelical = false;
//				_chainResidues[i][j].isStrictlyHelical = false;
//			} else {
//				if (!_chainResidues[i][j].rHbond.isHbondingStrictlyHelical) {
//					_chainResidues[i][j].isStrictlyHelical = false;
//				}
//			}
//			// remove by Z limits
//			if (!_chainResidues[i][j].isInMembrane) {
//				_chainResidues[i][j].isLooselyHelical = false;
//				_chainResidues[i][j].isStrictlyHelical = false;
//			}
//		}
//	}

	// WE NEED A SECTION FOR FLAGGING NON HELICAL HELANALS
					
	// lets create the helical segments
	for (uint i=0; i<_chainResidues.size(); i++) {
		cout << "UUU adding segments for chain " << i << endl;
		bool helixFound = false;
		uint segStart = 0;
		uint segEnd = 0;
		uint segLength = 0;
		bool addHelix = false;
		for (uint j=0; j<_chainResidues[i].size(); j++) {
			// select residues that are helical based on phi-psi
			if (_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical) {
				_chainResidues[i][j].isLooselyHelical = true;
				_chainResidues[i][j].isStrictlyHelical = true;
			} else {
				if (_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical) {
					_chainResidues[i][j].isLooselyHelical = true;
				}
			}
			// remove by Z limits
			if (!_chainResidues[i][j].isInMembrane) {
				_chainResidues[i][j].isLooselyHelical = false;
				_chainResidues[i][j].isStrictlyHelical = false;
			}
			// check for H-bonding only if we are at the fifth helical residue
			// since the first four won't have a partner
			if (helixFound && segLength > 4) {
				if (!_chainResidues[i][j].rHbond.isHbondingLooselyHelical) {
					_chainResidues[i][j].isLooselyHelical = false;
					_chainResidues[i][j].isStrictlyHelical = false;
				} else {
					if (!_chainResidues[i][j].rHbond.isHbondingStrictlyHelical) {
						_chainResidues[i][j].isStrictlyHelical = false;
					}
				}
			}

			cout << "   UUU scanning chain " << i << ", position" << j <<  endl;
			if (!helixFound) {
				if (_chainResidues[i][j].isLooselyHelical) {
					// the start of a helix
					helixFound = true;
					segStart = j;
					segEnd = j;
					segLength = 1;
					cout << "UUU found helix start at " << segStart << endl;
				}
				// if this residue is not bonded to the next, break the helix
				if (!_chainResidues[i][j].isBondedPosition) { 
					helixFound = false;
				}
			} else {
				// we are within a helix
				if (_chainResidues[i][j].isLooselyHelical) {
					// increment the helix
					segEnd = j;
					segLength = segEnd - segStart + 1;
					cout << "UUU increasing helix lenght, start " << segStart << ", end " << segEnd << endl;
					if (segLength >= _minSegmentSize) {
						// if the helix is long enough, set to add it to the list
						addHelix = true;
					}
					// if this residue is not bonded to the next, break the helix
					if (!_chainResidues[i][j].isBondedPosition) { 
						cout << "UUU breaking helix at " << j << endl;
						helixFound = false;
					}
				} else {
					// we found the end of the helix
					helixFound = false;
				}
			}
			// the below was previously j==_chainResidues.size()-1, which led to it not finding some chains. It works and finds chains for designed sequences now
			// the axialRotation prime is a bit off, but I'm going to try to run a test with it and see if it matters too much 
			if ((!helixFound || j==_chainResidues[i].size()-1) && addHelix) { 
				cout << "UUUU adding segment: chain " << i << ", start " << segStart << ", end " << segEnd << endl;
				// if we found the end of the helix or we are at the end of the chain
				// add the helix, if the addHelix flag is on
				//_segmentList.push_back(vector<uint>(3, 0));
				//_segmentList.back()[0] = i; // the index of the chain;
				//_segmentList.back()[1] = segStart; // the index of the start position;
				//_segmentList.back()[2] = segEnd; // the index of the end position;
				_helicalSegmentResidues.push_back(vector<ResidueProperties>());
				// let's record that these position belong to the last segment
				for (uint k=segStart; k <= segEnd; k++) {
					//_chainResidues[i][k].belongsToSegment = _segmentList.size()-1;
					_chainResidues[i][k].belongsToSegment = _helicalSegmentResidues.size()-1;
					_helicalSegmentResidues.back().push_back(_chainResidues[i][k]);
					if (k>segEnd-3) {
						// mark the last 3 helanal of the helix as non existing
						_helicalSegmentResidues.back().back().rHelanal.hasHelanal = false;
					}
				}
				j = segEnd; // move forward the position past the last one we included (because of the +3 addition)
				segStart = 0;
				segEnd = 0;
				addHelix = false;
				helixFound = false;
			}
		}
	}
}

void printProteinReport(string _PDBcode, vector<vector<ResidueProperties> > & _chainResidues, stringstream & _ss) {
	/*********************************************************************
	 * Here we print all the position with their various associated values
	 * in a tabbed CSV format.
	 * Field:
	 *                * * * U P D A T E    L I S T * * *
	 *  0 Chain index
	 *  1 Position index
	 *  2 Residue ID
	 *  3 Flag bonded positions (0 means a backbone break)
	 *  4 Flag positions with a CA atom (0 means no CA)
	 *  5 Flag positions with valid helanal (0 means no helanal)
	 *  6 Helanal height value (blank if no helanal)
	 *  7 Helanal twist value (blank if no helanal)
	 *  8 Helanal radius value (blank if no helanal)
	 *  9 Flag helanal within strict helical range (0 means it is not helical)
	 * 10 Flag helanal within loose helical range (0 means it is not helical)
	 * 11 Flag position with valid phi-psi
	 * 12 Phi
	 * 13 Psi
	 * 14 Index of the helical segment this position belongs to (blank if it does not belong to a segment)
	 *********************************************************************/

	 _ss.str(string());
	
	//_ss << "PDB Id\tPos Id\tChain #\tPos #\tIs bonded?\tHas CA?\tValid helanal?\tHel Height\tHel Twist\tHel Radius\tHelanal strictly helical?\tHelanal loosely helical?\tHas phi?\tHas psi?\tPhi\tPsi\tPhi-Psi strictly helical?\tPhi-Psi loosely helical?\tHas i-4 h-bonding?\tHas i+4 h-bonding?\tH-bond d i-4\tH-bond d i+4\tHbonding strictly helical?\tHbonding loosely helical?\tHelical segment #" << endl;
	_ss << "PDB Id\tPos Id\tChain #\tPos #\tIs bonded?\tHas CA?\tValid helanal?\tHel Height\tHel Twist\tHel Radius\tHelanal strictly helical?\tHelanal loosely helical?\tHas phi?\tHas psi?\tPhi\tPsi\tPhi-Psi strictly helical?\tPhi-Psi loosely helical?\tHas i-4 h-bonding?\tH-bond d i-4\tHbonding strictly helical?\tHbonding loosely helical?\tHelical segment #" << endl;

	for (uint i=0; i<_chainResidues.size(); i++) {
		for (uint j=0; j<_chainResidues[i].size(); j++) {
			_ss << _PDBcode << "\t" << _chainResidues[i][j].pPosition->getCurrentIdentity().getIdentityId() << "\t" << i << "\t" << j << "\t" << _chainResidues[i][j].isBondedPosition << "\t" << _chainResidues[i][j].hasCA << "\t" << _chainResidues[i][j].rHelanal.hasHelanal <<  "\t";
			if (_chainResidues[i][j].rHelanal.hasHelanal) {
				_ss << _chainResidues[i][j].rHelanal.h.getHeight() << "\t" << _chainResidues[i][j].rHelanal.h.getTwist() << "\t" << _chainResidues[i][j].rHelanal.h.getRadius() << "\t";
			} else {
				_ss << "\t\t\t";
			}
			_ss << _chainResidues[i][j].rHelanal.isHelanalStrictlyHelical << "\t" << _chainResidues[i][j].rHelanal.isHelanalLooselyHelical << "\t" << _chainResidues[i][j].rPhiPsi.hasPhiPsi[0] << "\t" << _chainResidues[i][j].rPhiPsi.hasPhiPsi[1] << "\t";
			if(_chainResidues[i][j].rPhiPsi.hasPhiPsi[0]) {
				_ss << _chainResidues[i][j].rPhiPsi.angle[0] << "\t";
			} else {
				_ss << "\t";
			}
			if(_chainResidues[i][j].rPhiPsi.hasPhiPsi[1]) {
				_ss << _chainResidues[i][j].rPhiPsi.angle[1] << "\t";
			} else {
				_ss << "\t";
			}
			_ss << _chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical << "\t" << _chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical << "\t" << _chainResidues[i][j].rHbond.hasHbonding[0] << "\t";
			if(_chainResidues[i][j].rHbond.hasHbonding[0]) {
				_ss << _chainResidues[i][j].rHbond.distance[0] << "\t";
			} else {
				_ss << "\t";
			}
			_ss << _chainResidues[i][j].rHbond.isHbondingStrictlyHelical << "\t" << _chainResidues[i][j].rHbond.isHbondingLooselyHelical << "\t";
			if (_chainResidues[i][j].belongsToSegment > -1) {
				_ss << _chainResidues[i][j].belongsToSegment << endl;
			} else {
				_ss << endl;
			}
		}
	}
}

void printSegmentListReport(string _PDBcode, vector<vector<ResidueProperties> > & _helicalSegmentResidues, stringstream & _ss) {
	 _ss.str(string());
	//_ss << "Segment #\tPos # in segment\tPos Id\tChain #\tPos#\tHel Height\tHel Twist\tHel Radius\tValid helanal?\tHel Height\tHel Twist\tHel Radius\tHelanal strictly helical?\tHelanal loosely helical?\tHas phi?\tHas psi?\tPhi\tPsi\tPhi-Psi strictly helical?\tPhi-Psi loosely helical?\tHas i-4 h-bonding?\tHas i+4 h-bonding?\tH-bond d i-4\tH-bond d i+4\tHbonding strictly helical?\tHbonding loosely helical?" << endl;
	_ss << "PDB Id\tSegment #\tPos Id start\tPos Id end\tChain #\tPos start #\tPos end #\tLength" << endl;
	for (uint i=0; i<_helicalSegmentResidues.size(); i++) {
		_ss << _PDBcode << "\t" << i << "\t" << _helicalSegmentResidues[i][0].pPosition->getCurrentIdentity().getIdentityId() << "\t" << _helicalSegmentResidues[i].back().pPosition->getCurrentIdentity().getIdentityId() << "\t" << _helicalSegmentResidues[i][0].chainIndex << "\t" << _helicalSegmentResidues[i][0].posIndexInChain  << "\t" << _helicalSegmentResidues[i].back().posIndexInChain << "\t" << _helicalSegmentResidues[i].back().posIndexInChain - _helicalSegmentResidues[i][0].posIndexInChain + 1 << endl;
	}
}

void findInteractingSegments(vector<vector<ResidueProperties> > & _helicalSegmentResidues, double _minCAdistance, uint _minNumOfContacts, uint _excludeEndCAs, vector<InteractingHelicalSegments> & _interactingSegments) {
	/*****************************************************************
	 * find segments with a minimum of short distance CA contacts
	 *  - _minNumOfContacts   : the minimal number of CA-CA contacts between
	 *                          segments
	 *  - _minCAdistance      : the distance under which we define a contact
	 *  - _excludeEndCAs      : do not consider the first or last N CA in a
	 *****************************************************************/
	_interactingSegments.clear();

	for (uint i=0; i<_helicalSegmentResidues.size(); i++) {
		for (uint j=i+1; j<_helicalSegmentResidues.size(); j++) {
			uint contactCounter = 0;
			//uint firstContactI = 0;
			//uint firstContactJ = 0;
			//uint lastContactI = 0;
			//uint lastContactJ = 0;
			for (uint ii=0; ii<_helicalSegmentResidues[i].size(); ii++) {
				if (ii<_excludeEndCAs || ii >= _helicalSegmentResidues[i].size()-_excludeEndCAs) {
					continue;
				}
				for (uint jj=0; jj<_helicalSegmentResidues[j].size(); jj++) {
					if (jj<_excludeEndCAs || jj >= _helicalSegmentResidues[j].size()-_excludeEndCAs) {
						continue;
					}
					if (_helicalSegmentResidues[i][ii].pCA->distance(*_helicalSegmentResidues[j][jj].pCA) < _minCAdistance) {
					//	cout << i << "/" << j << " Distance " << *_helicalSegmentResidues[i][ii].pCA << " - " << *_helicalSegmentResidues[j][jj] << ": " << _helicalSegmentResidues[i][ii].pCA->distance(*_helicalSegmentResidues[j][jj].pCA << endl;
						//if (contactCounter == 0) {
						//	firstContactI = i;
						//	firstContactJ = j;
						//}
						//lastContactI = i;
						//lastContactJ = j;
						contactCounter++;
					}
				}
			}
			//cout << i << "/" << j << " contacts: " << contactCounter << endl;

			if (contactCounter >= _minNumOfContacts) {
				_interactingSegments.push_back(InteractingHelicalSegments());
				_interactingSegments.back().segment1Index = i;
				_interactingSegments.back().segment2Index = j;
				_interactingSegments.back().pSegment1 = &_helicalSegmentResidues[i];
				_interactingSegments.back().pSegment2 = &_helicalSegmentResidues[j];
				_interactingSegments.back().totalContacts = contactCounter;
			}
		}
	}
}

void calculateInterhelicalGeometry(vector<InteractingHelicalSegments> & _interactingSegments) {
	/*****************************************************************
	 * This function compute the interhelical geometry stored in the InteractingHelicalSegments structures
	 *****************************************************************/

	for (uint k=0; k<_interactingSegments.size(); k++) {
		double minD = 99999.9999;

		for (uint i=0; i<_interactingSegments[k].pSegment1->size(); i++) {
			if (!(*_interactingSegments[k].pSegment1)[i].rHelanal.hasHelanal) {
				continue;
			}
			Line localAxis1((*_interactingSegments[k].pSegment1)[i].rHelanal.h.getNpoint(), (*_interactingSegments[k].pSegment1)[i].rHelanal.h.getCpoint()-(*_interactingSegments[k].pSegment1)[i].rHelanal.h.getNpoint());
			for (uint j=0; j<_interactingSegments[k].pSegment2->size(); j++) {
				if (!(*_interactingSegments[k].pSegment2)[j].rHelanal.hasHelanal) {
					continue;
				}
				Line localAxis2((*_interactingSegments[k].pSegment2)[j].rHelanal.h.getNpoint(), (*_interactingSegments[k].pSegment2)[j].rHelanal.h.getCpoint()-(*_interactingSegments[k].pSegment2)[j].rHelanal.h.getNpoint());
				double d = localAxis1.segmentDistance(localAxis2);
				//cout << "  Helanals " << i << "/" << j << " min D = " << d << endl;
				if (d<minD) {
					_interactingSegments[k].axisOfSegment1 = localAxis1;
					_interactingSegments[k].axisOfSegment2 = localAxis2;
					_interactingSegments[k].axisIndex1 = i;
					_interactingSegments[k].axisIndex2 = j;
					vector<CartesianPoint> PoCA = localAxis1.segmentsClosestPoints(localAxis2); // the points of closest approach on the helanal segments
					_interactingSegments[k].PoCA1 = PoCA[0];
					_interactingSegments[k].PoCA2 = PoCA[1];
					_interactingSegments[k].axialDistance = localAxis1.distance(localAxis2);;
					_interactingSegments[k].segmentDistance = d;
					// to calculate the dihedral we need the points of closest approach on the entire line
					CartesianPoint pointOfMinD1 = localAxis1.pointOfMinDistanceFromLine(localAxis2);
					CartesianPoint pointOfMinD2 = localAxis2.pointOfMinDistanceFromLine(localAxis1);
					_interactingSegments[k].angle = (PoCA[0]+localAxis1.getDirection()).dihedral(PoCA[0], PoCA[1], (PoCA[1]+localAxis2.getDirection()));
					_interactingSegments[k].angle = (pointOfMinD1+localAxis1.getDirection()).dihedral(pointOfMinD1, pointOfMinD2, (pointOfMinD2+localAxis2.getDirection()));
				//	cout << "UUU dihe NEW " << i << endl;
				//	cout << "   " << (pointOfMinD1+localAxis1.getDirection()) << endl;
				//	cout << "   " << pointOfMinD1 << endl;
				//	cout << "   " << pointOfMinD2 << endl;
				//	cout << "   " << (pointOfMinD2+localAxis2.getDirection()) << endl;
				//	cout << "   " << _interactingSegments[k].axialDistance << " " << _interactingSegments[k].segmentDistance << " " << _interactingSegments[k].angle << endl;
					minD = d;
					//foundMinD = true;
					//cout << "    Found better segment min D = " << d << " angle = " << _interhelicalGeometry.back()[1] << endl;
				}

			}
		}

		//cout << "  **** Found best segment " << _interhelicalAxesIndex.back()[0] << "-" << _interhelicalAxesIndex.back()[1] << " min D = " << _interhelicalGeometry.back()[0] << " angle = " << _interhelicalGeometry.back()[1] << endl;
	}
}

void printInterhelicalListReport(string _PDBcode, vector<InteractingHelicalSegments> & _interactingSegments, stringstream & _ss) {
	_ss.str(string());
	_ss << "PDB Id\tSegment 1 #\tSeg 1 Chain #\tSeg 1 Pos start #\tSeg 1 Pos end #\tSeg 1 size\tSeg 1 Pos start Id\tSeg 1 Pos end Id\tSegment 2 #\tSeg 2 Chain #\tSeg 2 Pos start #\tSeg 2 Pos end #\tSeg 2 size\tSeg 2 Pos start Id\tSeg 2 Pos end Id\tHelanal 1 #\tHelanal 2 #\tAxial distance\tLocal distance\tAngle" << endl;

	for (uint i=0; i<_interactingSegments.size(); i++) {
		_ss << _PDBcode << "\t" << _interactingSegments[i].segment1Index  << "\t" << (*_interactingSegments[i].pSegment1)[0].chainIndex << "\t" << (*_interactingSegments[i].pSegment1)[0].posIndexInChain << "\t" << (*_interactingSegments[i].pSegment1).back().posIndexInChain << "\t" << (*_interactingSegments[i].pSegment1).back().posIndexInChain - (*_interactingSegments[i].pSegment1)[0].posIndexInChain + 1 << "\t" << (*_interactingSegments[i].pSegment1)[0].pPosition->getCurrentIdentity().getIdentityId() << "\t" << (*_interactingSegments[i].pSegment1).back().pPosition->getCurrentIdentity().getIdentityId() << "\t" << _interactingSegments[i].segment2Index << "\t" << (*_interactingSegments[i].pSegment2)[0].chainIndex << "\t" << (*_interactingSegments[i].pSegment2)[0].posIndexInChain << "\t" << (*_interactingSegments[i].pSegment2).back().posIndexInChain << "\t" << (*_interactingSegments[i].pSegment2).back().posIndexInChain - (*_interactingSegments[i].pSegment2)[0].posIndexInChain + 1 << "\t" << (*_interactingSegments[i].pSegment2)[0].pPosition->getCurrentIdentity().getIdentityId() << "\t" << (*_interactingSegments[i].pSegment2).back().pPosition->getCurrentIdentity().getIdentityId() << "\t" << _interactingSegments[i].axisIndex1 << "\t" << _interactingSegments[i].axisIndex2 << "\t" << _interactingSegments[i].axialDistance << "\t" << _interactingSegments[i].segmentDistance << "\t" << _interactingSegments[i].angle << endl;
	}
}


void printInterhelicalStandardListReport(string _PDBcode, vector<InteractingHelicalSegments> & _interactingSegments, stringstream & _ss) {
	_ss.str(string());
	_ss << "PDB Id\tSegment 1 #\tSeg 1 Chain #\tSeg 1 Pos start #\tSeg 1 Pos end #\tSeg 1 size\tSeg 1 Pos start Id\tSeg 1 Pos end Id\tSegment 2 #\tSeg 2 Chain #\tSeg 2 Pos start #\tSeg 2 Pos end #\tSeg 2 size\tSeg 2 Pos start Id\tSeg 2 Pos end Id\tHoCA 1 #\tHoCA 2 #\tAxial distance\tLocal distance\tAngle\tHas quadrant 1\tZ 1\tω 1\tZ' 1\tω' 1\tQuad start 1 #\tQuad end 1 #\tHas quadrant 2\tZ 2\tω 2\tZ' 2\tω' 2\tQuad start 2 #\tQuad end 2 #\tSequence 1\tSequence 2\tHelical mask 1\tHelical mask 2\tFit RMSD 1\tFit RMSD 2\tPoCA shift 1\tPoCA shift 2" << endl;

	for (uint i=0; i<_interactingSegments.size(); i++) {
		_ss << _PDBcode << "\t" << _interactingSegments[i].segment1Index  << "\t" << (*_interactingSegments[i].pSegment1)[0].chainIndex << "\t" << (*_interactingSegments[i].pSegment1)[0].posIndexInChain << "\t" << (*_interactingSegments[i].pSegment1).back().posIndexInChain << "\t" << (*_interactingSegments[i].pSegment1).back().posIndexInChain - (*_interactingSegments[i].pSegment1)[0].posIndexInChain + 1 << "\t" << (*_interactingSegments[i].pSegment1)[0].pPosition->getCurrentIdentity().getIdentityId() << "\t" << (*_interactingSegments[i].pSegment1).back().pPosition->getCurrentIdentity().getIdentityId() << "\t" << _interactingSegments[i].segment2Index << "\t" << (*_interactingSegments[i].pSegment2)[0].chainIndex << "\t" << (*_interactingSegments[i].pSegment2)[0].posIndexInChain << "\t" << (*_interactingSegments[i].pSegment2).back().posIndexInChain << "\t" << (*_interactingSegments[i].pSegment2).back().posIndexInChain - (*_interactingSegments[i].pSegment2)[0].posIndexInChain + 1 << "\t" << (*_interactingSegments[i].pSegment2)[0].pPosition->getCurrentIdentity().getIdentityId() << "\t" << (*_interactingSegments[i].pSegment2).back().pPosition->getCurrentIdentity().getIdentityId() << "\t" << _interactingSegments[i].standardHelices.axisIndex1 << "\t" << _interactingSegments[i].standardHelices.axisIndex2 << "\t" << _interactingSegments[i].standardHelices.axialDistance << "\t" << _interactingSegments[i].standardHelices.segmentDistance << "\t" << _interactingSegments[i].standardHelices.angle << "\t";
		if (_interactingSegments[i].standardHelices.PoCA_is_inner1) {
			_ss << "1\t" << _interactingSegments[i].standardHelices.Z1 << "\t" << _interactingSegments[i].standardHelices.omega1 << "\t" << _interactingSegments[i].standardHelices.Zp1 << "\t" << _interactingSegments[i].standardHelices.omegaP1 << "\t" << _interactingSegments[i].standardHelices.interhelicalQuadrant1-5 << "\t" << _interactingSegments[i].standardHelices.interhelicalQuadrant1 << "\t";  
		} else {
			_ss << "0\t\t\t\t\t\t\t";
		}
		if (_interactingSegments[i].standardHelices.PoCA_is_inner2) {
			_ss << "1\t" << _interactingSegments[i].standardHelices.Z2 << "\t" << _interactingSegments[i].standardHelices.omega2 << "\t" << _interactingSegments[i].standardHelices.Zp2 << "\t" << _interactingSegments[i].standardHelices.omegaP2 << "\t" << _interactingSegments[i].standardHelices.interhelicalQuadrant2-5 << "\t" << _interactingSegments[i].standardHelices.interhelicalQuadrant2 << "\t";  
		} else {
			_ss << "0\t\t\t\t\t\t\t";
		}
		_ss << _interactingSegments[i].sequence1 << "\t" << _interactingSegments[i].sequence2 << "\t" << _interactingSegments[i].helicalString1 << "\t" << _interactingSegments[i].helicalString2 << "\t" << _interactingSegments[i].standardHelices.fitRMSD1 << "\t" << _interactingSegments[i].standardHelices.fitRMSD2 << "\t";
		// check how much did the Poing of Closest Approach shift from original to standard helices
		double d1 = _interactingSegments[i].standardHelices.PoCA1.distance(_interactingSegments[i].PoCA1);
		double d2 = _interactingSegments[i].standardHelices.PoCA2.distance(_interactingSegments[i].PoCA2);
		_ss << d1 << "\t" << d2 << endl;
	}
}

bool writeSegmentPDBs(vector<vector<ResidueProperties> > & _helicalSegmentResidues, string _outputdir) {
	bool writeFlag = true;
	for (uint i=0; i<_helicalSegmentResidues.size(); i++) {	
		AtomPointerVector segmentAtoms;
		for (uint j=0; j<_helicalSegmentResidues[i].size(); j++) {	
			AtomPointerVector posAtoms = _helicalSegmentResidues[i][j].pPosition->getAtomPointers();
			segmentAtoms.insert(segmentAtoms.end(), posAtoms.begin(), posAtoms.end() );
		}

		// create the helanal dummy atoms
		AtomPointerVector helanalAtoms;
		for (uint j=0; j<_helicalSegmentResidues[i].size(); j++) {	
			// add the helanal dummy atoms to the atoms
			if (_helicalSegmentResidues[i][j].rHelanal.hasHelanal) {
				stringstream name;
				name << "Z," << i+1 << ",DUM,NPT"; 
				helanalAtoms.push_back(new Atom(name.str(), _helicalSegmentResidues[i][j].rHelanal.h.getNpoint(), "H"));
				name.str("");
				name << "Z," << i+1 << ",DUM,CEN"; 
				helanalAtoms.push_back(new Atom(name.str(), _helicalSegmentResidues[i][j].rHelanal.h.getCenter(), "H"));
				name.str("");
				name << "Z," << i+1 << ",DUM,CPT"; 
				helanalAtoms.push_back(new Atom(name.str(), _helicalSegmentResidues[i][j].rHelanal.h.getCpoint(), "H"));
			}
		}
		segmentAtoms.insert(segmentAtoms.end(), helanalAtoms.begin(), helanalAtoms.end() );

		char filename[1000];
		sprintf(filename, "%s/%s_%03d.pdb", _outputdir.c_str(), "helix" , i);
		cout << "Writing helical segment: " << filename << endl;
		//cout << segmentAtoms;
		PDBWriter writer;
		if (writer.open(filename)) {
			writer.write(segmentAtoms);
			writer.close();
		} else {
			writeFlag = false;
		}

		// release allocated memory
		for (AtomPointerVector::iterator k=helanalAtoms.begin(); k!=helanalAtoms.end(); k++) {
			delete *k;
			*k = NULL;
		}
	}
	return writeFlag;
}

bool writePairsPDBs(vector<InteractingHelicalSegments> & _interactingSegments, string _outputdir) {
	bool writeFlag = true;
	for (uint k=0; k<_interactingSegments.size(); k++) {
		uint helanal1 = _interactingSegments[k].axisIndex1;
		uint helanal2 = _interactingSegments[k].axisIndex2;

		AtomPointerVector segmentAtoms;
		for (uint j=0; j<_interactingSegments[k].pSegment1->size(); j++) {	
			AtomPointerVector posAtoms = (*_interactingSegments[k].pSegment1)[j].pPosition->getAtomPointers();
			segmentAtoms.insert(segmentAtoms.end(), posAtoms.begin(), posAtoms.end() );
		}
		for (uint j=0; j<_interactingSegments[k].pSegment2->size(); j++) {	
			AtomPointerVector posAtoms = (*_interactingSegments[k].pSegment2)[j].pPosition->getAtomPointers();
			segmentAtoms.insert(segmentAtoms.end(), posAtoms.begin(), posAtoms.end() );
		}

		// create the helanal dummy atoms
		AtomPointerVector helanalAtoms;
		stringstream name;
		if ((*_interactingSegments[k].pSegment1)[helanal1].rHelanal.hasHelanal) {
			name.str("");
			name << "Z,1,LAX,NPT"; 
			helanalAtoms.push_back(new Atom(name.str(), (*_interactingSegments[k].pSegment1)[helanal1].rHelanal.h.getNpoint(), "H"));
			name.str("");
			name << "Z,1,LAX,CEN"; 
			helanalAtoms.push_back(new Atom(name.str(), (*_interactingSegments[k].pSegment1)[helanal1].rHelanal.h.getCenter(), "H"));
			name.str("");
			name << "Z,1,LAX,CPT"; 
			helanalAtoms.push_back(new Atom(name.str(), (*_interactingSegments[k].pSegment1)[helanal1].rHelanal.h.getCpoint(), "H"));
		}
		if ((*_interactingSegments[k].pSegment2)[helanal2].rHelanal.hasHelanal) {
			name.str("");
			name << "Z,2,LAX,NPT"; 
			helanalAtoms.push_back(new Atom(name.str(), (*_interactingSegments[k].pSegment2)[helanal2].rHelanal.h.getNpoint(), "H"));
			name.str("");
			name << "Z,2,LAX,CEN"; 
			helanalAtoms.push_back(new Atom(name.str(), (*_interactingSegments[k].pSegment2)[helanal2].rHelanal.h.getCenter(), "H"));
			name.str("");
			name << "Z,2,LAX,CPT"; 
			helanalAtoms.push_back(new Atom(name.str(), (*_interactingSegments[k].pSegment2)[helanal2].rHelanal.h.getCpoint(), "H"));
		}

		segmentAtoms.insert(segmentAtoms.end(), helanalAtoms.begin(), helanalAtoms.end() );

		char filename[1000];
		sprintf(filename, "%s/%s_%03d-%03d.pdb", _outputdir.c_str(),  "pair" , _interactingSegments[k].segment1Index,  _interactingSegments[k].segment2Index); // REMOVE THE +1
		cout << "Writing helical pair: " << filename << endl;
		
		PDBWriter writer;
		if (writer.open(filename)) {
			writer.write(segmentAtoms);
			writer.close();
		} else {
			writeFlag = false;
		}

		// release allocated memory
		for (AtomPointerVector::iterator k=helanalAtoms.begin(); k!=helanalAtoms.end(); k++) {
			delete *k;
			*k = NULL;
		}


	}
	return writeFlag;
}

bool writeStandardPairsPDBs(vector<InteractingHelicalSegments> & _interactingSegments, string _outputdir) {
	for (uint k=0; k<_interactingSegments.size(); k++) {

		AtomPointerVector segmentAtoms;
		// add the standardized CA segments
		for (uint i=0; i<_interactingSegments[k].standardHelices.standardHelixCAs1.size(); i++) {
			segmentAtoms.push_back(_interactingSegments[k].standardHelices.standardHelixCAs1[i]);
		}
		for (uint i=0; i<_interactingSegments[k].standardHelices.standardHelixCAs2.size(); i++) {
			segmentAtoms.push_back(_interactingSegments[k].standardHelices.standardHelixCAs2[i]);
		}

		// create the helanal dummy atoms
		AtomPointerVector helanalAtoms;
		stringstream name;
		name.str("");
		name << "Z,1,HAX,NPT"; 
		helanalAtoms.push_back(new Atom(name.str(), _interactingSegments[k].standardHelices.axis1.getCenter(), "H"));
		name.str("");
		name << "Z,1,HAX,CPT"; 
		helanalAtoms.push_back(new Atom(name.str(), _interactingSegments[k].standardHelices.axis1.getCenter() + _interactingSegments[k].standardHelices.axis1.getDirection(), "H"));
		name.str("");
		name << "Z,2,HAX,NPT"; 
		helanalAtoms.push_back(new Atom(name.str(), _interactingSegments[k].standardHelices.axis2.getCenter(), "H"));
		name.str("");
		name << "Z,2,HAX,CPT"; 
		helanalAtoms.push_back(new Atom(name.str(), _interactingSegments[k].standardHelices.axis2.getCenter() + _interactingSegments[k].standardHelices.axis2.getDirection(), "H"));
		name.str("");
		name << "Z,1,PCA,PCA"; 
		helanalAtoms.push_back(new Atom(name.str(), _interactingSegments[k].standardHelices.PoCA1, "F"));
		name.str("");
		name << "Z,2,PCA,PCA"; 
		helanalAtoms.push_back(new Atom(name.str(), _interactingSegments[k].standardHelices.PoCA2, "F"));

		segmentAtoms.insert(segmentAtoms.end(), helanalAtoms.begin(), helanalAtoms.end() );

		char filename[1000];
		sprintf(filename, "%s/%s_%03d-%03d.pdb", _outputdir.c_str(), "fit" , _interactingSegments[k].segment1Index,  _interactingSegments[k].segment2Index);
		cout << "Writing standard helix fit to pair: " << filename << endl;
		PDBWriter writer;
		writer.open(filename);
		writer.write(segmentAtoms);
		writer.close();

		// release allocated memory
		for (AtomPointerVector::iterator k=helanalAtoms.begin(); k!=helanalAtoms.end(); k++) {
			delete *k;
			*k = NULL;
		}
	}
	return true;

}

void generateStandardHelicalPairs(vector<InteractingHelicalSegments> & _interactingSegments, uint _alignNumber) {

	if (_alignNumber < 4) {
		cerr << "ERROR 8938: need more than four atoms for alignment in function generateStandardHelicalPairs" << endl;
		exit(1);
	}

	// calculate the standard coordinates for each helical pair
	for (uint k=0; k<_interactingSegments.size(); k++) {
		// delete the old standard atoms just in case
		deleteAtoms(_interactingSegments[k].standardHelices.standardHelixCAs1);
		deleteAtoms(_interactingSegments[k].standardHelices.standardHelixCAs2);

		// extract the CA atoms of the original segments for the alignments
		AtomPointerVector segmentCAs1;
		AtomPointerVector segmentCAs2;
		for (uint i=0; i<_interactingSegments[k].pSegment1->size(); i++) {
			segmentCAs1.push_back((*_interactingSegments[k].pSegment1)[i].pCA);
		}
		for (uint i=0; i<_interactingSegments[k].pSegment2->size(); i++) {
			segmentCAs2.push_back((*_interactingSegments[k].pSegment2)[i].pCA);
		}

		int alignStart1 = (int)_interactingSegments[k].axisIndex1 - (((int)_alignNumber - 4)/2);
		uint alignNumber1 = _alignNumber;
		if (alignStart1 < 0) {
			alignNumber1 += alignStart1;
			alignStart1 = 0;
		}
		if (alignStart1 + alignNumber1 >= segmentCAs1.size()) {
			alignNumber1 = segmentCAs1.size() - alignStart1;
		}
		int alignStart2 = (int)_interactingSegments[k].axisIndex2 - (((int)_alignNumber - 4)/2);
		uint alignNumber2 = _alignNumber;
		if (alignStart2 < 0) {
			alignNumber2 = _alignNumber + alignStart2;
			alignStart2 = 0;
		}
		if (alignStart2 + alignNumber2 >= segmentCAs2.size()) {
			alignNumber2 = segmentCAs2.size() - alignStart2;
		}
		//cout << "UUUU helalan of closest approach 1: " << _interactingSegments[k].axisIndex1 << ", helix size " << segmentCAs1.size() << ", align start " << alignStart1 << " and align " << alignNumber1 << " atoms" << endl;
		//cout << "UUUU helalan of closest approach 2: " << _interactingSegments[k].axisIndex2 << ", helix size " << segmentCAs2.size() << ", align start " << alignStart2 << " and align " << alignNumber2 << " atoms" << endl;

		//AtomPointerVector standardAtoms;
		generateStandardAtoms(segmentCAs1, _interactingSegments[k].standardHelices.standardHelixCAs1, alignStart1, alignNumber1, _interactingSegments[k].standardHelices.risePerResidue, _interactingSegments[k].standardHelices.twist, _interactingSegments[k].standardHelices.fitRMSD1);
		//cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
		//cout << segmentCAs1;
		//cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
		//cout << _standardHelicalPairs.back()[0];
		//cout << "--------------------------------------" << endl;
		//cout << "RMSD = " << RMSD << endl;
		//cout << "--------------------------------------" << endl;
		generateStandardAtoms(segmentCAs2, _interactingSegments[k].standardHelices.standardHelixCAs2, alignStart2, alignNumber2, _interactingSegments[k].standardHelices.risePerResidue, _interactingSegments[k].standardHelices.twist, _interactingSegments[k].standardHelices.fitRMSD2);
		//cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
		//cout << segmentCAs2;
		//cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
		//cout << _standardHelicalPairs.back()[1];
		//cout << "--------------------------------------" << endl;
		//cout << "RMSD = " << RMSD << endl;
		//cout << "--------------------------------------" << endl;

		// compute the helanals
		computeStandardHelixHelanals(_interactingSegments[k].standardHelices.standardHelixCAs1, _interactingSegments[k].standardHelices.helanalsHelix1, _interactingSegments[k].standardHelices.hasHelanal1);
		computeStandardHelixHelanals(_interactingSegments[k].standardHelices.standardHelixCAs2, _interactingSegments[k].standardHelices.helanalsHelix2, _interactingSegments[k].standardHelices.hasHelanal2);
	}}

void computeStandardHelixHelanals(AtomPointerVector & _standardHelixCAs, vector<Helanal>  &_helanals, vector<bool> & _hasHelanal) {
	/****************************************************
	 * Return variables:
	 * _helanals: the helanal for each quadruplet of CAs that are bonded
	 * _hasHelanal: a set of flags to indicate positions that do not have helanals
	 ****************************************************/
	_helanals.clear();
	_hasHelanal.clear();

	for (uint i=0; i<_standardHelixCAs.size()-3; i++) {
		Helanal	h(_standardHelixCAs[i]->getCoor(), _standardHelixCAs[i+1]->getCoor(), _standardHelixCAs[i+2]->getCoor(), _standardHelixCAs[i+3]->getCoor());   
		_helanals.push_back(h);
		_hasHelanal.push_back(true);
	}
	for (uint i=0; i<3; i++) {
		// add three empty helanals at the end to have the same size of the positions
		_helanals.push_back(Helanal());
		_hasHelanal.push_back(false);
	}

}
void generateStandardAtoms(AtomPointerVector & _originalCAs, AtomPointerVector & _standardAtoms, uint _alignStart, uint _alignNumber, double & _risePerResidue, double & _twist, double & _RMSD) {
	_risePerResidue = 1.50687727;
	_twist = 99.89647544;

	deleteAtoms(_standardAtoms); // just in case it was full with something (it shouldn't)

	// coordinate of the CAs of a 69 residue standard helix
	vector<CartesianPoint> helix69;
	helix69.push_back(CartesianPoint(-2.08236, -0.90629,-51.23383));
	helix69.push_back(CartesianPoint( 1.25070, -1.89562,-49.72695));
	helix69.push_back(CartesianPoint( 1.65245,  1.55788,-48.22007));
	helix69.push_back(CartesianPoint(-1.81871,  1.36011,-46.71320));
	helix69.push_back(CartesianPoint(-1.02730, -2.02541,-45.20632));
	helix69.push_back(CartesianPoint( 2.17183, -0.66391,-43.69944));
	helix69.push_back(CartesianPoint( 0.28076,  2.25361,-42.19257));
	helix69.push_back(CartesianPoint(-2.26833, -0.11074,-40.68569));
	helix69.push_back(CartesianPoint( 0.49895, -2.21555,-39.17881));
	helix69.push_back(CartesianPoint( 2.09682,  0.87232,-37.67193));
	helix69.push_back(CartesianPoint(-1.21971,  1.91570,-36.16506));
	helix69.push_back(CartesianPoint(-1.67756, -1.53081,-34.65818));
	helix69.push_back(CartesianPoint( 1.79635, -1.38951,-33.15130));
	helix69.push_back(CartesianPoint( 1.06010,  2.00844,-31.64442));
	helix69.push_back(CartesianPoint(-2.16075,  0.69914,-30.13755));
	helix69.push_back(CartesianPoint(-0.31737, -2.24875,-28.63067));
	helix69.push_back(CartesianPoint( 2.26983,  0.07385,-27.12379));
	helix69.push_back(CartesianPoint(-0.46286,  2.22337,-25.61691));
	helix69.push_back(CartesianPoint(-2.11073, -0.83810,-24.11004));
	helix69.push_back(CartesianPoint( 1.18840, -1.93528,-22.60316));
	helix69.push_back(CartesianPoint( 1.70223,  1.50333,-21.09628));
	helix69.push_back(CartesianPoint(-1.77352,  1.41853,-19.58941));
	helix69.push_back(CartesianPoint(-1.09262, -1.99092,-18.08253));
	helix69.push_back(CartesianPoint( 2.14909, -0.73418,-16.57565));
	helix69.push_back(CartesianPoint( 0.35389,  2.24329,-15.06877));
	helix69.push_back(CartesianPoint(-2.27073, -0.03693,-13.56190));
	helix69.push_back(CartesianPoint( 0.42664, -2.23060,-12.05502));
	helix69.push_back(CartesianPoint( 2.12408,  0.80367,-10.54814));
	helix69.push_back(CartesianPoint(-1.15677,  1.95435, -9.04126));
	helix69.push_back(CartesianPoint(-1.72646, -1.47545, -7.53439));
	helix69.push_back(CartesianPoint( 1.75021, -1.44718, -6.02751));
	helix69.push_back(CartesianPoint( 1.12484,  1.97290, -4.52063));
	helix69.push_back(CartesianPoint(-2.13686,  0.76902, -3.01375));
	helix69.push_back(CartesianPoint(-0.39032, -2.23724, -1.50688));
	helix69.push_back(CartesianPoint( 2.27104, -0.00000,  0.00000));
	helix69.push_back(CartesianPoint(-0.39032,  2.23724,  1.50688));
	helix69.push_back(CartesianPoint(-2.13687, -0.76903,  3.01375));
	helix69.push_back(CartesianPoint( 1.12484, -1.97290,  4.52063));
	helix69.push_back(CartesianPoint( 1.75022,  1.44718,  6.02751));
	helix69.push_back(CartesianPoint(-1.72646,  1.47545,  7.53439));
	helix69.push_back(CartesianPoint(-1.15677, -1.95435,  9.04126));
	helix69.push_back(CartesianPoint( 2.12408, -0.80366, 10.54814));
	helix69.push_back(CartesianPoint( 0.42665,  2.23060, 12.05502));
	helix69.push_back(CartesianPoint(-2.27073,  0.03693, 13.56190));
	helix69.push_back(CartesianPoint( 0.35389, -2.24330, 15.06877));
	helix69.push_back(CartesianPoint( 2.14908,  0.73417, 16.57565));
	helix69.push_back(CartesianPoint(-1.09261,  1.99093, 18.08253));
	helix69.push_back(CartesianPoint(-1.77352, -1.41853, 19.58941));
	helix69.push_back(CartesianPoint( 1.70223, -1.50333, 21.09628));
	helix69.push_back(CartesianPoint( 1.18840,  1.93528, 22.60316));
	helix69.push_back(CartesianPoint(-2.11073,  0.83810, 24.11004));
	helix69.push_back(CartesianPoint(-0.46286, -2.22337, 25.61691));
	helix69.push_back(CartesianPoint( 2.26983, -0.07384, 27.12379));
	helix69.push_back(CartesianPoint(-0.31737,  2.24875, 28.63067));
	helix69.push_back(CartesianPoint(-2.16074, -0.69913, 30.13755));
	helix69.push_back(CartesianPoint( 1.06009, -2.00844, 31.64442));
	helix69.push_back(CartesianPoint( 1.79635,  1.38951, 33.15130));
	helix69.push_back(CartesianPoint(-1.67757,  1.53080, 34.65818));
	helix69.push_back(CartesianPoint(-1.21971, -1.91570, 36.16506));
	helix69.push_back(CartesianPoint( 2.09683, -0.87232, 37.67193));
	helix69.push_back(CartesianPoint( 0.49895,  2.21555, 39.17881));
	helix69.push_back(CartesianPoint(-2.26833,  0.11075, 40.68569));
	helix69.push_back(CartesianPoint( 0.28076, -2.25362, 42.19257));
	helix69.push_back(CartesianPoint( 2.17183,  0.66390, 43.69944));
	helix69.push_back(CartesianPoint(-1.02729,  2.02540, 45.20632));
	helix69.push_back(CartesianPoint(-1.81870, -1.36011, 46.71320));
	helix69.push_back(CartesianPoint( 1.65245, -1.55789, 48.22007));
	helix69.push_back(CartesianPoint( 1.25070,  1.89561, 49.72695));
	helix69.push_back(CartesianPoint(-2.08237,  0.90630, 51.23383));

	AtomPointerVector standardAtomsForRMSD;
	AtomPointerVector originalAtomsForRMSD;
	//stringstream name;
	for (uint i=0; i<_originalCAs.size(); i++) {
		//name.str("");
		//name << _originalCAs[i]->getAtomId(); 
		_standardAtoms.push_back(new Atom(_originalCAs[i]->getAtomOfIdentityId(), helix69[i], "C"));
		if (i >= _alignStart && i < _alignStart + _alignNumber) {
			// this is the subset of atoms used to align the standard helix to the original
			standardAtomsForRMSD.push_back(_standardAtoms[i]);
			originalAtomsForRMSD.push_back(_originalCAs[i]);
		}
	}

	Transforms trans;
	trans.rmsdAlignment(standardAtomsForRMSD, originalAtomsForRMSD, _standardAtoms);
	_RMSD = trans.getLastRMSD();


}

void calculateStandardInterhelicalGeometry(vector<InteractingHelicalSegments> & _interactingSegments) {

	for (uint k=0; k<_interactingSegments.size(); k++) {

		// identify the helanal of closest approach
		double minD = 99999.9999;
		//double angle = 0.0;
		uint minI = 0;
		uint minJ = 0;
		//CartesianPoint pointOfMinD1;
		//CartesianPoint pointOfMinD2;
		Line minAxis1;
		Line minAxis2;
		for (uint i=0; i<_interactingSegments[k].standardHelices.helanalsHelix1.size(); i++) {
			if (!_interactingSegments[k].standardHelices.hasHelanal1[i]) {
				continue;
			}
			Line localAxis1(_interactingSegments[k].standardHelices.helanalsHelix1[i].getNpoint(), _interactingSegments[k].standardHelices.helanalsHelix1[i].getCpoint()-_interactingSegments[k].standardHelices.helanalsHelix1[i].getNpoint());
			for (uint j=0; j<_interactingSegments[k].standardHelices.helanalsHelix2.size(); j++) {
				if (!_interactingSegments[k].standardHelices.hasHelanal2[j]) {
					continue;
				}
				Line localAxis2(_interactingSegments[k].standardHelices.helanalsHelix2[j].getNpoint(), _interactingSegments[k].standardHelices.helanalsHelix2[j].getCpoint()-_interactingSegments[k].standardHelices.helanalsHelix2[j].getNpoint());
				double d = localAxis1.segmentDistance(localAxis2);
				//cout << "  Helanals " << i << "/" << j << " min D = " << d << endl;
				if (d<minD) {
					minAxis1 = localAxis1;
					minAxis2 = localAxis2;
					minD = d;
					minI = i;
					minJ = j;
					//cout << "    Found better segment min D = " << d << " angle = " << _standardInterhelicalGeometry.back()[1] << endl;
				}

			}
		}
		_interactingSegments[k].standardHelices.axis1 = minAxis1;
		_interactingSegments[k].standardHelices.axis2 = minAxis2;
		CartesianPoint pointOfMinD1 = minAxis1.pointOfMinDistanceFromLine(minAxis2);
		CartesianPoint pointOfMinD2 = minAxis2.pointOfMinDistanceFromLine(minAxis1);
		_interactingSegments[k].standardHelices.PoCA1 = pointOfMinD1;
		_interactingSegments[k].standardHelices.PoCA2 = pointOfMinD2;

		_interactingSegments[k].standardHelices.axisIndex1 = minI;
		_interactingSegments[k].standardHelices.axisIndex2 = minJ;
		_interactingSegments[k].standardHelices.segmentDistance = minD; // segment distance (may be longer than the line distance)
		_interactingSegments[k].standardHelices.axialDistance = minAxis1.distance(minAxis2); // line distance
		_interactingSegments[k].standardHelices.angle = (pointOfMinD1+minAxis1.getDirection()).dihedral(pointOfMinD1, pointOfMinD2, (pointOfMinD2+minAxis2.getDirection()));
		vector<double> coordinates;
		int quadrant;
		bool PoCA_is_inner;
		//cout << "Segment " << k << " helix 1" << endl;
		calculateZ_OmegaCoordinates(pointOfMinD1, pointOfMinD2, _interactingSegments[k].standardHelices.standardHelixCAs1, minAxis1, minI, _interactingSegments[k].standardHelices.risePerResidue, _interactingSegments[k].standardHelices.twist, coordinates, quadrant, PoCA_is_inner);
		if (PoCA_is_inner) {
			_interactingSegments[k].standardHelices.PoCA_is_inner1 = true;
			_interactingSegments[k].standardHelices.Z1      = coordinates[0];
			_interactingSegments[k].standardHelices.omega1  = coordinates[1];
			_interactingSegments[k].standardHelices.Zp1     = coordinates[2];
			_interactingSegments[k].standardHelices.omegaP1 = coordinates[3];
			_interactingSegments[k].standardHelices.interhelicalQuadrant1 = quadrant;
		}
		//cout << "Added params for helix 1" << endl;
		//cout << "Segment " << k << " helix 2" << endl;
		coordinates.clear();
		calculateZ_OmegaCoordinates(pointOfMinD2, pointOfMinD1, _interactingSegments[k].standardHelices.standardHelixCAs2, minAxis2, minJ, _interactingSegments[k].standardHelices.risePerResidue, _interactingSegments[k].standardHelices.twist, coordinates, quadrant, PoCA_is_inner);
		if (PoCA_is_inner) {
			_interactingSegments[k].standardHelices.PoCA_is_inner2 = true;
			_interactingSegments[k].standardHelices.Z2      = coordinates[0];
			_interactingSegments[k].standardHelices.omega2  = coordinates[1];
			_interactingSegments[k].standardHelices.Zp2     = coordinates[2];
			_interactingSegments[k].standardHelices.omegaP2 = coordinates[3];
			_interactingSegments[k].standardHelices.interhelicalQuadrant2 = quadrant;
		}
		//cout << "Added params for helix 1" << endl;

	}
}


void deleteStandardAtomFits(vector<InteractingHelicalSegments> & _interactingSegments) {
	// release allocated memory
	for (uint i=0; i<_interactingSegments.size(); i++) {
		deleteAtoms(_interactingSegments[i].standardHelices.standardHelixCAs1);
		deleteAtoms(_interactingSegments[i].standardHelices.standardHelixCAs2);
	}
}

void deleteAtoms(AtomPointerVector & _standardAtoms) {
	for (AtomPointerVector::iterator k=_standardAtoms.begin(); k!=_standardAtoms.end(); k++) {
		delete *k;
		*k = NULL;
	}

}


void calculateZ_OmegaCoordinates(CartesianPoint & _pointOfMinDThisHelix, CartesianPoint & _pointOfMinDOtherHelix, AtomPointerVector & _standardHelix, Line & _helanalAxis, uint & _helanalIndex, double _risePerResidue, double _twist, vector<double> & _coordinates, int & _quadrant, bool & _PoCA_is_inner) {

	_PoCA_is_inner = false;
	_quadrant = 0;
	_coordinates = vector<double>(8, 0.0);

	// let's check it the point of closest approach of the helical axis is within the helices, otherwise the geometry
	// does not make sense
	CartesianPoint N_CA_projection_on_axis = _helanalAxis.projection(_standardHelix[0]->getCoor());
	CartesianPoint C_CA_projection_on_axis = _helanalAxis.projection(_standardHelix.back()->getCoor());
	double dX = abs(N_CA_projection_on_axis.getX() - C_CA_projection_on_axis.getX());
	double dY = abs(N_CA_projection_on_axis.getY() - C_CA_projection_on_axis.getY());
	double dZ = abs(N_CA_projection_on_axis.getZ() - C_CA_projection_on_axis.getZ());

	// we use if the dimension (X or Y or Z) with the coordinate with the broader range to check if the point is within range
	if (dX >= dY && dX >= dZ) {
		if (N_CA_projection_on_axis.getX() > C_CA_projection_on_axis.getX()) {
			if (_pointOfMinDThisHelix.getX() <= N_CA_projection_on_axis.getX() && _pointOfMinDThisHelix.getX() >= C_CA_projection_on_axis.getX()) {
				_PoCA_is_inner = true;
			}
		} else {
			if (_pointOfMinDThisHelix.getX() <= C_CA_projection_on_axis.getX() && _pointOfMinDThisHelix.getX() >= N_CA_projection_on_axis.getX()) {
				_PoCA_is_inner = true;
			}
		}
	} else {
		if (dY >= dX && dY >= dZ) {
			if (N_CA_projection_on_axis.getY() > C_CA_projection_on_axis.getY()) {
				if (_pointOfMinDThisHelix.getY() <= N_CA_projection_on_axis.getY() && _pointOfMinDThisHelix.getY() >= C_CA_projection_on_axis.getY()) {
					_PoCA_is_inner = true;
				}
			} else {
				if (_pointOfMinDThisHelix.getY() <= C_CA_projection_on_axis.getY() && _pointOfMinDThisHelix.getY() >= N_CA_projection_on_axis.getY()) {
					_PoCA_is_inner = true;
				}
			}
		} else {
			if (N_CA_projection_on_axis.getZ() > C_CA_projection_on_axis.getZ()) {
				if (_pointOfMinDThisHelix.getZ() <= N_CA_projection_on_axis.getZ() && _pointOfMinDThisHelix.getZ() >= C_CA_projection_on_axis.getZ()) {
					_PoCA_is_inner = true;
				}
			} else {
				if (_pointOfMinDThisHelix.getZ() <= C_CA_projection_on_axis.getZ() && _pointOfMinDThisHelix.getZ() >= N_CA_projection_on_axis.getZ()) {
					_PoCA_is_inner = true;
				}
			}
		}
	}

	if(!_PoCA_is_inner) {
		// makes no sense to calculate Z' and W' when the point of closest approach is not within the helix
		//cout << "Out of range" << endl;
		return;
	}

	int ref_atom_index = (int)_helanalIndex + 3;
	int relativeQuadrant;
	//cout << "Starting ref atom " << ref_atom_index << endl;
	Atom * refCA = _standardHelix[ref_atom_index];
	quadrantCoordinates(_pointOfMinDOtherHelix, _helanalAxis, refCA, _risePerResidue, _twist, _coordinates, relativeQuadrant);
	//cout << "     after TRY 1, adjust quadrant1 by " << relativeQuadrant << ", ref atom " << (int)_helanalIndex + 3 - relativeQuadrant << endl;
	//for (uint i=0; i<_coordinates.size(); i++) {
	//	cout << "   " << i << " " << _coordinates[i] << endl;
	//}
	if (relativeQuadrant != 0) {
		int ref_atom_index = (int)_helanalIndex + 3 - relativeQuadrant;
		if (ref_atom_index >= 5 && ref_atom_index < _standardHelix.size()) {
			refCA = _standardHelix[ref_atom_index];
			quadrantCoordinates(_pointOfMinDOtherHelix, _helanalAxis, refCA, _risePerResidue, _twist, _coordinates, relativeQuadrant);
			//cout << "        after TRY 2, adjust quadrant1 by " << relativeQuadrant << ", ref atom " << ref_atom_index - relativeQuadrant<< endl;
			//for (uint i=0; i<_coordinates.size(); i++) {
			//	cout << "   " << i << " " << _coordinates[i] << endl;
			//}
			if (relativeQuadrant == 0) {
				_quadrant = ref_atom_index;
				_PoCA_is_inner = true;
			} else {
				_PoCA_is_inner = false;
			//	cout << "UUU quadrant not valid " << relativeQuadrant << ", ref atom " << ref_atom_index << " (size " << _standardHelix.size() << ")" << endl;
			}
		} else {
			_PoCA_is_inner = false;
			//cout << "UUU ref atom index not valid " << ref_atom_index << " (size " << _standardHelix.size() << ")" << endl;
		}
	} else {
		_quadrant = ref_atom_index;
		_PoCA_is_inner = true;
	}

}

void quadrantCoordinates(CartesianPoint & _pointOnOtherHelixAxis, Line & _helAxis, Atom * _refCAs, double _risePerResidue, double _twist, vector<double> & _coordinates, int & _quadrant) {

	/*******************************************
	 * This function calculates the Z and omega (W) coordinates
	 * of an interaction given the following parameters
	 *   - the point of closest approach on the other helix
	 *   - the helical axis obtained with helanal
	 *   - the reference CA corresponding to the C2 position (see scheme)
	 *   - the rise per residue
	 *   - the twist (degrees per residue)
	 * 
	 * It returns the coordinates in a 8 member vector of double _coordinates
	 *   1: the Z coordinate with respect to the reference CA
	 *   2: the W coordinate   "    "     "   "   "        "
	 *   3: the Z' coordinate  "    "     "   "   "        "
	 *   4: the W' coordinate  "    "     "   "   "        "
	 *   5: the Z coordinate with respect to the CA of the quadrant of closest approach
	 *   6: the W coordinate   "      "    "  "   "  "  "      "     "    "        "
	 *   7: the Z' coordinate  "      "    "  "   "  "  "      "     "    "        "
	 *   8: the W' coordinate  "      "    "  "   "  "  "      "     "    "        "
	 * 
	 * It also returns the relative index of the closest approach quadrant _quadrant
	 *
	 * In the Z W coordinate system, the coordinates of the quadrants
	 * are expressed in rise-per-residue (R) and twist (T)
	 *
	 *                           -----N1 (5R,5T)
	 *                      -----     /
	 *       (4R,4T) N2-----         /  
	 *               /              /   
	 *              /              /    
	 *             /              /     
	 *            /              / 
	 *           /              /       
	 *          /              /   
	 *         /         -----C1 (R,T)
	 *        /     -----     
	 * (0,0) C2-----         
	 *
	 * In the Z' W' coordinate system, the coordinates of the quadrants
	 * are precisely
	 *
	 *                           -----N1 (6,100)
	 *                      -----     /
	 *         (6,0) N2-----         /  
	 *               /              /   
	 *              /              /    
	 *             /              /     
	 *            /              / 
	 *           /              /       
	 *          /              /   
	 *         /         -----C1 (0,100)
	 *        /     -----     
	 * (0,0) C2-----         
	 *
	 * Z' = alpha * Z + beta * W
	 * W' = gamma * Z + delta * W
	 *******************************************/
	// the conversion factors
	double alpha = _twist / (60.0 * _risePerResidue);
	double beta = -1.0 / 60.0;
	double gamma = (900.0 - 10.0 * _twist) / (9.0 * _risePerResidue);
	double delta = 10.0 / 9.0;

	// reset the coordinates
	_coordinates.clear();
	_quadrant = 0;


	CartesianPoint CA_projection = _helAxis.projection(_refCAs->getCoor());
	CartesianPoint pointProjection = _helAxis.projection(_pointOnOtherHelixAxis);

	vector<double> out(8, 0.0);
	out[0] = (CA_projection - pointProjection) * _helAxis.getDirection().getUnit(); // point rise
	out[1] = _pointOnOtherHelixAxis.dihedral(_helAxis.getCenter(), _helAxis.getCenter()+_helAxis.getDirection(), _refCAs->getCoor()); // dihedral
	//cout << "RISE: (" << CA_projection << " - " << pointProjection << ") * " << _helAxis.getDirection().getUnit() << " = " << out[0] << endl;
	//cout << "DIHE: " << _pointOnOtherHelixAxis << ", " << _helAxis.getCenter() << ", " << _helAxis.getCenter()+_helAxis.getDirection() << ", " << _refCAs->getCoor() << " = " << out[1] << endl;
	if (out[1]<0.0) {
		out[1] += 360.0;
	}
	out[2] = alpha * out[0] + beta * out[1];
	out[3] = gamma * out[0] + delta * out[1];
	out[4] = out[0];
	out[5] = out[1];
	out[6] = out[2];
	out[7] = out[3];
	bool done = false;
	while (!done) {
		done = true;
		if (out[6] > 6.0) {
			out[4] -= 4.0 * _risePerResidue;
			out[5] -= 4.0 * _twist;
			out[6] -= 6.0;
			_quadrant += 4;
			done = false;
		}
		if (out[6] < 0.0) {
			out[4] += 4.0 * _risePerResidue;
			out[5] += 4.0 * _twist;
			out[6] += 6.0;
			_quadrant -= 4;
			done = false;
		}
		if (out[7] > 100.0) {
			out[4] -= _risePerResidue;
			out[5] -= _twist;
			out[7] -= 100.0;
			_quadrant += 1;
			done = false;
		}
		if (out[7] < 0.0) {
			out[4] += _risePerResidue;
			out[5] += _twist;
			out[7] += 100.0;
			_quadrant -= 1;
			done = false;
		}
	}
	while (out[5] >= 360.0) {
		out[5] -= 360.0;
	}
	while (out[5] < 0.0) {
		out[5] += 360.0;
	}
	_coordinates = out;
	
}

void classifyHelicalPhiPsi(vector<vector<ResidueProperties> > & _chainResidues) {
	vector<vector<int > > allowed;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -180;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -170;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -170;
	allowed.back()[1] = -20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -170;
	allowed.back()[1] = -60;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -170;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -170;
	allowed.back()[1] = 0;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -170;
	allowed.back()[1] = 10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -170;
	allowed.back()[1] = 20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -170;
	allowed.back()[1] = 30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -160;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -160;
	allowed.back()[1] = -20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -160;
	allowed.back()[1] = -30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -160;
	allowed.back()[1] = -40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -160;
	allowed.back()[1] = -60;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -160;
	allowed.back()[1] = 0;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -160;
	allowed.back()[1] = 10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -160;
	allowed.back()[1] = 20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -160;
	allowed.back()[1] = 30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -160;
	allowed.back()[1] = 40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = -20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = -30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = -40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = -50;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = -60;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = 0;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = 10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = 20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = 30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -150;
	allowed.back()[1] = 40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = -20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = -30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = -40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = -50;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = -60;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = 0;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = 10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = 20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = 30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -140;
	allowed.back()[1] = 40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = -20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = -30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = -40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = -50;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = -60;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = 0;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = 10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = 20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = 30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -130;
	allowed.back()[1] = 40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -120;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -120;
	allowed.back()[1] = -40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -120;
	allowed.back()[1] = -50;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -120;
	allowed.back()[1] = -60;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -120;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -120;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -120;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -120;
	allowed.back()[1] = 20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -120;
	allowed.back()[1] = 30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -120;
	allowed.back()[1] = 40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -110;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -110;
	allowed.back()[1] = -50;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -110;
	allowed.back()[1] = -60;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -110;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -110;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -110;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -110;
	allowed.back()[1] = 30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -110;
	allowed.back()[1] = 40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -100;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -100;
	allowed.back()[1] = -60;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -100;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -100;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -100;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -100;
	allowed.back()[1] = 20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -100;
	allowed.back()[1] = 30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -100;
	allowed.back()[1] = 40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -90;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -90;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -90;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -90;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -90;
	allowed.back()[1] = 20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -90;
	allowed.back()[1] = 30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -90;
	allowed.back()[1] = 40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -80;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -80;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -80;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -80;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -80;
	allowed.back()[1] = 10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -80;
	allowed.back()[1] = 20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -80;
	allowed.back()[1] = 30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -80;
	allowed.back()[1] = 40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -70;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -70;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -70;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -70;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -70;
	allowed.back()[1] = 0;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -70;
	allowed.back()[1] = 10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -70;
	allowed.back()[1] = 20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -70;
	allowed.back()[1] = 30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -70;
	allowed.back()[1] = 40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -60;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -60;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -60;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -60;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -60;
	allowed.back()[1] = 0;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -60;
	allowed.back()[1] = 10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -60;
	allowed.back()[1] = 20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -50;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -50;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -50;
	allowed.back()[1] = -20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -50;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -50;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -50;
	allowed.back()[1] = 0;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -40;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -40;
	allowed.back()[1] = -10;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -40;
	allowed.back()[1] = -20;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -40;
	allowed.back()[1] = -30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -40;
	allowed.back()[1] = -40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -40;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -40;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -30;
	allowed.back()[1] = -30;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -30;
	allowed.back()[1] = -40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -30;
	allowed.back()[1] = -50;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -30;
	allowed.back()[1] = -60;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -30;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -30;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -30;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -20;
	allowed.back()[1] = -40;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -20;
	allowed.back()[1] = -50;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -20;
	allowed.back()[1] = -60;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -20;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -20;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -20;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -10;
	allowed.back()[1] = -60;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -10;
	allowed.back()[1] = -70;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -10;
	allowed.back()[1] = -80;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = -10;
	allowed.back()[1] = -90;
	allowed.push_back(vector<int>(2, 0.0));
	allowed.back()[0] = 0;
	allowed.back()[1] = -80;

	vector<vector<int> > favorable;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -130;
	favorable.back()[1] = -10;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -120;
	favorable.back()[1] = -10;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -120;
	favorable.back()[1] = -20;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -120;
	favorable.back()[1] = -30;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -120;
	favorable.back()[1] = 0;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -120;
	favorable.back()[1] = 10;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -110;
	favorable.back()[1] = -10;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -110;
	favorable.back()[1] = -20;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -110;
	favorable.back()[1] = -30;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -110;
	favorable.back()[1] = -40;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -110;
	favorable.back()[1] = 0;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -110;
	favorable.back()[1] = 10;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -110;
	favorable.back()[1] = 20;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -100;
	favorable.back()[1] = -10;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -100;
	favorable.back()[1] = -20;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -100;
	favorable.back()[1] = -30;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -100;
	favorable.back()[1] = -40;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -100;
	favorable.back()[1] = -50;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -100;
	favorable.back()[1] = 0;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -100;
	favorable.back()[1] = 10;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -90;
	favorable.back()[1] = -10;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -90;
	favorable.back()[1] = -20;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -90;
	favorable.back()[1] = -30;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -90;
	favorable.back()[1] = -40;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -90;
	favorable.back()[1] = -50;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -90;
	favorable.back()[1] = -60;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -90;
	favorable.back()[1] = 0;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -90;
	favorable.back()[1] = 10;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -80;
	favorable.back()[1] = -10;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -80;
	favorable.back()[1] = -20;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -80;
	favorable.back()[1] = -30;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -80;
	favorable.back()[1] = -40;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -80;
	favorable.back()[1] = -50;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -80;
	favorable.back()[1] = -60;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -80;
	favorable.back()[1] = 0;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -70;
	favorable.back()[1] = -10;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -70;
	favorable.back()[1] = -20;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -70;
	favorable.back()[1] = -30;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -70;
	favorable.back()[1] = -40;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -70;
	favorable.back()[1] = -50;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -70;
	favorable.back()[1] = -60;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -70;
	favorable.back()[1] = -70;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -60;
	favorable.back()[1] = -20;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -60;
	favorable.back()[1] = -30;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -60;
	favorable.back()[1] = -40;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -60;
	favorable.back()[1] = -50;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -60;
	favorable.back()[1] = -60;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -60;
	favorable.back()[1] = -70;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -50;
	favorable.back()[1] = -30;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -50;
	favorable.back()[1] = -40;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -50;
	favorable.back()[1] = -50;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -50;
	favorable.back()[1] = -60;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -50;
	favorable.back()[1] = -70;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -40;
	favorable.back()[1] = -50;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -40;
	favorable.back()[1] = -60;
	favorable.push_back(vector<int>(2, 0.0));
	favorable.back()[0] = -40;
	favorable.back()[1] = -70;
	
	for (uint i=0; i<_chainResidues.size(); i++) {
		for (uint j=0; j<_chainResidues[i].size(); j++) {
			if (!_chainResidues[i][j].rPhiPsi.hasPhiPsi[0] && !_chainResidues[i][j].rPhiPsi.hasPhiPsi[1]) {
				// no phi or psi
				///cout << "UUUU position has no phi/psi " << _chainResidues[i][j].rPhiPsi.angle[0] << "/" << _chainResidues[i][j].rPhiPsi.angle[1] << endl;
				_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical = false;
				_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical = false;
				continue;
			}
			if (_chainResidues[i][j].rPhiPsi.hasPhiPsi[0] && _chainResidues[i][j].rPhiPsi.hasPhiPsi[1]) {
				// has both phi and psi
				int floorPhi = int(floor(_chainResidues[i][j].rPhiPsi.angle[0]/10.0)*10.0);
				int floorPsi = int(floor(_chainResidues[i][j].rPhiPsi.angle[1]/10.0)*10.0);
				bool foundFavorable = false;
				bool foundAllowed = false;
				for (uint k=0; k<favorable.size(); k++) {
					if (floorPhi == favorable[k][0] && floorPsi == favorable[k][1]) {
						///cout << "UUUU phi/psi " << _chainResidues[i][j].rPhiPsi.angle[0] << "/" << _chainResidues[i][j].rPhiPsi.angle[1] << "(floors " << floorPhi << "/" << floorPsi << ") is favorable" << endl;
						foundFavorable = true;
						_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical = true;
						_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical = true;
						break;
					}
				}
				if (foundFavorable) {
					continue;
				} else {
					for (uint k=0; k<allowed.size(); k++) {
						if (floorPhi == allowed[k][0] && floorPsi == allowed[k][1]) {
							///cout << "UUUU phi/psi " << _resPhiPsis[i][j].angle[0] << "/" << _resPhiPsis[i][j].angle[1] << "(floors " << floorPhi << "/" << floorPsi << ") is allowed" << endl;
							foundAllowed = true;
							_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical = false;
							_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical = true;
							break;
						}
					}
				}
				if (!foundAllowed) {
					///cout << "UUUU phi/psi " << _resPhiPsis[i][j].angle[0] << "/" << _chainResidues[i][j].rPhiPsi.angle[1] << "(floors " << floorPhi << "/" << floorPsi << ") is not allowed" << endl;
					_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical = false;
					_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical = false;
				} 
			} else {
				if (_chainResidues[i][j].rPhiPsi.hasPhiPsi[0]) {
					// phi only
					int floorPhi = int(floor(_chainResidues[i][j].rPhiPsi.angle[0]/10.0)*10.0);
					bool foundFavorable = false;
					bool foundAllowed = false;
					for (uint k=0; k<favorable.size(); k++) {
						if (floorPhi == favorable[k][0]) {
							///cout << "UUUU position with phi only phi/psi " << _chainResidues[i][j].rPhiPsi.angle[0] << "/" << _chainResidues[i][j].rPhiPsi.angle[1] << "(floors " << floorPhi << ") is favorable" << endl;
							foundFavorable = true;
							_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical = true;
							_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical = true;
							break;
						}
					}
					if (foundFavorable) {
						continue;
					} else {
						for (uint k=0; k<allowed.size(); k++) {
							if (floorPhi == allowed[k][0]) {
								///cout << "UUUU position with phi only phi/psi " << _chainResidues[i][j].rPhiPsi.angle[0] << "/" << _chainResidues[i][j].rPhiPsi.angle[1] << "(floors " << floorPhi << ") is allowed" << endl;
								foundAllowed = true;
								_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical = false;
								_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical = true;
								break;
							}
						}
					}
					if (!foundAllowed) {
						///cout << "UUUU position with phi only phi/psi " << _chainResidues[i][j].rPhiPsi.angle[0] << "/" << _chainResidues[i][j].rPhiPsi.angle[1] << "(floors " << floorPhi << ") is not allowed" << endl;
						_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical = false;
						_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical = false;
					} 
				} else {
					// psi only
					int floorPsi = int(floor(_chainResidues[i][j].rPhiPsi.angle[1]/10.0)*10.0);
					bool foundFavorable = false;
					bool foundAllowed = false;
					for (uint k=0; k<favorable.size(); k++) {
						if (floorPsi == favorable[k][1]) {
							///cout << "UUUU position with psi only phi/psi " << _chainResidues[i][j].rPhiPsi.angle[0] << "/" << _chainResidues[i][j].rPhiPsi.angle[1] << "(floors " << floorPsi << ") is favorable" << endl;
							foundFavorable = true;
							_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical = true;
							_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical = true;
							break;
						}
					}
					if (foundFavorable) {
						continue;
					} else {
						for (uint k=0; k<allowed.size(); k++) {
							if (floorPsi == allowed[k][1]) {
								///cout << "UUUU position with psi only phi/psi " << _chainResidues[i][j].rPhiPsi.angle[0] << "/" << _chainResidues[i][j].rPhiPsi.angle[1] << "(floors " << floorPsi << ") is allowed" << endl;
								foundAllowed = true;
								_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical = false;
								_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical = true;
								break;
							}
						}
					}
					if (!foundAllowed) {
						///cout << "UUUU position with psi only phi/psi " << _chainResidues[i][j].rPhiPsi.angle[0] << "/" << _chainResidues[i][j].rPhiPsi.angle[1] << "(floors " << floorPsi << ") is not allowed" << endl;
						_chainResidues[i][j].rPhiPsi.isPhiPsiStrictlyHelical = false;
						_chainResidues[i][j].rPhiPsi.isPhiPsiLooselyHelical = false;
					} 
				}
			}
		}
	}
}

/*
void classifyHelicalHbonding(vector<vector<ResidueProperties> > & _chainResidues, double _minHBondDistanceStrict, double _minHBondDistanceLoose) {
	for (uint i=0; i<_chainResidues.size(); i++) {
		for (uint j=0; j<_chainResidues[i].size(); j++) {
			if (!_chainResidues[i][j].rHbond.hasHbonding[0] && !_chainResidues[i][j].rHbond.hasHbonding[1]) {
				// no Hbonding either way
				///cout << "UUUU position has no phi/psi " << _hBonding[i][j][0] << "/" << _hBonding[i][j][1] << endl;
				_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = false;
				_chainResidues[i][j].rHbond.isHbondingLooselyHelical = false;
				continue;
			}
			if (_chainResidues[i][j].rHbond.hasHbonding[0] && _chainResidues[i][j].rHbond.hasHbonding[1]) {
				// has both N and C terminal hbonding
				if (_chainResidues[i][j].rHbond.distance[0] <= _minHBondDistanceLoose && _chainResidues[i][j].rHbond.distance[1] <= _minHBondDistanceLoose) {
					_chainResidues[i][j].rHbond.isHbondingLooselyHelical = true;
					if (_chainResidues[i][j].rHbond.distance[0] <= _minHBondDistanceStrict && _chainResidues[i][j].rHbond.distance[1] <= _minHBondDistanceStrict) {
						_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = true;
					} else {
						_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = false;
					}
				} else {
					_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = false;
					_chainResidues[i][j].rHbond.isHbondingLooselyHelical = false;
				}
			} else {
				if (_chainResidues[i][j].rHbond.hasHbonding[0]) {
					// N terminal only
					if (_chainResidues[i][j].rHbond.distance[0] <= _minHBondDistanceLoose) {
						_chainResidues[i][j].rHbond.isHbondingLooselyHelical = true;
						if (_chainResidues[i][j].rHbond.distance[0] <= _minHBondDistanceStrict) {
							_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = true;
						} else {
							_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = false;
						}
					} else {
						_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = false;
						_chainResidues[i][j].rHbond.isHbondingLooselyHelical = false;
					}
				} else {
					// C terminal only
					if (_chainResidues[i][j].rHbond.distance[1] <= _minHBondDistanceLoose) {
						_chainResidues[i][j].rHbond.isHbondingLooselyHelical = true;
						if (_chainResidues[i][j].rHbond.distance[1] <= _minHBondDistanceStrict) {
							_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = true;
						} else {
							_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = false;
						}
					} else {
						_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = false;
						_chainResidues[i][j].rHbond.isHbondingLooselyHelical = false;
					}
				}
			}
		}
	}
}
*/

void classifyHelicalHbonding(vector<vector<ResidueProperties> > & _chainResidues, double _minHBondDistanceStrict, double _minHBondDistanceLoose) {
	// revised function that no longer uses the i,i+4 H-bond, only the i-4
	for (uint i=0; i<_chainResidues.size(); i++) {
		for (uint j=0; j<_chainResidues[i].size(); j++) {
			if (!_chainResidues[i][j].rHbond.hasHbonding[0]) {
				// no Hbonding
				_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = false;
				_chainResidues[i][j].rHbond.isHbondingLooselyHelical = false;
				continue;
			} else {
				// has both N and C terminal hbonding
				if (_chainResidues[i][j].rHbond.distance[0] <= _minHBondDistanceLoose) {
					_chainResidues[i][j].rHbond.isHbondingLooselyHelical = true;
					if (_chainResidues[i][j].rHbond.distance[0] <= _minHBondDistanceStrict) {
						_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = true;
					} else {
						_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = false;
					}
				} else {
					_chainResidues[i][j].rHbond.isHbondingStrictlyHelical = false;
					_chainResidues[i][j].rHbond.isHbondingLooselyHelical = false;
				}
			}
		}
	}
}

void classifyByZcoordinate(vector<vector<ResidueProperties> > & _chainResidues, double _minZ, double _maxZ) {
	for (uint i=0; i<_chainResidues.size(); i++) {
		for (uint j=0; j<_chainResidues[i].size(); j++) {
			if (!_chainResidues[i][j].hasCA || _chainResidues[i][j].pCA->getZ() < _minZ || _chainResidues[i][j].pCA->getZ() > _maxZ) {
				_chainResidues[i][j].isInMembrane = false;
			} else {
				_chainResidues[i][j].isInMembrane = true;
			}
		}
	}
}

ICOptions parseICOptions(int _argc, char * _argv[], ICOptions defaults) {

	/******************************************
	 *  Pass the array of argument and the name of
	 *  a configuration file to the ArgumentParser
	 *  object.  Then ask for the value of the argument
	 *  and collect error and warnings.
	 *
	 *  This function returns a ICOptions structure
	 *  defined at the head of this file 
	 ******************************************/
	
	ICOptions opt;

	/******************************************
	 *  Set the allowed and required options:
	 *
	 *  Example of configuartion file:
	 *  
	 ******************************************/
	vector<string> required;
	vector<string> allowed;

	opt.required.push_back("filename");

	opt.allowed.push_back("PDBcode");
	opt.allowed.push_back("logFile");
	opt.allowed.push_back("outputdir");

	opt.allowed.push_back("backboneBondMin");
	opt.allowed.push_back("minHelicalSegmentSize");

	opt.allowed.push_back("minHBondDistanceStrict");
	opt.allowed.push_back("minHBondDistanceLoose");
	opt.allowed.push_back("minCAdistance");
	opt.allowed.push_back("minNumOfContacts");
	opt.allowed.push_back("excludeEndCAs");
	opt.allowed.push_back("numberOfCAinAlignment");
	opt.allowed.push_back("imposeZlimits");
	opt.allowed.push_back("minZ");
	opt.allowed.push_back("maxZ");

	opt.defaultArgs.push_back("filename"); // the default arguments can be specified in command line without "--option"

	//opt.allowed.push_back("verbose");

	//OptionParser OP(_argc, _argv);
	OptionParser OP;
	OP.setShortOptionEquivalent(opt.equivalent);
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
	if (OP.countOptions() == 0) {
		usage();
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
	 *  VERSION AND HELP
	 *
	 *  --version or -v arguments print the version number
	 *  --help prints the usage and help
	 *****************************************/
	opt.version = OP.getBool("version");
	//if (OP.fail()) {
	//	opt.version = OP.getBool("v");
	//}

	if (opt.version) {
		version();
		exit(0);
	}

	opt.help = OP.getBool("help");
	if (OP.fail()) {
		opt.help = OP.getBool("h");
	}

	if (opt.help) {
		help(defaults);
		exit(0);
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
	opt.filename = OP.getString("filename");
	if (OP.fail()) {
		opt.errorMessages += "PDB file name not specified\n";
		opt.errorFlag = true;
	}

	opt.PDBcode = OP.getDouble("PDBcode");
	if (OP.fail()) {
		opt.PDBcode = MslTools::getFileName(opt.filename);
		opt.warningMessages += "PDB code not specified, using file name based code: " + opt.PDBcode + "\n";
		opt.warningFlag = true;
	}

	opt.backboneBondMin = OP.getDouble("backboneBondMin");
	if (OP.fail()) {
		opt.backboneBondMin = defaults.backboneBondMin;
	}

	opt.minHelicalSegmentSize = OP.getDouble("minHelicalSegmentSize");
	if (OP.fail()) {
		opt.minHelicalSegmentSize = defaults.minHelicalSegmentSize;
	}

	opt.minHBondDistanceStrict = OP.getDouble("minHBondDistanceStrict");
	if (OP.fail()) {
		opt.minHBondDistanceStrict = defaults.minHBondDistanceStrict;
	}

	opt.minHBondDistanceLoose = OP.getDouble("minHBondDistanceLoose");
	if (OP.fail()) {
		opt.minHBondDistanceLoose = defaults.minHBondDistanceLoose;
	}

	opt.minCAdistance = OP.getDouble("minCAdistance");
	if (OP.fail()) {
		opt.minCAdistance = defaults.minCAdistance;
	}

	opt.minNumOfContacts = OP.getDouble("minNumOfContacts");
	if (OP.fail()) {
		opt.minNumOfContacts = defaults.minNumOfContacts;
	}

	opt.excludeEndCAs = OP.getDouble("excludeEndCAs");
	if (OP.fail()) {
		opt.excludeEndCAs = defaults.excludeEndCAs;
	}

	opt.numberOfCAinAlignment = OP.getDouble("numberOfCAinAlignment");
	if (OP.fail()) {
		opt.numberOfCAinAlignment = defaults.numberOfCAinAlignment;
	}

	// let's figure out if we need to impose Z limits
	opt.imposeZlimits = OP.getBool("imposeZlimits");
	bool Z_fail_flag = false;
	if (OP.fail()) {
		Z_fail_flag = true;
	}
	opt.minZ = OP.getDouble("minZ");
	bool minZ_fail_flag = false;
	if (OP.fail()) {
		minZ_fail_flag = true;
		opt.minZ = defaults.minZ;
	}
	opt.maxZ = OP.getDouble("maxZ");
	bool maxZ_fail_flag = false;
	if (OP.fail()) {
		maxZ_fail_flag = true;
		opt.maxZ = defaults.maxZ;
	}
	if (Z_fail_flag && (!minZ_fail_flag || !maxZ_fail_flag)) {
		// the --imposeZlimits wasn't give but minZ or maxZ were given: we should impose limits
		opt.imposeZlimits = true;
	}

	opt.outputdir = OP.getString("outputdir");
	if (OP.fail()) {
		opt.outputdir = opt.PDBcode + "-" + MslTools::getRandomAlphaNumString(6);
		opt.warningMessages += "outputdir not specified using " + opt.outputdir + "\n";
		opt.warningFlag = true;
	}

	opt.logFile = OP.getString("logFile");
	if (OP.fail()) {
		opt.logFile = opt.outputdir + "/" + opt.PDBcode + ".log";
		opt.warningMessages += "logFile not specified using " + opt.logFile + "\n";
		opt.warningFlag = true;
	}

	opt.rerunConf = "########################################################\n";
	opt.rerunConf += "#  This configuration file was automatically generated,\n";
	opt.rerunConf += "#  it will rerun this job with the same options. Run as:\n";
	opt.rerunConf += "#\n";
	opt.rerunConf += "#  Run as:\n";
	opt.rerunConf += "#\n";
	opt.rerunConf += "#    % " + programName + " --configfile " + opt.rerunConfFile + "\n";
	opt.rerunConf += "#\n";
	opt.rerunConf += "#  Job started on " + (string)ctime(&start_time);
	opt.rerunConf += "#  on host " + opt.host + ", path " + opt.pwd + "\n";
	opt.rerunConf += "########################################################\n";
	opt.rerunConf += "\n";
	opt.rerunConf += OP.getConfFile();
	
//	opt.verbose = OP.getBool("verbose");
//	if (OP.fail()) {
//		opt.warningMessages += "verbose not specified using false\n";
//		opt.warningFlag = true;
//		opt.verbose = false;
//	}

	
	// return the ICOptions structure

	return opt;

}

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}

void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

void help(ICOptions defaults) {
	cout << "This program runs as:" << endl;
	cout << " % interhelicalCoordinates <file.pdb> " << endl;
	cout << endl;
	cout << "   Optional Parameters " << endl;
	cout << "   --PDBcode                 the PDB code (the default is the file name, i.e. 1c3w.pdb -> 1c3w" << endl;
	cout << "   --minHelicalSegmentSize   the minimal lenght for a helical segment to be considered (default " << defaults.minHelicalSegmentSize << ")" << endl;
	cout << "   --minHBondDistanceStrict  the minimal distance between backbone N and O atoms at i, i+4 for considering the segment strictly helical (default " << defaults.minHBondDistanceStrict << ")"  << endl;
	cout << "   --minHBondDistanceLoose   the minimal distance between backbone N and O atoms at i, i+4 for considering the segment loosely helical (default " << defaults.minHBondDistanceLoose << ")"  << endl;
	cout << "   --minNumOfContacts        the min number of CA contacts between two helices to be considered a pair (default " << defaults.minNumOfContacts << ")" << endl;
	cout << "   --minCAdistance           the disstance between CAs in the two helices to be considered a contact (default " << defaults.minCAdistance << ")" << endl;
	cout << "   --excludeEndCAs           the number of terminal CAs to be exluded on each end from the contact count (default " << defaults.excludeEndCAs << ")" << endl;
	cout << "   --numberOfCAinAlignment   the number of CA used for aligning the natural helices to the standard ones at the point of closest approach (default " << defaults.numberOfCAinAlignment << ")" << endl;
	cout << "   --backboneBondMin         the minimal distance between backbone atoms (N, CA, C) for considering the backbone not broken (default " << defaults.backboneBondMin << ")" << endl;
	cout << endl;
}

void generateInterhelicalString(vector<InteractingHelicalSegments> & _interactingSegments) {
	/*************************************************************
	 * This function crates a string (a mask) for the sequence indicating
	 * if a position is strictly helical.
	 *  H = strictly helical
	 *  h = non strictly helical (thus loosely helical)
	 * If a position is at the verteces of the quadrant of closest approach
	 * these are marked as 
	 *  P if strictly helical
	 *  p if not strictly helical
	 *************************************************************/
	for (uint i=0; i<_interactingSegments.size(); i++) {
		// create the helical code
		string hString1;
		string hString2;
		uint quadrant1 = 0;
		uint quadrant2 = 0;
		bool hasQuadrant1 = false;
		bool hasQuadrant2 = false;
		if (_interactingSegments[i].standardHelices.PoCA_is_inner1) {
			quadrant1 = _interactingSegments[i].standardHelices.interhelicalQuadrant1;;
			hasQuadrant1 = true;
		}
		if (_interactingSegments[i].standardHelices.PoCA_is_inner2) {
			quadrant2 =  _interactingSegments[i].standardHelices.interhelicalQuadrant2;
			hasQuadrant2 = true;
		}
		for (uint j=0; j<_interactingSegments[i].pSegment1->size(); j++) {
			if ((*_interactingSegments[i].pSegment1)[j].isStrictlyHelical) {
				if (hasQuadrant1 && (j==quadrant1-5 || j==quadrant1-4 || j==quadrant1-1 || j==quadrant1)) {
					hString1 += "P";
				} else {
					hString1 += "H";
				}
			} else {
				if (hasQuadrant1 && (j==quadrant1-5 || j==quadrant1-4 || j==quadrant1-1 || j==quadrant1)) {
					hString1 += "p";
				} else {
					hString1 += "h";
				}
			}
		}
		for (uint j=0; j<_interactingSegments[i].pSegment2->size(); j++) {
			if ((*_interactingSegments[i].pSegment2)[j].isStrictlyHelical) {
				if (hasQuadrant2 && (j==quadrant2-5 || j==quadrant2-4 || j==quadrant2-1 || j==quadrant2)) {
					hString2 += "P";
				} else {
					hString2 += "H";
				}
			} else {
				if (hasQuadrant2 && (j==quadrant2-5 || j==quadrant2-4 || j==quadrant2-1 || j==quadrant2)) {
					hString2 += "p";
				} else {
					hString2 += "h";
				}
			}
		}
		_interactingSegments[i].helicalString1 = hString1;
		_interactingSegments[i].helicalString2 = hString2;

		// get the sequence (1 letter string) i.e. LIVLGILLFNIAACS
		string sequence1;
		for (uint j=0; j<_interactingSegments[i].pSegment1->size(); j++) {
			string resname = (*_interactingSegments[i].pSegment1)[j].pPosition->getResidueName();
			sequence1 += MslTools::getOneLetterCode(resname);
		}
		string sequence2;
		for (uint j=0; j<_interactingSegments[i].pSegment2->size(); j++) {
			string resname = (*_interactingSegments[i].pSegment2)[j].pPosition->getResidueName();
			sequence2 += MslTools::getOneLetterCode(resname);
		}

		_interactingSegments[i].sequence1 = sequence1;
		_interactingSegments[i].sequence2 = sequence2;
	}
}

bool writeReport(string _filename, string _outputdir, stringstream & _ss) {

	ofstream fout;
	char filename[1000];
	sprintf(filename, "%s/%s.csv", _outputdir.c_str(), _filename.c_str());
	cout << "Writing report: " << filename << endl;
	fout.open(filename);
	if(!fout.is_open()) {
		cerr << "Unable to open " << filename << endl;
		return false;
	}
	fout << _ss.str();
	fout.close();
	return true;
}

