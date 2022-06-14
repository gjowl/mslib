#include "MslTools.h"
#include "OptionParser.h"
using namespace std;
using namespace MSL;

string programName = "pdbFilter";
string programDescription = "This program reads a pdb file and writes a filtered version based on user options (ex. only the peptide atoms for chain B of the first model in the file).  This filtered file can then be used for further analysis in MSL."
string programAuthor = "Samson Condon";
string programVersion = "0.0.1";
string programDate = "20 January 2014";
string mslVersion = MSLVERSION;
string mslDate = MSLDATE;
