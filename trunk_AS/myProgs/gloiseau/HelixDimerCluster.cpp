#include <sstream>
#include <iterator>
#include <unistd.h>
#include "HelixDimerCluster.h"

using namespace std;
using namespace MSL;

HelixDimerCluster::HelixDimerCluster(HelixDimer * _structure) {
	members.push_back(_structure);
}

void HelixDimerCluster::addHelixDimer(HelixDimer * _structure) {
	if(members.size() < 1) {
		members.push_back(_structure);
		return;
	}
	//if(_structure->getEnergy() < 0) {
	members.push_back(_structure);
	//}
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
	for (uint i=0; i<members.size(); i++){
		map<string,double>& deltaEByTerm = members[i]->getDeltaEnergyByTerm();
		eOut << "Model " << i << endl;
		for(map<string,double>::iterator it = deltaEByTerm.begin(); it != deltaEByTerm.end(); it++) {
			if (it->first.compare("CHARMM_ELEC") != 0){
				if (it->second != 0) {
					eOut << it->first << " " << it->second << endl;
				}
			}
		}
		string elec = "CHARMM_ELEC";
		if (deltaEByTerm.find(elec) != deltaEByTerm.end()){
			auto it = deltaEByTerm.find(elec);
			eOut << it->first << " " << it->second << endl;
		}
		eOut << endl;
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
		//TODO: when reorganizing I had to get rid of the below: can I put that info into a file that all of my code for this program inherits from?
		//writer.addRemark("MSL " + mslVersion + ", CATM " + programVersion);

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