/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A "Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)" (2012) J. Comput. Chem, 33, 1645-61 
 DOI: 10.1002/jcc.22968

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/

#include "MslTools.h"
#include "System.h"
#include "ChiStatistics.h"
#include "SasaCalculator.h"
#include "OptionParser.h"
#include "release.h"


// Define Input Options and Store them.

struct Options {

	std::string origList;
	std::string repackedList;
	std::string chiDefFile;
	double tolerance;
	std::string sasaThreshFile;

	bool version; // ask for the program version
	bool help; // ask for the program help


	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;
	vector<vector<string> > equivalent; // this links short options to long ones (for example -x can be given for --extended)
	std::map<std::string,double> refSasa;

	string errorMessages;
	string warningMessages;

	bool errorFlag;
	bool warningFlag;
};


Options parseOptions(int theArgc, char * theArgv[]);
void usage();
void version();
void help();

