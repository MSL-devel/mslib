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


#include "DegreeOfFreedomReader.h"

using namespace MSL;
using namespace std;

DegreeOfFreedomReader::DegreeOfFreedomReader() : Reader() {
}

DegreeOfFreedomReader::~DegreeOfFreedomReader() {
}

bool DegreeOfFreedomReader::read(string _defiFile) { 
	/**************************************************
	 *  This overwrites the definitions from file.
	 *
	 *    File format:
	 *
	 * ALA phi -C N CA C 
	 * ALA psi N CA C +N 
	 * ARG phi -C N CA C 
	 * ARG psi N CA C +N 
	 * ARG chi1 N CA CB CG 
	 * ARG chi2 CA CB CG CD 
	 * ARG chi3 CB CG CD NE 
	 * ARG chi4 CG CD NE CZ 
	 * ASN phi -C N CA C 
	 * ASN psi N CA C +N 
	 * ASN chi1 N CA CB CG 
	 * ASN chi2 CA CB CG OD1 
	 * ...
	 **************************************************/

	if (!open(_defiFile)) {
		cerr << "WARNING 34842: cannot read definition file " << _defiFile << endl;
		return false;
	}
	vector<string> lines = getAllLines();
	close();

	degOfFreedomlabels.clear();

	for (unsigned int i=0; i<lines.size(); i++) {
		string line = MslTools::uncomment(lines[i]);
		vector<string> tokens = MslTools::tokenizeAndTrim(line, " ");
		if (tokens.size() == 0) {
			continue;
		}
		//for (unsigned int j=0; j<tokens.size(); j++) {
		//	cout << tokens[j] << " ";
		//}
		//cout << endl;
		if (tokens.size() < 4) {
			// not enough elements
			cerr << "WARNING 34847: syntax error in definition line " << line << " in file " << _defiFile << endl;
			continue;
		}
		for (unsigned int j=2; j<tokens.size(); j++) {
			// 0: residue   1: label   2-5: atoms
			degOfFreedomlabels[tokens[0]][tokens[1]].push_back(tokens[j]);
		}
	}


	return true;

}

std::vector<std::string> DegreeOfFreedomReader::getSingleDegreeOfFreedom(std::string _resname, std::string _deegreOfFreedom) {
	std::vector<std::string> out;

	std::map<std::string, std::map<std::string, std::vector<std::string> > >::const_iterator found = degOfFreedomlabels.find(_resname);
	if (found == degOfFreedomlabels.end()) {
		// residue not found in labels
		std::cerr << "WARNING 44842: definition for residue name  " << _resname << " not found" << std::endl;
		return out;
	}

	std::map<std::string, std::vector<std::string> >::const_iterator found2 = found->second.find(_deegreOfFreedom);

	if (found2 == found->second.end()) {
		// residue not found in labels
		std::cerr << "WARNING 44847: degree of freedom " << _deegreOfFreedom << " not found for residue name  " << _resname << std::endl;
		return out;
	}

	return found2->second;

}

vector<vector<string> >  DegreeOfFreedomReader::getChiAtoms(string _resName) {
	if(degOfFreedomlabels.find(_resName) != degOfFreedomlabels.end()) {
		if(degOfFreedomlabels[_resName].find("chi1") != degOfFreedomlabels[_resName].end()) {
			vector<vector<string> > chis;
			chis.push_back(degOfFreedomlabels[_resName]["chi1"]);

			if(degOfFreedomlabels[_resName].find("chi2") != degOfFreedomlabels[_resName].end()) { 
				chis.push_back(degOfFreedomlabels[_resName]["chi2"]);

				if(degOfFreedomlabels[_resName].find("chi3") != degOfFreedomlabels[_resName].end()) {
					chis.push_back(degOfFreedomlabels[_resName]["chi3"]);

					if(degOfFreedomlabels[_resName].find("chi4") != degOfFreedomlabels[_resName].end()) {
						chis.push_back(degOfFreedomlabels[_resName]["chi4"]);
					}
				}
			}
			return chis;

		} else {
			return vector<vector<string> >();
		}
	} else {
		return vector<vector<string> >();
	}

}

