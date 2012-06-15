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

#ifndef DEGREEOFFREEDOMREADER_H
#define DEGREEOFFREEDOMREADER_H
/*

 */

//MSL Includes
#include "Reader.h"
#include "MslTools.h"

// STL Includes
#include <vector>

/**
  This class will provide an object which is able
  to read a file that defines degrees of freedom such as
  chi1, chi2, phi and psi

	     File format:
	 
ALA phi -C N CA C 
ALA psi N CA C +N 
ARG phi -C N CA C 
ARG psi N CA C +N 
ARG chi1 N CA CB CG 
ARG chi2 CA CB CG CD 
ARG chi3 CB CG CD NE 
ARG chi4 CG CD NE CZ 
ASN phi -C N CA C 
ASN psi N CA C +N 
ASN chi1 N CA CB CG 
ASN chi2 CA CB CG OD1 
...


 */
namespace MSL { 
class DegreeOfFreedomReader : public Reader {

	public:
		// Constructors/Destructors
		DegreeOfFreedomReader();
		~DegreeOfFreedomReader();

		bool read(std::string _defiFile);

		// the structure is map of a map such as ["LEU"]["chi1"] = (N, CA, CB, CG)
		std::map<std::string, std::map<std::string, std::vector<std::string> > > getDegreesOfFreedom() const; 

		// this returns the atoms for a certain residue (LEU) and label (chi1)
		std::vector<std::string> getSingleDegreeOfFreedom(std::string _resname, std::string _deegreOfFreedom); 
		// this returns the atoms for all chis of certain residue (LEU) in order
		std::vector< std::vector<std::string> >  getChiAtoms(std::string _resName);

	private:
		std::map<std::string, std::map<std::string, std::vector<std::string> > > degOfFreedomlabels;

};

inline std::map<std::string, std::map<std::string, std::vector<std::string> > > DegreeOfFreedomReader::getDegreesOfFreedom() const {
	return degOfFreedomlabels;
}

}

#endif
