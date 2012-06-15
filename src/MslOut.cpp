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

#include "MslOut.h"
using namespace MSL;


// STL maps do not seem to work as static variables on OSX..
//std::map<std::string, bool> MslOut::outputFlag = std::map<std::string,bool>();

//std::vector<std::string> MslOut::outputFlagOn = std::vector<std::string>();
//std::vector<std::string> MslOut::allFlags  = std::vector<std::string>();

//std::vector<std::string> MslOut::outputFlagOn;
//std::vector<std::string> MslOut::allFlags;


MslOut::MslOut() : std::ostream(NULL) {

	static std::vector<std::string> &allFlags = MslOut::getAllFlags();
	static std::vector<std::string> &outputOnFlags = MslOut::getOutputOnFlags();

	uniqueAdd(allFlags,"MSL");
	uniqueAdd(outputOnFlags,"MSL");
	name = "MSL";
	turnAllOff();
}
MslOut::MslOut(std::string _name): std::ostream(NULL) {

	static std::vector<std::string> &allFlags = MslOut::getAllFlags();
	static std::vector<std::string> &outputOnFlags = MslOut::getOutputOnFlags();

	uniqueAdd(allFlags,_name);
        uniqueAdd(outputOnFlags,_name);

	name = _name;
	turnAllOff();
}

std::vector<std::string>& MslOut::getAllFlags(){
  static std::vector<std::string>* allFlagsVar = new std::vector<std::string>();
  return *allFlagsVar;
}

std::vector<std::string>& MslOut::getOutputOnFlags(){
  static std::vector<std::string>* outputOnFlagsVar = new std::vector<std::string>();
  return *outputOnFlagsVar;
}

std::map<std::string, double>& MslOut::getStaticLookup(){
  static std::map<std::string,double>* aStaticMap = new std::map<std::string,double>();
  return *aStaticMap;
}
