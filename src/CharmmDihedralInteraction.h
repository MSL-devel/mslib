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

#ifndef CHARMMDIHEDRALINTERACTION_H
#define CHARMMDIHEDRALINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "FourBodyInteraction.h"
#include "CharmmEnergy.h"


namespace MSL { 
class CharmmDihedralInteraction: public FourBodyInteraction {

	/*******************************************************
	 *   Inherits from TwoBodyInteraction (a prototype object
	 *   for the interaction of two atoms)
	 *******************************************************/

	public:
		CharmmDihedralInteraction();
		CharmmDihedralInteraction(Atom & _a1, Atom & _a2, Atom & _a3, Atom & _a4, double _Kchi, double _N, double _deltaRadians); // Delta must be in Radians
		CharmmDihedralInteraction(Atom & _a1, Atom & _a2, Atom & _a3, Atom & _a4, std::vector<std::vector <double > > _multipleParams); // Delta must be in Radians

		// should implement an operator= as well 
		CharmmDihedralInteraction(const CharmmDihedralInteraction & _interaction);
		~CharmmDihedralInteraction();

		/* setting and getting the parameters */
		void setParams(std::vector<double> _params);
		void setParams(std::vector<std::vector<double> > _multipleParams);
		void setParams(double _Kchi, double _N, double _deltaRadians);
		std::vector<double> & getParams();
		std::vector<std::vector<double> > & getMultipleParams();

		double getEnergy();
		double getEnergy(double _angleRadians,std::vector<double> *_ad=NULL);
		double getEnergy(std::vector<double> *_ad);

		std::vector<double> getEnergyGrad();
		std::vector<double> getEnergyGrad(Atom& a1, Atom& a2, Atom& a3, Atom& a4, double Kchi, double N, double deltaRadians);


		friend std::ostream & operator<<(std::ostream &_os, CharmmDihedralInteraction & _term) {_os << _term.toString(); return _os;};
		std::string toString() ;

		//unsigned int getType() const;
		std::string getName() const;

		bool isSelected(std::string _selection1, std::string _selection2) const;
		std::pair<double,std::vector<double> > partialDerivative();
		
	private:
		void setup(Atom * _pA1, Atom * _pA2, Atom * _pA3, Atom * _pA4, std::vector<std::vector <double> >  _params);
		void copy(const CharmmDihedralInteraction & _interaction);
		//static const unsigned int type = 5;
		static const std::string typeName;

		std::vector<std::vector<double> > multipleParams; //dihedrals can have multiple etries with different multiplicity (N)
		

};
inline void CharmmDihedralInteraction::setParams(std::vector<double> _params) { std::vector<std::vector<double> > multi; multi.push_back(_params); setParams(multi);}
inline void CharmmDihedralInteraction::setParams(std::vector<std::vector<double> > _multipleParams) {
	if (_multipleParams.size() == 0) {
		multipleParams.push_back(std::vector<double>(3, 0.0));
		multipleParams[0][1] = 1.0;
	} else {
		for (unsigned int i=0; i<_multipleParams.size(); i++) {
			if (_multipleParams[i].size() != 3) {
				std::cerr << "ERROR 54119: invalid number of parameters in inline void CharmmDihedralInteraction::setParams(std::vector<std::vector<double> > _multipleParams)" << std::endl; exit(54119);
			}
			multipleParams = _multipleParams;
		}
	}
}
inline void CharmmDihedralInteraction::setParams(double _Kchi, double _N, double _deltaRadians) { multipleParams = std::vector<std::vector<double> >(1, std::vector<double>(3, 0.0)); multipleParams[0][0] = _Kchi; multipleParams[0][1] = _N; multipleParams[0][2] = _deltaRadians;}
inline std::vector<double> & CharmmDihedralInteraction::getParams() {if(multipleParams.size() > 1) {std::cerr << "WARNING 48199: dihedral might contain multiple parameters" << std::endl;} return multipleParams[0];}
inline std::vector<std::vector<double> > & CharmmDihedralInteraction::getMultipleParams() {return multipleParams;}
inline double CharmmDihedralInteraction::getEnergy() {
	double energy = 0.0;
	double angle = pAtoms[0]->dihedralRadians(*pAtoms[1], *pAtoms[2], *pAtoms[3]);
	for (unsigned int i=0; i<multipleParams.size(); i++) {
		energy += CharmmEnergy::instance()->dihedralEner(angle, multipleParams[i][0], multipleParams[i][1], multipleParams[i][2]);
	}
	return energy;
}
inline double CharmmDihedralInteraction::getEnergy(double _angleRadians,std::vector<double> *_ad) {
	double energy = 0.0;
	for (unsigned int i=0; i<multipleParams.size(); i++) {
		energy += CharmmEnergy::instance()->dihedralEner(_angleRadians, multipleParams[i][0], multipleParams[i][1], multipleParams[i][2],_ad);
	}
	return energy;
}
inline double CharmmDihedralInteraction::getEnergy(std::vector<double> *_ad) {
	if(_ad) {
		double angleRadians = CartesianGeometry::dihedralDerivative(pAtoms[0]->getCoor(),pAtoms[1]->getCoor(),pAtoms[2]->getCoor(),pAtoms[3]->getCoor(),_ad);
		return getEnergy(angleRadians,_ad);
	}
	return getEnergy();
}
inline std::string CharmmDihedralInteraction::toString() {
	std::string out;
	char c [1000];
	sprintf(c, "CHARMM DIHE %s %s %s %s", pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), pAtoms[2]->toString().c_str(), pAtoms[3]->toString().c_str());
	out += c;
	for (unsigned int i=0; i<multipleParams.size(); i++) {
		sprintf(c, " %9.4f %9.4f %9.4f", multipleParams[i][0], multipleParams[i][1], multipleParams[i][2] * 180.0 / M_PI);
		out += c;
		if (i != multipleParams.size() -1) {
			out += "/";
		} else {
			out += " ";
		}
	}
	sprintf(c, "%9.4f %20.6f", pAtoms[0]->dihedral(*pAtoms[1], *pAtoms[2], *pAtoms[3]) * 180.0 / M_PI, getEnergy());
	out += c;
	return out;
}
//inline unsigned int CharmmDihedralInteraction::getType() const {return type;}
inline std::string CharmmDihedralInteraction::getName() const {return typeName;}
inline bool CharmmDihedralInteraction::isSelected(std::string _selection1, std::string _selection2) const {
	if ( (pAtoms[1]->getSelectionFlag(_selection1) && pAtoms[2]->getSelectionFlag(_selection2)) || (pAtoms[1]->getSelectionFlag(_selection2) && pAtoms[2]->getSelectionFlag(_selection1)) ) {
		return true;
	} else {
		return false;
	}
}


inline std::vector<double> CharmmDihedralInteraction::getEnergyGrad(){
	
	std::vector<double> result(12,0.0);
	for (unsigned int i=0; i<multipleParams.size(); i++) {
		std::vector<double> grad = getEnergyGrad(*pAtoms[0],*pAtoms[1],*pAtoms[2],*pAtoms[3],multipleParams[i][0], multipleParams[i][1], multipleParams[i][2]);

		for (uint j = 0; j < 12;j++){
			result[j] += grad[j];
		}
	}

	return result;
}
inline std::pair<double,std::vector<double> > CharmmDihedralInteraction::partialDerivative() {
	std::pair<double, std::vector<double> > partials;
	partials.first = CartesianGeometry::dihedralDerivative(pAtoms[0]->getCoor(),pAtoms[1]->getCoor(),pAtoms[2]->getCoor(),pAtoms[3]->getCoor(),&(partials.second));
	return partials;
}

}

#endif

