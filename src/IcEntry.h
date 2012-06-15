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

#ifndef ICENTRY_H
#define ICENTRY_H

#include <vector>
#include <math.h>
#include <map>


#include "Atom.h"

	/****************************************************
	 *  
	 *           I M P O R T A N T ! ! ! 
	 *  
	 *  Test storage buffers functions!!!
	 *  
	 ****************************************************/



#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/map.hpp>
#endif

namespace MSL { 
class IcTable;

class IcEntry {
	/****************************************************
	 *  This class builds atom a1 based on a2 a3 a4 and
	 *  the a1-a2 distance, a1-a2-a3 angle and a1-a2-a3-a4
	 *  dihedral
	 *
	 *  a1-a2-a3-a4
	 ****************************************************/
	public:
		IcEntry();
		IcEntry(Atom & _atom1, Atom & _atom2, Atom & _atom3, Atom & _atom4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag=false);
		IcEntry(const IcEntry & _ic);
		~IcEntry();

		void operator=(const IcEntry & _ic);

		/*******************************************
		 *   Print
		 *******************************************/
		std::string toString() const;
		friend std::ostream & operator<<(std::ostream &_os, const IcEntry & _ic)  {_os << _ic.toString(); return _os;};

		/*******************************************
		 *   Setters and getters
		 *******************************************/
		void setAtom1(Atom & _atom);
		void setAtom2(Atom & _atom);
		void setAtom3(Atom & _atom);
		void setAtom4(Atom & _atom);
		bool removeAtom(Atom * _pAtom);

		void setDistance1(double _d);
		void setDistance2(double _d);
		void setAngle1(double _a);
		void setAngle2(double _a);
		void setDihedral(double _dihe);
		void setAngle1Radians(double _a);
		void setAngle2Radians(double _a);
		void setDihedralRadians(double _dihe);
		void setImproper(bool _isImproper);

		Atom * getAtom1() const;
		Atom * getAtom2() const;
		Atom * getAtom3() const;
		Atom * getAtom4() const;
		double getDistance1() const;
		double getDistance2() const;
		double getAngle1() const;
		double getAngle2() const;
		double getDihedral() const;
		double getAngle1Radians() const;
		double getAngle2Radians() const;
		double getDihedralRadians() const;
		bool isImproper() const;
		std::vector<double> & getValues();

		// functions to check what atoms make the ic entry or make specific degrees of freedom
		bool match(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3, Atom * _pAtom4) const;
		bool areD1Atoms(Atom * _pAtom1, Atom * _pAtom2) const;
		bool areD2Atoms(Atom * _pAtom1, Atom * _pAtom2) const;
		bool areA1Atoms(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3) const;
		bool areA2Atoms(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3) const;

		/********************************************************
		 *   Building functions
		 *
		 *   The ic can build either pAtom1 (build1) or pAtom4(build4), or 
		 *   the pointer can tell the function what atom we want to build
		 *
		 *   This is a recursive function.  In order, for exampe to build 1
		 *   the function needs 2 3 and 4 to be built.  The function calls 
		 *   atoms' 2 3 and 4 buildFromIc() function: if they do not have
		 *   coordinates they will call the build() function in the IcEntry 
		 *   that can build them.  That might call other atoms to call
		 *   their IcEntry build() function, and so on, until everything
		 *   that needs to be built to build atom 1 is built.
		 *
		 ********************************************************/
		//bool build1(std::map<Atom*, bool> & _exclude, bool _onlyActive=true);
		//bool build4(std::map<Atom*, bool> & _exclude, bool _onlyActive=true);
		//bool build(Atom * _pAtom, std::map<Atom*, bool> & _exclude, bool _onlyActive=true);

		bool build1(std::map<IcEntry*, bool> & _exclude, bool _onlyActive=true);
		bool build4(std::map<IcEntry*, bool> & _exclude, bool _onlyActive=true);
		bool build(Atom * _pAtom, std::map<IcEntry*, bool> & _exclude, bool _onlyActive=true);

		/********************************************************
		 *  Set the internal coordinates from existing coordinates
		 ********************************************************/
		void fillFromCoor();

		/********************************************************
		 *  Save and restore IC entries to buffers
		 *
		 *  NEED TESTS!
		 ********************************************************/
		void saveBuffer(std::string _name);
		bool restoreFromBuffer(std::string _name);
		void clearAllBuffers();
		std::map<std::string, std::vector<double> > getStoredValues() const;
		void setStoredValues(std::map<std::string, std::vector<double> > _buffers);
		void setParentIcTable(IcTable * _table);
		IcTable * getParentTable() const;
		
		bool isValid() const; // valid if not null key positions
		void blankIC();

	private:
		void setup(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3, Atom * _pAtom4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag);
		void copy(const IcEntry & _ic);
		void updateParentMap(Atom * _pOldAtom, Atom * _pNewAtom);

		Atom * pAtom1;
		Atom * pAtom2;
		Atom * pAtom3;
		Atom * pAtom4;

		std::vector<double> vals;
		//double d1;
		//double a1;
		//double dihe;
		//double a2;
		//double d2;

		std::map<std::string, std::vector<double> > storedValues;
		
		/********************************************************
		 *   Atom 4 is built differently if it is an improper dihedral
		 *   or a
		 ********************************************************/
		bool improperFlag;

		IcTable * pParentTable ;


#ifdef __BOOST__
	public:
		void save_checkpoint(std::string filename) const{
			std::ofstream fout(filename.c_str());
			boost::archive::text_oarchive oa(fout);
			oa << (*this);
		}

		void load_checkpoint(std::string filename){
			std::ifstream fin(filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive ia(fin);
			ia >> (*this);
		}


		friend class boost::serialization::access;		

	private:
		template<class Archive> void serialize(Archive & ar, const unsigned int version){
			ar & pAtom1;
			ar & pAtom2;
			ar & pAtom3;
			ar & pAtom4;

			ar & vals;
			ar & storedValues; 
			ar & improperFlag;
		}
#endif

};

inline void IcEntry::setDistance1(double _d) {vals[0] = _d;}
inline void IcEntry::setDistance2(double _d) {vals[4] = _d;}
inline void IcEntry::setAngle1(double _a) {vals[1] = _a * M_PI / 180.0;}; // store as radians
inline void IcEntry::setAngle2(double _a) {vals[3] = _a * M_PI / 180.0;}
inline void IcEntry::setDihedral(double _dihe) {vals[2] = _dihe * M_PI / 180.0;}
inline void IcEntry::setAngle1Radians(double _a) {vals[1] = _a;}; // store as radians
inline void IcEntry::setAngle2Radians(double _a) {vals[3] = _a;}
inline void IcEntry::setDihedralRadians(double _dihe) {vals[2] = _dihe;}
inline void IcEntry::setImproper(bool _isImproper) {improperFlag = _isImproper;}
inline Atom * IcEntry::getAtom1() const {return pAtom1;}
inline Atom * IcEntry::getAtom2() const {return pAtom2;}
inline Atom * IcEntry::getAtom3() const {return pAtom3;}
inline Atom * IcEntry::getAtom4() const {return pAtom4;}
inline double IcEntry::getDistance1() const {return vals[0];}
inline double IcEntry::getDistance2() const {return vals[4];}
inline double IcEntry::getAngle1() const {return vals[1] * 180.0 / M_PI;}
inline double IcEntry::getAngle2() const {return vals[3] * 180.0 / M_PI;}
inline double IcEntry::getDihedral() const {return vals[2] * 180.0 / M_PI;}
inline double IcEntry::getAngle1Radians() const {return vals[1];}
inline double IcEntry::getAngle2Radians() const {return vals[3];}
inline double IcEntry::getDihedralRadians() const {return vals[2];}
inline bool IcEntry::isImproper() const {return improperFlag;}
inline std::vector<double> & IcEntry::getValues() {return vals;}
inline void IcEntry::saveBuffer(std::string _name) {storedValues[_name] = vals;}
inline bool IcEntry::restoreFromBuffer(std::string _name) {
	std::map<std::string, std::vector<double> >::iterator found=storedValues.find(_name);
	if (found==storedValues.end()) {
		return false;
	} else {
		for (unsigned int i=0; i<vals.size(); i++) {
			vals[i] = found->second[i];
		}
		return true;
	}
}
inline void IcEntry::clearAllBuffers() {storedValues.clear();}
inline std::map<std::string, std::vector<double> > IcEntry::getStoredValues() const {return storedValues;}
inline void IcEntry::setStoredValues(std::map<std::string, std::vector<double> > _buffers) {storedValues = _buffers;}
inline void IcEntry::setParentIcTable(IcTable * _table) { pParentTable = _table; }
inline IcTable * IcEntry:: getParentTable() const { return pParentTable;}

// fucntion to recognize if the IC contains certain atoms
inline bool IcEntry::match(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3, Atom * _pAtom4) const {
	if (pAtom1 == _pAtom1 && pAtom2 == _pAtom2 && pAtom3 == _pAtom3 && pAtom4 == _pAtom4) {
		// same order
		return true;
	} else {
		// reverse order is different for improper and dihe
		if (improperFlag) {
			if (pAtom1 == _pAtom4 && pAtom2 == _pAtom2 && pAtom3 == _pAtom3 && pAtom4 == _pAtom1) {
				return true;
			}
		} else {
			if (pAtom1 == _pAtom4 && pAtom2 == _pAtom3 && pAtom3 == _pAtom2 && pAtom4 == _pAtom1) {
				return true;
			}
		}
	}
	return false;
} 
inline bool IcEntry::areD1Atoms(Atom * _pAtom1, Atom * _pAtom2) const {
	// is distance 1 relative to these two atoms
	if (improperFlag) {
		if ((pAtom1 == _pAtom1 && pAtom3 == _pAtom2) || (pAtom1 == _pAtom2 && pAtom3 == _pAtom1)) {
			return true;
		}
	} else {
		if ((pAtom1 == _pAtom1 && pAtom2 == _pAtom2) || (pAtom1 == _pAtom2 && pAtom2 == _pAtom1)) {
			return true;
		}
	}
	return false;
}
inline bool IcEntry::areD2Atoms(Atom * _pAtom1, Atom * _pAtom2) const {
	// is distance 2 relative to these two atoms
	if ((pAtom4 == _pAtom1 && pAtom3 == _pAtom2) || (pAtom4 == _pAtom2 && pAtom3 == _pAtom1)) {
		return true;
	}
	return false;
}
inline bool IcEntry::areA1Atoms(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3) const {
	// is angle 1 relative to these two atoms
	if (improperFlag) {
		if (pAtom3 == _pAtom2 && ((pAtom1 == _pAtom1 && pAtom2 == _pAtom3) || (pAtom1 == _pAtom3 && pAtom2 == _pAtom1))) {
			return true;
		}
	} else {
		if (pAtom2 == _pAtom2 && ((pAtom1 == _pAtom1 && pAtom3 == _pAtom3) || (pAtom1 == _pAtom3 && pAtom3 == _pAtom1))) {
			return true;
		}
	}
	return false;
}
inline bool IcEntry::areA2Atoms(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3) const {
	// is angle 2 relative to these two atoms
	if (pAtom3 == _pAtom2 && ((pAtom4 == _pAtom1 && pAtom2 == _pAtom3) || (pAtom4 == _pAtom3 && pAtom2 == _pAtom1))) {
		return true;
	}
	return false;
}
inline bool IcEntry::isValid() const {
	/**********************************************
	 * valid if not null key positions
	 *  1: 1-2-X-X  valid not improper   
	 *  2: 1-X-3-X  valid improper       
	 *  3: 1-X-X-4  not valid            
	 *  4: X-2-3-X  not valid            
	 *  5: X-2-X-4  not valid      
	 *  6: X-X-3-4  valid
	 *  7: 1-2-3-X  valid
	 *  8: 1-2-X-4  valid not improper
	 *  9: 1-X-3-4  valid
	 * 10: X-2-3-4  valid
	 * 11: 1-2-3-4  valid
	 * 12: 1-X-X-X  not valid
	 * 13: X-2-X-X  not valid
	 * 14: X-X-3-X  not valid
	 * 15: X-X-X-4  not valid
	 * 16: X-X-X-X  not valid
	 **********************************************/
	 if (pAtom1 == NULL && pAtom4 == NULL) {
		//  4: X-2-3-X  not valid            
		// 13: X-2-X-X  not valid
		// 14: X-X-3-X  not valid
		// 16: X-X-X-X  not valid
		return false;
	 } else if (pAtom1 == NULL) {
	 	//  4: X-2-3-X  not valid            
	 	//  5: X-2-X-4  not valid      
	 	//  6: X-X-3-4  valid
	 	// 10: X-2-3-4  valid
		// 13: X-2-X-X  not valid
		// 14: X-X-3-X  not valid
		// 15: X-X-X-4  not valid
		// 16: X-X-X-X  not valid
		 if (pAtom3 == NULL || pAtom4 == NULL) {
			//  4: X-2-3-X  not valid            
			//  5: X-2-X-4  not valid      
			// 13: X-2-X-X  not valid
			// 14: X-X-3-X  not valid
			// 15: X-X-X-4  not valid
			// 16: X-X-X-X  not valid
			 return false;
		 }
	 } else if (pAtom4 == NULL) {
	 	//  1: 1-2-X-X  valid not improper   
	 	//  2: 1-X-3-X  valid improper       
	 	//  4: X-2-3-X  not valid            
	 	//  7: 1-2-3-X  valid
		// 12: 1-X-X-X  not valid
		// 13: X-2-X-X  not valid
		// 14: X-X-3-X  not valid
		// 16: X-X-X-X  not valid
		 if (improperFlag) {
			 if (pAtom1 == NULL || pAtom3 == NULL) {
				//  1: 1-2-X-X  valid not improper   
				//  4: X-2-3-X  not valid            
				// 12: 1-X-X-X  not valid
				// 13: X-2-X-X  not valid
				// 14: X-X-3-X  not valid
				// 16: X-X-X-X  not valid
				 return false;
			 }
		 } else {
			 if (pAtom1 == NULL || pAtom2 == NULL) {
				//  2: 1-X-3-X  valid improper       
				//  4: X-2-3-X  not valid            
				// 12: 1-X-X-X  not valid
				// 13: X-2-X-X  not valid
				// 14: X-X-3-X  not valid
				// 16: X-X-X-X  not valid
				 return false;
			 }
		 }
	 } else if (pAtom2 == NULL && pAtom3 == NULL) {
		//  3: 1-X-X-4  not valid            
		// 12: 1-X-X-X  not valid
		// 15: X-X-X-4  not valid
		// 16: X-X-X-X  not valid
		 return false;
	 } else if (improperFlag && pAtom3 == NULL) {
		//  1: 1-2-X-X  valid not improper   
		//  3: 1-X-X-4  not valid            
		//  5: X-2-X-4  not valid      
		//  8: 1-2-X-4  valid not improper
		// 12: 1-X-X-X  not valid
		// 13: X-2-X-X  not valid
		// 15: X-X-X-4  not valid
		// 16: X-X-X-X  not valid
		 return false;
	 }
	 return true;
}
inline void IcEntry::blankIC() {
	if (pAtom1 != NULL) {
		// remove from the old atom 1 the pointer to the IC
		pAtom1->removeIcEntry(this);
	}
	if (pAtom2 != NULL) {
		pAtom2->removeIcEntry(this);
	}
	if (pAtom3 != NULL) {
		pAtom3->removeIcEntry(this);
	}
	if (pAtom4 != NULL) {
		pAtom4->removeIcEntry(this);
	}
	vals = std::vector<double>(5, 0.0);
}
}

#endif

