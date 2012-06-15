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

#ifndef CHARMMTOPOLOGYRESIDUE_H
#define CHARMMTOPOLOGYRESIDUE_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <sys/types.h>


/********************************************************
 *  This class stores the information from a RESI or PATCH
 *  from a CHARMM topology file
 *
 *  It is filled by the CharmmTopologyReader and stored
 *  in the CharmmTopology object
 *
 *  RESI GLY          0.00
 *  GROUP   
 *  ATOM N    NH1    -0.47  !     |
 *  ATOM HN   H       0.31  !     N-H
 *  ATOM CA   CT2    -0.02  !     |  
 *  ATOM HA1  HB      0.09  !     |  
 *  ATOM HA2  HB      0.09  ! HA1-CA-HA2
 *  GROUP                   !     |  
 *  ATOM C    C       0.51  !     |  
 *  ATOM O    O      -0.51  !     C=O
 *  BOND N HN  N  CA  C CA   
 *  BOND C +N  CA HA1 CA HA2  
 *  DOUBLE O  C 
 *  IMPR N -C  CA HN  C CA   +N O   
 *  DONOR HN N   
 *  ACCEPTOR O C   
 *  IC -C   CA   *N   HN    1.3475 122.8200  180.0000 115.6200  0.9992
 *  IC -C   N    CA   C     1.3475 122.8200  180.0000 108.9400  1.4971
 *  IC N    CA   C    +N    1.4553 108.9400  180.0000 117.6000  1.3479
 *  IC +N   CA   *C   O     1.3479 117.6000  180.0000 120.8500  1.2289
 *  IC CA   C    +N   +CA   1.4971 117.6000  180.0000 124.0800  1.4560
 *  IC N    C    *CA  HA1   1.4553 108.9400  117.8600 108.0300  1.0814
 *  IC N    C    *CA  HA2   1.4553 108.9400 -118.1200 107.9500  1.0817
 *  PATCHING FIRS GLYP   
 *  
 ********************************************************/

namespace MSL { 
class CharmmTopologyResidue {
	public:
		CharmmTopologyResidue();
		CharmmTopologyResidue(std::string _name, bool _isPatch=false, double _charge=0.0, std::string _firstPatch="", std::string _lastPatch="");
		CharmmTopologyResidue(const CharmmTopologyResidue & _res);
		~CharmmTopologyResidue();

		void operator=(const CharmmTopologyResidue & _res);

		void setName(const std::string & _name);
		std::string getName() const;

		void setIsPatch(bool _flag);
		bool getIsPatch() const;

		void setFirstDefaultPatch(const std::string & _patch);
		std::string getFirstDefaultPatch() const;

		void setLastDefaultPatch(const std::string & _patch);
		std::string getLastDefaultPatch() const;

		void setCharge(const double & _charge);
		double getCharge() const;

		bool applyPatch(const CharmmTopologyResidue & _patch);

		void addTopolAtom(std::string _name, std::string _atomType, double _partialCharge, std::string _element, int group);
		std::string toString() const;

		unsigned int atomSize() const;
		void getTopolAtom(unsigned int _i, std::string & _name, std::string & _type, double & _partialCharge, std::string & _element, int & _group);
		bool deleteTopolAtom(std::string _name);
		void clearAllTopolAtoms();
		std::vector<std::string> getAllTopoAtomNames();

		void addIcLine(std::string _atom1, std::string _atom2, std::string _atom3, std::string _atom4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag);
		void getIcLine(unsigned int _i, std::vector<std::string> & _names, std::vector<double> & _values, bool & _improperFlag);
		unsigned int icSize() const;
		bool deleteIcLine(std::string _atom);
		bool deleteIcLine(std::string _atom1, std::string _atom2, std::string _atom3, std::string _atom4);
		void clearAllIcLines();

		void addBond(std::string _atom1, std::string _atom2, unsigned int _type=1);
		unsigned int bondSize() const;
		void getBond(unsigned int _i, std::string & _atom1, std::string & atom2, unsigned int & _type);
		bool deleteBond(std::string _atom);
		bool deleteBond(std::string _atom1, std::string _atom2);
		void clearAllBonds();

		void addAngle(std::string _atom1, std::string _atom2, std::string _atom3);
		unsigned int angleSize() const;
		void getAngle(unsigned int _i, std::string & _atom1, std::string & _atom2, std::string & _atom3);
		bool deleteAngle(std::string _atom);
		bool deleteAngle(std::string _atom1, std::string _atom2, std::string _atom3);
		void clearAllAngles();

		void addDihedral(std::string _atom1, std::string _atom2, std::string _atom3, std::string _atom4);
		unsigned int dihedralSize() const;
		void getDihedral(unsigned int _i, std::string & _atom1, std::string & _atom2, std::string & _atom3, std::string & _atom4);
		bool deleteDihedral(std::string _atom);
		bool deleteDihedral(std::string _atom1, std::string _atom2, std::string _atom3, std::string _atom4);
		void clearAllDihedrals();

		void addImproper(std::string _atom1, std::string _atom2, std::string _atom3, std::string _atom4);
		unsigned int improperSize() const;
		void getImproper(unsigned int _i, std::string & _atom1, std::string & _atom2, std::string & _atom3, std::string & _atom4);
		bool deleteImproper(std::string _atom);
		bool deleteImproper(std::string _atom1, std::string _atom2, std::string _atom3, std::string _atom4);
		void clearAllImpropers();

		void addDonor(std::string _hydrogen, std::string _heavy, std::string _antecendent1="", std::string _antecendent2="");
		unsigned int donorSize() const;
		void getDonor(unsigned int _i, std::string & _hydrogen, std::string & _heavy, std::string & _antecendent1, std::string & _antecendent2);
		bool deleteDonor(std::string _atom);
		bool deleteDonor(std::string _hydrogen, std::string _heavy, std::string _antecendent1="", std::string _antecendent2="");
		void clearAllDonors();

		void addAcceptor(std::string _acceptor, std::string _bonded="", std::string _bonded2="");
		unsigned int acceptorSize() const;
		void getAcceptor(unsigned int _i, std::string & _acceptor, std::string & _bonded, std::string & _bonded2);
	//	bool deleteAcceptor(std::string _atom);
		bool deleteAcceptor(std::string _acceptor, std::string _bonded="", std::string _bonded2="");
		void clearAllAcceptors();

		void addDelete(std::string _type, std::vector<std::string> atoms);
		unsigned int deleteSize() const;
		void getDelete(unsigned int _i, std::string & _type, std::vector<std::string> & _atoms);
		void clearAllDelete();

		void setPreviousTopologyResidue(CharmmTopologyResidue * _pPrevTopologyRes);
		CharmmTopologyResidue * getPreviousTopologyResidue() const;
		void setNextTopologyResidue(CharmmTopologyResidue * _pNextTopologyRes);
		CharmmTopologyResidue * getNextTopologyResidue() const;

		void linkAtoms();
		void unlinkAtoms();
		
	private:
		void setup(std::string _name, bool _isPatch, double _charge, std::string _firstPatch, std::string _lastPatch);
		void copy(const CharmmTopologyResidue & _res);
		void deletePointers();
		void reset();
		void sortAtomsByGroup();

		CharmmTopologyResidue * pPrevTopologyRes;
		CharmmTopologyResidue * pNextTopologyRes;


		std::string name;
		bool isPatch;
		double charge;
		std::string firstPatch;
		std::string lastPatch;

		struct TopolAtom {
			std::string name;
			std::string type;
			double partialCharge;
			int group;
			std::string element;
		};
		std::vector<TopolAtom*> atoms;
		std::map<std::string, TopolAtom*> atomMap;

		struct IcLine {
			std::vector<std::string> atoms;
			std::vector<double> values;
			bool improperFlag;
			std::vector<TopolAtom*> pAtoms;
		};
		std::vector<IcLine> IcTable;

		struct Bond {
			std::vector<std::string> atoms;
			unsigned int type; // 1 single, 2 double, 3 triple bond
			std::vector<TopolAtom*> pAtoms;
		};
		std::vector<Bond> bonds;

		struct Angle {
			std::vector<std::string> atoms;
			std::vector<TopolAtom*> pAtoms;
		};
		std::vector<Angle> angles;

		struct Dihedral {
			std::vector<std::string> atoms;
			std::vector<TopolAtom*> pAtoms;
			
		};
		std::vector<Dihedral> dihedrals;

		struct Improper {
			std::vector<std::string> atoms;
			std::vector<TopolAtom*> pAtoms;
		};
		std::vector<Improper> impropers;

		struct Donor {
			std::vector<std::string> atoms;
			std::vector<TopolAtom*> pAtoms;
		};
		std::vector<Donor> donors;

		struct Acceptor {
			std::vector<std::string> atoms;
			std::vector<TopolAtom*> pAtoms;
		};
		std::vector<Acceptor> acceptors;

		struct Delete {
			std::string type;
			std::vector<std::string> atoms;
		};
		std::vector<Delete> deletes;


};

inline void CharmmTopologyResidue::setName(const std::string & _name) {name = _name;}
inline std::string CharmmTopologyResidue::getName() const {return name;}
inline void CharmmTopologyResidue::setIsPatch(bool _flag) {isPatch = _flag;}
inline bool CharmmTopologyResidue::getIsPatch() const {return isPatch;}
inline void CharmmTopologyResidue::setFirstDefaultPatch(const std::string & _patch) {firstPatch = _patch;}
inline std::string CharmmTopologyResidue::getFirstDefaultPatch() const {return firstPatch;}
inline void CharmmTopologyResidue::setLastDefaultPatch(const std::string & _patch) {lastPatch = _patch;}
inline std::string CharmmTopologyResidue::getLastDefaultPatch() const {return lastPatch;}
inline void CharmmTopologyResidue::setCharge(const double & _charge) {charge = _charge;}
inline double CharmmTopologyResidue::getCharge() const {return charge;}

inline void CharmmTopologyResidue::addTopolAtom(std::string _name, std::string _atomType, double _partialCharge, std::string _element, int _group) {TopolAtom * tmp = new TopolAtom; tmp->name = _name; tmp->type = _atomType; tmp->partialCharge = _partialCharge; tmp->element = _element; tmp->group = _group; atoms.push_back(tmp); atomMap[_name] = atoms.back(); sortAtomsByGroup();}
inline unsigned int CharmmTopologyResidue::atomSize() const {return atoms.size();}
inline void CharmmTopologyResidue::getTopolAtom(unsigned int _i, std::string & _name, std::string & _type, double & _partialCharge, std::string & _element, int & _group) {_name = (atoms[_i])->name; _type = (atoms[_i])->type; _partialCharge = (atoms[_i])->partialCharge; _element = (atoms[_i])->element; _group = (atoms[_i])->group;}
inline void CharmmTopologyResidue::clearAllTopolAtoms() {deletePointers();}
inline bool CharmmTopologyResidue::deleteTopolAtom(std::string _name) {for (std::vector<TopolAtom*>::iterator k=atoms.begin(); k!=atoms.end(); k++) {if ((*k)->name == _name) {delete *k; atoms.erase(k); deleteIcLine(_name); deleteBond(_name); deleteAngle(_name); deleteDihedral(_name); deleteImproper(_name); deleteDonor(_name); deleteAcceptor(_name); sortAtomsByGroup(); return true;}} return false;}

inline void CharmmTopologyResidue::addIcLine(std::string _atom1, std::string _atom2, std::string _atom3, std::string _atom4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag) {IcLine tmp; tmp.atoms.push_back(_atom1); tmp.atoms.push_back(_atom2); tmp.atoms.push_back(_atom3); tmp.atoms.push_back(_atom4); tmp.values.push_back(_d1); tmp.values.push_back(_a1); tmp.values.push_back(_dihe); tmp.values.push_back(_a2); tmp.values.push_back(_d2); tmp.improperFlag = _improperFlag; IcTable.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::icSize() const {return IcTable.size();}
inline void CharmmTopologyResidue::getIcLine(unsigned int _i, std::vector<std::string> & _atoms, std::vector<double> & _values, bool & _improperFlag) {_atoms = IcTable[_i].atoms; _values = IcTable[_i].values; _improperFlag = IcTable[_i].improperFlag;}
inline void CharmmTopologyResidue::clearAllIcLines() {IcTable.clear();}
inline bool CharmmTopologyResidue::deleteIcLine(std::string _atom1, std::string _atom2, std::string _atom3, std::string _atom4) {for (std::vector<IcLine>::iterator k=IcTable.begin(); k!=IcTable.end(); k++) {if (((*k).atoms[0] == _atom1 && (*k).atoms[1] == _atom2 && (*k).atoms[2] == _atom3 && (*k).atoms[3] == _atom4) || ((*k).atoms[0] == _atom4 && (*k).atoms[1] == _atom3 && (*k).atoms[2] == _atom2 && (*k).atoms[3] == _atom1)) {IcTable.erase(k); return true;}} return false;}
inline bool CharmmTopologyResidue::deleteIcLine(std::string _atom) {for (std::vector<IcLine>::iterator k=IcTable.begin(); k!=IcTable.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom || (*k).atoms[2] == _atom || (*k).atoms[3] == _atom) {IcTable.erase(k); return true;}} return false;}

inline void CharmmTopologyResidue::addBond(std::string _atom1, std::string _atom2, unsigned int _type) {Bond tmp; tmp.atoms = std::vector<std::string>(2, ""); tmp.atoms[0] = _atom1; tmp.atoms[1] = _atom2; tmp.type = _type; bonds.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::bondSize() const {return bonds.size();}
inline void CharmmTopologyResidue::getBond(unsigned int _i, std::string & _atom1, std::string & _atom2, unsigned int & _type) {_atom1 = bonds[_i].atoms[0]; _atom2 = bonds[_i].atoms[1]; _type = bonds[_i].type;}
inline void CharmmTopologyResidue::clearAllBonds() {bonds.clear();}
inline bool CharmmTopologyResidue::deleteBond(std::string _atom1, std::string _atom2) {for (std::vector<Bond>::iterator k=bonds.begin(); k!=bonds.end(); k++) {if (((*k).atoms[0] == _atom1 && (*k).atoms[1] == _atom2) || ((*k).atoms[0] == _atom2 && (*k).atoms[1] == _atom1)) {bonds.erase(k); return true;}} return false;}
inline bool CharmmTopologyResidue::deleteBond(std::string _atom) {for (std::vector<Bond>::iterator k=bonds.begin(); k!=bonds.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom) {bonds.erase(k); return true;}} return false;}

inline void CharmmTopologyResidue::addAngle(std::string _atom1, std::string _atom2, std::string _atom3) {Angle tmp; tmp.atoms = std::vector<std::string>(3, ""); tmp.atoms[0] = _atom1; tmp.atoms[1] = _atom2; tmp.atoms[2] = _atom3; angles.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::angleSize() const {return angles.size();}
inline void CharmmTopologyResidue::getAngle(unsigned int _i, std::string & _atom1, std::string & _atom2, std::string & _atom3) {_atom1 = angles[_i].atoms[0]; _atom2 = angles[_i].atoms[1]; _atom3 = angles[_i].atoms[2];}
inline void CharmmTopologyResidue::clearAllAngles() {angles.clear();}
inline bool CharmmTopologyResidue::deleteAngle(std::string _atom1, std::string _atom2, std::string _atom3) {for (std::vector<Angle>::iterator k=angles.begin(); k!=angles.end(); k++) {if (((*k).atoms[0] == _atom1 && (*k).atoms[1] == _atom2 && (*k).atoms[2] == _atom3) || ((*k).atoms[0] == _atom3 && (*k).atoms[1] == _atom2 && (*k).atoms[2] == _atom1)) {angles.erase(k); return true;}} return false;}
inline bool CharmmTopologyResidue::deleteAngle(std::string _atom) {for (std::vector<Angle>::iterator k=angles.begin(); k!=angles.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom || (*k).atoms[2] == _atom) {angles.erase(k); return true;}} return false;}

inline void CharmmTopologyResidue::addDihedral(std::string _atom1, std::string _atom2, std::string _atom3, std::string _atom4) {Dihedral tmp; tmp.atoms = std::vector<std::string>(4, ""); tmp.atoms[0] = _atom1; tmp.atoms[1] = _atom2; tmp.atoms[2] = _atom3; tmp.atoms[3] = _atom4; dihedrals.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::dihedralSize() const {return dihedrals.size();}
inline void CharmmTopologyResidue::getDihedral(unsigned int _i, std::string & _atom1, std::string & _atom2, std::string & _atom3, std::string & _atom4) {_atom1 = dihedrals[_i].atoms[0]; _atom2 = dihedrals[_i].atoms[1]; _atom3 = dihedrals[_i].atoms[2]; _atom4 = dihedrals[_i].atoms[3];}
inline void CharmmTopologyResidue::clearAllDihedrals() {dihedrals.clear();}
inline bool CharmmTopologyResidue::deleteDihedral(std::string _atom1, std::string _atom2, std::string _atom3, std::string _atom4) {for (std::vector<Dihedral>::iterator k=dihedrals.begin(); k!=dihedrals.end(); k++) {if (((*k).atoms[0] == _atom1 && (*k).atoms[1] == _atom2 && (*k).atoms[2] == _atom3 && (*k).atoms[3] == _atom4) || ((*k).atoms[0] == _atom4 && (*k).atoms[1] == _atom3 && (*k).atoms[2] == _atom2 && (*k).atoms[3] == _atom1)) {dihedrals.erase(k); return true;}} return false;}
inline bool CharmmTopologyResidue::deleteDihedral(std::string _atom) {for (std::vector<Dihedral>::iterator k=dihedrals.begin(); k!=dihedrals.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom || (*k).atoms[2] == _atom || (*k).atoms[3] == _atom) {dihedrals.erase(k); return true;}} return false;}

inline void CharmmTopologyResidue::addImproper(std::string _atom1, std::string _atom2, std::string _atom3, std::string _atom4) {Improper tmp; tmp.atoms = std::vector<std::string>(4, ""); tmp.atoms[0] = _atom1; tmp.atoms[1] = _atom2; tmp.atoms[2] = _atom3; tmp.atoms[3] = _atom4; impropers.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::improperSize() const {return impropers.size();}
inline void CharmmTopologyResidue::getImproper(unsigned int _i, std::string & _atom1, std::string & _atom2, std::string & _atom3, std::string & _atom4) {_atom1 = impropers[_i].atoms[0]; _atom2 = impropers[_i].atoms[1]; _atom3 = impropers[_i].atoms[2]; _atom4 = impropers[_i].atoms[3];}
inline void CharmmTopologyResidue::clearAllImpropers() {impropers.clear();}
inline bool CharmmTopologyResidue::deleteImproper(std::string _atom1, std::string _atom2, std::string _atom3, std::string _atom4) {for (std::vector<Improper>::iterator k=impropers.begin(); k!=impropers.end(); k++) {if (((*k).atoms[0] == _atom1 && (*k).atoms[1] == _atom2 && (*k).atoms[2] == _atom3 && (*k).atoms[3] == _atom4) || ((*k).atoms[0] == _atom4 && (*k).atoms[1] == _atom3 && (*k).atoms[2] == _atom2 && (*k).atoms[3] == _atom1)) {impropers.erase(k); return true;}} return false;}
inline bool CharmmTopologyResidue::deleteImproper(std::string _atom) {for (std::vector<Improper>::iterator k=impropers.begin(); k!=impropers.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom || (*k).atoms[2] == _atom || (*k).atoms[3] == _atom) {impropers.erase(k); return true;}} return false;}

inline void CharmmTopologyResidue::addDonor(std::string _hydrogen, std::string _heavy, std::string _antecendent1, std::string _antecendent2) {Donor tmp; tmp.atoms = std::vector<std::string>(4, ""); tmp.atoms[0] = _hydrogen; tmp.atoms[1] = _heavy; tmp.atoms[1] = _antecendent1;   tmp.atoms[1] = _antecendent2; donors.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::donorSize() const {return donors.size();}
inline void CharmmTopologyResidue::getDonor(unsigned int _i, std::string & _hydrogen, std::string & _heavy, std::string & _antecendent1, std::string & _antecendent2) {_hydrogen = donors[_i].atoms[0]; _heavy = donors[_i].atoms[1]; _antecendent1 = donors[_i].atoms[2]; _antecendent2 = donors[_i].atoms[3];}
inline void CharmmTopologyResidue::clearAllDonors() {donors.clear();}
inline bool CharmmTopologyResidue::deleteDonor(std::string _hydrogen, std::string _heavy, std::string _antecendent1, std::string _antecendent2) {for (std::vector<Donor>::iterator k=donors.begin(); k!=donors.end(); k++) {if ((*k).atoms[0] == _hydrogen && (*k).atoms[1] == _heavy && ((*k).atoms[2] == _antecendent1 ||  _antecendent1 == "") && ((*k).atoms[3] == _antecendent2 || _antecendent2 == "")) {donors.erase(k); return true;}} return false;}
inline bool CharmmTopologyResidue::deleteDonor(std::string _atom) {for (std::vector<Donor>::iterator k=donors.begin(); k!=donors.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom || (*k).atoms[2] == _atom || (*k).atoms[3] == _atom) {donors.erase(k); return true;}} return false;}

inline void CharmmTopologyResidue::addAcceptor(std::string _acceptor, std::string _bonded, std::string _bonded2) {Acceptor tmp; tmp.atoms = std::vector<std::string>(3, ""); tmp.atoms[0] = _acceptor; tmp.atoms[1] = _bonded; tmp.atoms[2] = _bonded2; acceptors.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::acceptorSize() const {return acceptors.size();}
inline void CharmmTopologyResidue::getAcceptor(unsigned int _i, std::string & _acceptor, std::string & _bonded, std::string & _bonded2) {_acceptor = acceptors[_i].atoms[0]; _bonded = acceptors[_i].atoms[1]; _bonded2 = acceptors[_i].atoms[2];}
inline void CharmmTopologyResidue::clearAllAcceptors() {acceptors.clear();}
inline bool CharmmTopologyResidue::deleteAcceptor(std::string _acceptor, std::string _bonded, std::string _bonded2) {
	for (std::vector<Acceptor>::iterator k=acceptors.begin(); k!=acceptors.end(); k++) {
		if ((_bonded == "" && _bonded2 == "") && ((*k).atoms[0] == _acceptor || (*k).atoms[1] == _acceptor)) {
			acceptors.erase(k);
			return true;
		} else if ((*k).atoms[0] == _acceptor && ((*k).atoms[1] == _bonded || _bonded == "") && ((*k).atoms[2] == _bonded2 || _bonded2 == "")) {
			acceptors.erase(k);
			return true;
		}
	}
	return false;
}
//inline bool CharmmTopologyResidue::deleteAcceptor(std::string _atom) {for (std::vector<Acceptor>::iterator k=acceptors.begin(); k!=acceptors.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom) {acceptors.erase(k); return true;}} return false;}
inline void CharmmTopologyResidue::addDelete(std::string _type, std::vector<std::string> _atoms) {Delete tmp; tmp.type = _type; tmp.atoms = _atoms; deletes.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::deleteSize() const {return deletes.size();}
inline void CharmmTopologyResidue::getDelete(unsigned int _i, std::string & _type, std::vector<std::string> & _atoms) {_type = deletes[_i].type; _atoms = deletes[_i].atoms;}
inline void CharmmTopologyResidue::clearAllDelete() {deletes.clear();}

inline void CharmmTopologyResidue::setPreviousTopologyResidue(CharmmTopologyResidue * _pPrevTopologyRes) {pPrevTopologyRes = _pPrevTopologyRes;}
inline CharmmTopologyResidue * CharmmTopologyResidue::getPreviousTopologyResidue() const {return pPrevTopologyRes;}
inline void CharmmTopologyResidue::setNextTopologyResidue(CharmmTopologyResidue * _pNextTopologyRes) {_pNextTopologyRes = _pNextTopologyRes;}
inline CharmmTopologyResidue * CharmmTopologyResidue::getNextTopologyResidue() const {return pNextTopologyRes;}


}

#endif
