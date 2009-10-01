/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

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

using namespace std;

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

class CharmmTopologyResidue {
	public:
		CharmmTopologyResidue();
		CharmmTopologyResidue(string _name, bool _isPatch=false, double _charge=0.0, string _firstPatch="", string _lastPatch="");
		CharmmTopologyResidue(const CharmmTopologyResidue & _res);
		~CharmmTopologyResidue();

		void operator=(const CharmmTopologyResidue & _res);

		void setName(const string & _name);
		string getName() const;

		void setIsPatch(bool _flag);
		bool getIsPatch() const;

		void setFirstDefaultPatch(const string & _patch);
		string getFirstDefaultPatch() const;

		void setLastDefaultPatch(const string & _patch);
		string getLastDefaultPatch() const;

		void setCharge(const double & _charge);
		double getCharge() const;

		bool applyPatch(const CharmmTopologyResidue & _patch);

		void addTopolAtom(string _name, string _atomType, double _partialCharge, string _element, int group);
		string toString() const;

		unsigned int atomSize() const;
		void getTopolAtom(unsigned int _i, string & _name, string & _type, double & _partialCharge, string & _element, int & _group);
		bool deleteTopolAtom(string _name);
		void clearAllTopolAtoms();
		vector<string> getAllTopoAtomNames();

		void addIcLine(string _atom1, string _atom2, string _atom3, string _atom4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag);
		void getIcLine(unsigned int _i, vector<string> & _names, vector<double> & _values, bool & _improperFlag);
		unsigned int icSize() const;
		bool deleteIcLine(string _atom);
		bool deleteIcLine(string _atom1, string _atom2, string _atom3, string _atom4);
		void clearAllIcLines();

		void addBond(string _atom1, string _atom2, unsigned int _type=1);
		unsigned int bondSize() const;
		void getBond(unsigned int _i, string & _atom1, string & atom2, unsigned int & _type);
		bool deleteBond(string _atom);
		bool deleteBond(string _atom1, string _atom2);
		void clearAllBonds();

		void addAngle(string _atom1, string _atom2, string _atom3);
		unsigned int angleSize() const;
		void getAngle(unsigned int _i, string & _atom1, string & _atom2, string & _atom3);
		bool deleteAngle(string _atom);
		bool deleteAngle(string _atom1, string _atom2, string _atom3);
		void clearAllAngles();

		void addDihedral(string _atom1, string _atom2, string _atom3, string _atom4);
		unsigned int dihedralSize() const;
		void getDihedral(unsigned int _i, string & _atom1, string & _atom2, string & _atom3, string & _atom4);
		bool deleteDihedral(string _atom);
		bool deleteDihedral(string _atom1, string _atom2, string _atom3, string _atom4);
		void clearAllDihedrals();

		void addImproper(string _atom1, string _atom2, string _atom3, string _atom4);
		unsigned int improperSize() const;
		void getImproper(unsigned int _i, string & _atom1, string & _atom2, string & _atom3, string & _atom4);
		bool deleteImproper(string _atom);
		bool deleteImproper(string _atom1, string _atom2, string _atom3, string _atom4);
		void clearAllImpropers();

		void addDonor(string _hydrogen, string _heavy, string _antecendent1="", string _antecendent2="");
		unsigned int donorSize() const;
		void getDonor(unsigned int _i, string & _hydrogen, string & _heavy, string & _antecendent1, string & _antecendent2);
		bool deleteDonor(string _atom);
		bool deleteDonor(string _hydrogen, string _heavy, string _antecendent1="", string _antecendent2="");
		void clearAllDonors();

		void addAcceptor(string _acceptor, string _bonded="", string _bonded2="");
		unsigned int acceptorSize() const;
		void getAcceptor(unsigned int _i, string & _acceptor, string & _bonded, string & _bonded2);
	//	bool deleteAcceptor(string _atom);
		bool deleteAcceptor(string _acceptor, string _bonded="", string _bonded2="");
		void clearAllAcceptors();

		void addDelete(string _type, vector<string> atoms);
		unsigned int deleteSize() const;
		void getDelete(unsigned int _i, string & _type, vector<string> & _atoms);
		void clearAllDelete();

		void setPreviousTopologyResidue(CharmmTopologyResidue * _pPrevTopologyRes);
		CharmmTopologyResidue * getPreviousTopologyResidue() const;
		void setNextTopologyResidue(CharmmTopologyResidue * _pNextTopologyRes);
		CharmmTopologyResidue * getNextTopologyResidue() const;

		void linkAtoms();
		void unlinkAtoms();
		
	private:
		void setup(string _name, bool _isPatch, double _charge, string _firstPatch, string _lastPatch);
		void copy(const CharmmTopologyResidue & _res);
		void deletePointers();
		void reset();
		void sortAtomsByGroup();

		CharmmTopologyResidue * pPrevTopologyRes;
		CharmmTopologyResidue * pNextTopologyRes;


		string name;
		bool isPatch;
		double charge;
		string firstPatch;
		string lastPatch;

		struct TopolAtom {
			string name;
			string type;
			double partialCharge;
			int group;
			string element;
		};
		vector<TopolAtom*> atoms;
		map<string, TopolAtom*> atomMap;

		struct IcLine {
			vector<string> atoms;
			vector<double> values;
			bool improperFlag;
			vector<TopolAtom*> pAtoms;
		};
		vector<IcLine> IcTable;

		struct Bond {
			vector<string> atoms;
			unsigned int type; // 1 single, 2 double, 3 triple bond
			vector<TopolAtom*> pAtoms;
		};
		vector<Bond> bonds;

		struct Angle {
			vector<string> atoms;
			vector<TopolAtom*> pAtoms;
		};
		vector<Angle> angles;

		struct Dihedral {
			vector<string> atoms;
			vector<TopolAtom*> pAtoms;
			
		};
		vector<Dihedral> dihedrals;

		struct Improper {
			vector<string> atoms;
			vector<TopolAtom*> pAtoms;
		};
		vector<Improper> impropers;

		struct Donor {
			vector<string> atoms;
			vector<TopolAtom*> pAtoms;
		};
		vector<Donor> donors;

		struct Acceptor {
			vector<string> atoms;
			vector<TopolAtom*> pAtoms;
		};
		vector<Acceptor> acceptors;

		struct Delete {
			string type;
			vector<string> atoms;
		};
		vector<Delete> deletes;


};

inline void CharmmTopologyResidue::setName(const string & _name) {name = _name;}
inline string CharmmTopologyResidue::getName() const {return name;}
inline void CharmmTopologyResidue::setIsPatch(bool _flag) {isPatch = _flag;}
inline bool CharmmTopologyResidue::getIsPatch() const {return isPatch;}
inline void CharmmTopologyResidue::setFirstDefaultPatch(const string & _patch) {firstPatch = _patch;}
inline string CharmmTopologyResidue::getFirstDefaultPatch() const {return firstPatch;}
inline void CharmmTopologyResidue::setLastDefaultPatch(const string & _patch) {lastPatch = _patch;}
inline string CharmmTopologyResidue::getLastDefaultPatch() const {return lastPatch;}
inline void CharmmTopologyResidue::setCharge(const double & _charge) {charge = _charge;}
inline double CharmmTopologyResidue::getCharge() const {return charge;}

inline void CharmmTopologyResidue::addTopolAtom(string _name, string _atomType, double _partialCharge, string _element, int _group) {TopolAtom * tmp = new TopolAtom; tmp->name = _name; tmp->type = _atomType; tmp->partialCharge = _partialCharge; tmp->element = _element; tmp->group = _group; atoms.push_back(tmp); atomMap[_name] = atoms.back(); sortAtomsByGroup();}
inline unsigned int CharmmTopologyResidue::atomSize() const {return atoms.size();}
inline void CharmmTopologyResidue::getTopolAtom(unsigned int _i, string & _name, string & _type, double & _partialCharge, string & _element, int & _group) {_name = (atoms[_i])->name; _type = (atoms[_i])->type; _partialCharge = (atoms[_i])->partialCharge; _element = (atoms[_i])->element; _group = (atoms[_i])->group;}
inline void CharmmTopologyResidue::clearAllTopolAtoms() {deletePointers();}
inline bool CharmmTopologyResidue::deleteTopolAtom(string _name) {for (vector<TopolAtom*>::iterator k=atoms.begin(); k!=atoms.end(); k++) {if ((*k)->name == _name) {delete *k; atoms.erase(k); deleteIcLine(_name); deleteBond(_name); deleteAngle(_name); deleteDihedral(_name); deleteImproper(_name); deleteDonor(_name); deleteAcceptor(_name); sortAtomsByGroup(); return true;}} return false;}

inline void CharmmTopologyResidue::addIcLine(string _atom1, string _atom2, string _atom3, string _atom4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag) {IcLine tmp; tmp.atoms.push_back(_atom1); tmp.atoms.push_back(_atom2); tmp.atoms.push_back(_atom3); tmp.atoms.push_back(_atom4); tmp.values.push_back(_d1); tmp.values.push_back(_a1); tmp.values.push_back(_dihe); tmp.values.push_back(_a2); tmp.values.push_back(_d2); tmp.improperFlag = _improperFlag; IcTable.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::icSize() const {return IcTable.size();}
inline void CharmmTopologyResidue::getIcLine(unsigned int _i, vector<string> & _atoms, vector<double> & _values, bool & _improperFlag) {_atoms = IcTable[_i].atoms; _values = IcTable[_i].values; _improperFlag = IcTable[_i].improperFlag;}
inline void CharmmTopologyResidue::clearAllIcLines() {IcTable.clear();}
inline bool CharmmTopologyResidue::deleteIcLine(string _atom1, string _atom2, string _atom3, string _atom4) {for (vector<IcLine>::iterator k=IcTable.begin(); k!=IcTable.end(); k++) {if (((*k).atoms[0] == _atom1 && (*k).atoms[1] == _atom2 && (*k).atoms[2] == _atom3 && (*k).atoms[3] == _atom4) || ((*k).atoms[0] == _atom4 && (*k).atoms[1] == _atom3 && (*k).atoms[2] == _atom2 && (*k).atoms[3] == _atom1)) {IcTable.erase(k); return true;}} return false;}
inline bool CharmmTopologyResidue::deleteIcLine(string _atom) {for (vector<IcLine>::iterator k=IcTable.begin(); k!=IcTable.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom || (*k).atoms[2] == _atom || (*k).atoms[3] == _atom) {IcTable.erase(k); return true;}} return false;}

inline void CharmmTopologyResidue::addBond(string _atom1, string _atom2, unsigned int _type) {Bond tmp; tmp.atoms = vector<string>(2, ""); tmp.atoms[0] = _atom1; tmp.atoms[1] = _atom2; tmp.type = _type; bonds.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::bondSize() const {return bonds.size();}
inline void CharmmTopologyResidue::getBond(unsigned int _i, string & _atom1, string & _atom2, unsigned int & _type) {_atom1 = bonds[_i].atoms[0]; _atom2 = bonds[_i].atoms[1]; _type = bonds[_i].type;}
inline void CharmmTopologyResidue::clearAllBonds() {bonds.clear();}
inline bool CharmmTopologyResidue::deleteBond(string _atom1, string _atom2) {for (vector<Bond>::iterator k=bonds.begin(); k!=bonds.end(); k++) {if (((*k).atoms[0] == _atom1 && (*k).atoms[1] == _atom2) || ((*k).atoms[0] == _atom2 && (*k).atoms[1] == _atom1)) {bonds.erase(k); return true;}} return false;}
inline bool CharmmTopologyResidue::deleteBond(string _atom) {for (vector<Bond>::iterator k=bonds.begin(); k!=bonds.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom) {bonds.erase(k); return true;}} return false;}

inline void CharmmTopologyResidue::addAngle(string _atom1, string _atom2, string _atom3) {Angle tmp; tmp.atoms = vector<string>(3, ""); tmp.atoms[0] = _atom1; tmp.atoms[1] = _atom2; tmp.atoms[2] = _atom3; angles.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::angleSize() const {return angles.size();}
inline void CharmmTopologyResidue::getAngle(unsigned int _i, string & _atom1, string & _atom2, string & _atom3) {_atom1 = angles[_i].atoms[0]; _atom2 = angles[_i].atoms[1]; _atom3 = angles[_i].atoms[2];}
inline void CharmmTopologyResidue::clearAllAngles() {angles.clear();}
inline bool CharmmTopologyResidue::deleteAngle(string _atom1, string _atom2, string _atom3) {for (vector<Angle>::iterator k=angles.begin(); k!=angles.end(); k++) {if (((*k).atoms[0] == _atom1 && (*k).atoms[1] == _atom2 && (*k).atoms[2] == _atom3) || ((*k).atoms[0] == _atom3 && (*k).atoms[1] == _atom2 && (*k).atoms[2] == _atom1)) {angles.erase(k); return true;}} return false;}
inline bool CharmmTopologyResidue::deleteAngle(string _atom) {for (vector<Angle>::iterator k=angles.begin(); k!=angles.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom || (*k).atoms[2] == _atom) {angles.erase(k); return true;}} return false;}

inline void CharmmTopologyResidue::addDihedral(string _atom1, string _atom2, string _atom3, string _atom4) {Dihedral tmp; tmp.atoms = vector<string>(4, ""); tmp.atoms[0] = _atom1; tmp.atoms[1] = _atom2; tmp.atoms[2] = _atom3; tmp.atoms[3] = _atom4; dihedrals.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::dihedralSize() const {return dihedrals.size();}
inline void CharmmTopologyResidue::getDihedral(unsigned int _i, string & _atom1, string & _atom2, string & _atom3, string & _atom4) {_atom1 = dihedrals[_i].atoms[0]; _atom2 = dihedrals[_i].atoms[1]; _atom3 = dihedrals[_i].atoms[2]; _atom4 = dihedrals[_i].atoms[3];}
inline void CharmmTopologyResidue::clearAllDihedrals() {dihedrals.clear();}
inline bool CharmmTopologyResidue::deleteDihedral(string _atom1, string _atom2, string _atom3, string _atom4) {for (vector<Dihedral>::iterator k=dihedrals.begin(); k!=dihedrals.end(); k++) {if (((*k).atoms[0] == _atom1 && (*k).atoms[1] == _atom2 && (*k).atoms[2] == _atom3 && (*k).atoms[3] == _atom4) || ((*k).atoms[0] == _atom4 && (*k).atoms[1] == _atom3 && (*k).atoms[2] == _atom2 && (*k).atoms[3] == _atom1)) {dihedrals.erase(k); return true;}} return false;}
inline bool CharmmTopologyResidue::deleteDihedral(string _atom) {for (vector<Dihedral>::iterator k=dihedrals.begin(); k!=dihedrals.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom || (*k).atoms[2] == _atom || (*k).atoms[3] == _atom) {dihedrals.erase(k); return true;}} return false;}

inline void CharmmTopologyResidue::addImproper(string _atom1, string _atom2, string _atom3, string _atom4) {Improper tmp; tmp.atoms = vector<string>(4, ""); tmp.atoms[0] = _atom1; tmp.atoms[1] = _atom2; tmp.atoms[2] = _atom3; tmp.atoms[3] = _atom4; impropers.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::improperSize() const {return impropers.size();}
inline void CharmmTopologyResidue::getImproper(unsigned int _i, string & _atom1, string & _atom2, string & _atom3, string & _atom4) {_atom1 = impropers[_i].atoms[0]; _atom2 = impropers[_i].atoms[1]; _atom3 = impropers[_i].atoms[2]; _atom4 = impropers[_i].atoms[3];}
inline void CharmmTopologyResidue::clearAllImpropers() {impropers.clear();}
inline bool CharmmTopologyResidue::deleteImproper(string _atom1, string _atom2, string _atom3, string _atom4) {for (vector<Improper>::iterator k=impropers.begin(); k!=impropers.end(); k++) {if (((*k).atoms[0] == _atom1 && (*k).atoms[1] == _atom2 && (*k).atoms[2] == _atom3 && (*k).atoms[3] == _atom4) || ((*k).atoms[0] == _atom4 && (*k).atoms[1] == _atom3 && (*k).atoms[2] == _atom2 && (*k).atoms[3] == _atom1)) {impropers.erase(k); return true;}} return false;}
inline bool CharmmTopologyResidue::deleteImproper(string _atom) {for (vector<Improper>::iterator k=impropers.begin(); k!=impropers.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom || (*k).atoms[2] == _atom || (*k).atoms[3] == _atom) {impropers.erase(k); return true;}} return false;}

inline void CharmmTopologyResidue::addDonor(string _hydrogen, string _heavy, string _antecendent1, string _antecendent2) {Donor tmp; tmp.atoms = vector<string>(4, ""); tmp.atoms[0] = _hydrogen; tmp.atoms[1] = _heavy; tmp.atoms[1] = _antecendent1;   tmp.atoms[1] = _antecendent2; donors.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::donorSize() const {return donors.size();}
inline void CharmmTopologyResidue::getDonor(unsigned int _i, string & _hydrogen, string & _heavy, string & _antecendent1, string & _antecendent2) {_hydrogen = donors[_i].atoms[0]; _heavy = donors[_i].atoms[1]; _antecendent1 = donors[_i].atoms[2]; _antecendent2 = donors[_i].atoms[3];}
inline void CharmmTopologyResidue::clearAllDonors() {donors.clear();}
inline bool CharmmTopologyResidue::deleteDonor(string _hydrogen, string _heavy, string _antecendent1, string _antecendent2) {for (vector<Donor>::iterator k=donors.begin(); k!=donors.end(); k++) {if ((*k).atoms[0] == _hydrogen && (*k).atoms[1] == _heavy && ((*k).atoms[2] == _antecendent1 ||  _antecendent1 == "") && ((*k).atoms[3] == _antecendent2 || _antecendent2 == "")) {donors.erase(k); return true;}} return false;}
inline bool CharmmTopologyResidue::deleteDonor(string _atom) {for (vector<Donor>::iterator k=donors.begin(); k!=donors.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom || (*k).atoms[2] == _atom || (*k).atoms[3] == _atom) {donors.erase(k); return true;}} return false;}

inline void CharmmTopologyResidue::addAcceptor(string _acceptor, string _bonded, string _bonded2) {Acceptor tmp; tmp.atoms = vector<string>(3, ""); tmp.atoms[0] = _acceptor; tmp.atoms[1] = _bonded; tmp.atoms[2] = _bonded2; acceptors.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::acceptorSize() const {return acceptors.size();}
inline void CharmmTopologyResidue::getAcceptor(unsigned int _i, string & _acceptor, string & _bonded, string & _bonded2) {_acceptor = acceptors[_i].atoms[0]; _bonded = acceptors[_i].atoms[1]; _bonded2 = acceptors[_i].atoms[2];}
inline void CharmmTopologyResidue::clearAllAcceptors() {acceptors.clear();}
inline bool CharmmTopologyResidue::deleteAcceptor(string _acceptor, string _bonded, string _bonded2) {
	for (vector<Acceptor>::iterator k=acceptors.begin(); k!=acceptors.end(); k++) {
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
//inline bool CharmmTopologyResidue::deleteAcceptor(string _atom) {for (vector<Acceptor>::iterator k=acceptors.begin(); k!=acceptors.end(); k++) {if ((*k).atoms[0] == _atom || (*k).atoms[1] == _atom) {acceptors.erase(k); return true;}} return false;}
inline void CharmmTopologyResidue::addDelete(string _type, vector<string> _atoms) {Delete tmp; tmp.type = _type; tmp.atoms = _atoms; deletes.push_back(tmp);}
inline unsigned int CharmmTopologyResidue::deleteSize() const {return deletes.size();}
inline void CharmmTopologyResidue::getDelete(unsigned int _i, string & _type, vector<string> & _atoms) {_type = deletes[_i].type; _atoms = deletes[_i].atoms;}
inline void CharmmTopologyResidue::clearAllDelete() {deletes.clear();}

inline void CharmmTopologyResidue::setPreviousTopologyResidue(CharmmTopologyResidue * _pPrevTopologyRes) {pPrevTopologyRes = _pPrevTopologyRes;}
inline CharmmTopologyResidue * CharmmTopologyResidue::getPreviousTopologyResidue() const {return pPrevTopologyRes;}
inline void CharmmTopologyResidue::setNextTopologyResidue(CharmmTopologyResidue * _pNextTopologyRes) {_pNextTopologyRes = _pNextTopologyRes;}
inline CharmmTopologyResidue * CharmmTopologyResidue::getNextTopologyResidue() const {return pNextTopologyRes;}


#endif
