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

#ifndef ATOM_H
#define ATOM_H


// STL Includes
#include <iostream>
#include <vector>
#include <string>


// MSL Includes
#include "Real.h"
#include "Selectable.h"
#include "CartesianPoint.h"
#include "CartesianGeometry.h"


// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#endif


// Forward Declarations
class IcEntry;
class AtomGroup;
class AtomContainer;
class Residue;
class Position;
class Chain;
class System;



// Namespaces
using namespace std;



/*********************************************************
 *   NOTES
 *
 *   TODO:
 *    - add a "clipboard" to save coordinates
 *        - save only the current state or take a snapshot
 *          of all conformation?  Or allow to store multiple
 *          shapshots, indipendently on the number of alt conf,
 *          perhaps saving them with a name?
 *********************************************************/


class Atom : public Selectable<Atom> {
	public:
		Atom();
		Atom(const string _name, string _element="");
		Atom(string _name, Real _x, Real _y, Real _z, string _element=""); 
		Atom(string name, const CartesianPoint & _p, string _element="");
		Atom(const Atom & _atom);
		virtual ~Atom();

		void operator=(const Atom & _atom); // assignment

		void setName(string _name);
		string getName() const;
		void setResidueName(string _resname);
		string getResidueName() const;
		void setResidueNumber(int _resnum);
		int getResidueNumber() const;
		void setResidueIcode(string _icode);
		string getResidueIcode() const;
		void setChainId(string _chainId);
		string getChainId() const;
		void setElement(string _element);
		string getElement() const;
		void setType(string _type);
		string & getType();
		void setCharge(double _charge);
		double & getCharge();
		void setGroupNumber(unsigned int _groupNumber);
		unsigned int getGroupNumber() const;
		unsigned int getIdentityIndex(); // return the index of its parent identity in the position
		bool isInAlternativeIdentity(Atom * _pAtom) const; // checks if the two atoms happen to be in the same position but different residue types (cannot coexist)

		void setTempFactor(double _bfactor);
		double getTempFactor() const;

		void setSegID(string _segid);
		string getSegID() const;

		void setNameSpace(string _nameSpace);
		string getNameSpace() const;

		void setParentGroup(AtomGroup * _group);
		AtomGroup * getParentGroup() const;
		void setParentContainer(AtomContainer * _container);
		AtomContainer * getParentContainer() const;
		Residue * getParentResidue() const;
		Position * getParentPosition() const;
		Chain * getParentChain() const;
		System * getParentSystem() const;


		// set and get the coordinates
		void setCoor(CartesianPoint _p);
		void setCoor(Real _x, Real _y, Real _z);
		CartesianPoint & getCoor();
		vector<CartesianPoint *> & getAllCoor();
		Real getX() const;
		Real getY() const;
		Real getZ() const;
		Real operator[](size_t _n); // return X Y Z as atom[0], [1], [2] operators

		// print atom information
		string toString() const;
		friend ostream & operator<<(ostream &_os, const Atom & _atom)  {_os << _atom.toString(); return _os;};

		// measure relatiohships
		double distance(const Atom & _atom) const;
		double distance2(const Atom & _atom) const;
		double angle(const Atom & _atom) const;
		double angle(const Atom & _center, const Atom & _third) const;
		double angleRadians(const Atom & _atom) const;
		double angleRadians(const Atom & _center, const Atom & _third) const;
		double dihedral(const Atom & _second, const Atom & _third, const Atom & _fourth) const;
		double dihedralRadians(const Atom & _second, const Atom & _third, const Atom & _fourth) const;


		// virtual functions from Selectable class, allows selection of atoms
		virtual void addSelectableFunctions(); 

		/***************************************************
		 *  Building from internal coordinates
		 ***************************************************/
		bool hasCoor() const;
		void wipeCoordinates();
		void setHasCoordinates();
		void addIcEntry(IcEntry * _ic); // add an ic entry
		bool buildFromIc(bool _onlyFromActive=true); // try to build from the atom's ic entries (icEntries table)
		//bool buildFromIc(const map<IcEntry*, bool> & _exclude); // try to build from the atom's ic entries (icEntries table)
		bool buildFromIc(map<Atom*, bool> & _exclude, bool _onlyFromActive=true); // try to build from the atom's ic entries (icEntries table)
		void removeIcEntry(IcEntry * _ic); // remove the ic pointer with this address from icEntries
		vector<IcEntry*> & getIcEntries();

		/***************************************************
		 *  Alternate conformations
		 *
		 *  Alternative conformations are stored by adding
		 *  CartesianPoint ponters to the pCoorVec vector;
		 *
		 *  The current conformation is obtained using an
		 *  iterator (currentCoorIterator) that points to
		 *  the address of the current coordinates (NOTE:
		 *  the iterator is actually a pointer to a pointer
		 *  of CartesianPoint)
		 *
		 *  Using an interator for the current coordinates,
		 *  instead of a CartesianPoint pointer, allow to
		 *  keep track of the current coordinate without
		 *  having another variable (but makes understanding
		 *  the code a bit harder).
		 *
		 *  The current cartesian point is
		 *  *(*currentCoorIterator)
		 *
		 *  IMPORTANT NOTE: iterators are no longer valid after
		 *  an internal resize of the vector.  Always repoint
		 *  the iterators when element are inserted or deleted
		 *
		 ***************************************************/
		void setActiveConformation(size_t _i);
		unsigned int getActiveConformation() const;
		unsigned int getNumberOfAltConformations() const;
		void addAltConformation(); //default, same as current conformation
		void addAltConformation(const CartesianPoint & _point);  // it is important to regenerate the iterator to the current coor in case the vector is resized
		void removeAltConformation(size_t _i); // remove all alternate conformations (keeping the first conformation
		void removeAllAltConformations(); // remove all alternate conformations, keeping the first conformation

		/***************************************************
		 * As atoms have alternate conformations, residues can
		 * have alternate identities.
		 *
		 * Atoms can become not active when they belonging to 
		 * and identity that is not currently active (i.e. a
		 * VAL atom when LEU is active at the position)
		 ***************************************************/
		bool getActive() const;

		/***************************************************
		 *  Saved sets:
		 *
		 *  coordinates can be saved to named buffers, and copied back
		 *  from them
		 *
		 *  The difference between saved coord and alt coor is that
		 *  the saved coord are never active, they can only be
		 *  used to store and copy back coordinates
		 ***************************************************/
		void saveCoor(string _coordName);
		bool applySavedCoor(string _coordName);
		void clearSavedCoor();

		/***************************************************
		 *  Bonding information
		 *
		 *  The atom can store a list of pointers to atoms
		 *  it is bound to, is 1-3 or 1-4
		 *
		 *  - - - - - - - - - - - - - - - - - - - - - - - - - - 
		 *  |   I M P O R T A N T :   N E E D   T O   C O D E  |
		 *  |   M E C H A N I S M   F O R   D E L E T I N G    |
		 *  |   B O N D S                                      |
		 *  - - - - - - - - - - - - - - - - - - - - - - - - - - 
		 *
		 ***************************************************/
		//void setBondedTo(Atom * _pAtom, bool _bound=true);
		void setBoundTo(Atom * _pAtom);
		//void setOneThree(Atom * _pAtom, bool _bound=true);
		//void setOneFour(Atom * _pAtom, bool _bound=true);
		vector<vector<Atom*> > getBoundAtoms() const;
		bool isBoundTo(Atom * _pAtom) const;
		bool isOneThree(Atom * _pAtom) const;
		//bool isOneThree(Atom * _pAtom, const Atom * _14caller=NULL) const;
		bool isOneFour(Atom * _pAtom) const;
		vector<Atom*> getOneThreeMiddleAtoms(Atom * _pAtom) const; // returns a vector with the middle atoms for all 1-3 relationships with a given atom _pAtom (generally there is just one unless it is a 4 members ring)
		vector<vector<Atom*> > getOneFourMiddleAtoms(Atom * _pAtom) const; // returns a vector with the second and third atoms for all 1-4 relationships with a given atom _pAtom (can be multiple in 6 member rings)

	private:
		void setup(CartesianPoint _point, string _name, string _element);
		void copy(const Atom & _atom);
		void reset();
		void deletePointers();
		void updateResidueMap();
		void updateContainerMap();
		//void removeBonds();

		string name;
		string residueName;
		int residueNumber;
		string residueIcode;
		string chainId;
		string element;
		double charge;
		string type;
		unsigned int groupNumber;
		double tempFactor;
		string segId;
		// we can put the charge in the electrostatic interactions like the vdw radii etc
		
		string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd

		bool hasCoordinates;
		vector<IcEntry*> icEntries;

		//CartesianPoint * coor;
		
		vector<CartesianPoint*> pCoorVec;
		vector<CartesianPoint*>::iterator currentCoorIterator;

		map<string, CartesianPoint*> savedCoor;

		// pointer to parent electrostatic group
		AtomGroup * pParentGroup;
		// pointer to parent atom container
		AtomContainer * pParentContainer;

		map<Atom*, bool> bonds;
		// first level 1-2, second level 1-3, 3rd is 1-4
		map<Atom*, map<Atom*, map<Atom*, bool> > > boundAtoms;
		map<Atom*, map<Atom*, bool> > oneThreeAtoms; // oneThreeAtoms[X][Y] corresponds to this-Y-X
		map<Atom*, map<Atom*, map<Atom*, bool> > > oneFourAtoms; // oneFourAtoms[X][Y][Z] corresponds to this-Y-Z-X
		//map<Atom*, bool> oneThreeAtoms;
		//map<Atom*, bool> oneFourAtoms;


		// BOOST-RELATED FUNCTIONS , keep them away from main class def.
#ifdef __BOOST__
	public:

		void save_checkpoint(string filename) const{
			std::ofstream fout(filename.c_str());
			boost::archive::text_oarchive oa(fout);
			oa << (*this);
		}

		void load_checkpoint(string filename){
			std::ifstream fin(filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive ia(fin);
			ia >> (*this);
		}



	private:
		friend class boost::serialization::access;		


		template<class Archive> void serialize(Archive & ar, const unsigned int version){

			using boost::serialization::make_nvp;
			ar & make_nvp("name",name);
			ar & make_nvp("residueName",residueName);
			ar & make_nvp("residueNumber",residueNumber);
			ar & make_nvp("residueIcode",residueIcode);
			ar & make_nvp("chainId",chainId);
			ar & make_nvp("element",element);
			ar & make_nvp("segId",segId);
			ar & make_nvp("pCoorVec",pCoorVec);

			/*
			ar & name;
			ar & residueName;
			ar & residueNumber;
			ar & residueIcode;
			ar & chainId;
			ar & element;
			ar & segId;

			ar & nameSpace;
			
			ar & hasCoordinates;


			
 			ar & pCoorVec; 
			ar & bonds;
			*/
 			//ar & currentCoorIterator;
			//ar & savedCoor; 

 			//ar & icEntries; 
			//ar & pParentGroup;
			//ar & pParentContainer;
		}
#endif
		
};

// INLINED FUNCTIONS
inline void Atom::setName(string _name) {name = _name; updateResidueMap(); updateContainerMap();};
inline string Atom::getName() const {return name;};
inline void Atom::setElement(string _element) {element = _element;};
inline string Atom::getElement() const {return element;};
inline void Atom::setType(string _type) {type = _type;};
inline string & Atom::getType() {return type;};
inline void Atom::setCharge(double _charge) {charge = _charge;};
inline double & Atom::getCharge() {return charge;};
inline void Atom::setGroupNumber(unsigned int _groupNumber) {groupNumber = _groupNumber;};
inline void Atom::setTempFactor(double _bfactor) {tempFactor = _bfactor;};
inline double Atom::getTempFactor() const {return tempFactor;};
inline void Atom::setSegID(string _segid) { segId = _segid; }
inline string Atom::getSegID() const { return segId; }
inline void Atom::setParentGroup(AtomGroup * _group) {pParentGroup = _group; pParentContainer = NULL;};
inline AtomGroup * Atom::getParentGroup() const {return pParentGroup;};
inline void Atom::setParentContainer(AtomContainer * _container) {pParentContainer = _container; pParentGroup = NULL;};
inline AtomContainer * Atom::getParentContainer() const {return pParentContainer;};
inline void Atom::setCoor(CartesianPoint _p) {(*currentCoorIterator)->setCoor(_p); hasCoordinates = true;};
inline void Atom::setCoor(Real _x, Real _y, Real _z) {(*currentCoorIterator)->setCoor(_x, _y, _z); hasCoordinates = true;};
inline CartesianPoint & Atom::getCoor() { return *(*currentCoorIterator); };
inline vector<CartesianPoint *> & Atom::getAllCoor() { return pCoorVec; };
inline Real Atom::getX() const { return (*currentCoorIterator)->getX(); };
inline Real Atom::getY() const { return (*currentCoorIterator)->getY(); };
inline Real Atom::getZ() const { return (*currentCoorIterator)->getZ(); };
inline Real Atom::operator[](size_t _n) { return (*(*currentCoorIterator))[_n]; }; // return X Y Z as atom[0], [1], [2] operators
inline string Atom::toString() const { string qm = " "; if (!hasCoordinates) {qm = "?";}; string act = "+"; if (!getActive()) {act="-";} char c [100]; sprintf(c, "%-4s %-3s %4u%1s %1s [%10.3f %10.3f %10.3f]%1s(conf %3u/%3u) %1s", name.c_str(), getResidueName().c_str(), getResidueNumber(), getResidueIcode().c_str(), getChainId().c_str(), (*currentCoorIterator)->getX(), (*currentCoorIterator)->getY(), (*currentCoorIterator)->getZ(), qm.c_str(), getActiveConformation()+1, getNumberOfAltConformations(), act.c_str()); return (string)c; };
inline double Atom::distance(const Atom & _atom) const {return CartesianGeometry::instance()->distance(*(*currentCoorIterator), *(*(_atom.currentCoorIterator)));};
inline double Atom::distance2(const Atom & _atom) const {return CartesianGeometry::instance()->distance2(*(*currentCoorIterator), *(*(_atom.currentCoorIterator)));};
inline double Atom::angle(const Atom & _atom) const {return CartesianGeometry::instance()->angle(*(*currentCoorIterator), *(*(_atom.currentCoorIterator)));};
inline double Atom::angle(const Atom & _center, const Atom & _third) const {return CartesianGeometry::instance()->angle(*(*currentCoorIterator), *(*(_center.currentCoorIterator)), *(*(_third.currentCoorIterator)));};
inline double Atom::angleRadians(const Atom & _atom) const {return CartesianGeometry::instance()->angleRadians(*(*currentCoorIterator), *(*(_atom.currentCoorIterator)));};
inline double Atom::angleRadians(const Atom & _center, const Atom & _third) const {return CartesianGeometry::instance()->angleRadians(*(*currentCoorIterator), *(*(_center.currentCoorIterator)), *(*(_third.currentCoorIterator)));};
inline double Atom::dihedral(const Atom & _second, const Atom & _third, const Atom & _fourth) const {return CartesianGeometry::instance()->dihedral(*(*currentCoorIterator), *(*(_second.currentCoorIterator)), *(*(_third.currentCoorIterator)), *(*(_fourth.currentCoorIterator)));};
inline double Atom::dihedralRadians(const Atom & _second, const Atom & _third, const Atom & _fourth) const {return CartesianGeometry::instance()->dihedralRadians(*(*currentCoorIterator), *(*(_second.currentCoorIterator)), *(*(_third.currentCoorIterator)), *(*(_fourth.currentCoorIterator)));};
inline bool Atom::hasCoor() const {return hasCoordinates;};
inline void Atom::wipeCoordinates() {(*currentCoorIterator)->setCoor(0.0, 0.0, 0.0); hasCoordinates = false;};
inline void Atom::setHasCoordinates() {hasCoordinates = true;};
inline vector<IcEntry*> & Atom::getIcEntries() {return icEntries;}
inline void Atom::setActiveConformation(size_t _i) {currentCoorIterator = pCoorVec.begin() + _i;};
inline unsigned int Atom::getActiveConformation() const {return currentCoorIterator - pCoorVec.begin();};
inline unsigned int Atom::getNumberOfAltConformations() const {return pCoorVec.size();};
inline void Atom::addAltConformation() {addAltConformation(getCoor());}; //default, same as current conformation
inline void Atom::addAltConformation(const CartesianPoint & _point) {unsigned int curr = currentCoorIterator - pCoorVec.begin(); pCoorVec.push_back(new CartesianPoint(_point)); currentCoorIterator = pCoorVec.begin() + curr;};  // it is important to regenerate the iterator to the current coor in case the vector is resized
inline void Atom::removeAllAltConformations() {for (vector<CartesianPoint*>::iterator k=pCoorVec.begin()+1; k!=pCoorVec.end(); k++) {delete *k; pCoorVec.erase(k);}; currentCoorIterator = pCoorVec.begin();}; // remove all alternate conformations, keeping the first conformation
inline void Atom::saveCoor(string _coordName) { savedCoor[_coordName] = new CartesianPoint(**currentCoorIterator); }
inline bool Atom::applySavedCoor(string _coordName) { map<string, CartesianPoint*>::iterator found = savedCoor.find(_coordName); if (found != savedCoor.end()) { (*currentCoorIterator)->setCoor(*(found->second)); return true; } return false; }
//inline void Atom::setBondedTo(Atom * _pAtom, bool _bound) {if (_bound) {bonds[_pAtom] = true;} else {map<Atom*, bool>::iterator found=bonds.find(_pAtom); if (found!=bonds.end()) {bonds.erase(found);}}}
//inline vector<Atom*> Atom::getBoundAtoms() const {vector<Atom*> bonded; for (map<Atom*, bool>::const_iterator k=bonds.begin(); k!=bonds.end(); k++) {bonded.push_back(k->first);} return bonded;}
//inline vector<Atom*> Atom::getBoundAtoms() const {vector<Atom*> bonded; for (map<Atom*, map<Atom*, map<Atom*, bool> > >::const_iterator k=boundAtoms.begin(); k!=boundAtoms.end(); k++) {bonded.push_back(k->first);} return bonded;}
inline vector<vector<Atom*> > Atom::getBoundAtoms() const {
	vector<vector<Atom*> > bonded;
	for (map<Atom*, map<Atom*, map<Atom*, bool> > >::const_iterator k=boundAtoms.begin(); k!=boundAtoms.end(); k++) {
		bonded.push_back(vector<Atom*>(1, k->first));
		for (map<Atom*, map<Atom*, bool> >::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {
			bonded.back().push_back(l->first);
			for (map<Atom*, bool>::const_iterator m=l->second.begin(); m!=l->second.end(); m++) {
				bonded.back().push_back(m->first);
			}
		}
	}
	return bonded;
}
inline vector<Atom*> Atom::getOneThreeMiddleAtoms(Atom *_pAtom) const {
	vector<Atom*> middle;
	map<Atom*, map<Atom*, bool> >::const_iterator found=oneThreeAtoms.find(_pAtom);
	if (found != oneThreeAtoms.end()) {
		for (map<Atom*, bool>::const_iterator k=found->second.begin(); k!=found->second.end(); k++) {
			middle.push_back(k->first);
		}
	}
	return middle;
}
inline vector<vector<Atom*> > Atom::getOneFourMiddleAtoms(Atom *_pAtom) const {
	vector<vector<Atom*> > middle;
	map<Atom*, map<Atom*, map<Atom*, bool> > >::const_iterator found=oneFourAtoms.find(_pAtom);
	if (found != oneFourAtoms.end()) {
		for (map<Atom*, map<Atom*, bool> >::const_iterator k=found->second.begin(); k!=found->second.end(); k++) {
			middle.push_back(vector<Atom*>());
			middle.back().push_back(k->first); // second atom of the 1-4 relationship
			for (map<Atom*, bool>::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {
				middle.back().push_back(l->first); // 3rd atom
			}
		}
	}
	return middle;
}

//inline bool Atom::isBoundTo(Atom * _pAtom) const {map<Atom*, bool>::const_iterator found=bonds.find(_pAtom); return found!=bonds.end();}
//inline bool Atom::isOneThree(Atom * _pAtom, const Atom * _14caller) const {if (_pAtom == this) {return false;} for (map<Atom*, bool>::const_iterator k=bonds.begin(); k!=bonds.end(); k++) {if (k->first == _pAtom) {return false;} else if (k->first->isBoundTo(_pAtom) && k->first != _14caller) {return true;}} return false;}
//inline bool Atom::isOneFour(Atom * _pAtom) const {if (_pAtom == this) {return false;} for (map<Atom*, bool>::const_iterator k=bonds.begin(); k!=bonds.end(); k++) {if (k->first == _pAtom) {return false;} else if (k->first->isOneThree(_pAtom, this)) {return true;}} return false;}

inline bool Atom::isBoundTo(Atom * _pAtom) const { return boundAtoms.find(_pAtom) != boundAtoms.end(); }
inline bool Atom::isOneThree(Atom * _pAtom) const { return oneThreeAtoms.find(_pAtom) != oneThreeAtoms.end(); }
inline bool Atom::isOneFour(Atom * _pAtom) const { return oneFourAtoms.find(_pAtom) != oneFourAtoms.end(); }
inline bool Atom::isInAlternativeIdentity(Atom * _pAtom) const {return getParentPosition() == _pAtom->getParentPosition() && getParentResidue() != _pAtom->getParentResidue();}



#endif
