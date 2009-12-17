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

#ifndef RESIDUE_H
#define RESIDUE_H

#include <string>
#include <map>

#include "AtomGroup.h"

using namespace std;

class Position;
class Chain;
class System;


class Residue : public Selectable<Residue> {
	public:
		Residue();
		Residue(string _resName, int _resNum, string _icode="");
		Residue(const AtomVector & _atoms, string _resName, int _resNum, string _icode="");
		Residue(const Residue & _residue);
		~Residue();

		void operator=(const Residue & _residue); // assignment
	
		void setResidueName(string _resname);
		string getResidueName() const;

		void setResidueNumber(int _resnum);
		int getResidueNumber() const;

		void setResidueIcode(string _icode);
		string getResidueIcode() const;

		void setChainId(string _chainId);
		string getChainId() const;

		unsigned int getIdentityIndex(); // return the index of this identity in the position

		void setNameSpace(string _nameSpace);
		string getNameSpace() const;

		void setParentPosition(Position * _position);
		Position * getParentPosition() const;
		Chain * getParentChain() const;
		System * getParentSystem() const;

		double distance(Residue &_residue, string _atomOfInterest="CENTROID");

		unsigned int getGroupNumber(const AtomGroup * _pGroup) const;

		/* ADD ATOMS */
		void addAtom(const Atom & _atom);
		void addAtom(string _name, const CartesianPoint & _coor=CartesianPoint(0.0, 0.0, 0.0), size_t _group=0);
		void addAtoms(const AtomVector & _atoms);
		void addAltConformationToAtom(string _name, const CartesianPoint & _coor=CartesianPoint(0.0, 0.0, 0.0));
		/* REMOVE ATOMS */
		bool removeAtom(string _name);
		void removeAllAtoms();
		/* UPDATES REQUESTED BY ATOMS */
		void updateAtomMap(Atom * _atom);

		/* GET ATOMS */
		unsigned int size() const;
		Atom & operator[](size_t _n);
		Atom & operator()(string _name);
		Atom & getAtom(size_t _n);
		Atom & getAtom(string _name);
		AtomVector & getAtoms();
		map<string, Atom*> & getAtomMap();

		CartesianPoint getCentroid();

		bool getActive() const;  // is this the active identity of the position?
		// check the existance of atoms
		bool exists(string _name);
		Atom & getLastFoundAtom();


		/*
		  A method for quick finding all residue neighbors (within parent System object) within some geometric criteria...
		  right now it is just a spherical cutoff, but can be smarter in future... dot prodocts...ask bretth.
		  
		  will return vector of neighboring residue indices into parent System object
		 */
		vector<int> findNeighbors(double _distance);
		vector<int> findNeighbors(double _distance, string _atomInThisResidue, string _atomInOtherResidue);


		/***************************************************
		 *  Manage the alternative conformations of the atoms
		 ***************************************************/
		void setActiveConformation(size_t _i);
		//int getActiveConformation() const;
		unsigned int getNumberOfAltConformations() const;
		void addAltConformation();
		void addAltConformation(const vector<CartesianPoint> & _points);
		void removeAltConformation(size_t _i);
		void removeAllAltConformations();

		void wipeAllCoordinates(); // flag all active and inactive atoms as not having cartesian coordinates

		friend ostream & operator<<(ostream &_os, const Residue & _res)  {_os << _res.toString(); return _os;};
		string toString() const;

		double getSasa() const;

		// virtual functions from Selectable class, allows selection of residues
		virtual void addSelectableFunctions(); 
	private:

		void deletePointers();
		void setup(string _resName, int _resNum, string _insertionCode, string _chainId);
		void copy(const Residue & _residue);
		void updatePositionMap();

		Position * pParentPosition;
		
		string residueName;
		int residueNumber;
		string residueIcode;
		string chainId;
		string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd

		AtomVector atoms;
		vector<AtomGroup*> electrostaticGroups;
		map<string, Atom*> atomMap;

		map<string, Atom*>::iterator foundAtom;
		
};

// INLINE FUNCTIONS
inline void Residue::setResidueName(string _resname) {residueName = _resname; updatePositionMap();}
inline string Residue::getResidueName() const {return residueName;}
inline void Residue::setParentPosition(Position * _position) {pParentPosition = _position;}
inline Position * Residue::getParentPosition() const {return pParentPosition;}
inline unsigned int Residue::size() const {return atoms.size();}
inline Atom & Residue::operator[](size_t _n) {return *atoms[_n];}
inline Atom & Residue::operator()(string _name) {return *atomMap[_name];}
inline Atom & Residue::getAtom(size_t _n) {return *atoms[_n];}
inline Atom & Residue::getAtom(string _name) {
	if (exists(_name)) {
		return *(foundAtom->second);
	} else {
		cerr << "ERROR 3812: atom " << _name << " does not exist in residue at inline Atom & Residue::getAtom(string _name)" << endl;
		exit(3812);
	}
}
inline AtomVector & Residue::getAtoms() {return atoms;}
inline bool Residue::exists(string _name) {foundAtom = atomMap.find(_name); return foundAtom != atomMap.end();}

/**
 * This method will calculate the centorid of this residue
 * by calculating the average of the x, y, and z coordinates
 * for each atom in the residue.
 *
 * @return Returns a CartesianPoint representing the centroid of this residue.
 */
inline CartesianPoint Residue::getCentroid() {
    CartesianPoint centroid(0.0, 0.0, 0.0);
    
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); ++it) {
        centroid += (*it)->getCoor();
    }

    if (size() != 0.0){
	    centroid /= (double)size();
    }else {
	    cerr << "Residue::getCentroid() -> Residue has no atoms: "<<toString()<<endl;
    }

    return centroid;
}

/**
 * This method calculates the distance between this residue and another given
 * residue.  By default, it will calculate the distance between the two
 * centroids of the residues.  If a particular atom of interest is specified,
 * it calculate the distance between these two atoms.  Note that if either
 * residue lacks an atom with the given name, the method will fail and
 * will call exit(1).
 *
 * @param _residue  The other residue to use in our distance calculation.
 * @param _atomOfInterest  An optional string indicating which atom in the
 *                         residue should be used to calculate the distance.
 *                         By default this method will use the CENTROID, not a 
 *                         specific atom.
 *
 * @return A double giving the distance between the two residues.
 */
inline double Residue::distance (Residue &_residue, string _atomOfInterest) {
    double dist = 0.0f;
    if(_atomOfInterest == "CENTROID") {
        dist = getCentroid().distance( _residue.getCentroid() );
    }
    else {
        if( exists(_atomOfInterest) && _residue.exists(_atomOfInterest) )
            dist = getAtom(_atomOfInterest).distance( _residue(_atomOfInterest) );
        else {
            cout << "ERROR 54222: Did not find atom " << _atomOfInterest << " in residue.  Can't calculate distance in inline double Residue::distance (Residue &_residue, string _atomOfInterest)" << endl;
            exit(1);
        }
            
    }

    return dist;
}

inline Atom & Residue::getLastFoundAtom() {return *(foundAtom->second);}
inline void Residue::wipeAllCoordinates() {for (AtomVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {(*k)->wipeCoordinates();}}
inline unsigned int Residue::getGroupNumber(const AtomGroup * _pGroup) const {
	for (vector<AtomGroup*>::const_iterator k=electrostaticGroups.begin(); k!=electrostaticGroups.end(); k++) {
		if (*k == _pGroup) {
			return k-electrostaticGroups.begin();
		}
	}
	cerr << "ERROR 19210: AtomGroup address " << _pGroup << " not found in the electrostatic groups in unsigned int Residue::getGroupNumber(AtomGroup * _pGroup)" << endl;
	exit(19210);
	return 0;
}
inline string Residue::toString() const {
	/*
	stringstream ss;
	ss << getChainId() << " " << getResidueNumber() << getResidueIcode() << " " << residueName;
	return ss.str();
	*/
	char tmp[100];
	//sprintf(tmp," [ %1s %5d %3s %1s ] " , getChainId().c_str(),getResidueNumber(),getResidueName().c_str(),getResidueIcode().c_str());
	sprintf(tmp,"[%1s %5d %1s %3s] " , getChainId().c_str(), getResidueNumber(), getResidueIcode().c_str(), getResidueName().c_str());

	return (string)tmp;
}

inline double Residue::getSasa() const {
	double sasa = 0.0;
	for (AtomVector::const_iterator k=atoms.begin(); k!=atoms.end(); k++) {
		sasa += (*k)->getSasa();
	}
	return sasa;
}
#endif
