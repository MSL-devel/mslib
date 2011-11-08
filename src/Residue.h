/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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


namespace MSL { 
class Position;
class Chain;
class System;

/****************************************************************
  TO DO
 	// INDEX GETTERS
	unsigned int getAtomIndex(string _id); // "CA"
	unsigned int getAtomIndex(const Atom * _pAtom) const;	
 ****************************************************************/


class Residue : public Selectable<Residue> {
	public:
		Residue();
		Residue(std::string _resName, int _resNum, std::string _icode="");
		Residue(const AtomPointerVector & _atoms, std::string _resName, int _resNum, std::string _icode="");
		Residue(const Residue & _residue);
		~Residue();

		void operator=(const Residue & _residue); // assignment
	
		std::string getIdentityId(unsigned int _skip=0) const;
		std::string getPositionId(unsigned int _skip=0) const;

		void setResidueName(std::string _resname);
		std::string getResidueName() const;

		void setResidueNumber(int _resnum);
		int getResidueNumber() const;

		void setResidueIcode(std::string _icode);
		std::string getResidueIcode() const;

		void setChainId(std::string _chainId);
		std::string getChainId() const;

		unsigned int getIdentityIndex(); // return the index of this identity in the position

		void setNameSpace(std::string _nameSpace);
		std::string getNameSpace() const;

		void setParentPosition(Position * _position);
		Position * getParentPosition() const;
		Chain * getParentChain() const;
		System * getParentSystem() const;

		double distance(Residue &_residue, std::string _atomOfInterest="CENTROID");

		unsigned int getGroupNumber(const AtomGroup * _pGroup) const;

		/* ADD ATOMS */
		void addAtom(const Atom & _atom);
		void addAtom(std::string _atomId, const CartesianPoint & _coor=CartesianPoint(0.0, 0.0, 0.0), unsigned int _group=0, std::string _element="");
		void addAtoms(const AtomPointerVector & _atoms);
		void addAltConformationToAtom(std::string _atomId, const CartesianPoint & _coor=CartesianPoint(0.0, 0.0, 0.0));
		/* REMOVE ATOMS */
		bool removeAtom(std::string _name);
		void removeAllAtoms();
		/* UPDATES REQUESTED BY ATOMS */
		void updateAtomMap(Atom * _atom);

		/* GET ATOMS */
		unsigned int size() const;
		unsigned int atomSize() const;
		Atom & operator[](unsigned int _index);
		Atom & operator[](std::string _atomId);
		Atom & operator()(unsigned int _index); // redundant to [] operator
		Atom & operator()(std::string _atomId); // redundant to [] operator
		Atom & getAtom(unsigned int _index);
		Atom & getAtom(std::string _atomId);
		AtomPointerVector & getAtomPointers();
		std::map<std::string, Atom*> & getAtomMap();

		CartesianPoint getCentroid();

		bool getActive() const;  // is this the active identity of the position?
		// check the existance of atoms
		bool atomExists(std::string _atomId);
		//DEPRECATED exists function
	//	bool exists(std::string _name);
		Atom & getLastFoundAtom();

		/*
		  A method for quick finding all residue neighbors (within parent System object) within some geometric criteria...
		  right now it is just a spherical cutoff, but can be smarter in future... dot prodocts...ask bretth.
		  
		  will return std::vector of neighboring residue indices into parent System object

		*** ALESSANDRO 02/28/2010 These function should belong to the system

		 */

		// By default this uses sidechain centroids to find distances
		std::vector<int> findNeighbors(double _distance);

		// By default we look for ANY atom in other residue ( meaning _atomInOtherResidue is equal to "")
		std::vector<int> findNeighbors(double _distance, std::string _atomInThisResidue, std::string _atomInOtherResidue="");

		void findNeighborsAllConformations(double _distance,std::string _atomInThisResidue, std::string _atomInOtherResidue, std::vector<int> & _resnums, std::vector<int> & _altConformations);


		/***************************************************
		 *  Manage the alternative conformations of the atoms
		 ***************************************************/
		void setActiveConformation(unsigned int _i);
		//int getActiveConformation() const;
		unsigned int getNumberOfAltConformations() const;
		void addAltConformation();
		void addAltConformation(const std::vector<CartesianPoint> & _points);
		void removeAltConformation(unsigned int _i);
		void removeAllAltConformations();

		void setMaxNumberOfRotamers(unsigned int _n); // limit the number of rots
		void resetMaxNumberOfRotamers(); // remove all limits
		unsigned int getNumberOfRotamers() const;
		//  define and set rotamer sampling levels (number of rotamers given a label)
		void defineRotamerSamplingLevel(std::string _label, unsigned int _n);
		void setRotamerSamplingLevel(std::string _label);

		void wipeAllCoordinates(); // flag all active and inactive atoms as not having cartesian coordinates

		/***************************************************
		 *  Saving coordinates to buffers:
		 *
		 *  coordinates can be saved to named buffers (string _coordName),
		 *  and copied back from them
		 *
		 *  The difference between save coordinates to a buffer, and 
		 *  having multiple alternate coor is that the saved coord 
		 *  are simply a buffer that can be restored
		 *
		 *  Coor can be saved to buffer with two different commands:
		 *    saveCoor:
		 *      - saveCoor saves ONLY the current coor
		 *      - when restored with applySavedCoor, a buffer created with
		 *        saveCoor will replace the CURRENT coorinate only
		 *    saveAltCoor:
		 *      - saveAltCoor saves ALL alternative coordinates and
		 *        also remembers what was the current coordinate
		 *      - when restored with the same applySavedCoor, a buffer
		 *        created with saveAltCoor will wipe off all alternative
		 *        cordinates and recreate the situation that was present
		 *        when the buffer was saved
		 *
		 *  More details in Atom.h
		 ***************************************************/
		void saveCoor(std::string _coordName);
		void saveAltCoor(std::string _coordName);
		bool applySavedCoor(std::string _coordName);
		void clearSavedCoor(std::string _coordName="");		

		friend std::ostream & operator<<(std::ostream &_os, const Residue & _res)  {_os << _res.toString(); return _os;};
		std::string toString() const;

		double getSasa() const;

		// virtual functions from Selectable class, allows selection of residues
		virtual void addSelectableFunctions(); 
	private:

		void deletePointers();
		void setup(std::string _resName, int _resNum, std::string _insertionCode, std::string _chainId);
		void copy(const Residue & _residue);
		void updatePositionMap();

		Position * pParentPosition;
		
		std::string residueName;
		int residueNumber;
		std::string residueIcode;
		std::string chainId;
		std::string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd

		AtomPointerVector atoms;
		std::vector<AtomGroup*> electrostaticGroups;
		std::map<std::string, Atom*> atomMap;

		std::map<std::string, Atom*>::iterator foundAtom;

		bool limitRotamers;
		unsigned int maxNumOfRotamers;
		std::map<std::string, unsigned int> rotamerSamplingLevels;
		
};

// INLINE FUNCTIONS
inline void Residue::setResidueName(std::string _resname) {residueName = _resname; updatePositionMap();}
inline std::string Residue::getResidueName() const {return residueName;}
inline void Residue::setParentPosition(Position * _position) {pParentPosition = _position;}
inline Position * Residue::getParentPosition() const {return pParentPosition;}
inline unsigned int Residue::size() const {return atoms.size();}
inline unsigned int Residue::atomSize() const {return atoms.size();}
inline Atom & Residue::operator[](unsigned int _index) {return *atoms[_index];}
inline Atom & Residue::operator[](std::string _atomId) {return getAtom(_atomId);}
inline Atom & Residue::operator()(unsigned int _index) {return *atoms[_index];}
inline Atom & Residue::operator()(std::string _atomId) {return getAtom(_atomId);}
inline Atom & Residue::getAtom(unsigned int _index) {return *atoms[_index];}
inline Atom & Residue::getAtom(std::string _atomId) {
	if (atomExists(_atomId)) {
		return *(foundAtom->second);
	} else {
		// we should add try... catch support here
		std::cerr << "ERROR 3812: atom " << _atomId << " does not exist in residue at inline Atom & Residue::getAtom(string _atomId)" << std::endl;
		exit(3812);
	}
}
inline AtomPointerVector & Residue::getAtomPointers() {return atoms;}
//inline bool Residue::atomExists(std::string _name) {foundAtom = atomMap.find(_name); return foundAtom != atomMap.end();}
inline bool Residue::atomExists(std::string _atomId) {
	foundAtom = atomMap.find(_atomId);
	if (foundAtom != atomMap.end()) {
		return true;
	} else {
		// try to parse an atomId
		std::string chainid;
		int resnum;
		std::string icode;
		std::string identity;
		std::string atomName;
		bool OK = MslTools::parseAtomOfIdentityId(_atomId, chainid, resnum, icode, identity, atomName, 3);
		// if "CA" was given, identity will be = "" and the next function will return the atom for the current
		// identity
		if(OK) {
			foundAtom = atomMap.find(atomName);
			if (foundAtom != atomMap.end()) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}

}
//inline bool Residue::exists(std::string _name) {std::cerr << "DEPRECATED: Residue::exists(string), use Residue::atomExist(string)" << std::endl; return atomExists(_name);}

/**
 * This method will calculate the centorid of this residue
 * by calculating the average of the x, y, and z coordinates
 * for each atom in the residue.
 *
 * @return Returns a CartesianPoint representing the centroid of this residue.
 */
inline CartesianPoint Residue::getCentroid() {
    CartesianPoint centroid(0.0, 0.0, 0.0);
    
    for(AtomPointerVector::iterator it = atoms.begin(); it != atoms.end(); ++it) {
        centroid += (*it)->getCoor();
    }

    if (size() != 0.0){
	    centroid /= (double)size();
    }else {
	    std::cerr << "Residue::getCentroid() -> Residue has no atoms: "<<toString()<<std::endl;
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
 * @param _atomOfInterest  An optional std::string indicating which atom in the
 *                         residue should be used to calculate the distance.
 *                         By default this method will use the CENTROID, not a 
 *                         specific atom.
 *
 * @return A double giving the distance between the two residues.
 */
inline double Residue::distance (Residue &_residue, std::string _atomOfInterest) {
    double dist = 0.0f;
    if(_atomOfInterest == "CENTROID") {
        dist = getCentroid().distance( _residue.getCentroid() );
    }
    else {
        if( atomExists(_atomOfInterest) && _residue.atomExists(_atomOfInterest) )
            dist = getAtom(_atomOfInterest).distance( _residue(_atomOfInterest) );
        else {
            std::cout << "ERROR 54222: Did not find atom " << _atomOfInterest << " in residue.  Can't calculate distance in inline double Residue::distance (Residue &_residue, std::string _atomOfInterest)" << std::endl;
            exit(1);
        }
            
    }

    return dist;
}

inline Atom & Residue::getLastFoundAtom() {return *(foundAtom->second);}
inline void Residue::wipeAllCoordinates() {for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {(*k)->wipeCoordinates();}}
inline unsigned int Residue::getGroupNumber(const AtomGroup * _pGroup) const {
	for (std::vector<AtomGroup*>::const_iterator k=electrostaticGroups.begin(); k!=electrostaticGroups.end(); k++) {
		if (*k == _pGroup) {
			return k-electrostaticGroups.begin();
		}
	}
	std::cerr << "ERROR 19210: AtomGroup address " << _pGroup << " not found in the electrostatic groups in unsigned int Residue::getGroupNumber(AtomGroup * _pGroup)" << std::endl;
	exit(19210);
	return 0;
}
inline std::string Residue::toString() const {
	/*
	std::stringstream ss;
	ss << getChainId() << " " << getResidueNumber() << getResidueIcode() << " " << residueName;
	return ss.str();
	*/
	/*
	char tmp[100];
	//sprintf(tmp," [ %1s %5d %3s %1s ] " , getChainId().c_str(),getResidueNumber(),getResidueName().c_str(),getResidueIcode().c_str());
	sprintf(tmp,"[%1s %5d %1s %3s] " , getChainId().c_str(), getResidueNumber(), getResidueIcode().c_str(), getResidueName().c_str());

	return (std::string)tmp;
	*/
	return getIdentityId();
}

inline double Residue::getSasa() const {
	double sasa = 0.0;
	for (AtomPointerVector::const_iterator k=atoms.begin(); k!=atoms.end(); k++) {
		sasa += (*k)->getSasa();
	}
	return sasa;
}

inline std::string Residue::getIdentityId(unsigned int _skip) const {
	return MslTools::getIdentityId(getChainId(), getResidueNumber(), getResidueIcode(), getResidueName(), _skip);
}
inline std::string Residue::getPositionId(unsigned int _skip) const {
	return MslTools::getPositionId(getChainId(), getResidueNumber(), getResidueIcode(), _skip);
}
inline void Residue::saveCoor(std::string _coordName) {atoms.saveCoor(_coordName);}
inline void Residue::saveAltCoor(std::string _coordName) {atoms.saveAltCoor(_coordName);}
inline bool Residue::applySavedCoor(std::string _coordName) {return atoms.applySavedCoor(_coordName);}
inline void Residue::clearSavedCoor(std::string _coordName) {atoms.clearSavedCoor(_coordName);}
/***************************************************
 *  Comment...
 ***************************************************/
inline unsigned int Residue::getNumberOfRotamers() const {
	unsigned int rots = getNumberOfAltConformations();
	if (limitRotamers && maxNumOfRotamers < rots) {
		return maxNumOfRotamers;
	} else {
		return rots;
	}
}
inline void Residue::setMaxNumberOfRotamers(unsigned int _n) {
	// limit the number of rots
	limitRotamers = true;
	maxNumOfRotamers = _n;
}

inline void Residue::resetMaxNumberOfRotamers() { // remove all limits
	// limit the number of rots
	limitRotamers = false;
	maxNumOfRotamers = 0;
}
inline void Residue::defineRotamerSamplingLevel(std::string _label, unsigned int _n) {
	rotamerSamplingLevels[_label] = _n;
}
inline void Residue::setRotamerSamplingLevel(std::string _label) {
	if (rotamerSamplingLevels.find(_label) != rotamerSamplingLevels.end()) {
		limitRotamers = true;
		maxNumOfRotamers = rotamerSamplingLevels[_label];
	} else {
		std::cerr << "Error 12834: sampling level " << _label << " not found in inline void Residue::setRotamerSamplingLevel(string _label" << std::endl;
		exit(12834);
	}

}
}
#endif
