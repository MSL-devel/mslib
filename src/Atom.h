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


#ifndef ATOM_H
#define ATOM_H


// STL Includes
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>


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
namespace MSL { 
class IcEntry;
class IcTable;
class AtomGroup;
class AtomContainer;
class Residue;
class Position;
class Chain;
class System;

class Atom : public Selectable<Atom> {
	public:
		Atom();
		Atom(const std::string _atomId, std::string _element=""); // atomId "A,37,ILE,CA"
		Atom(std::string _atomId, Real _x, Real _y, Real _z, std::string _element=""); 
		Atom(std::string _atomId, const CartesianPoint & _p, std::string _element="");
		Atom(const Atom & _atom);
		virtual ~Atom();

		void operator=(const Atom & _atom); // assignment
	
		std::string getAtomId(unsigned int _skip=0) const;
		std::string getAtomOfIdentityId(unsigned int _skip=0) const;
		std::string getIdentityId(unsigned int _skip=0) const;
		std::string getPositionId(unsigned int _skip=0) const;

		void setName(std::string _name);
		std::string getName() const;
		void setResidueName(std::string _resname);
		std::string getResidueName() const;
		void setResidueNumber(int _resnum);
		int getResidueNumber() const;
		void setResidueIcode(std::string _icode);
		std::string getResidueIcode() const;
		void setChainId(std::string _chainId);
		std::string getChainId() const;
		void setElement(std::string _element);
		std::string getElement() const;
		void setType(std::string _type);
		std::string & getType();
		void setRadius(double _radius);
		double getRadius() const;
		void setCharge(double _charge);
		double getCharge() const;
		void setGroupNumber(unsigned int _groupNumber);
		unsigned int getGroupNumber() const;
		CartesianPoint& getGroupGeometricCenter(unsigned int _stamp=0);
		unsigned int getIdentityIndex(); // return the index of its parent identity in the position
		bool isInAlternativeIdentity(Atom * _pAtom) const; // checks if the two atoms happen to be in the same position but different residue types (cannot coexist)

		void setTempFactor(double _bfactor);
		double getTempFactor() const;

		void setSasa(double _sasa);
		double getSasa() const;

		void setSegID(std::string _segid);
		std::string getSegID() const;

		void setNameSpace(std::string _nameSpace);
		std::string getNameSpace() const;

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
		void copyAllCoor(const Atom _a); // copy all coordinates from another atoms (including alt coors)
		CartesianPoint & getCoor();
		std::vector<CartesianPoint *> & getAllCoor();
		std::vector<CartesianPoint *> & getHiddenCoor();
		Real getX() const;
		Real getY() const;
		Real getZ() const;
		Real operator[](unsigned int _n); // return X Y Z as atom[0], [1], [2] operators

		/*************************************************************
		 * print atom information
		 * Formats:
		 *   0   default. Current coor only                         coordinates wiped (? if yes, otherwise blank)
		 *       name        atom number           x,y,z            |   conf of N   active atom (+ or -)
		 *       |           |     ------------------|------------- | -----|------  |
		 *       CA   LEU   37    [     0.000      0.000      0.000] (conf   1/  3) +
		 *
		 *   1   all alt confs, not including hidden
		 *       CA   LEU   37    [     0.000      0.000      0.000] (conf   1/  3) +
		 *       CA   LEU   37   *[     3.000      0.000      0.000] (conf   2/  3) +   active conf (*)
		 *       CA   LEU   37    [     4.000      0.000      0.000] (conf   3/  3) +
		 *
		 *   2   all alt confs, including hidden
		 *       CA   LEU   37    [     0.000      0.000      0.000] (conf   1/  3) +
		 *       CA   LEU   37    [     2.000      0.000      0.000] (conf   H/  1) +   hidden conf
		 *       CA   LEU   37   *[     3.000      0.000      0.000] (conf   2/  3) +
		 *       CA   LEU   37    [     4.000      0.000      0.000] (conf   3/  3) +
		 *************************************************************/
		std::string toString() const;
		void setToStringFormat(unsigned int _format);
		friend std::ostream & operator<<(std::ostream &_os, const Atom & _atom)  {_os << _atom.toString(); return _os;};

		// MEASURE RELATIOHSHIPS
		double distance(const Atom & _atom) const;
		double distance2(const Atom & _atom) const;
		double angle(const Atom & _atom) const;
		double angle(const Atom & _center, const Atom & _third) const;
		double angleRadians(const Atom & _atom) const;
		double angleRadians(const Atom & _center, const Atom & _third) const;
		double dihedral(const Atom & _second, const Atom & _third, const Atom & _fourth) const;
		double dihedralRadians(const Atom & _second, const Atom & _third, const Atom & _fourth) const;

		double groupDistance(Atom & _atom, unsigned int _stamp=0);

		// virtual functions from Selectable class, allows selection of atoms
		virtual void addSelectableFunctions(); 
		void clearAllFlags();

		// Minimization index is -1 for fixed atoms and >= 1 for atoms that will be minimized. It should not be 0 
		int getMinimizationIndex();
		void setMinimizationIndex(int _index);


		/***************************************************
		 *  Building from internal coordinates
		 ***************************************************/
		bool hasCoor() const;
		void wipeCoordinates();
		void setHasCoordinates(bool _flag=true);
		void addIcEntry(IcEntry * _ic); // add an ic entry
		bool buildFromIc(bool _onlyFromActive=true); // try to build from the atom's ic entries (icEntries table)
		//bool buildFromIc(const std::map<IcEntry*, bool> & _exclude); // try to build from the atom's ic entries (icEntries table)
		//bool buildFromIc(std::map<Atom*, bool> & _exclude, bool _onlyFromActive=true); // try to build from the atom's ic entries (icEntries table)
		bool buildFromIc(std::map<IcEntry*, bool> & _exclude, bool _onlyFromActive=true); // try to build from the atom's ic entries (icEntries table)
		void removeIcEntry(IcEntry * _ic); // remove the ic pointer with this address from icEntries
		std::vector<IcEntry*> & getIcEntries();
		void printIcEntries() const;

		/***************************************************
		 *  Alternate conformations
		 *
		 *  Alternative conformations are stored by adding
		 *  CartesianPoint ponters to the pCoorVec std::vector;
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
		 *  an internal resize of the std::vector.  Always repoint
		 *  the iterators when element are inserted or deleted
		 *
		 ***************************************************/
		bool setActiveConformation(unsigned int _i);
		unsigned int getActiveConformation() const;
		unsigned int getNumberOfAltConformations(bool _absolute=false) const; // if _absolute=true it counts the hidden
		void addAltConformation(); //default, same as current conformation
		void addAltConformation(const CartesianPoint & _point);  // it is important to regenerate the iterator to the current coor in case the std::vector is resized
		void addAltConformation(Real _x, Real _y, Real _z);  // it is important to regenerate the iterator to the current coor in case the std::vector is resized
		void removeAltConformation(unsigned int _i); // remove one specific alt conformation
		void removeAllAltConformations(); // remove all alternate conformations, keeping the current conformation

		/***************************************************
		 *  Alternate conformations can be temporarily hidden
		 *
		 *  If hidden the alternative coordinates will be still stored but 
		 *  it is like it is not present.  Useful for turning 
		 *  on and off rotamers, which are implemented using the
		 *  alternative coordinates
		 *
		 *  HOW TO OPERATE
		 *  For example, let's say there are 5 alt-coor loaded.
		 *  The absolute index is 0 to 4 independengly if a
		 *  coor is hidden. The relative index only considers 
		 *  unhidden conformations. When everything is active
		 *  they are the same:
		 *    Relative index: 0 1 2 3 4
		 *    Absative index: 0 1 2 3 4
		 * 
		 *  Hide alt coor 1 with absolute index (the atom will
		 *  behave like it had only 4 alt coors)
		 *    atm.hideAltCoorAbsIndex(1);
		 *    Rel: 0 1 2 3 4  >>  0 - 1 2 3  (4)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *           ^              ^
		 *
		 *    atm.hideAltCoorAbsIndex(3);
		 *    Rel: 0 - 1 2 3  >>  0 - 1 - 2  (3)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *               ^              ^
		 *
		 *  Hide by relative index
		 *    atm.hideAltCoorRelIndex(1);
		 *    Rel: 0 - 1 - 2  >>  0 - - - 1  (2)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *             ^              ^
		 *
		 *  Unhide a alt coor
		 *    atm.unhideAltCoorAbsIndex(3);
		 *    Rel: 0 - - - 1  >>  0 - - 1 2  (3)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *             ^              ^
		 *
		 *  Hide all alt coors but one (absolute undex)
		 *    atm.hideAllAltCoorsButOneAbsIndex(3)
		 *    Rel: 0 - - 1 2  >>  - - - 0 -  (1)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *               ^              ^
		 *
		 *  Unhide a alt coor
		 *    atm.unhideAltCoorAbsIndex(1);
		 *    Rel: - - - 0 -  >>  - 0 - 1 -  (2)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *           ^              ^
		 *
		 *  Hide all alt coors but one (relative undex)
		 *    atm.hideAllAltCoorsButOneRelIndex(0)
		 *    Rel: - 0 - 1 -  >>  - 0 - - -  (1)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *           ^              ^
		 *
		 *  Hide all alt coors except the first 3
		 *    atm.hideAllAltCoorsButFirstN(3);
		 *    Rel: - 0 - - -  >>  0 1 2 - -  (3)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *         ^ ^ ^          ^ ^ ^
		 *
		 *  Unhide all alt coors
		 *    atm.unhideAllAltCoors()
		 *    Rel: 0 1 2 - -  >>  0 1 2 3 4  (5)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *         ^ ^ ^ ^ ^      ^ ^ ^ ^ ^
		 *
		 ***************************************************/
		bool hideAltCoorRelIndex(unsigned int _relativeIndex); // hide a coor based on relative index
		bool hideAltCoorAbsIndex(unsigned int _absoluteIndex); // hide a coor based on absolute index
		bool hideAllAltCoorsButOneRelIndex(unsigned int _keepThisIndex); // turns all alt coor off except one, expressed as relative index
		bool hideAllAltCoorsButOneAbsIndex(unsigned int _keepThisIndex); // turns all alt coor off except one, expressed as absolute index
		bool hideAllAltCoorsButFirstN(unsigned int _numberToKeepAbsIndex); // turns all alt coor off except the first N, expressed as absolute index
		bool unhideAltCoorAbsIndex(unsigned int _absoluteIndex); // unhide a specific coor based on absolute index
		bool unhideAllAltCoors();
	
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
		 * As atoms have alternate conformations, residues can
		 * have alternate identities.
		 *
		 * Atoms can become hidden when the identity they belong 
		 * to is hidden
		 ***************************************************/


		bool getHidden() const;

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
		 *  Examples
		 *   saveCoor, simple case
		 *       Coor (acvive = *)        buffer
		 *       *[0.0  0.0  0.0]
		 *                       saveCoor
		 *       *[0.0  0.0  0.0]    ->  [0.0  0.0  0.0]
		 *                  change the position
		 *       *[3.7  4.1  2.3]    ->  [0.0  0.0  0.0]
		 *                  applySavedCoor
		 *       *[0.0  0.0  0.0]   <-   [0.0  0.0  0.0]
		 *
		 *   saveCoor, with alt coordinates
		 *       Coor (active = *)        buffer
		 *        [0.0  0.0  0.0]
		 *       *[1.0  0.0  0.0]
		 *        [2.0  0.0  0.0]
		 *                       saveCoor
		 *        [0.0  0.0  0.0]    ->  [1.0  0.0  0.0]
		 *       *[1.0  0.0  0.0]
		 *        [2.0  0.0  0.0]
		 *                  change active coordinate
		 *        [0.0  0.0  0.0]    ->  [1.0  0.0  0.0]
		 *        [1.0  0.0  0.0]
		 *       *[2.0  0.0  0.0]
		 *                  applySavedCoor
		 *        [0.0  0.0  0.0]   <-  [1.0  0.0  0.0]
		 *        [1.0  0.0  0.0]
		 *       *[1.0  0.0  0.0]
		 *
		 *   saveAltCoor
		 *       Coor (active = *)        buffer
		 *        [0.0  0.0  0.0]
		 *       *[1.0  0.0  0.0]
		 *        [2.0  0.0  0.0]
		 *                       saveAltCoor
		 *        [0.0  0.0  0.0]    ->  [0.0  0.0  0.0]
		 *       *[1.0  0.0  0.0]       *[1.0  0.0  0.0]
		 *        [2.0  0.0  0.0]        [2.0  0.0  0.0]
		 *                  remove all alt coors
		 *       *[2.0  0.0  0.0]    ->  [0.0  0.0  0.0]
		 *                              *[1.0  0.0  0.0]
		 *                               [2.0  0.0  0.0]
		 *                  applySavedCoor
		 *        [0.0  0.0  0.0]   <-   [0.0  0.0  0.0]
		 *       *[1.0  0.0  0.0]       *[1.0  0.0  0.0]
		 *        [2.0  0.0  0.0]        [2.0  0.0  0.0]
		 *
		 ***************************************************/
		void saveCoor(std::string _coordName);
		void saveAltCoor(std::string _coordName);
		bool applySavedCoor(std::string _coordName);
		void clearSavedCoor(std::string _coordName="");

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
		 *  |                      DONE                        |
		 *  - - - - - - - - - - - - - - - - - - - - - - - - - - 
		 *
		 ***************************************************/
		void setBoundTo(Atom * _pAtom);
		void setUnboundFrom(Atom * _pAtom, bool _propagate=true);
		void setUnboundFromAll(bool _propagate=true); // remove all bonds
		std::vector<Atom*> getBonds();
		std::vector<std::vector<Atom*> > getBoundAtoms() const;
		bool isBoundTo(Atom * _pAtom) const;
		bool isOneThree(Atom * _pAtom) const;
		//bool isOneThree(Atom * _pAtom, const Atom * _14caller=NULL) const;
		bool isOneFour(Atom * _pAtom) const;
		std::vector<Atom*> getOneThreeMiddleAtoms(Atom * _pAtom) const; // returns a std::vector with the middle atoms for all 1-3 relationships with a given atom _pAtom (generally there is just one unless it is a 4 members ring)
		std::vector<std::vector<Atom*> > getOneFourMiddleAtoms(Atom * _pAtom) const; // returns a std::vector with the second and third atoms for all 1-4 relationships with a given atom _pAtom (can be multiple in 6 member rings)
		std::set<Atom*> findLinkedAtoms(const std::set<Atom*> & _excluded);

		/***************************************************
		 *  Ask the chain if this position is N- or C-terminal
		 ***************************************************/
		bool isPositionNterminal() const;
		bool isPositionCterminal() const;

	private:
		void setup(CartesianPoint _point, std::string _name, std::string _element);
		void copy(const Atom & _atom);
		void reset();
		void deletePointers();
		void updateResidueMap();
		void updateContainerMap();
		void purge13(Atom * _pAtom2, Atom * _pAtom3);
		void purge14mid(Atom * _pAtom2, Atom * _pAtom3);
		void purge14end(Atom * _pAtom3, Atom * _pAtom4);
		void removeFromIc();
		//void removeBonds();

		bool hideAltCoors(unsigned int _absoluteIndex, unsigned int _relativeIndex, unsigned int _indexInHiddenn);
		void convertRelToAbs(unsigned int _relativeIndex, unsigned int & _absoluteIndex, unsigned int & _indexInHidden) const;
		void convertAbsToRel(unsigned int _absoluteIndex, unsigned int & _relativeIndex, unsigned int & _indexInHidden) const;

		std::string name;
		std::string residueName;
		int residueNumber;
		std::string residueIcode;
		std::string chainId;
		std::string element;
		double charge;
		std::string type;
		double radius;
		unsigned int groupNumber;
		double tempFactor;
		double sasa;
		std::string segId;
		// we can put the charge in the electrostatic interactions like the vdw radii etc
		
		std::string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd

		bool hasCoordinates;
		std::vector<IcEntry*> icEntries;
		unsigned int toStringFormat;

		int minIndex; //minimzation index 

		std::vector<CartesianPoint*> pCoorVec;
		std::vector<CartesianPoint*>::iterator currentCoorIterator;
		std::vector<CartesianPoint*> pHiddenCoorVec;
		std::vector<unsigned int> hiddenCoorIndeces;

		std::map<std::string, CartesianPoint*> savedCoor;
		std::map<std::string, std::vector<CartesianPoint*> > savedAltCoor;
		std::map<std::string, std::vector<CartesianPoint*> > savedHiddenCoor;
		std::map<std::string, std::vector<unsigned int> > savedHiddenCoorIndeces;
		std::map<std::string, unsigned int> savedAltCoorCurrent;

		// pointer to parent electrostatic group
		AtomGroup * pParentGroup;
		// pointer to parent atom container
		AtomContainer * pParentContainer;

		/*********************************************************
		 *  Structure that record what other atoms are bound (directly
		 *  or indirectly up to 1-4) to this atom
		 *
		 *      E   F          bound atoms of A (this)
		 *       \ /             boundAtoms[B]            A-B 1-2
		 *        B                boundAtoms[B][E]       A-B-E: 1-3
		 *        |                boundAtoms[B][F]
		 *       *A      K       boundAtoms[D]            A-D: 1-2
		 *       / \    /          boundAtoms[D][I]       A-D-I: 1-3
		 *   G--C   D--I             boundAtoms[D][I][K]  A-D-I-K: 1-4
		 *      |   |   \
		 *      H   J    L     1-3 atoms of A
		 *         /             oneThreeAtoms[E][B]      E is 1-3 through B
		 *        M              oneThreeAtoms[I][D]      D is 1-3 through D
		 *                         
		 *                     1-4 atoms of A (note the unexpected order 4->2->3)
		 *                       oneFourAtoms[K][D][I]    K is 1-4 through D-I
		 *                       oneFourAtoms[M][D][J]    M is 1-4 through D-J
		 *********************************************************/
		std::map<Atom*, std::map<Atom*, std::map<Atom*, bool> > > boundAtoms;
		std::map<Atom*, std::map<Atom*, bool> > oneThreeAtoms; // oneThreeAtoms[X][Y] corresponds to this-Y-X
		std::map<Atom*, std::map<Atom*, std::map<Atom*, bool> > > oneFourAtoms; // oneFourAtoms[X][Y][Z] corresponds to this-Y-Z-X

		// BOOST-RELATED FUNCTIONS , keep them away from main class def.
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
			// does not work ?
			//ar & make_nvp("currentCoorIterator",currentCoorIterator);

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
inline void Atom::setName(std::string _name) {name = _name; updateResidueMap(); updateContainerMap();};
inline std::string Atom::getName() const {return name;};
inline void Atom::setElement(std::string _element) {element = _element;};
inline std::string Atom::getElement() const {return element;};
inline void Atom::setType(std::string _type) {type = _type;};
inline std::string & Atom::getType() {return type;};
inline void Atom::setRadius(double _radius) {radius = _radius;};
inline double Atom::getRadius() const {return radius;};
inline void Atom::setCharge(double _charge) {charge = _charge;};
inline double Atom::getCharge () const {return charge;};
inline void Atom::setGroupNumber(unsigned int _groupNumber) {groupNumber = _groupNumber;};
inline void Atom::setTempFactor(double _bfactor) {tempFactor = _bfactor;};
inline double Atom::getTempFactor() const {return tempFactor;};
inline void Atom::setSasa(double _sasa) {sasa = _sasa;};
inline double Atom::getSasa() const {return sasa;};
inline void Atom::setSegID(std::string _segid) { segId = _segid; }
inline std::string Atom::getSegID() const { return segId; }
inline void Atom::setParentGroup(AtomGroup * _group) {pParentGroup = _group; pParentContainer = NULL;};
inline AtomGroup * Atom::getParentGroup() const {return pParentGroup;};
inline void Atom::setParentContainer(AtomContainer * _container) {pParentContainer = _container; pParentGroup = NULL;};
inline AtomContainer * Atom::getParentContainer() const {return pParentContainer;};
inline void Atom::setCoor(CartesianPoint _p) {(*currentCoorIterator)->setCoor(_p); hasCoordinates = true;};
inline void Atom::setCoor(Real _x, Real _y, Real _z) {(*currentCoorIterator)->setCoor(_x, _y, _z); hasCoordinates = true;};
inline CartesianPoint & Atom::getCoor() { return *(*currentCoorIterator); };
inline std::vector<CartesianPoint *> & Atom::getAllCoor() { return pCoorVec; };
inline std::vector<CartesianPoint *> & Atom::getHiddenCoor() { return pHiddenCoorVec; };
inline Real Atom::getX() const { return (*currentCoorIterator)->getX(); };
inline Real Atom::getY() const { return (*currentCoorIterator)->getY(); };
inline Real Atom::getZ() const { return (*currentCoorIterator)->getZ(); };
inline Real Atom::operator[](unsigned int _n) { return (*(*currentCoorIterator))[_n]; }; // return X Y Z as atom[0], [1], [2] operators
inline double Atom::distance(const Atom & _atom) const {return CartesianGeometry::distance(*(*currentCoorIterator), *(*(_atom.currentCoorIterator)));};
inline double Atom::distance2(const Atom & _atom) const {return CartesianGeometry::distance2(*(*currentCoorIterator), *(*(_atom.currentCoorIterator)));};
inline double Atom::angle(const Atom & _atom) const {return CartesianGeometry::angle(*(*currentCoorIterator), *(*(_atom.currentCoorIterator)));};
inline double Atom::angle(const Atom & _center, const Atom & _third) const {return CartesianGeometry::angle(*(*currentCoorIterator), *(*(_center.currentCoorIterator)), *(*(_third.currentCoorIterator)));};
inline double Atom::angleRadians(const Atom & _atom) const {return CartesianGeometry::angleRadians(*(*currentCoorIterator), *(*(_atom.currentCoorIterator)));};
inline double Atom::angleRadians(const Atom & _center, const Atom & _third) const {return CartesianGeometry::angleRadians(*(*currentCoorIterator), *(*(_center.currentCoorIterator)), *(*(_third.currentCoorIterator)));};
inline double Atom::dihedral(const Atom & _second, const Atom & _third, const Atom & _fourth) const {return CartesianGeometry::dihedral(*(*currentCoorIterator), *(*(_second.currentCoorIterator)), *(*(_third.currentCoorIterator)), *(*(_fourth.currentCoorIterator)));};
inline double Atom::dihedralRadians(const Atom & _second, const Atom & _third, const Atom & _fourth) const {return CartesianGeometry::dihedralRadians(*(*currentCoorIterator), *(*(_second.currentCoorIterator)), *(*(_third.currentCoorIterator)), *(*(_fourth.currentCoorIterator)));};
inline bool Atom::hasCoor() const {return hasCoordinates;};
inline void Atom::wipeCoordinates() {(*currentCoorIterator)->setCoor(0.0, 0.0, 0.0); hasCoordinates = false;};
inline void Atom::setHasCoordinates(bool _flag) {hasCoordinates = _flag;};
inline std::vector<IcEntry*> & Atom::getIcEntries() {return icEntries;}
inline bool Atom::setActiveConformation(unsigned int _i) {
	if (_i<pCoorVec.size()) {
		currentCoorIterator = pCoorVec.begin() + _i;
		return true;
	}
	return false;
};
inline unsigned int Atom::getActiveConformation() const {return currentCoorIterator - pCoorVec.begin();};
inline unsigned int Atom::getNumberOfAltConformations(bool _absolute) const {
	if(_absolute) {
		// count also any hidden conformations
		return pCoorVec.size() + pHiddenCoorVec.size();
	} else {
		// only those that are not hidden
		return pCoorVec.size();
	}
}
inline void Atom::addAltConformation() {addAltConformation(getCoor());}; //default, same as current conformation
inline void Atom::addAltConformation(const CartesianPoint & _point) {unsigned int curr = currentCoorIterator - pCoorVec.begin(); pCoorVec.push_back(new CartesianPoint(_point)); currentCoorIterator = pCoorVec.begin() + curr;};  // it is important to regenerate the iterator to the current coor in case the std::vector is resized
inline void Atom::addAltConformation(Real _x, Real _y, Real _z) {addAltConformation(CartesianPoint(_x, _y, _z));};  // it is important to regenerate the iterator to the current coor in case the std::vector is resized
inline void Atom::setToStringFormat(unsigned int _format) {
	if (_format < 3) {
		toStringFormat = _format;
	}
}
// remove all alternate conformations, keeping the first conformation
inline void Atom::removeAllAltConformations() {
	if (pCoorVec.size() > 1) {
		CartesianPoint tmpCoor(**currentCoorIterator);
		for (std::vector<CartesianPoint*>::iterator k=pCoorVec.begin()+1; k!=pCoorVec.end(); k++) {
			delete *k;
			*k = NULL;
		}
		pCoorVec.erase(pCoorVec.begin()+1, pCoorVec.end());
		currentCoorIterator = pCoorVec.begin();
		(*currentCoorIterator)->setCoor(tmpCoor);
	}
	// clear any hidden alt-coors
	for (std::vector<CartesianPoint*>::iterator k=pHiddenCoorVec.begin(); k!=pHiddenCoorVec.end(); k++) {
		delete *k;
	}
	pHiddenCoorVec.clear();
	hiddenCoorIndeces.clear();

}
inline void Atom::saveCoor(std::string _coordName) {
	if(_coordName == "") {
		return;
	}
	if (savedCoor.find(_coordName) != savedCoor.end()) {
		// already existing, assign coordinates
		*(savedCoor[_coordName]) = **currentCoorIterator;
	} else {
		if (savedAltCoor.find(_coordName) != savedAltCoor.end()) {
			// if the name exists in the alt coor array, remove it
			clearSavedCoor(_coordName);
		}
		// create a new alt coor entry
		savedCoor[_coordName] = new CartesianPoint(**currentCoorIterator);
	}
}
inline void Atom::saveAltCoor(std::string _coordName) {
	if(_coordName == "") {
		return;
	}
	// clear the entry if pre-existing
	clearSavedCoor(_coordName);

	// presize the vector
	savedAltCoor[_coordName] = std::vector<CartesianPoint*>(pCoorVec.size(), (CartesianPoint*)NULL);
	//for (std::vector<CartesianPoint*>::iterator k=pCoorVec.begin(); k!=pCoorVec.end(); k++) {
	for (unsigned int i=0; i<pCoorVec.size(); i++) {
		savedAltCoor[_coordName][i] = new CartesianPoint(*pCoorVec[i]);
	}
	// save the index of the current coor
	savedAltCoorCurrent[_coordName] = currentCoorIterator - pCoorVec.begin();

	// save any hidden coordinates
	savedHiddenCoor[_coordName] = std::vector<CartesianPoint*>(pHiddenCoorVec.size(), (CartesianPoint*)NULL);
	savedHiddenCoorIndeces[_coordName] = std::vector<unsigned int>(pHiddenCoorVec.size(), 0);
	for (unsigned int i=0; i<pHiddenCoorVec.size(); i++) {
		savedHiddenCoor[_coordName][i] = new CartesianPoint(*pHiddenCoorVec[i]);
		savedHiddenCoorIndeces[_coordName][i] = hiddenCoorIndeces[i];
	}

}
inline bool Atom::applySavedCoor(std::string _coordName) {
	std::map<std::string, CartesianPoint*>::iterator found = savedCoor.find(_coordName);

	if (found != savedCoor.end()) {
		(*currentCoorIterator)->setCoor(*(found->second));
		return true;
	} else {
		std::map<std::string, std::vector<CartesianPoint*> >::iterator found2 = savedAltCoor.find(_coordName);
		if (found2 != savedAltCoor.end()) {
			// make sure that pCoorVector is resized correctly by adding or removing
			while (pCoorVec.size() < found2->second.size()) {
				pCoorVec.push_back(new CartesianPoint());
			}
			if (pCoorVec.size() > found2->second.size()) {
				for (std::vector<CartesianPoint*>::iterator k=pCoorVec.begin()+found2->second.size(); k != pCoorVec.end(); k++) {
					delete *k;
				}
				pCoorVec.erase(pCoorVec.begin()+found2->second.size(), pCoorVec.end());
			}
			// assign values
			for (unsigned int i=0; i<pCoorVec.size(); i++) {
				*(pCoorVec[i]) = *((found2->second)[i]);
			}
			setActiveConformation(savedAltCoorCurrent[_coordName]);

			// copy any hidden alt coords
			found2 = savedHiddenCoor.find(_coordName);
			if (found2 != savedHiddenCoor.end()) {
				// make sure that pHiddenCoorVec is resized correctly by adding or removing
				while (pHiddenCoorVec.size() < found2->second.size()) {
					pHiddenCoorVec.push_back(new CartesianPoint());
				}
				if (pHiddenCoorVec.size() > found2->second.size()) {
					for (std::vector<CartesianPoint*>::iterator k=pHiddenCoorVec.begin()+found2->second.size(); k != pHiddenCoorVec.end(); k++) {
						delete *k;
					}
					pHiddenCoorVec.erase(pHiddenCoorVec.begin()+found2->second.size(), pHiddenCoorVec.end());
				}
				// resize the container of the indeces
				hiddenCoorIndeces.resize(savedHiddenCoorIndeces[_coordName].size(), 0);
				// assign values
				for (unsigned int i=0; i<pHiddenCoorVec.size(); i++) {
					*(pHiddenCoorVec[i]) = *((found2->second)[i]);
					hiddenCoorIndeces[i] = savedHiddenCoorIndeces[_coordName][i];
				}
			}
			return true;
		}
	}
	return false;
}
//inline void Atom::setBondedTo(Atom * _pAtom, bool _bound) {if (_bound) {bonds[_pAtom] = true;} else {std::map<Atom*, bool>::iterator found=bonds.find(_pAtom); if (found!=bonds.end()) {bonds.erase(found);}}}
//inline std::vector<Atom*> Atom::getBoundAtoms() const {std::vector<Atom*> bonded; for (std::map<Atom*, bool>::const_iterator k=bonds.begin(); k!=bonds.end(); k++) {bonded.push_back(k->first);} return bonded;}
//inline std::vector<Atom*> Atom::getBoundAtoms() const {std::vector<Atom*> bonded; for (std::map<Atom*, std::map<Atom*, std::map<Atom*, bool> > >::const_iterator k=boundAtoms.begin(); k!=boundAtoms.end(); k++) {bonded.push_back(k->first);} return bonded;}
//inline std::map<Atom*, bool> & Atom::getBonds() {return bonds;}
inline std::vector<Atom*> Atom::getBonds() {
	std::vector<Atom*> bonded;
	for (std::map<Atom*, std::map<Atom*, std::map<Atom*, bool> > >::const_iterator k=boundAtoms.begin(); k!=boundAtoms.end(); k++) {
		bonded.push_back(k->first);
	}
	return bonded;
}

inline std::vector<std::vector<Atom*> > Atom::getBoundAtoms() const {
	std::vector<std::vector<Atom*> > bonded;
	for (std::map<Atom*, std::map<Atom*, std::map<Atom*, bool> > >::const_iterator k=boundAtoms.begin(); k!=boundAtoms.end(); k++) {
		bonded.push_back(std::vector<Atom*>(1, k->first));
		for (std::map<Atom*, std::map<Atom*, bool> >::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {
			bonded.back().push_back(l->first);
			for (std::map<Atom*, bool>::const_iterator m=l->second.begin(); m!=l->second.end(); m++) {
				bonded.back().push_back(m->first);
			}
		}
	}
	return bonded;
}
inline std::vector<Atom*> Atom::getOneThreeMiddleAtoms(Atom *_pAtom) const {
	std::vector<Atom*> middle;
	std::map<Atom*, std::map<Atom*, bool> >::const_iterator found=oneThreeAtoms.find(_pAtom);
	if (found != oneThreeAtoms.end()) {
		for (std::map<Atom*, bool>::const_iterator k=found->second.begin(); k!=found->second.end(); k++) {
			middle.push_back(k->first);
		}
	}
	return middle;
}
inline std::vector<std::vector<Atom*> > Atom::getOneFourMiddleAtoms(Atom *_pAtom) const {
	std::vector<std::vector<Atom*> > middle;
	std::map<Atom*, std::map<Atom*, std::map<Atom*, bool> > >::const_iterator found=oneFourAtoms.find(_pAtom);
	if (found != oneFourAtoms.end()) {
		for (std::map<Atom*, std::map<Atom*, bool> >::const_iterator k=found->second.begin(); k!=found->second.end(); k++) {
			middle.push_back(std::vector<Atom*>());
			middle.back().push_back(k->first); // second atom of the 1-4 relationship
			for (std::map<Atom*, bool>::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {
				middle.back().push_back(l->first); // 3rd atom
			}
		}
	}
	return middle;
}

//inline bool Atom::isBoundTo(Atom * _pAtom) const {std::map<Atom*, bool>::const_iterator found=bonds.find(_pAtom); return found!=bonds.end();}
//inline bool Atom::isOneThree(Atom * _pAtom, const Atom * _14caller) const {if (_pAtom == this) {return false;} for (std::map<Atom*, bool>::const_iterator k=bonds.begin(); k!=bonds.end(); k++) {if (k->first == _pAtom) {return false;} else if (k->first->isBoundTo(_pAtom) && k->first != _14caller) {return true;}} return false;}
//inline bool Atom::isOneFour(Atom * _pAtom) const {if (_pAtom == this) {return false;} for (std::map<Atom*, bool>::const_iterator k=bonds.begin(); k!=bonds.end(); k++) {if (k->first == _pAtom) {return false;} else if (k->first->isOneThree(_pAtom, this)) {return true;}} return false;}

inline bool Atom::isBoundTo(Atom * _pAtom) const { return boundAtoms.find(_pAtom) != boundAtoms.end(); }
inline bool Atom::isOneThree(Atom * _pAtom) const { return oneThreeAtoms.find(_pAtom) != oneThreeAtoms.end(); }
inline bool Atom::isOneFour(Atom * _pAtom) const { return oneFourAtoms.find(_pAtom) != oneFourAtoms.end(); }
inline bool Atom::isInAlternativeIdentity(Atom * _pAtom) const {return getParentPosition() == _pAtom->getParentPosition() && getParentResidue() != _pAtom->getParentResidue();}
inline double Atom::groupDistance(Atom & _atom, unsigned int _stamp) {
	return MSL::CartesianGeometry::distance(getGroupGeometricCenter(_stamp), _atom.getGroupGeometricCenter(_stamp));
}

inline std::string Atom::getAtomId(unsigned int _skip) const {
	return MslTools::getAtomId(getChainId(), getResidueNumber(), getResidueIcode(), getName(), _skip);
}

inline std::string Atom::getAtomOfIdentityId(unsigned int _skip) const {
	return MslTools::getAtomOfIdentityId(getChainId(), getResidueNumber(), getResidueIcode(), getResidueName(), getName(), _skip);
}
inline std::string Atom::getIdentityId(unsigned int _skip) const {
	return MslTools::getIdentityId(getChainId(), getResidueNumber(), getResidueIcode(), getResidueName(), _skip);
}
inline std::string Atom::getPositionId(unsigned int _skip) const {
	return MslTools::getPositionId(getChainId(), getResidueNumber(), getResidueIcode(), _skip);
}
inline void Atom::clearAllFlags() {
	Selectable<Atom>::clearAllFlags();
	setSelectionFlag("all",true);
}

inline int  Atom::getMinimizationIndex() { return minIndex; }
inline void Atom::setMinimizationIndex(int _index) { minIndex = _index; }


}

//}

#endif
