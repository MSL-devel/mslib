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

#ifndef BBQ_TABLE_H
#define BBQ_TABLE_H

#include <map>
#include "CoordAxes.h"
#include "AtomPointerVector.h"
#include "CartesianPoint.h"
#include "Residue.h"
#include "Chain.h"

typedef std::pair<Residue *, Residue *> ResiduePtrPair;


/**
 * This class holds the information needed for the BBQ (Backbone Building from
 * Quadrilaterals) algorithm.  It is basically a 3-D LUT, where each dim
 * is a coordinate in the R-coordinate system for a given C-alpha atom
 * in a C-alpha backbone chain.  See "Backbone Building from Quadrilaterals:..."
 * by Dominik Gront, Sebastian Kmiecik, and Andrezej Kolinski in the
 * Journal of Compuatational Chemistry Vol 28: 1593-1597, 2007 for more details.
 */
class BBQTable : public std::map<CartesianPoint, AtomPointerVector *> {
public:
    BBQTable();
    BBQTable(string _bbqTableFileName);
    BBQTable(const BBQTable &_table);
    ~BBQTable();

    void setBinSizes(Real x, Real y, Real z);
    void getBinSizes(Real &x, Real &y, Real &z);

    void fillInMissingBBAtoms(std::vector<Residue *> &_rv);

    int fillInMissingBBAtoms(Chain &_ch);
    //void addQuadrilateralInfoFromResidues(std::vector<Residue *> &_rv, std::map<string, bool> &_atomsOfInterest);
    void addQuadrilateralInfoFromResidues(std::vector<Residue *> &_rv);

    
    void normalize();

    void setDebugFlag(bool _flag) { debugFlag = _flag; };

    void deleteTableEntries();
    void openReader(string _bbqTableFileName);
private:
    void addAtomsToResidue(Real _r02, Real _r03, Real _r13, CoordAxes &_axes, Residue &_res);
    void addAtomPointerVector(Real _r02, Real _r03, Real _r13, AtomPointerVector *_av, CoordAxes &_axes);
    
    void calcRDistances(std::vector<Residue *> &_rv, std::map<ResiduePtrPair, Real> &_rDistances);
    void calcRDistances(Chain &_ch, std::map<ResiduePtrPair, Real> &_rDistances);

    void calcLCoords(std::vector<Residue *> &_rv, std::map<Residue *, CoordAxes> &_lCoords);
    void calcLCoords(Chain &_ch, std::map<Residue *, CoordAxes> &_lCoords);


    void sumAtomPointerVectors(AtomPointerVector *_av1, AtomPointerVector *_av2);
    double calcCADihedral(Residue *pRes0, Residue *pRes1, Residue *pRes2, Residue *pRes3);
    bool isLegalQuad(Residue *pRes0, Residue *pRes1, Residue *pRes2, Residue *pRes3);

    void findClosestTableEntry(Residue &_res, CartesianPoint &_key);
    bool search3x3Neighborhood(Residue &_res, CartesianPoint &_key);
    void bruteForceFindClosestTableEntry(Residue &_res, CartesianPoint &_key);
    bool doCADistanceCheck(Residue *pRes0, Residue *pRes1);
    bool doAllFourResiduesHaveGivenAtom(Residue *pRes0, Residue *pRes1, Residue *pRes2, Residue *pRes3, string atomName);

    unsigned int dims[3];
    Real binSizes[3];
    bool debugFlag;
    std::map<CartesianPoint, unsigned int> counts;
};

// inlines
inline BBQTable::BBQTable() {
    debugFlag = false;
}
/**
 * This function will delete all AtomPointerVectors in our table.
 **/
inline void BBQTable::deleteTableEntries() {
    // The table owns the AtomPointerVectors it contains, as well as the Atoms in the atom
    // vectors, so delete them when the table is deleted.
    for(std::map<CartesianPoint, AtomPointerVector *>::iterator currAVIter = begin(); currAVIter != end(); ++currAVIter) {
        AtomPointerVector *currAv = currAVIter->second;

        for(AtomPointerVector::iterator currAtomIter = currAv->begin(); currAtomIter != currAv->end(); ++currAtomIter) {
            delete *currAtomIter;
        }

        delete currAv;
    }
}
/**
 * This is the destructor for the BBQTable object.
 */
inline BBQTable::~BBQTable() {
    deleteTableEntries();
    clear();
}

/**
 * This function will set the size of the bins
 * of our table (in Angstroms).
 */
inline void BBQTable::setBinSizes(Real x, Real y, Real z) {
    binSizes[0] = x;
    binSizes[1] = y;
    binSizes[2] = z;
};

inline void BBQTable::getBinSizes(Real &x, Real &y, Real &z) {
    x = binSizes[0];
    y = binSizes[1];
    z = binSizes[2];
}

#endif // BBQ_TABLE_H
