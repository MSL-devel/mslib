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


#include "BBQTable.h"
#include "Frame.h"
#include "Matrix.h"
#include "MslTools.h"
#include "CartesianGeometry.h"
#include "BBQTableReader.h"
#include "System.h"

using namespace std;

using namespace MSL;


#define CA_CA_DISTANCE          ((Real)3.78f)
#define CA_CA_DISTANCE_PROLINE  ((Real)3.2f)
#define CA_CA_TOLERANCE         ((Real)0.16f)
//#define CA_CA_TOLERANCE         ((Real)0.50f)


BBQTable::BBQTable(string _bbqTableFileName){
	openReader(_bbqTableFileName);


}
void BBQTable::openReader(string _bbqTableFileName){

    BBQTableReader bbqTableReader(_bbqTableFileName);
    bbqTableReader.open();
    bbqTableReader.read(*this);
    bbqTableReader.close();
}



/**
 * This function will fill in the missing N, C, and O atoms for
 * all Residues in the given vector.
 *
 * @param _rv  A vector of Residues with only C-alpha's given.  This function
 *             will update the Residues to give coord positions of N, C, and O.
 **/

void BBQTable::fillInMissingBBAtoms(vector<Residue *> &_rv) {
    map<ResiduePtrPair, Real> rDistances;
    map<Residue *, CoordAxes> lCoords;

    calcRDistances(_rv, rDistances);
    calcLCoords(_rv, lCoords);

    for(int currIndex = 0; currIndex <= (int)_rv.size()-4; ++currIndex) {
        // Skip this if it is not a legal quadrilateral
        if(isLegalQuad( _rv[currIndex], _rv[currIndex+1], _rv[currIndex+2], _rv[currIndex+3]) == false)
            continue;
        
        Real distances[3];
        CoordAxes axes = lCoords[_rv[currIndex]];

        distances[0] = rDistances[ResiduePtrPair(_rv[currIndex], _rv[currIndex+2])];
        // We store the information regarding whether the dihedral is positive
        // or negative on the distance metric between the current atom and current atom + 2.
        // ignore the signs on other distances.
        distances[1] = fabs(rDistances[ResiduePtrPair(_rv[currIndex], _rv[currIndex+3])]);
        distances[2] = fabs(rDistances[ResiduePtrPair(_rv[currIndex+1], _rv[currIndex+3])]);

        Residue newRes;
        addAtomsToResidue(distances[0], distances[1], distances[2], axes, newRes);
        // The position of the C, O, and N atoms are relative to their C-alpha atoms.
        newRes.getAtom("C").getCoor() += _rv[currIndex+1]->getAtom("CA").getCoor();
        newRes.getAtom("O").getCoor() += _rv[currIndex+1]->getAtom("CA").getCoor();
        newRes.getAtom("N").getCoor() += _rv[currIndex+2]->getAtom("CA").getCoor();

        // The C and O atoms belong to residue 1, while the N atom belongs to residue 2.
        _rv[currIndex+1]->addAtom( newRes.getAtom("C") );
        _rv[currIndex+1]->addAtom( newRes.getAtom("O") );
        _rv[currIndex+2]->addAtom( newRes.getAtom("N") );
    }
}

int BBQTable::fillInMissingBBAtoms(Chain &_chain) {

	if (_chain.positionSize()== 0){
		cerr << "WARNING chain has no residues to fillInBBAtoms for. doing nothing.\n";
		return -1;
	}

	map<ResiduePtrPair, Real> rDistances;
	map<Residue *, CoordAxes> lCoords;

	calcRDistances(_chain, rDistances);
	calcLCoords(_chain, lCoords);

	int illegalQuads = 0;
	for(int currIndex = 0; currIndex < (int)_chain.positionSize()-3; currIndex++) {


		AtomPointerVector ats;
		if (isLegalQuad( &_chain.getResidue(currIndex), &_chain.getResidue(currIndex+1),&_chain.getResidue(currIndex+2),  &_chain.getResidue(currIndex+3)) == false){
			fprintf(stdout,"ILLEGAL QUAD %d\n",currIndex);
			illegalQuads++;
			continue;
		} else {

			Real distances[3];
			CoordAxes axes = lCoords[&_chain.getResidue(currIndex)];
			CartesianPoint cAlphaCoord = _chain.getResidue(currIndex).getAtom("CA").getCoor();

			distances[0] = rDistances[ResiduePtrPair(&_chain.getResidue(currIndex), &_chain.getResidue(currIndex+2))];
			// We store the information regarding whether the dihedral is positive
			// or negative on the distance metric between the current atom and current atom + 2.
			// ignore the signs on other distances.
			distances[1] = fabs(rDistances[ResiduePtrPair(&_chain.getResidue(currIndex), &_chain.getResidue(currIndex+3))]);
			distances[2] = fabs(rDistances[ResiduePtrPair(&_chain.getResidue(currIndex+1), &_chain.getResidue(currIndex+3))]);

			Residue newRes;

			//AtomPointerVector av = getAtomPointerVector(distances[0], distances[1], distances[2], axes, cAlphaCoord);
			addAtomsToResidue(distances[0], distances[1], distances[2], axes, newRes);

			newRes.getAtom("C").getCoor() += _chain.getResidue(currIndex+1).getAtom("CA").getCoor();
			newRes.getAtom("O").getCoor() += _chain.getResidue(currIndex+1).getAtom("CA").getCoor();
			newRes.getAtom("N").getCoor() += _chain.getResidue(currIndex+2).getAtom("CA").getCoor();

			/*
			  newRes.setChainId(_chain.getResidue(currIndex).getChainId());
			  newRes.setResidueName(_chain.getResidue(currIndex).getResidueName());
			  newRes.setResidueNumber(_chain.getResidue(currIndex).getResidueNumber());
			  newRes.setResidueIcode(_chain.getResidue(currIndex).getResidueIcode());
			*/


			ats.push_back(new Atom(newRes.getAtom("C")));
			ats.back()->setChainId(_chain.getResidue(currIndex+1).getChainId());
			ats.back()->setResidueName(_chain.getResidue(currIndex+1).getResidueName());
			ats.back()->setResidueNumber(_chain.getResidue(currIndex+1).getResidueNumber());
			ats.back()->setResidueIcode(_chain.getResidue(currIndex+1).getResidueIcode());
			ats.back()->setElement("C");

			ats.push_back(new Atom(newRes.getAtom("O")));
			ats.back()->setChainId(_chain.getResidue(currIndex+1).getChainId());
			ats.back()->setResidueName(_chain.getResidue(currIndex+1).getResidueName());
			ats.back()->setResidueNumber(_chain.getResidue(currIndex+1).getResidueNumber());
			ats.back()->setResidueIcode(_chain.getResidue(currIndex+1).getResidueIcode());
			ats.back()->setElement("O");

			ats.push_back(new Atom(newRes.getAtom("N")));
			ats.back()->setChainId(_chain.getResidue(currIndex+2).getChainId());
			ats.back()->setResidueName(_chain.getResidue(currIndex+2).getResidueName());
			ats.back()->setResidueNumber(_chain.getResidue(currIndex+2).getResidueNumber());
			ats.back()->setResidueIcode(_chain.getResidue(currIndex+2).getResidueIcode());
			ats.back()->setElement("N");
		}
		
		if (_chain.getParentSystem() != NULL){
			_chain.getParentSystem()->addAtoms(ats);
		} else {
			_chain.addAtoms(ats);
		}

		ats.deletePointers();

		// Force updating of atomvectors in chain/system
		_chain.getPosition(currIndex).setActiveIdentity(_chain.getResidue(currIndex).getResidueName());

	}
        
        for(int currIndex = (int)_chain.positionSize()-3; currIndex < (int)_chain.positionSize(); currIndex++)
            _chain.getPosition(currIndex).setActiveIdentity(_chain.getResidue(currIndex).getResidueName());


	return illegalQuads;
}

/**
 * This function is used when building a BBQTable (not simply reading from file.)
 * The user will pass in the given r coordinates, and the vector of atoms associated
 * with it (typicall the N, C, and O atoms).  It will also pass in the offset of the
 * current C-alpha atom.  This will be subtracted off all of atoms so that they
 * are now relative to the C-alpha.  Finally, the user will also supply the axes of
 * this quadrilateral.  The atoms will be rotated to a standard reference frame.
 *
 * @param r02           The distance between the 0th and 2nd c-alpha.
 * @param r03           The distance between the 0th and 3rd c-alpha.
 * @param r13           The distance between the 1st and 3rd c-alpha.
 * @param av            The vector of atoms to add.  Typically this vector includes C, O, and N.
 * @param axes          The axes relative to standard x, y, and z.
 */
void BBQTable::addAtomPointerVector(Real _r02, Real _r03, Real _r13, AtomPointerVector *_av, CoordAxes &_axes) {
    CartesianPoint key;
    AtomPointerVector *oldAv;
    Frame newFrame;

    newFrame.computeFrameFromAxes(_axes);

    // Transform to global basis.
    newFrame.transformToGlobalBasis(*_av);

    try {
        // First see if we already have an entry for this combination of r coords.
        map<CartesianPoint, AtomPointerVector *>::iterator searchResult;

        key.setX( MslTools::round( _r02 / binSizes[0] ) );
        key.setY( MslTools::round( _r03 / binSizes[1] ) );
        key.setZ( MslTools::round( _r13 / binSizes[2] ) );

        searchResult = find(key);
        if(searchResult != end()) {
            oldAv = searchResult->second;
            sumAtomPointerVectors(oldAv, _av);
            counts[key] += 1;
        }
        else {
            // Create new atom vector, copying over the atoms just passed in.
            // The BBQTable owns all elements in the table and will be responsible
            // for deleting them later on.
            // Create a new atom vector for this key.  Create new Atoms too
            // so that BBQTable owns all objects that it holds.
            AtomPointerVector *avCopy = new AtomPointerVector();
            for(AtomPointerVector::iterator currIter = _av->begin(); currIter != _av->end(); ++currIter) {
                avCopy->push_back( new Atom(**currIter) );
            }

            insert( pair<CartesianPoint, AtomPointerVector *>(key, avCopy) );
            counts[key] = 1;
        }
    }
    catch (...) {
        cout << "Error in BBQTable::addAtomPointerVector.\n";
        exit(1);
    }
}


/**
 * This function will get a vector of atoms from the bin
 * closest to the given r coordinates (the distance between
 * c-alpha's for residues 0 and 2, 0 and 3, and 1 and 3).  These
 * atoms will be added to the residue provided.
 *
 * @param r02           The distance between the 0th and 2nd c-alpha.
 * @param r03           The distance between the 0th and 3rd c-alpha.
 * @param r13           The distance between the 1st and 3rd c-alpha.
 * @param axes          The axes relative to standard x, y, and z.
 * @param res           The residue to which we should add the atoms.  Note: this should be an empty residue.
 *                      We will place the C and O atoms for residue 1 in our quadrilateral, and the N atom for
 *                      residue 2 in our quadrilateral in this residue, and return it.
 */
void BBQTable::addAtomsToResidue(Real _r02, Real _r03, Real _r13, CoordAxes &_axes, Residue &_res) {
    CartesianPoint key;
    AtomPointerVector av;
    Frame newFrame;
    
    newFrame.computeFrameFromAxes(_axes);
    try {
        map<CartesianPoint, AtomPointerVector *>::iterator searchResult;

        key.setX( MslTools::round( _r02 / binSizes[0] ) );
        key.setY( MslTools::round( _r03 / binSizes[1] ) );
        key.setZ( MslTools::round( _r13 / binSizes[2] ) );

        searchResult = find(key);
        if(searchResult != end()) 
            _res.addAtoms(*(searchResult->second));
        else 
            findClosestTableEntry(_res, key);

        // Get the new atom vector from the residue.
        AtomPointerVector &resAV = _res.getAtomPointers();
        // Rotate atoms around the given axes.
        newFrame.transformFromGlobalBasis(resAV);
    }
    catch (...) {
        cout << "Error in BBQTable::getAtomPointerVector.\n";
        exit(1);
    }
}

/**
 * This function will normalize the atom vectors that have been 
 * added to our table.
 */
void BBQTable::normalize() {
    // Loop over all of the counts that we have.
    for(map<CartesianPoint, unsigned int>::iterator currIter = counts.begin(); currIter != counts.end(); ++currIter) {
        // Find the atom vector that corresponds to this count entry.
        AtomPointerVector *av = find(currIter->first)->second;
        
        // Normalize this atom vector.
        for( AtomPointerVector::iterator currAVIter = av->begin(); currAVIter != av->end(); ++currAVIter ) {
            Atom *currAtom = *currAVIter;
            currAtom->setCoor( currAtom->getCoor() / (Real) currIter->second );
        }
    }
}

/**
 * This function loops over all quads of Residues in the supplied vector, and
 * calculate the distances between C-alpha 0 and C-alpha 2, C-alpha 0, and
 * C-alpha 3, and C-alpha 1 and C-alpha 3.  These three distances are called
 * the R-Distances, and they will be used to index our BBQ Table.
 * Note that for each distance, the sign indicates whether the dihedral angle
 * was positive or negative for the set of four points starting with the first residue
 * in that key.
 *
 * @param _rv  A vector of Residues with only C-alpha values given.
 * @param _rDistances A map of Residue pointer pairs to reals.  This is where
 *                    we will put our output.
 */
void BBQTable::calcRDistances(vector<Residue *> &_rv, map<ResiduePtrPair, Real> &_rDistances) {
    // We need at least 4 C-alpha's...
    if( _rv.size() < 4 )
        return;
    int currIndex = -1;
    Residue * pRes0 = NULL;
    Residue * pRes1 = NULL;
    Residue * pRes2 = NULL;
    Residue * pRes3 = NULL;
    bool success = false;

    // Find the first legal quadrilateral.
    while( (success == false) && (currIndex < (int)_rv.size()-4) )  {
        ++currIndex;
        pRes0 = _rv[currIndex];
        pRes1 = _rv[currIndex+1];
        pRes2 = _rv[currIndex+2];
        pRes3 = _rv[currIndex+3];

        success = isLegalQuad(pRes0, pRes1, pRes2, pRes3);
    }
    // If we reached the end of our residue vector before finding 
    // a legal quadrilateral, get out of here.
    if(currIndex > ((int)_rv.size()-4) )
        return;

    ResiduePtrPair key(pRes0, pRes2);
    _rDistances[key] = pRes0->distance(*pRes2, "CA");

    // Loop over the residues in the residue vector calculating the r distances.
    for(; currIndex <= (int)_rv.size()-4; ++currIndex){
        pRes0 = _rv[currIndex];
        pRes1 = _rv[currIndex + 1];
        pRes2 = _rv[currIndex + 2];
        pRes3 = _rv[currIndex+ 3];
        
        // Skip if this is not a legal quadrilateral.
        if( isLegalQuad(pRes0, pRes1, pRes2, pRes3) == false )
            continue;

        key.first = pRes0;
        key.second = pRes3;
        _rDistances[key] = pRes0->distance(*pRes3, "CA");
        
        key.first = pRes1;
        key.second = pRes3;
        _rDistances[key] = pRes1->distance(*pRes3, "CA");

        // If the dihedral angle is positive, then leave
        // the distance between 0 and 2 positive.  If the
        // dihedral is negative, then make it negative.
        key.first = pRes0;
        key.second = pRes2;
        if(calcCADihedral(pRes0, pRes1, pRes2, pRes3) < 0)
            _rDistances[key] = -_rDistances[key];
    }
}

void BBQTable::calcRDistances(Chain &_ch, map<ResiduePtrPair, Real> &_rDistances) {
    // We need at least 4 C-alpha's...
    if( _ch.positionSize() < 4 )
        return;

    Residue *pRes0 = &_ch.getResidue(0), *pRes1 = NULL, *pRes2 = &_ch.getResidue(2), *pRes3 = NULL;

    ResiduePtrPair key(pRes0, pRes2);
    _rDistances[key] = pRes0->distance(*pRes2, "CA");

    // Loop over the residues in the residue vector calculating the r distances.
    for(int currIndex = 0; currIndex <= (int)_ch.positionSize()-4; ++currIndex){
        pRes0 = &_ch.getResidue(currIndex);
        pRes1 = &_ch.getResidue(currIndex + 1);
        pRes2 = &_ch.getResidue(currIndex + 2);
        pRes3 = &_ch.getResidue(currIndex + 3);

        key.first = pRes0;
        key.second = pRes3;
        _rDistances[key] = pRes0->distance(*pRes3, "CA");
        
        key.first = pRes1;
        key.second = pRes3;
        _rDistances[key] = pRes1->distance(*pRes3, "CA");

        // If the dihedral angle is positive, then leave
        // the distance between 0 and 2 positive.  If the
        // dihedral is negative, then make it negative.
        key.first = pRes0;
        key.second = pRes2;
        if(calcCADihedral(pRes0, pRes1, pRes2, pRes3) < 0)
            _rDistances[key] = -_rDistances[key];
    }
}

/**
 * This function will loop through all of the residues in the
 * given residue vector and form quadrilaterals.  It will then add entries
 * into our BBQTable for each of the quadrilaterals it sees.
 *
 * @param _rv               A vector of Residues with all backbone atoms specified.
 **/
void BBQTable::addQuadrilateralInfoFromResidues(vector<Residue *> &_rv) {
    map<ResiduePtrPair, Real> rDistances;
    map<Residue *, CoordAxes> lCoords;

    calcRDistances(_rv, rDistances);
    calcLCoords(_rv, lCoords);

    for(int currIndex = 0; currIndex <= (int)_rv.size()-4; ++currIndex) {
        // Skip this if it is not a legal quadrilateral, or if any of the 4 residues
        // don't have the "C", "O", and "N" backbone atoms.
        if(isLegalQuad( _rv[currIndex], _rv[currIndex+1], _rv[currIndex+2], _rv[currIndex+3]) == false)
            continue;
        if( (doAllFourResiduesHaveGivenAtom(_rv[currIndex], _rv[currIndex+1], _rv[currIndex+2], _rv[currIndex+3], "C") &&
             doAllFourResiduesHaveGivenAtom(_rv[currIndex], _rv[currIndex+1], _rv[currIndex+2], _rv[currIndex+3], "O") &&
             doAllFourResiduesHaveGivenAtom(_rv[currIndex], _rv[currIndex+1], _rv[currIndex+2], _rv[currIndex+3], "N") ) == false )
            continue;

        Real distances[3];

        distances[0] = rDistances[ResiduePtrPair(_rv[currIndex], _rv[currIndex+2])];
        // We store the information regarding whether the dihedral is positive
        // or negative on the distance metric between the current atom and current atom + 2.
        // ignore the signs on other distances.
        distances[1] = fabs(rDistances[ResiduePtrPair(_rv[currIndex], _rv[currIndex+3])]);
        distances[2] = fabs(rDistances[ResiduePtrPair(_rv[currIndex+1], _rv[currIndex+3])]);

        // Copy over backbone atoms in this residue.
        Residue resCopy;
        resCopy.addAtom( _rv[currIndex+1]->getAtom("C") );
        resCopy.addAtom( _rv[currIndex+1]->getAtom("O") );
        resCopy.addAtom( _rv[currIndex+2]->getAtom("N") );

        CoordAxes axes = lCoords[_rv[currIndex]];
        // Get the position of the C, O, and N atoms relative to their C-alpha atom.
        resCopy.getAtom("C").getCoor() -= _rv[currIndex+1]->getAtom("CA").getCoor();
        resCopy.getAtom("O").getCoor() -= _rv[currIndex+1]->getAtom("CA").getCoor();
        resCopy.getAtom("N").getCoor() -= _rv[currIndex+2]->getAtom("CA").getCoor();
        addAtomPointerVector(distances[0], distances[1], distances[2], &resCopy.getAtomPointers(), axes);
    }
}

/**
 * This function loops over all quads of Residues in the supplied vector, and
 * calculates a coordinate axes for each Residue.  This coordinate axis
 * is defined as the unit vector created by summing unit vectors pointing from
 * C-alpha 0 to C-alpha 2 and C-alpha 0 to C-alpha 3.  The second unit vector
 * is created by taking the difference of the unit vectors pointing from these
 * same pairs.  The 3rd unit vector is simply the cross product of these first
 * two unit vectors.
 *
 * @param _rv  A vector of Residues with only C-alpha values given.
 * @param _lCoords A map of residue pointers to coordinate axes.  This is
 *                 where we'll put the output.
 */
void BBQTable::calcLCoords(vector<Residue *> &_rv, map<Residue *, CoordAxes> &_lCoords) {
    // We need at least 3 C-alpha's...
    if( _rv.size() < 3 )
        return;
    
    for(int currIndex = 0; currIndex <= (int)_rv.size()-3; ++currIndex) {
        CartesianPoint points[3], versor01, versor02;
        CoordAxes axes;

        // Make sure all three residues have a C-alpha. Otherwise, move on.
        if( _rv[currIndex]->atomExists("CA") && _rv[currIndex+1]->atomExists("CA") && _rv[currIndex+2]->atomExists("CA")) {
            points[0] = _rv[currIndex]->getAtom("CA").getCoor();
            points[1] = _rv[currIndex+1]->getAtom("CA").getCoor();
            points[2] = _rv[currIndex+2]->getAtom("CA").getCoor();
        }
        else
            continue;
        
        versor01 = points[1] - points[0];
        versor01 = versor01.getUnit();

        versor02 = points[2] - points[0];
        versor02 = versor02.getUnit();

        axes.first = versor01 + versor02;
        axes.first = axes.first.getUnit();
        axes.second = versor01 - versor02;
        axes.second = axes.second.getUnit();
        axes.third = axes.first.cross(axes.second);

        _lCoords[_rv[currIndex]] = axes;
    }
}

void BBQTable::calcLCoords(Chain &_ch, map<Residue *, CoordAxes> &_lCoords) {
    // We need at least 3 C-alpha's...
    if( _ch.positionSize() < 3 )
        return;
    
    for(int currIndex = 0; currIndex <= (int)_ch.positionSize()-3; ++currIndex) {
        CartesianPoint points[3], versor01, versor02;
        CoordAxes axes;

        points[0] = _ch.getResidue(currIndex).getAtom("CA").getCoor();
        points[1] = _ch.getResidue(currIndex+1).getAtom("CA").getCoor();
        points[2] = _ch.getResidue(currIndex+2).getAtom("CA").getCoor();
        
        versor01 = points[1] - points[0];
        versor01 = versor01.getUnit();

        versor02 = points[2] - points[0];
        versor02 = versor02.getUnit();

        axes.first = versor01 + versor02;
        axes.first = axes.first.getUnit();
        axes.second = versor01 - versor02;
        axes.second = axes.second.getUnit();
        axes.third = axes.first.cross(axes.second);

        _lCoords[&_ch.getResidue(currIndex)] = axes;
    }
}
/**
 * This function will take the coordinates from all of the atoms in
 * _av2 and add them to the corresponding atom's coords in _av1.
 * It expects there to be a 1:1 correspondence between the atoms in
 * _av1 and _av2, although the atoms can appear in different orders.
 *
 * @param _av1     The first atom vector, into which the sums will be placed.
 * @param _av2     The second atom vector, which we will add to the first.
 */
void BBQTable::sumAtomPointerVectors(AtomPointerVector *_av1, AtomPointerVector *_av2) {
    map<string, int> av1Map;
    // Look through all of the atoms in the first atom vector and make
    // a mapping of the atom name to the index in the atom vector.
    for(int currIndex = 0; currIndex < (int)_av1->size(); ++currIndex) {
        av1Map[(*_av1)[currIndex]->getName()] = currIndex;
    }

    // Now loop over all of the atoms in the second atom vector, see
    for(AtomPointerVector::iterator currIter = _av2->begin(); currIter != _av2->end(); ++currIter) {
        Atom *currAtom = *currIter;
        map<string, int>::iterator searchResult = av1Map.find( currAtom->getName() );
        
        // Check to see if we have this atom in av1 already.  
        // I would think that we always would have the atom, but you never know.
        if( searchResult != av1Map.end() ) {
            Atom *origAtom = (*_av1)[ searchResult->second ];
            origAtom->setCoor( origAtom->getCoor() + currAtom->getCoor() );
        }
        else {
            // See, I would just tack this atom onto the atom vector, because why not?
            // Well, because we keep one global counter for the number of times we have
            // added to this atom vector, and then we divide by that number to get the average
            // value.  If some atoms are seen different number of times, this won't work.
            cout << "For some reason, when adding the current AtomPointerVector,";
            cout << " we saw an atom that we hadn't seen before - " << currAtom->getName() << ".  Sorry, but I'm bailing out.\n";
            exit(1);
        }
    }
}

/**
 * This function will attempt to find the best atom vector for a set of R distances
 * that don't have an entry in our table.  Basically, since each entry is a key into
 * a 3-D LUT, this function will first look in the 26 neighboring entries to the given
 * key, and return the first entry that it finds.  If none of the 26 neighboring keys
 * have entries either, it will then do a brute force calculation, looking to see what
 * key that is in the table minimizes the distance between keys.
 *
 * @param _res  This is the Residue where we will place the atom vector that best matches the key.
 * @param _key This is the key that we wanted to use, but doesn't have an entry.
 */
void BBQTable::findClosestTableEntry(Residue &_res, CartesianPoint &_key) {
    if( !search3x3Neighborhood(_res, _key) )
        bruteForceFindClosestTableEntry(_res, _key);
}

/**
 * This function will look at all entries in the 3x3 cube surrounding the
 * given key.  It is called by findClosestTableEntry, which is in turn
 * called when there is no entry for the given key.  In searching the 3x3
 * area around this key, we are searching bins that are have r distances
 * very similar to the r distances we desire.  At the first entry this
 * function finds in the 3x3 neighborhood, it will set the atom vector
 * to the atom vector in this entry, and return true.  In the future we could
 * look into taking an average of all atom vectors found in the 3x3 neighborhood,
 * but since this function should be called relatively infrequently, its probably
 * not worth the effort.
 *
 * @param _res This is the Residue where we will place the atom vector that best matches the key.
 * @param _key This is the key that we wanted to use, but doesn't have an entry.
 * @return We return true if we found an entry in the 3x3 neighborhood, false otherwise.
 */
bool BBQTable::search3x3Neighborhood(Residue &_res, CartesianPoint &_key) {
    CartesianPoint newKey;
    map<CartesianPoint, AtomPointerVector *>::iterator searchResult;

    // Loop over the 3x3 region surrounding the key which did not
    // have an entry.  If any of these indices has an entry,
    // return that atom vector.  (Should we instead average
    // the results of all neighbors that have an entry?)
    // If none fo the 3x3 neighbors has an entry, return false.
    for(int x = int(_key.getX())-1; x <= int(_key.getX())+1; ++x) {
        for(int y = int(_key.getY())-1; y <= int(_key.getY())+1; ++y) {
            for(int z = int(_key.getZ())-1; z <= int(_key.getZ())+1; ++z){
                newKey.setCoor(x, y, z);
                searchResult = find(newKey);

                if( searchResult != end() ) {
                    _res.addAtoms( *(searchResult->second) );
                    return true;
                }
            }
        }
    }

    return false;
}

/**
 * This function is called as a last resort to find the best entry in our
 * LUT.  Basically, first we look for an entry with the given key.  If that doesn't
 * work, then we call search3x3Neighborhood and look for an entry in any of the
 * 26 neighboring locations (3x3 neighborhood).  If none of those has an entry,
 * then we iterate through all entries in our table, and calculate the cartesian
 * distance between the key for that entry and our desired key.  We return the value
 * in the entry which minimizes this distance.
 *
 * @param _res This is the Residue where we will place the atom vector that best matches the key.
 * @param _key This is the key that we wanted to use, but doesn't have an entry.
 */
void BBQTable::bruteForceFindClosestTableEntry(Residue &_res, CartesianPoint &_key) {
    Real minDistance = MslTools::floatMax;
    map<CartesianPoint, AtomPointerVector *>::iterator searchResult;

    // Loop over all entries in our table and find the entry that minimizes
    // the squared distance between the given key and the key of that entry.
    for(BBQTable::iterator currIter = begin(); currIter != end(); ++currIter) {
        CartesianPoint currKey = currIter->first;
        Real newDistance = currKey.distance(_key);

        if(newDistance < minDistance) {
            minDistance = newDistance;
            searchResult = currIter;
        }
    }
    
    _res.addAtoms( *(searchResult->second) );
}

/**
 * This helper function will calculate the dihedral angle between the C-alpha
 * atoms found in the four given residues.
 *
 * @param pRes0  The first residue in the chain.
 * @param pRes1  The second residue in the chain.
 * @param pRes2  The third residue in the chain.
 * @param pRes3  The fourth residue in the chain.
 *
 * @return The dihedral angle formed from connecting the 4 C-alpha atoms in the
 *         given residues.
 */
double BBQTable::calcCADihedral(Residue *pRes0, Residue *pRes1, Residue *pRes2, Residue *pRes3) {
    Atom &atom0 = pRes0->getAtom("CA"), &atom1 = pRes1->getAtom("CA"), &atom2 = pRes2->getAtom("CA"), &atom3 = pRes3->getAtom("CA");

    return(CartesianGeometry::dihedral(atom0.getCoor(), atom1.getCoor(), atom2.getCoor(), atom3.getCoor()));
}

/**
 * This function will look at the four residues and decide if it is a legal
 * quadrilateral or not.  For instance, all four residues must have a C-alpha,
 * and the distance between consecutive C-alpha's must be within some set tolerance.
 *
 * @param pRes0  The first residue in the chain.
 * @param pRes1  The second residue in the chain.
 * @param pRes2  The third residue in the chain.
 * @param pRes3  The fourth residue in the chain.
 *
 * @return A bool indicating whether this is a legal quadrilateral or not.
 */
bool BBQTable::isLegalQuad(Residue *pRes0, Residue *pRes1, Residue *pRes2, Residue *pRes3) {
    bool allResHaveCA, caDistancesAcceptable = true;

    // Do all 4 residues have a C-alpha?
    allResHaveCA = doAllFourResiduesHaveGivenAtom(pRes0, pRes1, pRes2, pRes3, "CA");
    // If they don't all have a C-alpha, then this is not a valid quad.
    if(!allResHaveCA) {
	    cout << "CANT FIND ALL CAAAAAAAAAAAAAs\n";
        return false;
    }


    // Make sure that distance between C-alphas is copacetic.
    caDistancesAcceptable = caDistancesAcceptable && doCADistanceCheck(pRes0, pRes1);
    if (!caDistancesAcceptable){
	    cout << "RES1 - RES2: "<<pRes0->getResidueNumber()<<" "<<pRes1->getResidueNumber()<<endl;
    }
    caDistancesAcceptable = caDistancesAcceptable && doCADistanceCheck(pRes1, pRes2);
    if (!caDistancesAcceptable){
	    cout << "RES2 - RES3: "<<pRes1->getResidueNumber()<<" "<<pRes2->getResidueNumber()<<endl;
    }
    caDistancesAcceptable = caDistancesAcceptable && doCADistanceCheck(pRes2, pRes3);
    if (!caDistancesAcceptable){
	    cout << "RES3 - RES4: "<<pRes2->getResidueNumber()<<" "<<pRes3->getResidueNumber()<<endl;
    }
    return caDistancesAcceptable;
}

/**
 * This simple method just checks each of the four residues for the atom of the given name.
 *
 * @param pRes0, pRes1, pRes2, pRes3  The four residues.
 * @param atomName The name of the atom we are interested in.
 * @return True if all four residues have the given atom, false otherwise.
 */
bool BBQTable::doAllFourResiduesHaveGivenAtom(Residue *pRes0, Residue *pRes1, Residue *pRes2, Residue *pRes3, string atomName) {
    bool allFourResiduesHaveGivenAtom = (pRes0->atomExists(atomName)) && (pRes1->atomExists(atomName)) && (pRes2->atomExists(atomName)) && (pRes3->atomExists(atomName));

    return allFourResiduesHaveGivenAtom;
}

/**
 * This little function will check that the distance between
 * the C-alpha atoms in the two residues is within a set
 * tolerance of the appropriate distance for consecutive
 * residues in a protein.
 *
 * @param pRes0 The first residue.
 * @param pRes1 The second residue.
 */
bool BBQTable::doCADistanceCheck(Residue *pRes0, Residue *pRes1) {
    double dist, diff, diffPro;

    dist = pRes0->distance(*pRes1, "CA");
    diff = fabs(dist - CA_CA_DISTANCE);
    // Note: cis-Proline has a different C-alpha to C-alpha distance.
    // Interestingly, according to th e paper by Gront et al., the
    // quadrilateral r-dims for Cis-Pro end up in their own bin most of the
    // time anyway, so we should allow both C-alpha distances.
    diffPro = fabs(dist - CA_CA_DISTANCE_PROLINE);

    if( (diff > CA_CA_TOLERANCE) && (diffPro > CA_CA_TOLERANCE) ) {
        if(debugFlag) {
            cout << "Failed CA distance check.  Dist:  " << dist << " Diff:" << diff << " ";
            cout << "Residues: " << pRes0->getChainId() << "-" << pRes0->getResidueNumber() << " " << pRes0->getResidueName() << " ; ";
            cout << pRes1->getChainId() << " " << pRes1->getResidueNumber() << "-" << pRes1->getResidueName() << "\n";
        }
        return false;
    }

    return true;
}
