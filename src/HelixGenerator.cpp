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

#include "HelixGenerator.h"
#include "Transforms.h"

using namespace MSL;
using namespace std;

//#define PI 3.14159265358979323846264338327950288

/**
 * This is the default constructor.  By default,
 * it will set up the Helix Generator to use
 * parameters that are appropriate for an Alpha-Helix.
 */
HelixGenerator::HelixGenerator() {
    // Default parameters are for an alpha helix.
    setHelixParameters(3.78, 91.6 / 180.0 * M_PI, 51.2 / 180.0 * M_PI);
}

/**
 * This constructor allows the user to set the helix
 * parameters at initialization.  They can also be
 * changed later using setHelixParameters.
 */
HelixGenerator::HelixGenerator(double _cAlphaDistance, double _cAlphaAngle, double _cAlphaDihedral) {
    setHelixParameters(_cAlphaDistance,_cAlphaAngle, _cAlphaDihedral);
}

/**
 * This constructor allows the user to set the BBQ filename
 * at initialization.  This can also be done later using
 * the setBBQTableFileName method.  The parameters
 * will be the default alpha helix parameters.
 */
HelixGenerator::HelixGenerator(string _bbqTableFilename) {
    // Default parameters are for an alpha helix.
    setHelixParameters(3.78, 91.6 / 180.0 * M_PI, 51.2 / 180.0 * M_PI);
    setBBQTableFileName(_bbqTableFilename);
}

/**
 * This constructor will allow the user to set the helix parameters and
 * BBQ table filename all at once.
 */
HelixGenerator::HelixGenerator(double _cAlphaDistance, double _cAlphaAngle, double _cAlphaDihedral, string _bbqTableFilename) {
    setHelixParameters(_cAlphaDistance, _cAlphaAngle, _cAlphaDihedral);
    setBBQTableFileName(_bbqTableFilename);
}

HelixGenerator::~HelixGenerator() {
    
}

void HelixGenerator::setHelixParameters(double _cAlphaDistance, double _cAlphaAngle, double _cAlphaDihedral) {
    cAlphaDistance = _cAlphaDistance;
    cAlphaAngle = _cAlphaAngle;
    cAlphaDihedral = _cAlphaDihedral;

    calcHelixParameters();
}

void HelixGenerator::generateHelix(AtomPointerVector &_av, int numCAlphas, bool fillInMissingBackboneAtoms, bool center) {
    CartesianPoint currCAlphaLocation;
    int start, end;
    fillInMissingBackboneAtoms = fillInMissingBackboneAtoms && (bbq.size() > 0);

    // Generate an extra C-Alpha at each end so that BBQ can work its magic.
    // We'll take off these atoms in the end.
    if(fillInMissingBackboneAtoms) {
        start = -2;
        end = numCAlphas + 2;
    }
    else {
        start = 0;
        end = numCAlphas;
    }
    
    for(int currCAlpha = start; currCAlpha < end; ++currCAlpha) {
        currCAlphaLocation.setX(radius * cos(period * currCAlpha));
        currCAlphaLocation.setY(radius * sin(period * currCAlpha));
        currCAlphaLocation.setZ(rise * currCAlpha);
        Atom *newAtom = new Atom("CA", currCAlphaLocation);
        newAtom->setResidueNumber(currCAlpha);
        newAtom->setElement("C");
        newAtom->setChainId("A");

        _av.push_back(newAtom);
    }

    if(fillInMissingBackboneAtoms) {
        generateMissingBackboneAtoms(_av, numCAlphas);
    }
    
    // Center at origin.
    if(center){
        Atom *firstAtom = _av.front();
        Atom *lastAtom = _av.back();
        double midZ = (firstAtom->getCoor().getZ() + lastAtom->getCoor().getZ()) / 2.0;
        Transforms t;
        t.translate(_av, CartesianPoint(0.0, 0.0, -midZ));
    }
}

void HelixGenerator::generateMissingBackboneAtoms(AtomPointerVector &_av, int numCAlphas) {
    // Copy the atom vector into a chain.
    Chain tempChain(_av, "A");
    // Delete the old atoms.
    _av.deletePointers();
    // Now fill in the missing atoms from the chain.
    bbq.fillInMissingBBAtoms(tempChain);

    // Now fill in CB atoms, and move atoms from the chain over to the atom vector.
    for(int i=0; i < tempChain.positionSize(); i++) {
        Residue &tempRes = tempChain.getResidue(i);

        if( tempRes.atomExists("CA") && tempRes.atomExists("N") && tempRes.atomExists("C") && tempRes.atomExists("O") ) {
            Atom *tempAtom = new Atom();
            tempAtom->setName("CB");
            tempAtom->setElement("C");
            tempAtom->setChainId("A");
            tempAtom->setResidueNumber(tempRes.getResidueNumber());
            tempAtom->setResidueIcode(tempRes.getResidueIcode());
            tempAtom->setCoor(CartesianGeometry::build( tempRes.getAtom("CA").getCoor(), tempRes.getAtom("N").getCoor(), tempRes.getAtom("C").getCoor(), 1.521, 110.5, -122.5));

            _av.push_back(new Atom(tempRes.getAtom("N")));
            _av.push_back(new Atom(tempRes.getAtom("CA")));
            _av.push_back(tempAtom);
            _av.push_back(new Atom(tempRes.getAtom("C")));
            _av.push_back(new Atom(tempRes.getAtom("O")));
        }
    }
}

void HelixGenerator::calcHelixParameters() {
    /****************************************************
     * Compute the radius based on internal coordinates *
     ****************************************************/
    double cAlphaDistanceSqrd = cAlphaDistance * cAlphaDistance;
    double cosCAlphaAngle = cos(cAlphaAngle);
    double cosCAlphaDihedral = cos(cAlphaDihedral);
    double numerator = sqrt(2.0 * cAlphaDistanceSqrd * (1.0 + cosCAlphaAngle));
    double denominator = 3.0 + cosCAlphaAngle - cosCAlphaDihedral * (1.0 - cosCAlphaAngle);

    radius = numerator / denominator;

    /*********************************************************
     * Compute the periodicity based on internal coordinates *
     *********************************************************/
    period = 0.5 * (-cosCAlphaAngle + cosCAlphaDihedral * (1.0 - cosCAlphaAngle) - 1.0);
    period = acos(period);

    /*****************************
     * Compute the z-translation *
     *****************************/
    numerator = cAlphaDistanceSqrd * (1.0 - cosCAlphaAngle) * (1.0 - cosCAlphaDihedral);
    rise = sqrt(numerator / denominator);
}

#ifdef __GSL__
Line getHelixAxis(AtomPointerVector &_av) {
    AtomPointerVector cAlphas, idealHelix;
    HelixGenerator hg;
    Transforms trans;
    CartesianPoint helixAxis(0, 0, 1);
    Line returnVal;

    // We'll only align the C-alpha atoms.
    for(int i = 0; i < _av.size(); ++i) {
        if(_av[i]->getName() == "CA")
            cAlphas.push_back(new Atom(*(_av[i])));
    }

    // Generate an ideal helix of the same size,
    // and align it to the given helix.
    hg.generateHelix(idealHelix, cAlphas.size(), false);
    trans.rmsdAlignment(cAlphas, idealHelix);

    // Set the center and direction.
    returnVal.setCenter( trans.getLastTranslation() );
    returnVal.setDirection( helixAxis * trans.getLastRotationMatrix() );

    // Delete the pointers.
    cAlphas.deletePointers();
    idealHelix.deletePointers();

    return returnVal;
}
#endif


