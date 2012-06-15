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
#ifndef _HELIX_GENERATOR_H
#define _HELIX_GENERATOR_H

#include "AtomPointerVector.h"
#include "BBQTable.h"
#include "Line.h"

/**
 * This class will be used to generate atom vectors of helices.
 */
namespace MSL { 
class HelixGenerator {
public:
    HelixGenerator();
    HelixGenerator(double _cAlphaDistance, double _cAlphaAngle, double _cAlphaDihedral);
    HelixGenerator(std::string _bbqTableFilename);
    HelixGenerator(double _cAlphaDistance, double _cAlphaAngle, double _cAlphaDihedral, std::string _bbqTableFilename);
    ~HelixGenerator();

    inline AtomPointerVector *generateHelix(int numCAlphas, bool fillInMissingBackboneAtoms = true, bool center = true) {
        AtomPointerVector *av = new AtomPointerVector();
        generateHelix(*av, numCAlphas, fillInMissingBackboneAtoms, center);
        return av;
    }

    void setHelixParameters(double _cAlphaDistance, double _cAlphaAngle, double _cAlphaDihedral);
    void generateHelix(AtomPointerVector &_av, int numCAlphas, bool fillInMissingBackboneAtoms = true, bool center = true);
    void setCAlphaParameters(double _distance, double _angle, double _dihedral) {
        cAlphaDistance = _distance;
        cAlphaAngle = _angle;
        cAlphaDihedral = _dihedral;
        calcHelixParameters();
    }
    void setBBQTableFileName(std::string _fileName) {
        bbq.openReader(_fileName);
    };

private:
    void calcHelixParameters();
    void generateMissingBackboneAtoms(AtomPointerVector &_av, int numCAlphas);

    BBQTable bbq;
    double cAlphaDistance;
    double cAlphaAngle;
    double cAlphaDihedral;
    double radius;
    double period;
    double rise;
};

#ifdef __GSL__
// Global function??? Is it ever used??? REMOVE?
Line getHelixAxis(AtomPointerVector &_av);
#endif


}

#endif // _HELIX_GENERATOR_H
