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
#ifndef _HELIX_GENERATOR_H
#define _HELIX_GENERATOR_H

#include "AtomVector.h"
#include "BBQTable.h"

/**
 * This class will be used to generate atom vectors of helices.
 */
class HelixGenerator {
public:
    HelixGenerator();
    HelixGenerator(double _cAlphaDistance, double _cAlphaAngle, double _cAlphaDihedral);
    HelixGenerator(string _bbqTableFilename);
    HelixGenerator(double _cAlphaDistance, double _cAlphaAngle, double _cAlphaDihedral, string _bbqTableFilename);
    ~HelixGenerator();

    inline AtomVector *generateHelix(int numCAlphas, bool fillInMissingBackboneAtoms = true) {
        AtomVector *av = new AtomVector();
        generateHelix(*av, numCAlphas, fillInMissingBackboneAtoms);
        return av;
    }

    void setHelixParameters(double _cAlphaDistance, double _cAlphaAngle, double _cAlphaDihedral);
    void generateHelix(AtomVector &_av, int numCAlphas, bool fillInMissingBackboneAtoms = true);
    void setCAlphaParameters(double _distance, double _angle, double _dihedral) {
        cAlphaDistance = _distance;
        cAlphaAngle = _angle;
        cAlphaDihedral = _dihedral;
        calcHelixParameters();
    }
    void setBBQTableFileName(string _fileName) {
        bbq.openReader(_fileName);
    };

private:
    void calcHelixParameters();
    void generateMissingBackboneAtoms(AtomVector &_av, int numCAlphas);

    BBQTable bbq;
    double cAlphaDistance;
    double cAlphaAngle;
    double cAlphaDihedral;
    double radius;
    double period;
    double rise;
};



#endif // _HELIX_GENERATOR_H
