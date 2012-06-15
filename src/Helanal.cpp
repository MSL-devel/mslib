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

#include "Helanal.h"

using namespace MSL;
using namespace std;


Helanal::Helanal() {
	pCA1 = NULL;
        pCA2 = NULL;
        pCA3 = NULL;
        pCA4 = NULL;
	setup();
}

Helanal::Helanal(CartesianPoint & _p1, CartesianPoint & _p2, CartesianPoint & _p3, CartesianPoint & _p4) {
	update(_p1, _p2, _p3, _p4);
}

Helanal::Helanal(const Helanal& _helanal) {
	axis = _helanal.axis;		
	center = _helanal.center;
	Npoint = _helanal.Npoint;
	Cpoint = _helanal.Cpoint;
	twist = _helanal.twist;
	height = _helanal.height;
	radius = _helanal.radius;
	resPerTurn = _helanal.resPerTurn;

	errorFlag = _helanal.errorFlag;
	pCA1 = _helanal.pCA1;
	pCA2 = _helanal.pCA2;
	pCA3 = _helanal.pCA3;
	pCA4 = _helanal.pCA4;
}

Helanal::~Helanal() {
}

void Helanal::setup() {
	axis.setCoor(0.0,0.0,0.0);		
	center.setCoor(0.0,0.0,0.0);		
	Npoint.setCoor(0.0,0.0,0.0);		
	Cpoint.setCoor(0.0,0.0,0.0);		
	CA1projection.setCoor(0.0,0.0,0.0);
	CA2projection.setCoor(0.0,0.0,0.0);
	CA3projection.setCoor(0.0,0.0,0.0);
	CA4projection.setCoor(0.0,0.0,0.0);
	twist = 0.0;
	height = 0.0;
	radius = 0.0;
	resPerTurn = 0.0;
	errorFlag = true;
}

void Helanal::update(CartesianPoint & _p1, CartesianPoint & _p2, CartesianPoint & _p3, CartesianPoint & _p4) {
	setup();
	pCA1 = &_p1;
	pCA2 = &_p2;
	pCA3 = &_p3;
	pCA4 = &_p4;
	update();
}

CartesianPoint & Helanal::getCA1() {
	return *pCA1;
}

CartesianPoint & Helanal::getCA2() {
	return *pCA2;
}

CartesianPoint & Helanal::getCA3() {
	return *pCA3;
}

CartesianPoint & Helanal::getCA4() {
	return *pCA4;
}

CartesianPoint & Helanal::getAxis() {
	return axis;
}

CartesianPoint & Helanal::getCenter() {
	return center;
}

double Helanal::getTwist() const {
	return twist;
}

double Helanal::getHeight() const {
	return height;
}

double Helanal::getRadius() const {
	return radius;
}

double Helanal::getResPerTurn() const {
	return resPerTurn;
}

CartesianPoint & Helanal::getNpoint() {
	return Npoint;
}

CartesianPoint & Helanal::getCpoint() {
	return Cpoint;
}

CartesianPoint & Helanal::getCA1projection() {
	return CA1projection;
}

CartesianPoint & Helanal::getCA2projection() {
	return CA2projection;
}

CartesianPoint & Helanal::getCA3projection() {
	return CA3projection;
}

CartesianPoint & Helanal::getCA4projection() {
	return CA4projection;
}

bool Helanal::fail() const {
	return errorFlag;
}

void Helanal::update() {
	//initialize();

/*
	if (pCA1 == NULL || pCA2 == NULL || pCA3 == NULL || pCA4 == NULL) {
		cerr << "WARNING 9410: cannot calculate helanal on initialized object in void Helanal::update()" << endl;
		initialize();
		return;
	}	
*/

	// calculate some difference vectors
	CartesianPoint AB = *pCA2 - *pCA1;
	CartesianPoint BC = *pCA3 - *pCA2;
	CartesianPoint CD = *pCA4 - *pCA3;

	CartesianPoint ABC = AB - BC;
	CartesianPoint BCD = BC - CD;

	if (ABC.length() == 0.0 || BCD.length() == 0.0) {
		errorFlag = true;
		return;
	}

	// the axis is the normalized cross product of ABC x BCD
	axis.setCoor((ABC.cross(BCD)).getUnit());

	// the twist (deg/residue) is the angle between ABC and BCD
	twist = ABC.angle(BCD);
	if (twist == 0.0) {
		errorFlag = true;
		return;
	}

	// residues per turn
	resPerTurn = 360 / twist;

	// dot product of ABC * BCD divided the product of their lengths
	double cosTheta = (ABC * BCD)/(ABC.length() * BCD.length());
	double cosTheta1 =  1 - cosTheta;

	if (cosTheta1 == 0.0) {
		errorFlag = true;
		return;
	}
	// radius of the helix
	radius = (sqrt(ABC.length() * BCD.length())) / (2 * cosTheta1);
	// the height is the helical pitch, the dot product of BC * axis
	height = BC * axis; 
	
	CartesianPoint origin2 = *pCA2 - (ABC.getUnit() * radius);
	CartesianPoint origin3 = *pCA3 - (BCD.getUnit() * radius);
	// center of the vector, Npoint is half the height from the center toward the N-terminus
	// C-point is the same toward the C-terminus
	//center = (origin2 + origin3) / 2;
	//Npoint = center - (axis * height / 2);
	//Cpoint = center + (axis * height / 2);
	center.setCoor((origin2 + origin3) / 2);
	Npoint.setCoor(center - (axis * height / 2));
	Cpoint.setCoor(center + (axis * height / 2));

	// the CA projections are the projections of the four CA to the helical axis
	CA1projection.setCoor(CartesianGeometry::projection(*pCA1, center, Npoint));
	CA2projection.setCoor(CartesianGeometry::projection(*pCA2, center, Npoint));
	CA3projection.setCoor(CartesianGeometry::projection(*pCA3, center, Npoint));
	CA4projection.setCoor(CartesianGeometry::projection(*pCA4, center, Npoint));

	errorFlag = false;
}
