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



#include "VectorPair.h"
#include "CartesianGeometry.h"

using namespace MSL;
using namespace MSL::CartesianGeometry;

#include "MslOut.h"
static MslOut MSLOUT("VectorPair");

VectorPair::VectorPair(){
	vectorAid = "";
	vectorBid = "";
	a1.setCoor(0,0,0);
	a2.setCoor(0,0,0);

	b1.setCoor(0,0,0);
	b2.setCoor(0,0,0);

	distance1 = 0.0;
	distance2 = 0.0;
	angle1 = 0.0;
	angle2 = 0.0;
	angle3 = 0.0;
	angle4 = 0.0;
	torsion1 = 0.0;
	torsion2 = 0.0;

	archiveType = "binary";

}

VectorPair::VectorPair(CartesianPoint &_a1,CartesianPoint &_a2,CartesianPoint &_b1,CartesianPoint &_b2,string _vectorAid, string _vectorBid) {


	a1 = _a1;
	a2 = _a2;

	b1 = _b1;
	b2 = _b2;

	vectorAid = _vectorAid;
	vectorBid = _vectorBid;

	distance1 = 0.0;
	distance2 = 0.0;
	angle1 = 0.0;
	angle2 = 0.0;
	angle3 = 0.0;
	angle4 = 0.0;
	torsion1 = 0.0;
	torsion2 = 0.0;


}

VectorPair::VectorPair(const VectorPair &_vp){
	copy(_vp);
}

VectorPair& VectorPair::operator=(const VectorPair &_vp){
	copy(_vp);
	return *this;
}

void VectorPair::copy(const VectorPair &_vp){
	a1  = _vp.a1;
	a2  = _vp.a2;

	b1  = _vp.b1;
	b2  = _vp.b2;

	vectorAid = _vp.vectorAid;
	vectorBid = _vp.vectorBid;

	calcAll();
}

VectorPair::~VectorPair(){
}

void VectorPair::calcAll(){
	calcDistance1();
	calcDistance2();
	calcAngle1();
	calcAngle2();
	calcAngle3();
	calcAngle4();
	calcTorsion1();
	calcTorsion2();
}
double VectorPair::getDistance1(){	 
	return distance1;
}
double VectorPair::calcDistance1(){
	distance1 = CartesianGeometry::distance(a1,b1);
	return distance1;
}
double VectorPair::getDistance2(){	 
	return distance2;
}

double VectorPair::calcDistance2(){
	distance2 = CartesianGeometry::distance(a2,b2);
	return distance2;
}

double VectorPair::getAngle1(){
	return angle1;
}
double VectorPair::calcAngle1(){
	angle1 = angle(a2,a1,b1);
	return angle1;
}

double VectorPair::getAngle2(){
	return angle2;
}
double VectorPair::calcAngle2(){
	angle2 = angle(b2,b1,a1);
	return angle2;
}
double VectorPair::getAngle3(){
	return angle3;
}
double VectorPair::calcAngle3(){
	angle3 = angle(a1,a2,b2);
	return angle3;
}

double VectorPair::getAngle4(){
	return angle4;
}
double VectorPair::calcAngle4(){
        angle4 = angle(b1,b2,a2);
	return angle4;
}

double VectorPair::getTorsion1(){
	return torsion1;
}

double VectorPair::calcTorsion1(){
	torsion1 = dihedral(a2,a1,b1,b2);
	return torsion1;
}

double VectorPair::getTorsion2(){
	return torsion2;
}

double VectorPair::calcTorsion2(){
        torsion2 = dihedral(a1,a2,b2,b1);
	return torsion2;
}

bool VectorPair::operator< ( const VectorPair &rhs ) const {
	return (distance1 < rhs.distance1);
}


string VectorPair::getVectorAid(){
	return vectorAid;
}

string VectorPair::getVectorBid(){
	return vectorBid;
}


string VectorPair::toString() const{

	char tmp[100];
	sprintf(tmp,"%-10s %-10s %8.3f %8.3f %8.3f %8.3f",vectorAid.c_str(),vectorBid.c_str(),distance1,angle1,angle2,torsion1);

	return (string)tmp;
}

string VectorPair::getVectorPairId(){
	stringstream ss;
	ss << vectorAid<<"-"<<vectorBid;
	return (ss.str());
}

string VectorPair::getVectorPairIdReverse(){
	stringstream ss;
	ss << vectorBid<<"-"<<vectorAid;
	return (ss.str());
}
