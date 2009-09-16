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

#include "Frame.h"
#include "PDBFormat.h"

Frame::Frame() {
    name = "";
}

Frame::Frame(const Frame &_frame) {
    copy(_frame);
}

Frame::~Frame() {
}

void Frame::operator=(const Frame & _frame) {
    copy(_frame);
}

map<string, Line> Frame::getLines() const {
    return lines;
}

PrincipleComponentAnalysis Frame::getPCA() const {
    return pca;
}

void Frame::computeFrameFromPCA(AtomVector &_atoms) {
    pca.computePrincipleComponents(_atoms);
    vector<Line> pcaLines = pca.getLines();
    lines["Z"] = pcaLines[0];
    lines["X"] = pcaLines[1];
    lines["Y"] = pcaLines[2];

    _atoms.updateGeometricCenter();
    center = _atoms.getGeometricCenter();
}

void Frame::computeFrameFrom3Atoms(Atom &_at1, Atom &_at2, Atom &_at3) {
    /*
      Define Frame by 3 lines:
      1. _at2 -> _at1    Line[0]
      2. _at2 -> _at3    Line[1]
      3. (1) cross (2)   Line[2]
     */

    //cout << "Ats: "<<endl<<_at1<<endl<<_at2<<endl<<_at3<<endl;

    CartesianPoint dir1 = _at2.getCoor() - _at1.getCoor();
    CartesianPoint dir2 = _at2.getCoor() - _at3.getCoor();

    dir1 = dir1.getUnit();
    dir2 = dir2.getUnit();

    CartesianPoint dir3 = dir1.cross(dir2);
    dir3 = dir3.getUnit(); // Needed?

    // dir1,2 and 3 are not all normal. dir3 is normal to dir1,2.
    CartesianPoint dir4 = dir3.cross(dir2);
    dir4 = dir4.getUnit();

    Line l1;
    l1.setCenter(_at2.getCoor());
    l1.setDirection(dir4);

    Line l2;
    l2.setCenter(_at2.getCoor());
    l2.setDirection(dir2);

    Line l3;
    l3.setCenter(_at2.getCoor());
    l3.setDirection(dir3);

    lines["X"] = l1;
    lines["Y"] = l2;
    lines["Z"] = l3;

    center = _at2.getCoor();
}

/**
 * This function will set the frame given 3 points on the coordinate axes.
 * It is assumed that these points are orthogonal to each other through the
 * origin.  It does not check this though.  They do not have to be unit vectors.
 * Basically we are defining our frame by 3 lines:
 *     1. (0,0,0) -> _point1    Line[0]
 *     2. (0,0,0) -> _point2    Line[1]
 *     3. (0,0,0) -> _point3    Line[2]
 */
void Frame::computeFrameFromAxes(CoordAxes &_axes) {
    /*
      Define Frame by 3 lines:
      1. (0,0,0) -> _point1    Line[0]
      2. (0,0,0) -> _point2    Line[1]
      3. (0,0,0) -> _point3    Line[2]
     */

    //cout << "Pts: "<<endl<<_point1<<endl<<_point2<<endl<<_point3<<endl;
    CartesianPoint origin((Real) 0.0f, (Real) 0.0f, (Real) 0.0f);

    Line l1;
    l1.setCenter(origin);
    l1.setDirection(_axes.first.getUnit());

    Line l2;
    l2.setCenter(origin);
    l2.setDirection(_axes.second.getUnit());

    Line l3;
    l3.setCenter(origin);
    l3.setDirection(_axes.third.getUnit());

    lines["X"] = l1;
    lines["Y"] = l2;
    lines["Z"] = l3;

    center = origin;
}

Matrix Frame::getBasisTransformMatrix(Frame &_fromFrame, Frame &_toFrame) {
    Matrix result(3, 3, 0);

    CartesianPoint i, j, k, ip, jp, kp;
    i.setCoor(_fromFrame.lines["X"].getDirection());
    j.setCoor(_fromFrame.lines["Y"].getDirection());
    k.setCoor(_fromFrame.lines["Z"].getDirection());

    ip.setCoor(_toFrame.lines["X"].getDirection());
    jp.setCoor(_toFrame.lines["Y"].getDirection());
    kp.setCoor(_toFrame.lines["Z"].getDirection());

    //  Transformation matrix going from this frame to input frame.

    // Dot products
    result[0][0] = i * ip;
    result[0][1] = i * jp;
    result[0][2] = i * kp;

    result[1][0] = j * ip;
    result[1][1] = j * jp;
    result[1][2] = j * kp;

    result[2][0] = k * ip;
    result[2][1] = k * jp;
    result[2][2] = k * kp;

    return result;
}

Matrix Frame::getBasisTransformMatrix(Frame &_frame) {
    return( getBasisTransformMatrix(*this, _frame) );
}

void Frame::copy(const Frame & _frame) {
    name = _frame.getName();
    map<string, Line> in;
    in = _frame.getLines();
    for (map<string, Line>::iterator currIter = in.begin(); currIter != in.end(); ++currIter) {
        string key = currIter->first;
        lines[key] = in[key];
    }

    pca = _frame.getPCA();
    center = _frame.getCenter();
}

string Frame::toString() {
    stringstream ss;
    for (map<string, Line>::iterator currIter = lines.begin(); currIter != lines.end(); ++currIter) {
        string key = currIter->first;
        lines[key].setOutputFormat("pymolAtom");

        if(key == "X")
            lines[key].setColor(1.0, 0.0, 0.0, 0.5, 0.0, 0.0);
        else if(key == "Y")
            lines[key].setColor(0.0, 1.0, 0.0, 0.0, 0.5, 0.0);
        else if(key == "Z")
            lines[key].setColor(0.0, 0.0, 1.0, 0.0, 0.0, 0.5);

        stringstream n;
        n << name << "Axis" << key;
        lines[key].setName(n.str().c_str());
        ss << lines[key].toString() << endl;
    }

    Atom a("DUM", center, "H");
    a.setResidueName("DUM");

    string pdbline = PDBFormat::createAtomLine(PDBFormat::createAtomData(a));

    string pdbName = name;
    if (pdbName == "") {
        pdbName = "Center";
    }
    ss << "cmd.read_pdbstr(\"" << pdbline << "\",\"" << pdbName << "\")" << endl;

    return ss.str();
}

void Frame::setName(string _name) {
    name = _name;
}

double Frame::distanceToFrame(Frame &_frame) {
    return center.distance(_frame.getCenter());
}

Matrix Frame::anglesBetweenFrame(Frame &_frame){
	
	Matrix m = getBasisTransformMatrix(_frame);	
	Matrix angles(3,3);
	
	for (uint ii = 0; ii < angles.getRows();ii++){
		for (uint jj = 0; jj < angles.getCols();jj++){
			angles[ii][jj] = acos(m[ii][jj]) * 180/M_PI;
		}
	}
	
	return angles;
}
void Frame::translate(const CartesianPoint & _vec) {
    for (map<string, Line>::iterator currIter = lines.begin(); currIter != lines.end(); ++currIter) {
        string key = currIter->first;
        lines[key].translate(_vec);
    }
    center += _vec;
}

Line& Frame::operator[](string _n) {
    if (lines.find(_n) != lines.end()) {
        return lines[_n];
    }

    cerr << "ERROR Frame::operator[n], element " << _n << " not found in map." << endl;
    exit(1);

}


void Frame::transformAtoms(AtomVector &_atoms, Frame &_fromFrame, Frame &_toFrame) {
    // Create Transformation Matrix
    Matrix result = getBasisTransformMatrix(_fromFrame, _toFrame);


    // Transform atoms.
    for (uint i = 0; i < _atoms.size(); i++) {
        _atoms(i).setCoor(_atoms(i).getCoor() - _fromFrame.center);
	_atoms(i).setCoor(_atoms(i).getCoor() * result);
    }
}

void Frame::transformToFromGlobalBasis(AtomVector &_atoms, bool bToGlobal) {
    Frame globalFrame;
    CartesianPoint origin((Real)0.0f, (Real)0.0f, (Real)0.0f);
    CartesianPoint xAxis((Real)1.0f, (Real)0.0f, (Real)0.0f);
    CartesianPoint yAxis((Real)0.0f, (Real)1.0f, (Real)0.0f);
    CartesianPoint zAxis((Real)0.0f, (Real)0.0f, (Real)1.0f);


    globalFrame.lines["X"] = Line( origin, xAxis);
    globalFrame.lines["Y"] = Line( origin, yAxis);
    globalFrame.lines["Z"] = Line( origin, zAxis);

    if(bToGlobal) {

	    // Create Transformation Matrix
	    Matrix result = getBasisTransformMatrix(*this,globalFrame);


	    // Transform atoms.
	    for (uint i = 0; i < _atoms.size(); i++) {
		    _atoms(i).setCoor(_atoms(i).getCoor() - center);
		    _atoms(i).setCoor(_atoms(i).getCoor() * result);
	    }
    //        transformAtoms(_atoms, *this, globalFrame);
    } else {
	    

	    // Create Transformation Matrix
	    Matrix result = getBasisTransformMatrix(globalFrame,*this);

	    // Transform atoms.
	    for (uint i = 0; i < _atoms.size(); i++) {
		    _atoms(i).setCoor(_atoms(i).getCoor() * result);
		    _atoms(i).setCoor(_atoms(i).getCoor() + center);
	    }

	    //transformAtoms(_atoms, globalFrame, *this);
    }
}

void Frame::transformToGlobalBasis(AtomVector &_atoms) {
    transformToFromGlobalBasis(_atoms, true);
}

void Frame::transformFromGlobalBasis(AtomVector &_atoms) {
    transformToFromGlobalBasis(_atoms, false);
}


bool Frame::computeFrameFromFunctionalGroup(Residue &_res){
	

	bool result = false;
	if (_res.getResidueName() == "ARG" &&
	    _res.exists("NE") &&
	    _res.exists("CZ") &&
	    _res.exists("NH1")){

		computeFrameFrom3Atoms(_res("NE"),_res("CZ"),_res("NH1"));

		result = true;
		
	}


	if (_res.getResidueName() == "LYS" &&
	    _res.exists("CD") &&
	    _res.exists("CE") &&
	    _res.exists("NZ")){

		computeFrameFrom3Atoms(_res("CD"),_res("CE"),_res("NZ"));
		translate(_res("NZ").getCoor() - _res("CE").getCoor());

		result = true;

	}
	if (_res.getResidueName() == "MLZ" &&
	    _res.exists("CD") &&
	    _res.exists("CE") &&
	    _res.exists("NZ")){

		computeFrameFrom3Atoms(_res("CD"),_res("CE"),_res("NZ"));

		result = true;

	}
	if (_res.getResidueName() == "M3L" &&
	    _res.exists("CD") &&
	    _res.exists("CE") &&
	    _res.exists("NZ")){

		computeFrameFrom3Atoms(_res("CD"),_res("CE"),_res("NZ"));

		result = true;

	}

	

	return result;
}

// Closest atom last
bool Frame::computeFrameFrom3AtomNames(Residue &_res, string & atom1, string & atom2, string & atom3){
	
	bool result = false;
	if (_res.exists(atom1) &&
	    _res.exists(atom2) &&
	    _res.exists(atom3)){

		computeFrameFrom3Atoms(_res(atom1),_res(atom2),_res(atom3));

		result = true;
		
	}

	return result;
}
