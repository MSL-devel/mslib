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

#include <iostream>

#include "PDBReader.h"
#include "PDBWriter.h"
#include "Transforms.h"

using namespace std;

using namespace MSL;


int main() {





	string pdbtext = "\
ATOM      1 1H   ALA A   1      -0.079  -2.341  -3.811  1.00  0.00           H  \n\
ATOM      2 2H   ALA A   1       0.350  -0.712  -3.601  1.00  0.00           H  \n\
ATOM      3 3H   ALA A   1      -0.809  -1.422  -2.584  1.00  0.00           H  \n\
ATOM      4  N   ALA A   1       0.072  -1.586  -3.112  1.00  0.00           N  \n\
ATOM      5  CA  ALA A   1       1.121  -1.981  -2.192  1.00  0.00           C  \n\
ATOM      6  HA  ALA A   1       0.834  -2.871  -1.644  1.00  0.00           H  \n\
ATOM      7  CB  ALA A   1       2.417  -2.319  -2.965  1.00  0.00           C  \n\
ATOM      8 1HB  ALA A   1       2.218  -3.144  -3.682  1.00  0.00           H  \n\
ATOM      9 2HB  ALA A   1       2.775  -1.441  -3.546  1.00  0.00           H  \n\
ATOM     10 3HB  ALA A   1       3.228  -2.646  -2.281  1.00  0.00           H  \n\
ATOM     11  C   ALA A   1       1.370  -0.901  -1.154  1.00  0.00           C  \n\
ATOM     12  O   ALA A   1       1.549  -1.196   0.026  1.00  0.00           O  \n\
TER      13                                                                     \n\
END                                                                             \n";

	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|      Write a pdb file and read it into an atom vector     |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	
	// write the input test PDB file
	ofstream pdb_fs;
	pdb_fs.open("testPdb.pdb");
	pdb_fs << pdbtext;
	if (pdb_fs.fail()) {
		cerr << "Cannot write test input pdb file testPdb.pdb" << endl;
		exit(1);
	} else {
		cout << "Written test input pdb file testPdb.pdb" << endl;
	}

	pdb_fs.close();

	PDBReader rAv;
	//rAv.open(argv[1]);
	rAv.open("testPdb.pdb");
	rAv.read();
	AtomPointerVector av = rAv.getAtomPointers();
	cout << "Read atom vector with size " << av.size() << endl;
	rAv.close();
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		cout << *(*k) << endl;
	}

	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|                  Copy the first atom and run              |" << endl;
	cout << "|                the atom transform tests on it             |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	Transforms tr;
	tr.setStoreTransformHistory(true);
	if (av.size() == 0) {
		cerr << "Empty atom vector! Exit" << endl;
		exit(1);
	}
	Atom a(*(av[0]));

	cout << "* Initial coordinates:" << endl;
	cout << a << endl;
	cout << "-----" << endl;
	cout << endl;
	a.addAltConformation();

	// translation
	CartesianPoint trans(1.0, 1.0, 1.0);
	tr.translate(a, trans);
	cout << "* Translation by " << trans << ", length = " << trans.length() << endl;
	a.setActiveConformation(1);
	CartesianPoint c1 = a.getCoor();
	a.setActiveConformation(0);
	CartesianPoint c2 = a.getCoor();
	a.addAltConformation();

	cout << a << endl;
	cout << "Distance to previous = " << (c2-c1).length() << endl;
	cout << "-----" << endl;
	cout << endl;

	// X rotation
	double rot = 73.45;
	tr.Xrotate(a, rot);
	cout << "* X rotation by " << rot << " degrees" << endl;
	a.setActiveConformation(2);
	c1 = a.getCoor();
	a.setActiveConformation(0);
	c2 = a.getCoor();
	CartesianPoint xO(c2.getX(), 0.0, 0.0);
	a.addAltConformation();

	cout << a << endl;
	cout << "Angle between old and new position with respect to X axis = " << c2.angle(xO, c1) << endl;
	cout << "-----" << endl;
	cout << endl;

	// Y rotation
	rot = -47.10;
	tr.Yrotate(a, rot);
	cout << "* Y rotation by " << rot << " degrees" << endl;
	a.setActiveConformation(3);
	c1 = a.getCoor();
	a.setActiveConformation(0);
	c2 = a.getCoor();
	CartesianPoint yO(0.0, c2.getY(), 0.0);
	a.addAltConformation();

	cout << a << endl;
	cout << "Angle between old and new position with respect to Y axis = " << c2.angle(yO, c1) << endl;
	cout << "-----" << endl;
	cout << endl;

	// Z rotation
	rot = -164.53;
	tr.Zrotate(a, rot);
	cout << "* Z rotation by " << rot << " degrees" << endl;
	a.setActiveConformation(4);
	c1 = a.getCoor();
	a.setActiveConformation(0);
	c2 = a.getCoor();
	CartesianPoint zO(0.0, 0.0, c2.getZ());
	a.addAltConformation();

	cout << a << endl;
	cout << "Angle between old and new position with respect to Z axis = " << c2.angle(zO, c1) << endl;
	cout << "-----" << endl;
	cout << endl;

	// axial rotation, axis from origin
	rot = 63.42;
	CartesianPoint axisFromCenter(12.3, 11.56, -2.45);
	CartesianPoint center(0.0, 0.0, 0.0);
	tr.rotate(a, rot, axisFromCenter);
	cout << "* Rotation around an arbitrary axis " << axisFromCenter << " by " << rot << " degrees" << endl;
	a.setActiveConformation(5);
	c1 = a.getCoor();
	a.setActiveConformation(0);
	c2 = a.getCoor();
	CartesianPoint prO = CartesianGeometry::projection(c2, axisFromCenter, center);
	a.addAltConformation();
	
	cout << a << endl;
	cout << "Angle between old and new position with respect to axis = " << c2.angle(prO, c1) << endl;
	cout << "-----" << endl;
	cout << endl;

	// axial rotation, axis not from origin
	rot = -23.67;
	axisFromCenter.setCoor(-12.5, 11.2, 9.76);
	center.setCoor(-4.6, 7.2, 9.1);
	tr.rotate(a, rot, axisFromCenter, center);
	cout << "* Rotation around arbitrary axis " << axisFromCenter << " and center " << center << " by " << rot << " degrees" << endl;
	a.setActiveConformation(6);
	c1 = a.getCoor();
	a.setActiveConformation(0);
	c2 = a.getCoor();
	prO = CartesianGeometry::projection(c2, axisFromCenter, center);
	a.addAltConformation();
	
	cout << a << endl;
	cout << "Angle between old and new position with respect to axis = " << c2.angle(prO, c1) << endl;
	cout << "-----" << endl;
	cout << endl;

	// align, target from origin
	CartesianPoint target(-8.45, 3.24, -12.32);
	center.setCoor(0.0, 0.0, 0.0);
	tr.align(a, target, center);
	cout << "* Alignment with target vector " << target << endl;
	a.setActiveConformation(7);
	c1 = a.getCoor();
	a.setActiveConformation(0);
	c2 = a.getCoor();
	a.addAltConformation();
	
	cout << "Angle with target before aligning = " << c1.angle(center, target) << endl;
	cout << a << endl;
	cout << "Angle with target after aligning = " << c2.angle(center, target) << endl;
	cout << "-----" << endl;
	cout << endl;

	// align, target from origin
	target.setCoor(8.76, -3.21, 14.32);
	center.setCoor(6.2, -4.2, 7.6);
	tr.align(a, target, center);
	cout << "* Alignment with target vector " << target << " and center " << center << endl;
	a.setActiveConformation(8);
	c1 = a.getCoor();
	a.setActiveConformation(0);
	c2 = a.getCoor();
	a.addAltConformation();

	cout << "Angle with target before aligning = " << c1.angle(center, target) << endl;
	cout << a << endl;
	cout << "Angle with target after aligning = " << c2.angle(center, target) << endl;
	cout << "-----" << endl;
	cout << endl;

	// axially orient
	target.setCoor(22.13, -9.34, 21.5);
	CartesianPoint axis1(-8.97, -3.21, -5.77);
	CartesianPoint axis2(5.66, 8.23, -7.21);
	tr.orient(a, target, axis1, axis2);
	cout << "* Orient as target vector " << target << " with respect to axis " << axis1 << axis2 << endl;
	a.setActiveConformation(9);
	c1 = a.getCoor();
	a.setActiveConformation(0);
	c2 = a.getCoor();
	a.addAltConformation();

	cout << "Dihedral angle before orienting = " << c1.dihedral(axis1, axis2, target) << endl;
	cout << a << endl;
	cout << "Dihedral angle after orienting = " << c2.dihedral(axis1, axis2, target) << endl;
	cout << "-----" << endl;
	cout << endl;

	cout << "* Apply all the history at once on the starting position" << endl;
	a.setActiveConformation(1);
	cout << "Start orientation" <<endl;
	cout << a << endl;
	tr.applyHistory(a);
	cout << "After history" <<endl;
	cout << a << endl;

	cout << endl;


	cout << endl;
	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|              Apply the same set of transformations        |" << endl;
	cout << "|                     on the whole atom vector              |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	cout << endl;
	tr.resetHistory();

	cout << "* Initial coordinates:" << endl;
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		cout << *(*k) << endl;
		(*k)->addAltConformation();
	}
	cout << "-----" << endl;
	cout << endl;

	// translation
	trans.setCoor(1.0, 1.0, 1.0);
	tr.translate(av, trans);
	cout << "* Translation by " << trans << ", length = " << trans.length() << endl;
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		(*k)->setActiveConformation(1);
		c1 = (*k)->getCoor();
		(*k)->setActiveConformation(0);
		c2 = (*k)->getCoor();
		(*k)->addAltConformation();

		cout << *(*k) << endl;
		cout << "Distance to previous = " << (c2-c1).length() << endl;
	}
	cout << "-----" << endl;
	cout << endl;

	// X rotation
	rot = 73.45;
	tr.Xrotate(av, rot);
	cout << "* X rotation by " << rot << " degrees" << endl;
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		(*k)->setActiveConformation(2);
		c1 = (*k)->getCoor();
		(*k)->setActiveConformation(0);
		c2 = (*k)->getCoor();
		(*k)->addAltConformation();

		cout << *(*k) << endl;
		CartesianPoint xO(c2.getX(), 0.0, 0.0);
		cout << "Angle between old and new position with respect to X axis = " << c2.angle(xO, c1) << endl;
	}
	cout << "-----" << endl;
	cout << endl;

	// Y rotation
	rot = -47.10;
	tr.Yrotate(av, rot);
	cout << "* Y rotation by " << rot << " degrees" << endl;
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		(*k)->setActiveConformation(3);
		c1 = (*k)->getCoor();
		(*k)->setActiveConformation(0);
		c2 = (*k)->getCoor();
		(*k)->addAltConformation();

		cout << *(*k) << endl;
		CartesianPoint yO(0.0, c2.getY(), 0.0);
		cout << "Angle between old and new position with respect to Y axis = " << c2.angle(yO, c1) << endl;
	}
	cout << "-----" << endl;
	cout << endl;

	// Z rotation
	rot = -164.53;
	tr.Zrotate(av, rot);
	cout << "* Z rotation by " << rot << " degrees" << endl;
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		(*k)->setActiveConformation(4);
		c1 = (*k)->getCoor();
		(*k)->setActiveConformation(0);
		c2 = (*k)->getCoor();
		(*k)->addAltConformation();

		cout << *(*k) << endl;
		CartesianPoint zO(0.0, 0.0, c2.getZ());
		cout << "Angle between old and new position with respect to Z axis = " << c2.angle(zO, c1) << endl;
	}
	cout << "-----" << endl;
	cout << endl;

	// axial rotation, axis from origin
	rot = 63.42;
	axisFromCenter.setCoor(12.3, 11.56, -2.45);
	center.setCoor(0.0, 0.0, 0.0);
	tr.rotate(av, rot, axisFromCenter);
	cout << "* Rotation around an arbitrary axis " << axisFromCenter << " by " << rot << " degrees" << endl;
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		(*k)->setActiveConformation(5);
		c1 = (*k)->getCoor();
		(*k)->setActiveConformation(0);
		c2 = (*k)->getCoor();
		(*k)->addAltConformation();

		cout << *(*k) << endl;
		CartesianPoint prO = CartesianGeometry::projection(c2, axisFromCenter, center);
		cout << "Angle between old and new position with respect to axis = " << c2.angle(prO, c1) << endl;
	}
	cout << "-----" << endl;
	cout << endl;

	// axial rotation, axis not from origin
	rot = -23.67;
	axisFromCenter.setCoor(-12.5, 11.2, 9.76);
	center.setCoor(-4.6, 7.2, 9.1);
	tr.rotate(av, rot, axisFromCenter, center);
	cout << "* Rotation around arbitrary axis " << axisFromCenter << " and center " << center << " by " << rot << " degrees" << endl;
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		(*k)->setActiveConformation(6);
		c1 = (*k)->getCoor();
		(*k)->setActiveConformation(0);
		c2 = (*k)->getCoor();
		(*k)->addAltConformation();

		cout << *(*k) << endl;
		prO = CartesianGeometry::projection(c2, axisFromCenter, center);
		cout << "Angle between old and new position with respect to axis = " << c2.angle(prO, c1) << endl;
	}
	cout << "-----" << endl;
	cout << endl;

	// align, target from origin
	target.setCoor(-8.45, 3.24, -12.32);
	center.setCoor(0.0, 0.0, 0.0);
	tr.align(av, av[0]->getCoor(), target, center);
	cout << "* Alignment of the first atom with target vector " << target << endl;
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		(*k)->setActiveConformation(7);
		c1 = (*k)->getCoor();
		(*k)->setActiveConformation(0);
		c2 = (*k)->getCoor();
		(*k)->addAltConformation();

		cout << "Angle with target before aligning = " << c1.angle(center, target) << endl;
		cout << *(*k) << endl;
		cout << "Angle with target after aligning = " << c2.angle(center, target) << endl;
	}
	cout << "-----" << endl;
	cout << endl;

	// align, target from origin
	target.setCoor(8.76, -3.21, 14.32);
	center.setCoor(6.2, -4.2, 7.6);
	tr.align(av, av[0]->getCoor(), target, center);
	cout << "* Alignment of the first atom with target vector " << target << " and center " << center << endl;
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		(*k)->setActiveConformation(8);
		c1 = (*k)->getCoor();
		(*k)->setActiveConformation(0);
		c2 = (*k)->getCoor();
		(*k)->addAltConformation();

		cout << "Angle with target before aligning = " << c1.angle(center, target) << endl;
		cout << *(*k) << endl;
		cout << "Angle with target after aligning = " << c2.angle(center, target) << endl;
	}
	cout << "-----" << endl;
	cout << endl;

	// axially orient
	target.setCoor(22.13, -9.34, 21.5);
	axis1.setCoor(-8.97, -3.21, -5.77);
	axis2.setCoor(5.66, 8.23, -7.21);
	tr.orient(av, av[0]->getCoor(), target, axis1, axis2);
	cout << "* Orient the first atom as target vector " << target << " with respect to axis " << axis1 << axis2 << endl;
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		(*k)->setActiveConformation(9);
		c1 = (*k)->getCoor();
		(*k)->setActiveConformation(0);
		c2 = (*k)->getCoor();
		(*k)->addAltConformation();

		cout << "Dihedral angle before orienting = " << c1.dihedral(axis1, axis2, target) << endl;
		cout << *(*k) << endl;
		cout << "Dihedral angle after orienting = " << c2.dihedral(axis1, axis2, target) << endl;
	}

	cout << "-----" << endl;
	cout << endl;

	cout << "* Apply all the history at once on the starting position" << endl;
	cout << "Start orientation" <<endl;
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		(*k)->setActiveConformation(1);
		cout << *(*k) << endl;
	}
	tr.applyHistory(av);
	cout << "After history" <<endl;
	for (AtomPointerVector::iterator k = av.begin(); k != av.end() ; k++){
		cout << *(*k) << endl;
	}

	cout << endl;

	return 0;



}
