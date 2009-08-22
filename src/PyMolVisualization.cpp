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

#include "PyMolVisualization.h"
#include "PDBFormat.h"


PyMolVisualization::PyMolVisualization(){
	pymolObjectStrings.clear();
	pymolAtomStrings.clear();
}


PyMolVisualization::~PyMolVisualization(){}


bool PyMolVisualization::createSphere(CartesianPoint &_cp, string _name, double _radius,int rgbR,int rgbG,int rgbB){

	stringstream ss;

	// Change to get a random name not in our map later...
	if (_name == ""){
		_name = "foo";
	}

	ss << _name << " = [ COLOR, "<<rgbR<<","<<rgbG<<","<<rgbB<<",";
	ss << "SPHERE, "<<_cp.getX()<<","<<_cp.getY()<<","<<_cp.getZ()<<","<<_radius<<" ]";

	pymolObjectStrings[_name] = ss.str();

	return true;
}

bool PyMolVisualization::createAtom(CartesianPoint &_a, string _name, double _radiusScaling){

	Atom a;
	a.setCoor(_a);
	a.setResidueNumber(1);
	a.setChainId("Z");
	a.setName("FOO");
	a.setResidueName("FOO");

	return createAtom(a,_name,_radiusScaling);
}
bool PyMolVisualization::createAtom(Atom &_a, string _name, double _radiusScaling){


	stringstream ss;

	if (_name== ""){
		_name = "foo";
	}

	ss << "cmd.read_pdbstr(\""<<PDBFormat::createAtomLine(PDBFormat::createAtomData(_a))<<"\",\""<<_name<<"\");"<<endl;
	
	if (_radiusScaling != 1.0){
		ss << "cmd.do(\"set sphere_scale,"<<_radiusScaling<<","<<_name<<"\");"<<endl;
	}

	pymolAtomStrings[_name] = ss.str();

	return true;
}

bool PyMolVisualization::existPyMolObject(string _name){

	map<string,string>::iterator it;
	it = pymolObjectStrings.find(_name);

	return (it != pymolObjectStrings.end());
}

string PyMolVisualization::getPyMolObject(string _name){

	map<string,string>::iterator it;
	it = pymolObjectStrings.find(_name);


	if (it != pymolObjectStrings.end()){
		return it->second;
	}

	return "";
}


string PyMolVisualization::toString() {
	
	stringstream ss;

	ss << "# MSL PyMOLVisualization Object"<<endl;
	ss << "from pymol.cgo import *"<<endl;
	ss << "from pymol import cmd"<<endl;

	for (map<string,string>::iterator it = pymolObjectStrings.begin(); it != pymolObjectStrings.end();it++){
		ss << it->second<<endl;
		ss << "cmd.load_cgo("<<it->first<<",'"<<it->first<<"')"<<endl;
	}

	ss <<endl;

	for (map<string,string>::iterator it = pymolAtomStrings.begin(); it != pymolAtomStrings.end();it++){
		ss << it->second<<endl;
	}
	ss <<endl;
	ss << "# End "<<endl;

	return ss.str();
}
