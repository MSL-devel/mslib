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

#include "PyMolVisualization.h"
#include "PDBFormat.h"


using namespace MSL;
using namespace std;



PyMolVisualization::PyMolVisualization(){
	pymolObjectStrings.clear();
	pymolAtomStrings.clear();
	rng.setSeed(49452);
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

                                     
bool PyMolVisualization::createArrow(CartesianPoint &_start, CartesianPoint &_vector, string _name,double _cylinderRadius, double _coneRadius, int rgbR,int rgbG,int rgbB){

	// Assign a random name, but for now...
	if (_name == ""){
		_name = "arrow";
	}
	stringstream ss;
	ss << _name << " = [ CYLINDER ,";
	ss <<_start[0]           <<","<<_start[1]           <<","<<_start[2]<<",";
	ss <<_start[0]+_vector[0]<<","<<_start[1]+_vector[1]<<","<<_start[2]+_vector[2]<<",";
	ss << _cylinderRadius <<","<<rgbR<<","<<rgbG<<","<<rgbB<<","<<rgbR<<","<<rgbG<<","<<rgbB<<",";

	ss << " CONE ,";
	ss <<_start[0]+_vector[0]<<","<<_start[1]+_vector[1]<<","<<_start[2]+_vector[2]<<",";
	ss <<(_start[0]+_vector[0]+0.2*_vector[0])<<","<<(_start[1]+_vector[1]+0.2*_vector[1])<<","<<(_start[2]+_vector[2]+0.2*_vector[2])<<",";
	ss << _coneRadius<<",0.0,"<<rgbR<<","<<rgbG+0.5<<","<<rgbB-0.5<<","<<rgbR<<","<<rgbG+0.5<<","<<rgbB-0.5<<",1.0,1.0 ]";

	pymolObjectStrings[_name] = ss.str();

	return true;
}

bool PyMolVisualization::createCylinder(CartesianPoint &_start,CartesianPoint &_end, string _name,double _radius,double rgbR,double rgbG,double rgbB){
	
	// Assign a random name, but for now...
	if (_name == ""){
		_name = "cylinder";
	}
	if (_name == "random"){


	  int random = rng(1000);
	  char c[80];
	  sprintf(c,"cylinder_%04d",random);
	  _name = (string) c;

	}

	stringstream ss;
	ss << _name << " = [ CYLINDER ,";
	ss <<_start[0]           <<","<<_start[1]           <<","<<_start[2]<<",";
	ss <<_end[0]<<","<<_end[1]<<","<<_end[2]<<",";
	ss << _radius <<","<<rgbR<<","<<rgbG<<","<<rgbB<<","<<rgbR<<","<<rgbG<<","<<rgbB<<" ] ";

	pymolObjectStrings[_name] = ss.str();

	return true;

}

bool PyMolVisualization::createCone(CartesianPoint &_start,CartesianPoint &_end, string _name,double _radius,int rgbR,int rgbG,int rgbB){
	
	// Assign a random name, but for now...
	if (_name == ""){
		_name = "cone";
	}

	stringstream ss;
	ss << _name << " = [ CONE ,";
	ss <<_start[0]           <<","<<_start[1]           <<","<<_start[2]<<",";
	ss <<_end[0]<<","<<_end[1]<<","<<_end[2]<<",";
	ss << _radius<<",0.0,"<<rgbR<<","<<rgbG<<","<<rgbB<<","<<rgbR<<","<<rgbG<<","<<rgbB<<",1.0,1.0 ]";


	pymolObjectStrings[_name] = ss.str();

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
	for (map<string,string>::iterator it = pymolSelectionStrings.begin(); it != pymolSelectionStrings.end();it++){
		ss << it->second<<endl;
	}
	ss <<endl;
	ss << "# End "<<endl;

	return ss.str();
}

bool PyMolVisualization::createSelection(string &_selName, string &_sel){

  stringstream ss;
  ss << "cmd.do(\"select "<<_selName<<","<<_sel<<"\")"<<endl;
  ss << "cmd.do(\"show sticks, "<<_selName<<"\")"<<endl;
  pymolSelectionStrings[_selName] = ss.str();
  return true;
}
