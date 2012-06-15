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

#include "EnvironmentDescriptor.h"
#include "MslTools.h"
#include "AtomSelection.h"

using namespace MSL;
using namespace std;


EnvironmentDescriptor::EnvironmentDescriptor(){
	core  = new AtomPointerVector();
	frame = new Frame();
}

EnvironmentDescriptor::EnvironmentDescriptor(EnvironmentDescriptor & _ed){
	copy(_ed);
}

void EnvironmentDescriptor::operator=(EnvironmentDescriptor & _ed){
	copy(_ed);
}

EnvironmentDescriptor::~EnvironmentDescriptor(){
	for (uint i = 0; i < core->size();i++){
		delete((*core)[i]);
	}
	delete(core);
	delete(frame);

	map<string, AtomPointerVector *>::iterator eMapIt;
	for (eMapIt = environmentMap.begin(); eMapIt != environmentMap.end();eMapIt++){
		AtomPointerVector *t = eMapIt->second;
		string type   = eMapIt->first;

		for (uint i =0; i < t->size();i++){
			delete((*t)[i]);
		}

		t->clear();		
	}
	environmentMap.clear();


	map<string, Frame *>::iterator fMapit;
	for (fMapit = frameMap.begin(); fMapit != frameMap.end();fMapit++){
		Frame *ftmp = fMapit->second;

		delete(ftmp);
	}

	frameMap.clear();


}


void EnvironmentDescriptor::copy(EnvironmentDescriptor & _ed){

	// Make new Atoms from _ed.getCore()'s atoms.
	core   = new AtomPointerVector();
	AtomPointerVector tmp = _ed.getCore();
	for (uint i = 0; i < tmp.size();i++){
		core->push_back(new Atom(tmp(i)));
	}

	// Copy the frame
	frame  = new Frame(_ed.getReferenceFrame());


	// Copy all the environment atomvectors
	map<string, AtomPointerVector *> envMap = _ed.getEnvironmentMap();
	map<string, AtomPointerVector *>::iterator eMapIt;
	for (eMapIt = envMap.begin(); eMapIt != envMap.end();eMapIt++){
		AtomPointerVector *t = eMapIt->second;
		string type   = eMapIt->first;

		AtomPointerVector *tNew = new AtomPointerVector();
		for (uint i =0; i < t->size();i++){
			tNew->push_back(new Atom((*t)(i)));
		}
		
		environmentMap[type] = tNew;
		tNew = NULL;

	}


	// Copy all the environment Frames
	map<string, Frame *> fMap = _ed.getFrameMap();
	map<string, Frame *>::iterator fMapit;
	for (fMapit = fMap.begin(); fMapit != fMap.end();fMapit++){
		Frame *ftmp = fMapit->second;
		string type = fMapit->first;

		frameMap[type] = new Frame(*ftmp);
	}
	

	
}


AtomPointerVector & EnvironmentDescriptor::getCore(){
	return *core;
}


void EnvironmentDescriptor::setCore(AtomPointerVector &_atoms){

	for (uint i = 0; i < _atoms.size();i++){
		core->push_back(new Atom(_atoms(i)));
	}
	
}

Frame & EnvironmentDescriptor::getReferenceFrame(){
	return *frame;
}


void EnvironmentDescriptor::setReferenceFrame(Frame &_frame){

	frame = new Frame(_frame);
}




Frame & EnvironmentDescriptor::getEnvironmentFrame(string _environmentType){

	// TODO test if it exsists
	return *frameMap[_environmentType];
	
}
AtomPointerVector & EnvironmentDescriptor::getEnvironment(string _environmentType){

	// TODO test if it exsists
	return *environmentMap[_environmentType];
}
void EnvironmentDescriptor::setEnvironment(string _environmentType, AtomPointerVector &_atoms){

	// Create copy of the atoms in environment..
	AtomPointerVector *tmp = new AtomPointerVector();
	
	for (uint i = 0; i < _atoms.size();i++){
		tmp->push_back(new Atom(_atoms(i)));
	}

	environmentMap[_environmentType] = tmp;

	// Create Copy of Environment Frame
	Frame *tmpFrame = new Frame();
	tmpFrame->computeFrameFromPCA(*tmp);

	frameMap[_environmentType] = tmpFrame;
}


map<string,AtomPointerVector*> & EnvironmentDescriptor::getEnvironmentMap(){
	return environmentMap;
}
map<string,Frame*> & EnvironmentDescriptor::getFrameMap() {
	return frameMap;
}




string EnvironmentDescriptor::generateLookupKey(string _envType){
	/*
	  Key lookup: 
	  DistanceBin:EigenVector1Bin:EigenVector2Bin:EigenVector3Bin
	*/

	// We could store this rather than computing here..
	Frame *f = frameMap[_envType];
	Matrix transform = f->getBasisTransformMatrix(*frame);
	
	Matrix angles(3,3);
	for (uint i = 0; i < angles.getRows();i++){
		for (uint j = 0; j < angles.getCols();j++){
			angles[i][j] = acos(transform[i][j]) * 180/M_PI;
		}
	}

	double dist    = MslTools::smartRound(f->distanceToFrame(*frame),2);
	double angles1 = MslTools::smartRound(angles[0][0],10); 
	double angles2 = MslTools::smartRound(angles[1][1],10); 
	double angles3 = MslTools::smartRound(angles[2][2],10); 

	stringstream ss;
	ss << dist<<":"<<angles1<<":"<<angles2<<":"<<angles3;

	return ss.str();
}


bool EnvironmentDescriptor::setupDescriptor(Residue  &_res, System &_sys, string type){


	// Add refFrame
	Frame refFrame;
	refFrame.setName("ref");
	if (_res.getResidueName() != "GLY"){
		refFrame.computeFrameFrom3Atoms(_res("N"), _res("CA"), _res("CB"));
	} else {
		CartesianPoint CBcoor = CartesianGeometry::build(_res("CA").getCoor(), _res("N").getCoor(), _res("C").getCoor(), 1.521, 110.5, -122.5);
		Atom CB;
		CB.setCoor(CBcoor);
		refFrame.computeFrameFrom3Atoms(_res("N"), _res("CA"), CB);
	}




	// Add environment..
	int startRes = _res.getResidueNumber() - 7;
	if (startRes < 0){
		startRes = 0;
	}
	int endRes = _res.getResidueNumber() + 7;
	if (endRes > _sys(_res.getChainId()).positionSize()){
		endRes = _sys(_res.getChainId()).positionSize();
	}

	char a[200];


	if (type == "POLAR"){
		sprintf(a,"res, ((resn CYS+TYR+SER+THR+ASP+GLU+ARG+LYS+HIS+ASN+GLN and name CA) and not resi %d-%d) WITHIN %d OF ((chain %s and resi %d) and name CA)",startRes,endRes,10,_res.getChainId().c_str(),_res.getResidueNumber());
	} else if (type == "HYDROPHOBIC"){
		sprintf(a,"res, ((resn VAL+ILE+LEU+PHE+TRP+MET+ALA+PRO and name CA) and not resi %d-%d) WITHIN %d OF ((chain %s and resi %d) and name CA)",startRes,endRes,10,_res.getChainId().c_str(),_res.getResidueNumber());
	} else if (type == "SMALL"){
		sprintf(a,"res, ((resn ALA+GLY+SER+THR+CYS and name CA) and not resi %d-%d) WITHIN %d OF ((chain %s and resi %d) and name CA)",startRes,endRes,10,_res.getChainId().c_str(),_res.getResidueNumber());
	} else if (type == "BB"){
		sprintf(a,"res, name N+CA+C+O+CB and not resi %d-%d WITHIN %d OF ((chain %s and resi %d) and name CA)",startRes,endRes,10,_res.getChainId().c_str(),_res.getResidueNumber());
	} else if (type == "POLARATOMS"){
		// Maybe CE2 this is an aromtic specific atmo type (only in PHE,TRP,TYR) ... however its not symmetric so may get weird results
		//   Also you don't want to include ALL atoms of aromatic rings because you will over weight them in the PCA.
		sprintf(a,"res, (name OD1+OD2+OE1+OE2+OH+NH1+NH2+NZ+ND2+NE2+ND1+NE2+OG+OG1+NE1 and not resi %d-%d) WITHIN %d OF ((chain %s and resi %d) and name OD1+OD2+OE1+OE2+OH+NH1+NH2+NZ+ND2+NE2+ND1+NE2+OG+OG1+NE1)",startRes,endRes,10,_res.getChainId().c_str(),_res.getResidueNumber());
	} else{
		sprintf(a,"res, name CA and not resi %d-%d WITHIN %d OF ((chain %s and resi %d) and name CA)",startRes,endRes,10,_res.getChainId().c_str(),_res.getResidueNumber());
	}

	//cout << "\tUsing environment selection: "<<(string)a<<endl;
	AtomSelection sel(_sys.getAtomPointers());
	AtomPointerVector env = sel.select(string(a));
	
	// Bail out and don't add if less than 3 residues in environment, most likely a surface residue.
	if (env.size() < 3) return false;


	
	setCore(_res.getAtomPointers());			
	setReferenceFrame(refFrame);
 	setEnvironment(type, env);

	return true;

}

void EnvironmentDescriptor::setName(string _name){
	name = _name;
}

string EnvironmentDescriptor::getName(){
	return name;
}


