#include "RotamerLibraryBuilder.h"

using namespace MSL;
using namespace std;


RotamerLibraryBuilder::RotamerLibraryBuilder() {
	pRotlib = NULL;
}

RotamerLibraryBuilder::RotamerLibraryBuilder(RotamerLibrary * _pRotlib) {
	pRotlib = _pRotlib;
}


RotamerLibraryBuilder::RotamerLibraryBuilder(const RotamerLibraryBuilder & _rotLibbuild) {
	pRotlib = _rotLibbuild.pRotlib;
}

RotamerLibraryBuilder::~RotamerLibraryBuilder() {
}


bool RotamerLibraryBuilder::addRotamer(Residue & _res, string _libName) {

	if(_libName == "") {
		_libName = (pRotlib->getLibraryNames())[0];
	}
// getter/setter for defaultLibraryName in RotamerLibrary
	string resName = _res.getResidueName();


	if(!pRotlib->libraryExists(_libName) || !pRotlib->residueExists(_libName, resName)) {
		cerr << "Set up the Rotamer Library object with the residue and library name before passing to the RotamerLibraryBuilder. Library " << _libName << " or Residue " << resName << " is not present" << endl;
		return false;
	}

//	cout << "UUU libName " << _libName << " resName " << resName << endl;

	vector<RotamerLibrary::InternalCoorDefi> defis = pRotlib->getInternalCoorDefinition(_libName,resName);	
	vector <double> values(defis.size()); 
	
	for (vector<RotamerLibrary::InternalCoorDefi>::iterator l = defis.begin(); l != defis.end(); l++) {
		vector<Atom*> atoms((*l).atomNames.size());
		for (unsigned int i=0; i < (*l).atomNames.size(); i++) {
			//cout << "UUUUU Atom: " << (*l).atomNames[i] << endl;
			if (_res.exists((*l).atomNames[i])) {
				atoms[i] = &_res.getLastFoundAtom();
			} else {
				cerr << "WARNING 38913: atom " << (*l).atomNames[i] << " not found in bool RotamerLibraryBuilder::addRotamer(Residue & _res, string _libName)" << endl;
				return false;
			}
		}
	
		switch((*l).type) {
			case 0:
				//bond

				values[l-defis.begin()] = atoms[0]->distance(*atoms[1]);
				//values[l-defis.begin()] = CartesianGeometry::instance()->distance(atoms[0]->getCoor(),atoms[1]->getCoor());
				
				break;
			case 1:
				//angle

				//values[l-defis.begin()] = CartesianGeometry::instance()->angle(atoms[0]->getCoor(),atoms[1]->getCoor(),atoms[2]->getCoor());
				values[l-defis.begin()] = atoms[0]->angle(*atoms[1],*atoms[2]);
				break;
			
			case 2:
				// dihedral uses same math of improper
			case 3:
				// improper

				//values[l-defis.begin()] = CartesianGeometry::instance()->dihedral(atoms[0]->getCoor(),atoms[1]->getCoor(),atoms[2]->getCoor(),atoms[3]->getCoor());

				values[l-defis.begin()] = atoms[0]->dihedral(*atoms[1],*atoms[2],*atoms[3]);
				break;

			default:
				cerr << "Warning: 31391 Wrong Defi-type. Storing 0s" << endl;
			
		}		

	}
	

	pRotlib->addConformation(_libName, resName, values);
	return true;
}


