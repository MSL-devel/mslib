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


bool RotamerLibraryBuilder::addRotamer(Residue & _res, std::string _copyLibName, std::string _newLibName, unsigned int _bin) {

  // Add new library if needed
  if (!pRotlib->libraryExists(_newLibName)){
    pRotlib->addLibrary(_newLibName);
  }

  // Add new residue name, mobileAtoms and defi lines if needed
  if (!pRotlib->residueExists(_newLibName,_res.getResidueName())){
      pRotlib->addResidue(_newLibName,_res.getResidueName());
      pRotlib->addMobileAtoms(_newLibName,_res.getResidueName(),pRotlib->getMobileAtoms(_copyLibName,_res.getResidueName()));



      // Add the Coor lines... can't believe I have to re-parse the original rotamer library lines!
      std::vector<std::string> lines = pRotlib->getInternalCoorDefinitionLines(_copyLibName,_res.getResidueName());

      for (uint l = 0; l < lines.size();l++){
	std::string line = MslTools::uncomment(lines[l]);
	vector<string> tokens = MslTools::tokenize(line);

	// found a bond, angle, improper or dihedral definiton
	if (tokens[0] == "DEFI") {

	  tokens.erase(tokens.begin()); // remove the DEFI token
	  while (tokens.size() < 4) {
	    // add one or two blanks if it was a bond or angle defintion
	    tokens.push_back("");
	  }

	  //cout << "Add internal coor defi: "<<tokens[0]<<" "<<tokens[1]<<" "<<tokens[2]<<" "<<tokens[3]<<endl;
	  pRotlib->addInternalCoorDefinition(_newLibName,_res.getResidueName(), tokens);
	}
      }
  }

  return addRotamer(_res,_newLibName,_bin);
}
bool RotamerLibraryBuilder::addRotamer(Residue & _res, string _libName, unsigned int _bin) {

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
			if (_res.atomExists((*l).atomNames[i])) {
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
				//values[l-defis.begin()] = CartesianGeometry::distance(atoms[0]->getCoor(),atoms[1]->getCoor());
				
				break;
			case 1:
				//angle

				//values[l-defis.begin()] = CartesianGeometry::angle(atoms[0]->getCoor(),atoms[1]->getCoor(),atoms[2]->getCoor());
				values[l-defis.begin()] = atoms[0]->angle(*atoms[1],*atoms[2]);
				break;
			
			case 2:
				// dihedral uses same math of improper
			case 3:
				// improper

				//values[l-defis.begin()] = CartesianGeometry::dihedral(atoms[0]->getCoor(),atoms[1]->getCoor(),atoms[2]->getCoor(),atoms[3]->getCoor());

				values[l-defis.begin()] = atoms[0]->dihedral(*atoms[1],*atoms[2],*atoms[3]);
				break;

			default:
				cerr << "Warning: 31391 Wrong Defi-type. Storing 0s" << endl;
			
		}		

	}
	

	pRotlib->addConformation(_libName, resName, values,_bin);
	return true;
}


