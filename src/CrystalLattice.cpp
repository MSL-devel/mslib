/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan
 Sabareesh Subramaniam, Ben Mueller

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
#include "CrystalLattice.h"

using namespace MSL;
using namespace std;


CrystalLattice::CrystalLattice(){
	pdbFileRead = false;
	pdbFile = "";
}
CrystalLattice::CrystalLattice(string _pdbFile){

	pdbFileRead = false;
	pdbFile = _pdbFile;
	
}

CrystalLattice::~CrystalLattice(){
	map<string,AtomPointerVector *>::iterator it;
	for (it = crystalUnits.begin();it != crystalUnits.end();it++){
		for (AtomPointerVector::iterator k=it->second->begin(); k!=it->second->end(); k++) {
			delete (*k);
		}
		delete(it->second);
	}
}


void CrystalLattice::generateCrystal(){

	// Read PDB file if haven't already.
	if (pdbFile != "" && !pdbFileRead){
		readPdb();
	}

	// Get everything relevant from PDBReader

	// Scale matrix gives ortho-to-fractional coordinates (SCALE) (we will invert to get fractional-to-ortho later in the code)
	Matrix &scaleMat                    = pin.getScaleRotation();
	//CartesianPoint &scaleTrans          = pin.getScaleTranslation();

	// Symmetry matrices gives symmetry related molecules (REMARK 290)
	vector<Matrix  *> &symMats          = pin.getSymmetryRotations();
	vector<CartesianPoint  *> &symTrans = pin.getSymmetryTranslations();

	// Bounding-box of PDB coordinates
	//map<string,double> &bounds          = pin.getBoundingCoordinates();

	// Unit cell parameters, redundant with SCALE, but sometimes doesn't work (3DVH)
	//vector<double>& unitCellParams      = pin.getUnitCellParameters();

	if (symMats.size() == 0){
		cerr << "ERROR 1915 CrystalLattice::generateCrystal() PDB file missing SCALE and/or REMARK 290, which are required in MSL for crystal lattice generation\n";
		exit(1915);
	}


	// Store starting structure
	AtomPointerVector *ats = new AtomPointerVector(pin.getAtomPointers());
	//ats->updateGeometricCenter();
	crystalUnits["orig"] = ats;

	// Invert the Scale Matrix to go from fractional to ortho..
	Matrix scaleMatInv(3,3,0.0);
	double det=scaleMat[0][0]*(scaleMat[1][1]*scaleMat[2][2]-scaleMat[2][1]*scaleMat[1][2])-scaleMat[0][1]*(scaleMat[1][0]*scaleMat[2][2]-scaleMat[1][2]*scaleMat[2][0])+scaleMat[0][2]*(scaleMat[1][0]*scaleMat[2][1]-scaleMat[1][1]*scaleMat[2][0]);//adjoin
    	scaleMatInv[0][0]=(scaleMat[1][1]*scaleMat[2][2]-scaleMat[2][1]*scaleMat[1][2])/det;
    	scaleMatInv[0][1]=-(scaleMat[1][0]*scaleMat[2][2]-scaleMat[1][2]*scaleMat[2][0])/det;
    	scaleMatInv[0][2]=(scaleMat[1][0]*scaleMat[2][1]-scaleMat[2][0]*scaleMat[1][1])/det;
    	scaleMatInv[1][0]=-(scaleMat[0][1]*scaleMat[2][2]-scaleMat[0][2]*scaleMat[2][1])/det;
    	scaleMatInv[1][1]=(scaleMat[0][0]*scaleMat[2][2]-scaleMat[0][2]*scaleMat[2][0])/det;
    	scaleMatInv[1][2]=-(scaleMat[0][0]*scaleMat[2][1]-scaleMat[2][0]*scaleMat[0][1])/det;
    	scaleMatInv[2][0]=(scaleMat[0][1]*scaleMat[1][2]-scaleMat[0][2]*scaleMat[1][1])/det;
    	scaleMatInv[2][1]=-(scaleMat[0][0]*scaleMat[1][2]-scaleMat[1][0]*scaleMat[0][2])/det;
    	scaleMatInv[2][2]=(scaleMat[0][0]*scaleMat[1][1]-scaleMat[1][0]*scaleMat[0][1])/det;

	Transforms tr;
		
	// For each symmetry matrix generate +/- 1 Unit Cell.
	for (uint i = 0; i < symMats.size();i++){
		for (int a = -1 ; a < 2 ;a++){
			for (int b = -1; b < 2; b++){	
				for (int c = -1; c < 2; c++){	

					// Store atoms under approriate key
					char key[80];
					sprintf(key,"%02d%02d%02d%02d",i,a,b,c);
					//cout << "\tKEY: "<<key<<endl;

					// Copy initial atoms to new vector
					AtomPointerVector *newAts = new AtomPointerVector();
					copyAtoms(ats,newAts);
					
					// Do Symmetry Related transformations
				//	newAts->translate( (*symTrans[i] * -1) );
				//	newAts->rotate(*symMats[i]);
					tr.translate(*newAts, (*symTrans[i] * -1));
					tr.rotate(*newAts, *symMats[i]);

					// Do Unit Cell Transformations
					CartesianPoint p(a,b,c);
					CartesianPoint translateP = CartesianGeometry::matrixTransposeTimesCartesianPoint(p,scaleMatInv);					
				//	newAts->translate(translateP);
					tr.translate(*newAts,translateP);
					//newAts->updateGeometricCenter();

					// Store atoms under approriate key
					crystalUnits[key] = newAts;

				}
			}
		}
	}

	



	// Unit Cell parameter usage...
	/*
	double a      = unitCellParams[0];
	double b      = unitCellParams[1];
	double c      = unitCellParams[2];
	double alpha  = unitCellParams[3] * M_PI / 180;
	double beta   = unitCellParams[4] * M_PI / 180;
	double gamma  = unitCellParams[5] * M_PI / 180;
	
	Matrix fractToOrtho (3,3,0.0);

	double cosAlpha  = ( ( cos(beta) * cos(gamma) - cos(alpha) ) / (sin(beta)*sin(gamma)));
	double sinAlpha  = sqrt(1 - cosAlpha*cosAlpha) ;
	fractToOrtho[0][0] = a;
	fractToOrtho[0][1] = b * cos(gamma);
	fractToOrtho[0][2] = c * cos(beta);
	fractToOrtho[1][0] = 0.0;
	fractToOrtho[1][1] = b *sin(gamma);
	fractToOrtho[1][2] = -c*sin(beta)* cosAlpha;
	fractToOrtho[2][0] = 0.0;
	fractToOrtho[2][1] = 0.0;
	fractToOrtho[2][2] = c*sin(beta) * sinAlpha;
	*/
	
		
}

void CrystalLattice::writeCrystalUnits(string _pathAndPrefix,bool _closeContactsOnly,bool _singleFile, string _renameChainsExcept,bool _nmrStyleFile){

	// sqrt3 * largest dimension should be length of diagnol (worse case)?
	double sqrt3 =  1.732050807568877;


	map<string,AtomPointerVector *>::iterator it;
	map<string,AtomPointerVector *>::iterator it2;

	// Check that we have the "original" unit
	it = crystalUnits.find("orig");
	if (it == crystalUnits.end()){
		cerr << "ERROR 1914 CrystalLattice::writeCrystalUnits...no original crystal unit found.\n";
		exit(1914);
	}

	AtomPointerVector *ats = it->second;

	// For each crystal unit we have
	int numUnitsPrinted = 1;

	// Print original Unit
	char nameOrig[180];
	sprintf(nameOrig,"%s%s-lattice.pdb",_pathAndPrefix.c_str(),_renameChainsExcept.c_str());

	if (_renameChainsExcept != ""){
		for (uint i = 0; i < ats->size();i++){
			if ((*ats)(i).getChainId() != _renameChainsExcept){
				(*ats)(i).setChainId("W");
			}
		}
	}

	pout.open(nameOrig);
	pout.write(*ats,true,false,_nmrStyleFile);  // addTer = true, removeHydrogens=false, nmrFile = _nmrStyleFile

	if (!_singleFile){
		pout.close();
	}

	/**********************************************************************
	 * the stamp is a random number that is used to recall the center of each AtomPointerVector
	 * group to avoid to calculate it multiple times (essentially the center is calculated
	 * the first time the getGeometricCenter is called and the value is cached and returned directly
	 * the next time.  The assumption is that the atoms do not move during the loop.
	 **********************************************************************/
	//unsigned int stamp = MslTools::getRandomInt(1000000);
	RandomNumberGenerator rng;
	unsigned int stamp = rng.getRandomInt(1000000);

	for (it = crystalUnits.begin();it != crystalUnits.end();it++){

		// Remove redundant units...slow but should work.
		bool redundant = false;
		for (it2 = crystalUnits.begin();it2 != crystalUnits.end();it2++){
			
			if (it != it2 && it->second->getGeometricCenter(stamp) == it2->second->getGeometricCenter(stamp)){
				redundant = true;
				break;
			}
		}
		if (redundant){
			continue;
		}

		// Only print when is not "orig" (dist == 0) and dist is less than sqrt(3)*largestDimensionOfProtein
		double dist = ats->getGeometricCenter(stamp).distance(it->second->getGeometricCenter(stamp));
		if ( dist != 0  && (!_closeContactsOnly || dist < sqrt3 * pin.getBoundingCoordinates()["maxDelta"])){

			if (!_singleFile){
				char name[180];
				sprintf(name,"%s-%s.pdb",_pathAndPrefix.c_str(),it->first.c_str());
				pout.open(name);
			}

			if (_renameChainsExcept != ""){
				for (uint i = 0; i < it->second->size();i++){	
					(*it->second)(i).setChainId("Z");
				}
			}
			pout.write((*it->second),true,false,_nmrStyleFile);

			if (!_singleFile){
				pout.close();
			}

			numUnitsPrinted++;
		}
	}

	fprintf(stdout, "%10d Crystal Units printed.\n",numUnitsPrinted);
}
void CrystalLattice::readPdb(){	

	if (pdbFile == ""){
		cerr << "ERROR 1931 CrystalLattice::readPdb(), no pdb file to read,use CrystalLattice::setPdbFile() or equivalent constructor."<<endl;
		exit(1931);
	}

	pin.open(pdbFile);
	pin.read();
	pin.close();
	pdbFileRead = true;
}


void CrystalLattice::copyAtoms(AtomPointerVector * _atoms, AtomPointerVector *newAts) {
	for (uint i = 0; i < _atoms->size();i++){
		newAts->push_back(new Atom(*(*_atoms)[i]));
	}
}


