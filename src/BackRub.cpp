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
#include "BackRub.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "PDBWriter.h"
#include "PyMolVisualization.h"


BackRub::BackRub(){
}

BackRub::~BackRub(){
}

string BackRub::localSample(Chain &_ch, int _startResIndex, int _endResIndex, int _numFragments){

	Residue &stem1 = _ch.getResidueByIndex(_startResIndex);
	Residue &stem2 = _ch.getResidueByIndex(_endResIndex);

	if (!stem1.exists("CA") || !stem2.exists("CA")){
		cerr << "ERROR BackRub::localSample() one or both stem residues does not contain a 'CA' atom."<<endl;	
		exit(1);
	}

	

	
	RandomNumberGenerator rng;
	rng.setRNGTimeBasedSeed();
	//double seed = rng.getRNGSeed();

	stringstream ss;
	Transforms t;
	for (uint f = 0; f < _numFragments;f++){

		// Randomly pick a 3 residue rotation vector, within our bounds of stem1,stem2
		int randRange=(_endResIndex - _startResIndex -3 );
		int randInt  = rng.getRandomIntLimit( randRange);
		int startRes = _startResIndex + randInt;
		int endRes   = startRes + 2;
		
		//CartesianPoint mainRotVector = stem2("CA").getCoor() - stem1("CA").getCoor();		
		Residue &res1 = _ch.getResidueByIndex(startRes);
		Residue &res2 = _ch.getResidueByIndex(startRes+1);
		Residue &res3 = _ch.getResidueByIndex(endRes);
		
		if (!res1.exists("CA")  || !res2.exists("CA") || !res3.exists("CA")){
			cerr << "ERROR BackRub::localSample() res1,res2,res3 for defining mainRotationAxis does not have CA atom.\n";
			exit(1);
		}

		if (!res2.exists("O")  || !res3.exists("O")){
			cerr << "ERROR BackRub::localSample() res2,res3 for defining mainRotationAxis does not have O atom.\n";
			exit(1);
		}


		/*
		  Store some atom coordinates before major rotation.  We will use carboxyl oxygens for the minor rotations.
		     That is:     O1-Calpha1-Calpha2-O2

		     This seemed to work better than: O1-Calpha1-Calpha2-N2
		 */
		CartesianPoint preO1 = res1("O").getCoor();
		CartesianPoint preO2 = res2("O").getCoor();


		// For second minor rotation we will try to use Oxygen of endRes+1 residue, if it doesn't exist use Nitrogen of endRes.
		CartesianPoint preAt2 = res2("N").getCoor();
		if (endRes+1 < _ch.size() && _ch.getResidueByIndex(endRes+1).exists("O")){
			
			preAt2 = _ch.getResidueByIndex(endRes+1)("O").getCoor();
		}


		CartesianPoint mainRotVector = res3("CA").getCoor();
		double omega = 0.0;
		while (omega == 0.0){
			omega = rng.getRandomIntLimit(20);
		}

		// assign  +/- randomly.
		if (rng.getRandomDouble() > 0.5){
			omega *= -1;
		}


		fprintf(stdout, "\nBackRub%04d: %s %3d %3s, %3d %3s : omega: %8.3f",f, res1.getChainId().c_str(), res1.getResidueNumber(), res1.getResidueName().c_str(), res3.getResidueNumber(), res3.getResidueName().c_str(),omega);
		// Rotate by omega

		// Rotate C=0, of res1 residue, by omega
		t.rotate(res1("C"),omega,mainRotVector, res1("CA").getCoor());
		t.rotate(res1("O"),omega,mainRotVector, res1("CA").getCoor());


		// Rotate all atoms of residues in between stems by omega  (for loop if larger than 3 residue fragment rotations...)
		t.rotate(res2.getAtoms(), omega, mainRotVector, res1("CA").getCoor());


		// Rotate N-H, of stem2 residue, by omega
		t.rotate(res3("N"),omega,mainRotVector, res1("CA").getCoor());
		if (res3.exists("H")){
			t.rotate(res3("H"),omega,mainRotVector, res1("CA").getCoor());
		}

		/*
		char pname[80];
		sprintf(pname,"/tmp/frag-%04d.pdb",f);
		PDBWriter pout;
		pout.open(pname);
		pout.write(_ch.getAtoms());
		pout.close();
		*/

		// ROTATE PEPTIDE BOND ATOMS BACK AS CLOSE TO ORIGINAL POSITION AS POSSIBLE

		// Do compensatory rotation along res1->res2 vector, try to put "O" atom back onto preO1. using O1-CA1-CA2-N1 dihedral angle
		//   this also should reduce strain placed on N-Calpha-C' bond angle due to initial rotation
		fprintf(stdout, " minor1 ");
		doMinorRotation(res1,res2,preO1,preO2);

		/*
		char pname2[80];
		sprintf(pname2,"/tmp/minor1-%04d.pdb",f);
		pout.open(pname2);
		pout.write(_ch.getAtoms());
		pout.close();
		*/

		// Do compensatory rotation along res2->res3 vector, try to put "O" atom back onto preO2. using O2-CA2-CA3-N2 dihedral angle
		//   this also should reduce strain placed on C'-Calpha-N bond angle due to initial rotation
		fprintf(stdout, " minor2 ");
		doMinorRotation(res2,res3,preO2,preAt2);
		fprintf(stdout,"\n");
		/*
		char pname3[80];
		sprintf(pname3,"/tmp/minor2-%04d.pdb",f);
		pout.open(pname3);
		pout.write(_ch.getAtoms());
		pout.close();
		*/

		// Write out to stringstream at this point...
		ss << "MODEL"<<endl;
		PDBWriter pss;
		pss.open(ss);
		pss.write(_ch.getAtoms());
		pss.close();		
		ss << "ENDMDL"<<endl;

	}


	return ss.str();
}

void BackRub::doMinorRotation(Residue &_r1, Residue &_r2, CartesianPoint &_targetAtom1, CartesianPoint &_targetAtom2){



	// Define minor rotation vector	
	CartesianPoint minorRotVector = _r2("CA").getCoor();

	// Get pre-dihedral
	double preDihedral = _targetAtom1.dihedral(_r1("CA").getCoor(),_r2("CA").getCoor(),_targetAtom2);

	// Get post-dihedral
	double postDihedral = _r1("O").getCoor().dihedral(_r1("CA").getCoor(),_r2("CA").getCoor(),_targetAtom2);

	// Rotate back by post-pre
	double alpha = postDihedral - preDihedral;
	fprintf(stdout, "preDihedral = %8.3f postDihedral = %8.3f alpha = %8.3f", preDihedral, postDihedral, alpha);




	// Create a local transform object to rotate atoms...
	Transforms t;
	t.rotate(_r1("O"),alpha, minorRotVector, _r1("CA").getCoor());
	t.rotate(_r1("C"),alpha, minorRotVector, _r1("CA").getCoor());
	t.rotate(_r2("N"),alpha, minorRotVector, _r1("CA").getCoor());
	if (_r2.exists("H")){
		t.rotate(_r2("H"),alpha, minorRotVector, _r1("CA").getCoor());
	}


}


