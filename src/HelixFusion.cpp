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

#include "HelixFusion.h"
#include "AtomVector.h"
#include "Transforms.h"
#include "MslTools.h"
#include "PyMolVisualization.h"
#include "System.h"

bool HelixFusion::fusionByAtomicAlignment(double _rmsdTolerance, string _newChainId){

	deleteFusedChains();

	nIndex.clear();
	cIndex.clear();
	vector<double> rmsds;
	// try 
	//   n = 0 --> residues 0-3
	//   n = 1 --> residues 1-4
	//   ...
	if (nterm == NULL){
		cerr << "NTERM chain is NULL.\n";	
		exit(3423);
	}

	if (cterm == NULL){
		cerr << "CTERM chain is NULL.\n";	
		exit(3424);
	}

	cterm->getAtoms().saveCoor("start");
	for (uint n = 0; n < 3 && n < nterm->size(); n++){
		AtomVector ngroup;
		ngroup.push_back( &(*nterm).getResidueByIndex(n)("CA") );
		ngroup.push_back( &(*nterm).getResidueByIndex(n+1)("CA") );
		ngroup.push_back( &(*nterm).getResidueByIndex(n+2)("CA") );
		ngroup.push_back( &(*nterm).getResidueByIndex(n+3)("CA") );
		
		// Save pre-alignment coordinates
		ngroup.saveCoor("pre");

		// try 
		//   c = size-8 --> residues size-8 - size-5
		//   c = size-7 --> residues size-7 - size-4
		//   ...
		for (uint c = cterm->size()-8; c < cterm->size()-3 && c > 0; c++){

			AtomVector cgroup;
			cgroup.push_back( &(*cterm).getResidueByIndex(c)("CA") );
			cgroup.push_back( &(*cterm).getResidueByIndex(c+1)("CA") );
			cgroup.push_back( &(*cterm).getResidueByIndex(c+2)("CA") );
			cgroup.push_back( &(*cterm).getResidueByIndex(c+3)("CA") );


		
			// Save pre-alignment coordinates
			cgroup.saveCoor("pre");

			Transforms t;
			if (!t.align(cgroup,ngroup)){
				cerr << "ERROR 3242 HelixFusion::fusionByAtomicAlignment alignment failed.\n";
				exit(3242);
			}

			
			double rmsd = cgroup.rmsd(ngroup);

			if (rmsd < _rmsdTolerance){
				nIndex.push_back(n);
				cIndex.push_back(c);
				rmsds.push_back(rmsd);
			}

			cgroup.applySavedCoor("pre");

		}

		
	}

	for (uint f = 0; f < nIndex.size();f++){

		AtomVector ngroup;
		ngroup.push_back( &(*nterm).getResidueByIndex(nIndex[f])("CA") );
		ngroup.push_back( &(*nterm).getResidueByIndex(nIndex[f]+1)("CA") );
		ngroup.push_back( &(*nterm).getResidueByIndex(nIndex[f]+2)("CA") );
		ngroup.push_back( &(*nterm).getResidueByIndex(nIndex[f]+3)("CA") );


		AtomVector cgroup;
		cgroup.push_back( &(*cterm).getResidueByIndex(cIndex[f])("CA") );
		cgroup.push_back( &(*cterm).getResidueByIndex(cIndex[f]+1)("CA") );
		cgroup.push_back( &(*cterm).getResidueByIndex(cIndex[f]+2)("CA") );
		cgroup.push_back( &(*cterm).getResidueByIndex(cIndex[f]+3)("CA") );



		cterm->getAtoms().saveCoor("pre");


		// Transform all atoms using best sets.
		Transforms t;
		if (!t.align(cgroup,ngroup,cterm->getAtoms())){
			cerr << "ERROR 3243 HelixFusion::fusionByAtomicAlignment alignment failed.\n";
			exit(3243);
		}


			
			
		// Now remove atoms that are not needed.
		//   Save nterm: all
		//   Save cterm: 0->cIndex-1
		//nterm->setChainId("A");
		Chain *fuse = new Chain();
		int newResNum = 1;
		for (uint i = 0; i <= cIndex[f]+1;i++){	

			AtomVector &tmp = cterm->getPositionByIndex(i).getAtoms();
			AtomVector newtmp;
			for (uint j = 0;j<tmp.size();j++){
				newtmp.push_back(new Atom());

				// New values
				newtmp.back()->setChainId(_newChainId);
				newtmp.back()->setResidueNumber(newResNum);

				// Take values from cterm atom
				newtmp.back()->setResidueName(tmp(j).getResidueName());
				newtmp.back()->setName(tmp(j).getName());
				newtmp.back()->setElement(tmp(j).getElement());
				newtmp.back()->setType(tmp(j).getType());
				newtmp.back()->setCharge(tmp(j).getCharge());
				newtmp.back()->setTempFactor(tmp(j).getTempFactor());
				newtmp.back()->setSegID(tmp(j).getSegID());
				newtmp.back()->setCoor(tmp(j).getCoor());
			}
			newResNum++;
			fuse->addAtoms(newtmp);
			newtmp.deletePointers();
		}

		//cout << "NINDEX: "<<nIndex[f]<<" for "<<f<<endl;
		for (uint i = nIndex[f]+2; i < nterm->size();i++){	


			AtomVector &tmp = nterm->getPositionByIndex(i).getAtoms();
			//cout << "Taking "<<nterm->getPositionByIndex(i).getResidueNumber()<<" "<<nterm->getPositionByIndex(i).getResidueName()<<" "<<nterm->getPositionByIndex(i).getChainId()<<" "<<tmp.size()<<endl;

			AtomVector newtmp;
			for (uint j = 0;j<tmp.size();j++){
				newtmp.push_back(new Atom());

				// New values
				newtmp.back()->setChainId(_newChainId);
				newtmp.back()->setResidueNumber(newResNum);

				// Take values from cterm atom
				newtmp.back()->setResidueName(tmp(j).getResidueName());
				newtmp.back()->setName(tmp(j).getName());
				newtmp.back()->setElement(tmp(j).getElement());
				newtmp.back()->setType(tmp(j).getType());
				newtmp.back()->setCharge(tmp(j).getCharge());
				newtmp.back()->setTempFactor(tmp(j).getTempFactor());
				newtmp.back()->setSegID(tmp(j).getSegID());
				newtmp.back()->setCoor(tmp(j).getCoor());
			}
			newResNum++;
			fuse->addAtoms(newtmp);
			newtmp.deletePointers();

		}
		
		fuse->setChainId(_newChainId);
		fused.push_back(fuse);
		/*

		PDBWriter pout;
		stringstream ss; 
		ss << "/tmp/f-"<<f<<".pdb";
		pout.open(ss.str());
		pout.write(ngroup);
		pout.write(cgroup);
		pout.close();

		stringstream ss2;
		ss2 <<"/tmp/Fuse-"<<f<<".pdb";
		pout.open(ss2.str());
		pout.write(fuse->getAtoms());
		pout.close();
		*/

		cterm->getAtoms().applySavedCoor("pre");


	}

	cterm->getAtoms().applySavedCoor("start");

	// No fusions  under rmsdTol...
	if (cIndex.size() == 0){
		return false;
	}


	return true;
	
}

bool HelixFusion::fusionByHelicalFrames(){


	deleteFusedChains();

	PyMolVisualization pviz;
	vector<CartesianPoint *> ntermHelixAxes;
	vector<int> nIndex;
	vector<int> cIndex;
	for (uint n = 0; n < 7 && n < nterm->size(); n++){
		
		Atom &a1 = (*nterm).getResidueByIndex(n)("CA");
		Atom &a2 = (*nterm).getResidueByIndex(n+1)("CA");
		Atom &a3 = (*nterm).getResidueByIndex(n+2)("CA");

		
		CartesianPoint v1 = (a1.getCoor() - a2.getCoor()).getUnit();
		CartesianPoint v2 = (a3.getCoor() - a2.getCoor()).getUnit();
		CartesianPoint *p = new CartesianPoint();
		
		*p = (v1 + v2).getUnit()*2.3 + a2.getCoor();

		ntermHelixAxes.push_back(p);
		nIndex.push_back(n);

		char name[30];
		sprintf(name,"nterm%04d",n);

		pviz.createAtom(*p,(string)name,0.5);

	}


	vector<CartesianPoint *> ctermHelixAxes;
	for (uint c = cterm->size()-8; c < cterm->size()-2 && c > 0; c++){

		AtomVector cgroup;
		Atom &a1 = (*cterm).getResidueByIndex(c)("CA");
		Atom &a2 = (*cterm).getResidueByIndex(c+1)("CA");
		Atom &a3 = (*cterm).getResidueByIndex(c+2)("CA");


		CartesianPoint v1 = (a1.getCoor() - a2.getCoor()).getUnit();
		CartesianPoint v2 = (a3.getCoor() - a2.getCoor()).getUnit();
		CartesianPoint *p = new CartesianPoint();
		
		*p = (v1 + v2).getUnit()*2.3 + a2.getCoor();

		ctermHelixAxes.push_back(p);
		cIndex.push_back(c);

		char name[30];
		sprintf(name,"cterm%04d",c);

		pviz.createAtom(*p,(string)name,0.5);
		
	}



	// Write out PyMol object with helix axes points
	ofstream fout;
	fout.open("axes.py");
	fout << pviz;
	fout.close();

	
	Atom *n1 = new Atom("N1");
	Atom *n2 = new Atom("N2");
	Atom *n3 = new Atom("N3");

	
	Atom *c1 = new Atom("C1");
	Atom *c2 = new Atom("C2");
	Atom *c3 = new Atom("C3");
	for (uint n = 0; n < ntermHelixAxes.size();n++){
		
		n1->setCoor(ntermHelixAxes[n]->getCoor());
		n2->setCoor(ntermHelixAxes[n+1]->getCoor());
		n3->setCoor(ntermHelixAxes[n+2]->getCoor());

		AtomVector nv;
		nv.push_back(n1);
		nv.push_back(n2);
		nv.push_back(n3);

		for (uint c = 0; c < ctermHelixAxes.size();c++){


			c1->setCoor(ctermHelixAxes[c]->getCoor());
			c2->setCoor(ctermHelixAxes[c+1]->getCoor());
			c3->setCoor(ctermHelixAxes[c+2]->getCoor());

			AtomVector cv;
			cv.push_back(c1);
			cv.push_back(c2);
			cv.push_back(c3);

			cterm->getAtoms().saveCoor("pre");

			Transforms t;
			if (!t.align(cv,nv,cterm->getAtoms())){
				
				// Calculate angle
				//  Nca2-n1-c1-Cca2 ?
				//double angle = nterm[nIndex[n]]->dihedral(n1,c1,*cterm[cIndex[c]]);

				// Rotate about n1-c1  , apply to cterm.
				//cterm->getAtoms().rotate(CartesianGeometry::instance()->getRotationMatrix(angle,(n1-n2).getUnit()));
				

				
				
			}

			cterm->getAtoms().applySavedCoor("pre");

			
		}
	}
			

	return true;
}

bool HelixFusion::fusionByHydrogenBonding(){
	return true;
}






