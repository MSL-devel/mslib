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

#include <queue>

// General Includes
#include "System.h"
#include "PDBReader.h"
#include "MslTools.h"

// Specific objects 
#include "ChiStatistics.h"
#include "Transforms.h"
#include "CartesianPoint.h"
#include "CartesianGeometry.h"
#include "BBQTable.h"
#include "CCD.h"
#include "RandomNumberGenerator.h"
#include "Quench.h"

using namespace std;
using namespace MslTools;

#include <Python.h>

// Global variables..

// Store backbone fragments as RMSD,PDBSTRING
priority_queue< pair<double,string>, vector< pair<double,string> >, less<pair<double,string> > > matchingFragments;

// FUNCTIONS
static PyObject* 
getChi(PyObject *self, PyObject *args) {


	char *pdbstr;
	if (!PyArg_ParseTuple(args,"s",&pdbstr))
		  return NULL;

	ChiStatistics chi;

	string pdb = (string)pdbstr;
	PDBReader rin;
	rin.read(pdb);
	rin.close();

	System sys;
	sys.addAtoms(rin.getAtoms());
	for (uint i = 0 ; i < sys.residueSize();i++){
		Residue &r = sys.getResidue(i);
		if (chi.getNumberChis(r) == -1) continue;

		fprintf(stdout, "%1s %3d %3s ",r.getChainId().c_str(),r.getResidueNumber(), r.getResidueName().c_str());
		for (uint c = 0; c < (uint)chi.getNumberChis(r);c++){
			if (!(chi.atomsExist(r,c+1))) {
				fprintf(stdout, " ---- MISSING ATOMS ---- ");
				break;
			}
			double angle = chi.getChi(r,c+1);
			if (angle != MslTools::doubleMax){
				fprintf(stdout,"%8.2f ",angle);
			}
		}
		fprintf(stdout,"\n");
	}

	return Py_BuildValue("s","");
} 

static PyObject* 
getFragments(PyObject *self, PyObject *args) {

	int numTopFragments;

	if (!PyArg_ParseTuple(args,"i",&numTopFragments))
		  return NULL;

	// Go through matches and make an NMR style file...
	stringstream ss;
	// Writer.writeln("MODEL") .... "ENDMDL"..
	int index = 1;
	while (!matchingFragments.empty()){
		char model[10];
		sprintf(model,"MODEL %5d",index++);
		cout << "WRITING MODEL : "<<index<<" ("<<model<<")"<<endl;
		ss << model<<endl;
		ss << matchingFragments.top().second<<endl;
		ss << "ENDMDL\n";
		matchingFragments.pop();

		if (numTopFragments != -1 && index == numTopFragments){
			break;
		}
		cout << "END LOOP"<<endl;
	}

	cout << "OUTSIDE LOOP"<<endl;
	cout << "STRINGSTREAM: "<<endl<<ss.str()<<endl;
	return Py_BuildValue("s",ss.str().c_str());

}
static PyObject* 
localSampling(PyObject *self, PyObject *args) {

	int numFragments;
	int maxDegree;
	char *fragment;
	char *bbqTable;

	if (!PyArg_ParseTuple(args,"siis",&fragment,&numFragments,&maxDegree,&bbqTable))
		  return NULL;

	// Convert PDB to string
	string pdb = (string)fragment;
	
	// Read in PDB to MSL
	PDBReader rin;
	rin.read(pdb);
	rin.close();


	// Get fragment to sample
	AtomVector &loop = rin.getAtoms();

	// Make CCD object with BBQTable .. (BBQTable is for other backbone atoms)
	CCD sampleCCD((string)bbqTable);

	// Do local sampling inside CCD object
	string pdbresult = sampleCCD.localSample(loop,numFragments,maxDegree);

	// Return PDB string of all atoms
	return Py_BuildValue("s",pdbresult.c_str());
	
}


static PyObject* 
searchForFragments(PyObject *self, PyObject *args) {

	int numResidues;
	char *stem1pdb,*stem2pdb,*database;

	if (!PyArg_ParseTuple(args,"sssi",&stem1pdb,&stem2pdb,&database,&numResidues))
		  return NULL;


	string pdb1 = (string)stem1pdb;
	string pdb2 = (string)stem2pdb;

	
	PDBReader rin1;
	rin1.read(pdb1);
	rin1.close();

	PDBReader rin2;
	rin2.read(pdb2);
	rin2.close();


	AtomVector &stem1 = rin1.getAtoms();
	AtomVector &stem2 = rin2.getAtoms();

	if (stem1.size() <= 1 || stem2.size() <= 1){
		cerr<< "Stems are too small."<<endl;
		return Py_BuildValue("s","");
	}

	if (numResidues == -1){
		numResidues = stem2(0).getResidueNumber() - stem1(stem1.size()-1).getResidueNumber() - 1;
	}
	//cout << "STEM1: "<<endl<<stem1.toString()<<endl;
	//cout << "STEM2: "<<endl<<stem2.toString()<<endl;
	AtomVector stems = stem1 + stem2;
	vector<double> stemDistanceSq;
	for (uint c = 0; c < stem1.size();c++){

		/*
		  fprintf(stdout,"C: %1s %4d %3s\n",
				stem1[c]->getChainId().c_str(),
				stem1[c]->getResidueNumber(),
			stem1[c]->getResidueName().c_str());
		*/
		for (uint n = 0; n < stem2.size();n++){

			double distSq = stem1(c).distance2(stem2(n));			
			/*
			fprintf(stdout,"\tN: %1s %4d %3s = %8.3f\n",
				stem2[n]->getChainId().c_str(),
				stem2[n]->getResidueNumber(),
				stem2[n]->getResidueName().c_str(),distSq);
			*/

			stemDistanceSq.push_back(distSq);
		}
	}

	AtomVector fragDB;
	fragDB.load_checkpoint((string)database);
	cout << "Loaded Fragment database: "<<database<<" total size: "<<fragDB.size()<<endl;

	//cout << "Stem sizes: "<<stem1.size()<<" "<<stem2.size()<<endl;

	// BBQ Table for backbone atoms
	BBQTable bbqT("/home/dwkulp/software/msl/tables/PiscesBBQTable.txt");

	double tol = 3;
	for (uint i = 0 ; i < fragDB.size()-(numResidues+stem1.size());i++){
		//cout << "I: "<<i<<endl;
		AtomVector ctermStem;
		for (uint n = 0; n < stem1.size();n++){
			ctermStem.push_back(fragDB[i+n]);
		}

		
		AtomVector ntermStem;
		for (uint n = 0; n < stem2.size();n++){
			ntermStem.push_back(fragDB[i+stem1.size()+numResidues+n]);
		}
		//cout << "Check gap"<<endl;

		if ( abs(ctermStem[1]->getResidueNumber() - ntermStem[0]->getResidueNumber()) != numResidues+1){
			fprintf(stdout," Gap found at residue %1s %4d %3s and %1s %4d %3s\n",
				ctermStem[0]->getChainId().c_str(),
				ctermStem[0]->getResidueNumber(),
				ctermStem[0]->getResidueName().c_str(),
				ntermStem[0]->getChainId().c_str(),
				ntermStem[0]->getResidueNumber(),
				ntermStem[0]->getResidueName().c_str());
			continue;
		}


		//cout << "Check distances"<<endl;
		bool passDistanceFilter = true;
		int index = 0;
		for (uint c = 0; c < ctermStem.size();c++){
			/*
				fprintf(stdout," Checking residue %1s %4d %3s\n",
					ctermStem[c]->getChainId().c_str(),
					ctermStem[c]->getResidueNumber(),
					ctermStem[c]->getResidueName().c_str());
			*/
			for (uint n = 0; n < ntermStem.size();n++){


				double distSq = ctermStem[c]->distance2(*ntermStem[n]);

				
				/*
					fprintf(stdout," \t%1s %4d %3s = %8.3f\n",
					ntermStem[n]->getChainId().c_str(),
					ntermStem[n]->getResidueNumber(),
					ntermStem[n]->getResidueName().c_str(),distSq);
				*/
				if (abs(stemDistanceSq[index++] - distSq) > tol){
					passDistanceFilter = false;
					break;
				}
			}

			if (!passDistanceFilter){
				break;
			}
		}

		// Continue if the distance filter not passed
		if (!passDistanceFilter) {
			continue;
		}

		//cout << "alignment.."<<endl;
		// align and print winning fragment
		AtomVector fragStem = ctermStem + ntermStem;
		AtomVector fragAll = ctermStem;
		for (int f = 0; f < numResidues+1;f++){
			fragAll.push_back(fragDB[i+ctermStem.size()+f]);
		} 		

		fragAll += ntermStem;

		fragAll.saveCoor("pre");

		Transforms tm;
		tm.align(fragStem,stems,fragAll);
		double rmsd = fragStem.rmsd(stems);

		fprintf(stdout,"%1s %3d %1s %3d  %8.3f\n",
				ctermStem[0]->getChainId().c_str(),
				ctermStem[0]->getResidueNumber(),
			        ntermStem[0]->getChainId().c_str(),
				ntermStem[0]->getResidueNumber(),
			        rmsd);



		System sys;
		sys.addAtoms(fragAll);
		bbqT.fillInMissingBBAtoms(sys.getChain(0));
		
		//cout << "write out: "<<endl<<fragAll.toString()<<endl;
		stringstream ss;
		PDBWriter pout;
		pout.open(ss);
		//pout.write(fragAll);
		pout.write(sys.getAtoms());
		pout.close();

		matchingFragments.push(pair<double,string>(rmsd,ss.str()));

		fragAll.applySavedCoor("pre");

		
	}

	fprintf(stdout, "Best RMSD: %8.3f\n",matchingFragments.top().first);
	return Py_BuildValue("s",matchingFragments.top().second.c_str());

} 

static PyObject* 
quickQuench(PyObject *self, PyObject *args) {


	char *str;
	int rounds;
	if (!PyArg_ParseTuple(args,"si",&str,&rounds))
		  return NULL;

	string pdb = (string)str;
	PDBReader rin;
	rin.read(pdb);
	rin.close();

	System initSys;
	initSys.addAtoms(rin.getAtoms());

	Quench quencher;
	quencher.setVariableNumberRotamers(10,100);
	System finalSys = quencher.runQuench(initSys,rounds);

	stringstream ss;
	PDBWriter pout;
	pout.open(ss);
	pout.write(finalSys.getAtoms());
	pout.close();

	return Py_BuildValue("s",ss.str().c_str());
} 

static PyObject* 
printHello(PyObject *self, PyObject *args) {


	char *str;
	if (!PyArg_ParseTuple(args,"s",&str))
		  return NULL;

	fprintf(stdout, "Hi: %s\n",str);

	return Py_BuildValue("s","");
} 



static char python_msl_doc[] = " commonMSL interface to python ";

static PyMethodDef msl_methods[] = {
	{"printHello", printHello, METH_VARARGS, python_msl_doc},
	{"getChi", getChi, METH_VARARGS, python_msl_doc},
	{"searchForFragments", searchForFragments, METH_VARARGS, python_msl_doc},
	{"getFragments",getFragments,METH_VARARGS,python_msl_doc},
	{"localSampling",localSampling,METH_VARARGS,python_msl_doc},
	{"quickQuench",quickQuench,METH_VARARGS,python_msl_doc},
        {NULL, NULL}
};

PyMODINIT_FUNC
initPythonMSL(void)
{
        Py_InitModule("PythonMSL", msl_methods);

}
