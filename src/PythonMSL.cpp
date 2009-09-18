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
#include "PDBFragments.h"
#include "BackRub.h"
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
localSamplingCCD(PyObject *self, PyObject *args) {

	int numFragments;
	int maxDegree;
	char *system;
	char *fragment;
	char *bbqTable;

	if (!PyArg_ParseTuple(args,"ssiis",&system,&fragment,&numFragments,&maxDegree,&bbqTable))
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
localSamplingPDB(PyObject *self, PyObject *args) {

	double rmsdTol;
	int numResidues;
	int numFragments;
	char *chain,*fragment,*database;

	if (!PyArg_ParseTuple(args,"sssidi",&chain,&fragment,&database,&numResidues,&rmsdTol,&numFragments))
		  return NULL;


	PDBFragments fragDB(database);
	fragDB.loadFragmentDatabase();

	string pdb = (string)fragment;
	PDBReader rin;
	rin.read(pdb);
	rin.close();

	AtomVector &ats = rin.getAtoms();
	if (ats.size() < 4){
		fprintf(stdout, "Fragment defines less than 4 atoms\n");
		return Py_BuildValue("s","");
	}

	for (uint i = 0; i < ats.size();i++){
		if (ats(i).getName() != "CA"){
			fprintf(stdout, "Fragment atom%d is not a CA -> %s\n",i,ats(0).toString().c_str());
			return Py_BuildValue("s","");
		}
	}

	string c = (string)chain;
	PDBReader rinSys;
	rinSys.read(c);
	rinSys.close();

	System sys;
	sys.addAtoms(rinSys.getAtoms());
	
	if (sys.size() != 1){
		fprintf(stdout,"System has %d chains, only give 1 chain\n",sys.size());
		return Py_BuildValue("s","");
	}

	vector<int> stems;
	stems.push_back( sys.getPositionIndex(ats(0).getChainId(),ats(0).getResidueNumber()));
	stems.push_back( sys.getPositionIndex(ats(1).getChainId(),ats(1).getResidueNumber()));
	stems.push_back( sys.getPositionIndex(ats(ats.size()-2).getChainId(),ats(ats.size()-2).getResidueNumber()));
	stems.push_back( sys.getPositionIndex(ats(ats.size()-1).getChainId(),ats(ats.size()-1).getResidueNumber()));

	int numMatchingFrags = fragDB.searchForMatchingFragments(sys.getChain(0),stems);

	if (rmsdTol == 0.0){


		stringstream ss;

		fprintf(stdout, "Number of Matching Fragments: %d\n",numMatchingFrags);
		
		// print out first 'numFragments'
		System &frags       = fragDB.getLastSearchResults();
		AtomVector &fragAts = frags.getAtoms();
		for (uint i = 0; i < numMatchingFrags;i++){

			if (i > numFragments){
				break;
			}

			for (uint a = 0; a < fragAts.size();a++){
				fragAts(a).setActiveConformation(i);
			}


			// Add SC here...
			
			ss << "MODEL"<<endl;
			PDBWriter pout;
			pout.open(ss);
			pout.write(fragAts);
			pout.close();
			ss << "ENDMDL"<<endl;

			/*
			char fname[100];
			sprintf(fname,"/tmp/frag-%06d.pdb",i);
			pout.open(fname);
			pout.write(fragAts);
			pout.close();		
			*/

		}

		return Py_BuildValue("s",ss.str().c_str());

	} else {
		// Go through each match, check all bb atom rmsd...
	}

	return Py_BuildValue("s","");

} 

static PyObject* 
localSamplingBR(PyObject *self, PyObject *args) {

	int numFragments;
	char *bbqTable;
	char *chain,*fragment,*database;

	if (!PyArg_ParseTuple(args,"ssi",&chain, &fragment,&numFragments))
		  return NULL;

	// Convert PDB to string
	string pdb = (string)chain;
	
	// Read in PDB to MSL
	PDBReader rin;
	rin.read(pdb);
	rin.close();

	System sys;
	sys.addAtoms(rin.getAtoms());

	string fragStr = (string)fragment;
	PDBReader rinFrag;
	rinFrag.read(fragStr);
	rinFrag.close();

	System frag;
	frag.addAtoms(rinFrag.getAtoms());

	int startIndex = sys.getPositionIndex(frag.getResidue(0).getChainId(),frag.getResidue(0).getResidueNumber());
	int endIndex   = sys.getPositionIndex(frag.getResidue(frag.residueSize()-1).getChainId(),frag.getResidue(frag.residueSize()-1).getResidueNumber());

	BackRub br;
	string brPdbs = br.localSample(sys.getChain(0),startIndex, endIndex, numFragments);

	return Py_BuildValue("s",brPdbs.c_str());
}

static PyObject* 
localSamplingMIN(PyObject *self, PyObject *args) {
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
	{"localSamplingPDB",localSamplingPDB, METH_VARARGS, python_msl_doc},
	{"localSamplingCCD",localSamplingCCD, METH_VARARGS,python_msl_doc},
	{"localSamplingBR" ,localSamplingBR,  METH_VARARGS,python_msl_doc},
	{"localSamplingMIN",localSamplingMIN, METH_VARARGS,python_msl_doc},
	{"quickQuench",quickQuench,METH_VARARGS,python_msl_doc},
        {NULL, NULL}
};

PyMODINIT_FUNC
initPythonMSL(void)
{
        Py_InitModule("PythonMSL", msl_methods);

}
