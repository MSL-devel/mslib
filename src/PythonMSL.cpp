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
#include "SasaCalculator.h"

using namespace std;

using namespace MSL;


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
	sys.addAtoms(rin.getAtomPointers());
	for (uint i = 0 ; i < sys.positionSize();i++){
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
	AtomPointerVector &loop = rin.getAtomPointers();

	// Make CCD object with BBQTable .. (BBQTable is for other backbone atoms)
	CCD sampleCCD((string)bbqTable);

	// Do local sampling inside CCD object
	sampleCCD.localSample(loop,numFragments,maxDegree);
	string pdbresult = sampleCCD.getNMRString();

	// Return PDB string of all atoms
	return Py_BuildValue("s",pdbresult.c_str());
	
}


static PyObject* 
localSamplingPDB(PyObject *self, PyObject *args) {

	double rmsdTol;
	int numResidues;
	int numFragments;
	char *chain,*fragment,*database,*bbqTable;

	if (!PyArg_ParseTuple(args,"sssidis",&chain,&fragment,&database,&numResidues,&rmsdTol,&numFragments,&bbqTable))
		  return NULL;


	PDBFragments fragDB(database,bbqTable);
	fragDB.loadFragmentDatabase();

	string pdb = (string)fragment;
	PDBReader rin;
	rin.read(pdb);
	rin.close();

	AtomPointerVector &ats = rin.getAtomPointers();
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
	sys.addAtoms(rinSys.getAtomPointers());
	
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
		AtomPointerVector &fragAts = frags.getAtomPointers();
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
	sys.addAtoms(rin.getAtomPointers());

	string fragStr = (string)fragment;
	PDBReader rinFrag;
	rinFrag.read(fragStr);
	rinFrag.close();

	System frag;
	frag.addAtoms(rinFrag.getAtomPointers());

	int startIndex = sys.getPositionIndex(frag.getResidue(0).getChainId(),frag.getResidue(0).getResidueNumber());
	int endIndex   = sys.getPositionIndex(frag.getResidue(frag.positionSize()-1).getChainId(),frag.getResidue(frag.positionSize()-1).getResidueNumber());

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
	initSys.addAtoms(rin.getAtomPointers());

	Quench quencher;
	quencher.setVariableNumberRotamers(10,100);
	System finalSys = quencher.runQuench(initSys,rounds);

	stringstream ss;
	PDBWriter pout;
	pout.open(ss);
	pout.write(finalSys.getAtomPointers());
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


static PyObject* 
getSasa(PyObject *self, PyObject *args) {

	char *str;
	int refSasaByRes;
	double probeRadius;
	if (!PyArg_ParseTuple(args,"sid",&str,&refSasaByRes,&probeRadius))
		  return NULL;

	string pdb = (string)str;
	PDBReader rin;
	rin.read(pdb);
	rin.close();

	System sys;
	sys.addAtoms(rin.getAtomPointers());


	SasaCalculator sas(sys.getAtomPointers());
	sas.setTempFactorWithSasa(true);
	sas.setProbeRadius(probeRadius);
	sas.calcSasa();
	sas.printSasaTable(false);
	cout << endl;

	/*  
	    Ref. Sasa
	    
Protein Engineering vol.15 no.8 pp.659–667, 2002
Quantifying the accessible surface area of protein residues in their local environment
Uttamkumar Samanta Ranjit P.Bahadur and  Pinak Chakrabarti
	*/
	if (refSasaByRes){
	  map<string,double> refSasa;
	  refSasa["G"] = 83.91;
	  refSasa["A"] = 116.40;
	  refSasa["S"] = 125.68;
	  refSasa["C"] = 141.48;
	  refSasa["P"] = 144.80;
	  refSasa["T"] = 148.06;
	  refSasa["D"] = 155.37;
	  refSasa["V"] = 162.24;
	  refSasa["N"] = 168.87;
	  refSasa["E"] = 187.16;
	  refSasa["Q"] = 189.17;
	  refSasa["I"] = 189.95;
	  refSasa["L"] = 197.99;
	  refSasa["H"] = 198.51;
	  refSasa["K"] = 207.49;
	  refSasa["M"] = 210.55;
	  refSasa["F"] = 223.29;
	  refSasa["Y"] = 238.30;
	  refSasa["R"] = 249.26;
	  refSasa["W"] = 265.42;
	
	  fprintf(stdout, "Normalized SASA:\n");
	  for (uint i = 0; i < sys.residueSize();i++){
	  
	    Residue &res = sys.getResidue(i);
	    double sasa = 0.0;
	    for (uint j = 0; j < res.size();j++){
	      sasa += res.getAtom(j).getTempFactor();
	    }

	    if (fabs(sasa) < 0.01) {
	      sasa = 0.0;
	    } else {
	      sasa = sasa / refSasa[MslTools::getOneLetterCode(res.getResidueName())];
	    }

	    if (sasa > 1.0){
	      sasa = 1.0;
	    }



	    fprintf(stdout, "%1s %3d %4.1f\n",res.getChainId().c_str(),res.getResidueNumber(),sasa);
	    for (uint j = 0; j < res.size();j++){
	      res.getAtom(j).setTempFactor(sasa);
	    }

	  }
	}


	stringstream ss;
	PDBWriter pout;
	pout.open(ss);
	pout.write(sys.getAtomPointers());
	pout.close();


	return Py_BuildValue("s",ss.str().c_str());
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
	{"getSasa",getSasa,METH_VARARGS,python_msl_doc},
        {NULL, NULL}
};

PyMODINIT_FUNC
initPythonMSL(void)
{
        Py_InitModule("PythonMSL", msl_methods);

}
