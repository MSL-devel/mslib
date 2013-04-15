/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2013 The MSL Developer Group (see README.TXT)
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

#include "EnergeticAnalysis.h"
#include "CharmmEnergyCalculator.h"

using namespace MSL;
using namespace std;


EnergeticAnalysis::EnergeticAnalysis(){
	paramFile = "";
	pymolOutput = true;
}
EnergeticAnalysis::~EnergeticAnalysis(){
}

void   EnergeticAnalysis::setParameterFile(string _paramFile) { 
	paramFile = _paramFile;
}
string EnergeticAnalysis::getParameterFile() { 
	return paramFile;
}

void EnergeticAnalysis::setPymolOutput(bool usePymolOutput) {
	pymolOutput = usePymolOutput;
}

bool EnergeticAnalysis::getPymolOutput() {
	return pymolOutput;
}

void EnergeticAnalysis::analyzePosition(System &_sys, int _position, int _rotamer){

	double pymolDistCutoff = 8.;
	double pymolVDWCutoff = 1.;
	double pymolElecCutoff = 5.;
	
	if (paramFile == "" || !MslTools::fileExists(paramFile)){
		cerr << "ERROR 4545 EnergeticAnalysis::analyzePosition(), paramFile does not exist: '"<<paramFile<<"'"<<endl;
		exit(4545);
	}

	// Compute energy by term and by grouping for this position
        CharmmEnergyCalculator calculator(paramFile);
	
	// Calculate Residue-Residue energies

	//Set active rotamer
	Position &pos1 = _sys.getPosition(_position);
	pos1.setActiveRotamer(_rotamer);

	// Get atoms of active rotamer
	AtomPointerVector &atoms1 = pos1.getAtomPointers();

	// Iterate over each other position, computing energies
	vector<map<string,double> > allResidueEnergies(_sys.positionSize(),map<string,double>());
	PyMolVisualization pymolViz;
	CartesianPoint centroid1 =pos1.getCurrentIdentity().getCentroid();
	for (uint i = 0; i < _sys.positionSize();i++){
		Position &pos = _sys.getPosition(i);

		// Skip Variable Position Side chains
		if (pos.getTotalNumberOfRotamers() > 1){
			continue;
		}

		// Skip self..
		if (i == _position) {
			continue;
		}

		allResidueEnergies[i] = calculator.calculatePairwiseNonBondedEnergy(atoms1, pos.getAtomPointers());

		string usedInPymol = "";
		if (pymolOutput) {
			// Generate a VDW energy
			CartesianPoint centroid2 =pos.getCurrentIdentity().getCentroid();
			CartesianPoint vect = centroid2 - centroid1;
			if (centroid1.distance(centroid2) < pymolDistCutoff && fabs(allResidueEnergies[i]["CHARMM_VDW"]) > pymolVDWCutoff){
				char name[80];
				sprintf(name,"VDW_%s%04d%3s_%s%04d%3s", 
					pos1.getChainId().c_str(), pos1.getResidueNumber(),pos1.getResidueName().c_str(), 
					 pos.getChainId().c_str(),  pos.getResidueNumber(), pos.getResidueName().c_str());
				
				vect = vect.getUnit();
				usedInPymol = "*";
		
				// # angstroms from centroid toward interacting residue.
				double angstromsFromCentroid = 4.0;
				CartesianPoint startVect =  centroid1 + (vect * angstromsFromCentroid); 
				CartesianPoint endVect   =  centroid1 + (vect * (angstromsFromCentroid+0.5));
			
				vector<double> blue;
				blue.push_back(0.0);
				blue.push_back(0.0);
				blue.push_back(1.0);
		
				vector<double> red;
				red.push_back(1.0);
				red.push_back(0.0);
				red.push_back(0.0);
		
				vector<double> rgb = MslTools::getRGB(blue,red,-5, 5,allResidueEnergies[i]["CHARMM_VDW"]);
				pymolViz.createCylinder(startVect, endVect, name,1.0,rgb[0],rgb[1],rgb[2]);
			}

			if (centroid1.distance(centroid2) < pymolDistCutoff && fabs(allResidueEnergies[i]["CHARMM_ELEC"]) > pymolElecCutoff){
				char name[80];
				sprintf(name,"ELE_%s%04d%3s_%s%04d%3s", 
					pos1.getChainId().c_str(), pos1.getResidueNumber(),pos1.getResidueName().c_str(), 
					 pos.getChainId().c_str(),  pos.getResidueNumber(), pos.getResidueName().c_str());
				
				vect = vect.getUnit();
				usedInPymol += "+";

				// # angstroms from centroid toward interacting residue.
				double angstromsFromCentroid = 5.0;
				CartesianPoint startVect =  centroid1 + (vect * angstromsFromCentroid); 
				CartesianPoint endVect   =  centroid1 + (vect * (angstromsFromCentroid+0.5));
			
				vector<double> blue;
				blue.push_back(0.0);
				blue.push_back(0.0);
				blue.push_back(1.0);

				vector<double> red;
				red.push_back(1.0);
				red.push_back(0.0);
				red.push_back(0.0);

				vector<double> rgb = MslTools::getRGB(blue,red,-10, 10,allResidueEnergies[i]["CHARMM_ELEC"]);
				pymolViz.createCone(startVect, endVect, name,1.0,int(rgb[0]),int(rgb[1]),int(rgb[2]));
			}
		}

		// Allow an energy cutoff like this where we can set energy cutoff.  Also change output so it looks better with this
		if (allResidueEnergies[i]["TOTAL"] > 10) {
			fprintf(stdout, "%1s %04d %3s  %8.3f  %8.3f %1s\n",	
				pos.getChainId().c_str(),  pos.getResidueNumber(), pos.getResidueName().c_str(),
				allResidueEnergies[i]["CHARMM_VDW"],allResidueEnergies[i]["CHARMM_ELEC"],usedInPymol.c_str());
		}
	}

	if (pymolOutput) {
		char name[80];
		sprintf(name,"%s_%04d_%3s.py",pos1.getChainId().c_str(),pos1.getResidueNumber(),pos1.getResidueName().c_str());
		ofstream fout(name);
		fout << pymolViz;
		fout <<endl;
		fout.close();
	}
}
