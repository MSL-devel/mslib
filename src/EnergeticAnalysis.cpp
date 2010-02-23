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
#include "EnergeticAnalysis.h"
#include "AtomicPairwiseEnergy.h"

using namespace MSL;
using namespace std;


EnergeticAnalysis::EnergeticAnalysis(){
}
EnergeticAnalysis::~EnergeticAnalysis(){
}


void EnergeticAnalysis::analyzePosition(System &_sys, int _position, int _rotamer){
	

	// Compute energy by term and by grouping for this position
	AtomicPairwiseEnergy ape("/library/charmmTopPar/par_all27_prot_lipid_extra.inp");
	
	// Calculate Residue-Residue energies

	//Set active rotamer
	Position &pos1 = _sys.getPosition(_position);
	pos1.setActiveRotamer(_rotamer);

	// Get atoms of active rotamer
	AtomPointerVector &atoms1 = pos1.getAtoms();


	// Iterate over each other position, computing energies
	vector<map<string,double> > allResidueEnergies(_sys.residueSize(),map<string,double>());
	double maxTotal = 0.0;
	double minTotal = 0.0;
	int minTotalIndex = 0;
	int maxTotalIndex = 0;
	PyMolVisualization pymolViz;
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

		//cout << "\t TEMPLATE POSITION IS "<<i<<" "<<pos.getResidueNumber()<<" "<<pos.getChainId()<<" "<<pos.getResidueName()<<endl;
		allResidueEnergies[i] = ape.calculatePairwiseEnergy(_sys, atoms1, pos.getAtoms());

		if (allResidueEnergies[i]["TOTAL"] > maxTotal){
			maxTotal = allResidueEnergies[i]["TOTAL"];
			maxTotalIndex = i;
		}

		if (allResidueEnergies[i]["TOTAL"] < minTotal){
			minTotal = allResidueEnergies[i]["TOTAL"];
			minTotalIndex = i;
		}


		// Generate a VDW thingy.
		CartesianPoint centroid1 =pos1.getCurrentIdentity().getCentroid();
		CartesianPoint centroid2 =pos.getCurrentIdentity().getCentroid();
		CartesianPoint vect = centroid2 - centroid1;
		string usedInPymol = "";
		if (centroid1.distance(centroid2) < 8 && abs(allResidueEnergies[i]["CHARMM_VDW"]) > 1.0 ){
			usedInPymol = "*";
			vect = vect.getUnit();

			// # angstroms from centroid toward interacting residue.
			double angstromsFromCentroid = 4.0;
			CartesianPoint startVect =  centroid1 + (vect * angstromsFromCentroid); 
			CartesianPoint endVect   =  centroid1 + (vect * (angstromsFromCentroid+0.5));
		
			char name[80];
			sprintf(name,"VDW_%s%04d%3s_%s%04d%3s", 
				pos1.getChainId().c_str(), pos1.getResidueNumber(),pos1.getResidueName().c_str(), 
				 pos.getChainId().c_str(),  pos.getResidueNumber(), pos.getResidueName().c_str());
			
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

		if (centroid1.distance(centroid2) < 8 && abs(allResidueEnergies[i]["CHARMM_ELEC"]) > 5.0 ){

			usedInPymol += "+";
			vect = vect.getUnit();

			// # angstroms from centroid toward interacting residue.
			double angstromsFromCentroid = 5;
			CartesianPoint startVect =  centroid1 + (vect * angstromsFromCentroid); 
			CartesianPoint endVect   =  centroid1 + (vect * (angstromsFromCentroid+0.5));
		
			char name[80];
			sprintf(name,"ELE_%s%04d%3s_%s%04d%3s", 
				pos1.getChainId().c_str(), pos1.getResidueNumber(),pos1.getResidueName().c_str(), 
				 pos.getChainId().c_str(),  pos.getResidueNumber(), pos.getResidueName().c_str());
			
			vector<double> blue;
			blue.push_back(0.0);
			blue.push_back(0.0);
			blue.push_back(1.0);

			vector<double> red;
			red.push_back(1.0);
			red.push_back(0.0);
			red.push_back(0.0);

			vector<double> rgb = MslTools::getRGB(blue,red,-10, 10,allResidueEnergies[i]["CHARMM_ELEC"]);
			pymolViz.createCone(startVect, endVect, name,1.0,rgb[0],rgb[1],rgb[2]);

		}

		fprintf(stdout, "%1s %04d %3s  %8.3f  %8.3f %1s\n",	
			pos.getChainId().c_str(),  pos.getResidueNumber(), pos.getResidueName().c_str(),
			allResidueEnergies[i]["CHARMM_VDW"],allResidueEnergies[i]["CHARMM_ELEC"],usedInPymol.c_str());


	}

	char name[80];
	sprintf(name,"%s_%04d_%3s.py",pos1.getChainId().c_str(),pos1.getResidueNumber(),pos1.getResidueName().c_str());

	ofstream fout(name);
	fout << pymolViz;
	fout <<endl;
	fout.close();

}
