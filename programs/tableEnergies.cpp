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

#include <string>
#include <sstream>
#include <vector>
#include <ostream>
#include <fstream>
#include <cmath>



#include "tableEnergies.h"
#include "Quench.h"
#include "OptionParser.h"
#include "PDBWriter.h"
#include "System.h"

using namespace std;

using namespace MSL;



int main(int argc, char *argv[]){

	
	// Option Parser
	Options opt = setupOptions(argc,argv);

	TwoBodyDistanceDependentPotentialTable tbd;
	tbd.readPotentialTable(opt.potfile);

	System initialSystem;
	initialSystem.readPdb(opt.pdb);

	OnTheFlyManager otfm(&initialSystem, opt.parfile);
	double energy = otfm.calculateTotalEnergy(initialSystem, tbd);

	cout << "Total: " << energy << endl;
	if (opt.deltaG) {
		vector<Chain *> theChains = initialSystem.getChains(); 
		double chainEnergies = 0.;
		for (uint i = 0; i < theChains.size(); i++) {
			System thisChain(*theChains[i]);
			double thisEnergy = otfm.calculateTotalEnergy(thisChain, tbd);
			chainEnergies += thisEnergy;
			cout << "Chain " << theChains[i]->getChainId() << ": " << thisEnergy << endl;
		}
		energy = energy - chainEnergies;
		cout << "Delta G energy: " << energy << endl;
	}

	return 0;
}


Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;
	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);
	OP.readArgv(theArgc, theArgv);


	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "tableEnergies --pdb PDB [--potfile KBFILE --parfile PARFILE --deltaG]\n";
		cout << "For deltaG option, energy is calculated for total, then individual chain energies is subtracted out" << endl;
		exit(0);
	}

	opt.pdb  = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.potfile = OP.getString("potfile");
	if (OP.fail()){
		cerr << "WARNING no potfile specified, using default /export/home/jedonald/membraneKB/dpher_beta_0.1.equated.msl.pmf\n";
		opt.potfile = "/export/home/jedonald/membraneKB/dpher_beta_0.1.equated.msl.pmf";
	}

	opt.parfile = OP.getString("parfile");
	if (OP.fail()){
		cerr << "WARNING no parfile specified, using default /library/charmmTopPar/par_all22_prot.inp\n";
		opt.parfile = "/library/charmmTopPar/par_all22_prot.inp";
	}

	opt.deltaG = OP.getBool("deltaG");

	return opt;
}
