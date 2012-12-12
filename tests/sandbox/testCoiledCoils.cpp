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

#include "CoiledCoils.h"
#include "CoiledCoilFitter.h"
#include "PDBWriter.h"
#include "PDBReader.h"
#include "Transforms.h"
#include "OptionParser.h"
#include "MslOut.h"
#include <cstdlib>
#include <string>

using namespace std;
using namespace MSL;

string programName = "testCoiledCoils";
string programDescription = "This test generates Coiled Coils based on a given set of parameters, then attempts to blindly fit them";
string programAuthor = "Dan Kulp";
string programVersion = "1.0.0";
string programDate = "19 September 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;


/******************************************
 *  
 *         =======  MAIN  =======
 *
 ******************************************/
int main(int argc, char *argv[]) {


        // Define Objects for the test
	CoiledCoils cc;
	CoiledCoilFitter ccf;
	PDBWriter pout;

	/******************************************************************************
	 *
	 *                   === GEVORG/CRICKS CC GENERATION ===
	 *
	 ******************************************************************************/





	//sys.writePdb(opt.fileName);
	cc.getCoiledCoilCricks(4.910, -4.027, -13.530, 2.267, 102.806, 149.984, 0.0, 26);
	AtomPointerVector cricks = cc.getAtomPointers();
	pout.open("/tmp/cricks.pdb");
	pout.write(cricks);
	pout.close();

	double shr = 4.910;
	double risePerRes = 1.51;
	double shp = 128.206;
	double ahr = 2.26;
	double ahp = 102.8;
	double ahphase = 149.984;
	double dZ  = 0.0; //meaningless parameter at the moment.
	cc.getCoiledCoil(shr, risePerRes, shp, ahr, ahp, ahphase, dZ, 26);
	AtomPointerVector norths = cc.getAtomPointers();
	pout.open("/tmp/norths.pdb");
	pout.write(norths);
	pout.close();

	Symmetry sym;
	sym.applyCN(norths,2);
	pout.open("/tmp/northC2.pdb");
	pout.write(sym.getAtomPointers());
	pout.close();

	System sys;
	sys.readPdb("/tmp/northC2.pdb");

	// Fitting procedure..

	// Add helix chain A
	ccf.addNextHelix(&sys.getChain("A").getAtomPointers());

	// Add helix chain B
	ccf.addNextHelix(&sys.getChain("B").getAtomPointers());	

	// Set symmetry
	ccf.setSymmetry("C");

	// Do fittin procedure
	ccf.fit();

	// Do something else..
	vector<double> params = ccf.getMinimizedParameters();

	if ( abs(params[0] - shr) > 0.1){
	  cerr << "ERROR Super-helical radius is off. Target = "<<shr<<" MinValue = "<<params[0]<<endl;
	} else if (abs(params[1] - risePerRes) > 0.1){
	  cerr << "ERROR RisePerRes is off. Target = "<<risePerRes<<" MinValue = "<<params[1]<<endl;
	} else if (abs(params[2] - shp) > 0.1){
	  cerr << "ERROR Super-helical pitch is off. Target = "<<shp<<" MinValue = "<<params[2]<<endl;
	} else if (abs(params[3] - ahr) > 0.1){
	  cerr << "ERROR Alpha-helical radius is off. Target = "<<ahr<<" MinValue = "<<params[3]<<endl;
	} else if (abs(params[4] - ahp) > 0.1){
	  cerr << "ERROR Alpha-helical pitch is off. Target = "<<ahp<<" MinValue = "<<params[4]<<endl;
	} else if (abs(params[5] - ahphase) > 0.1){
	  cerr << "ERROR Alpha-helical phase is off. Target = "<<ahphase<<" MinValue = "<<params[5]<<endl;
	} else if (abs(params[6] - dZ) > 0.1){
	  cerr << "ERROR deltaZ offset is off. Target = "<<dZ<<" MinValue = "<<params[6]<<endl;
	} else {
	  cout << "LEAD";
	}



}



