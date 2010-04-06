/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
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



// MSL Includes
#include "System.h"
#include "ChiStatistics.h"
#include "MslTools.h"
#include "OptionParser.h"
#include "PhiPsiReader.h"
#include "PhiPsiStatistics.h"
#include "getDihedrals.h"

// STL Includes
#include<iostream>

// R Includes
#ifdef __R__
   #include "RInside.h"
#endif


using namespace std;
using namespace MSL;




int main(int argc, char *argv[]){
    
    // Option Parser
    Options opt = setupOptions(argc,argv);

    // Read PDB
    System sys;
    sys.readPdb(opt.pdb);

    
    PhiPsiStatistics pps;

    if (opt.phiPsiTable != ""){
	    PhiPsiReader ppr(opt.phiPsiTable);
	    ppr.open();
	    ppr.read();
	    ppr.close();
	    pps = ppr.getPhiPsiStatistics();
    }


    string filename = MslTools::getFileName(opt.pdb);
    ChiStatistics chi;


    fprintf(stdout, "FILE CHAIN RESNUM RESICODE RESNAME PHI PSI ");
    if (opt.phiPsiTable != ""){
	    fprintf(stdout,"PP-COUNTS PP-PROB PP-PROBALL PP-PROP ");
    }
    fprintf(stdout, "CHI1 CHI2 CHI3 CHI4\n");
    vector<double> phiAngles;
    vector<double> psiAngles;

    for (uint i = 0 ; i < sys.positionSize();i++){

	    Residue & n   = sys.getResidue(i);

	    // Remove non-amino acids... bad way to do this, but should work.
	    string oneLetter = MslTools::getOneLetterCode(n.getResidueName());
	    if (oneLetter == "X") continue;

	    fprintf(stdout, "%s %1s %3d%1s %3s ",filename.c_str(),n.getChainId().c_str(),n.getResidueNumber(),n.getResidueIcode().c_str(),n.getResidueName().c_str());



	    double phi           = 0.0;
	    double psi    	 = 0.0;
	    int counts    	 = 0;
	    double prob   	 = 0.0;   
	    double probAll	 = 0.0;
	    double prop          = 0.0; 

	    if (i > 0 && (i < sys.positionSize()-1  &&  
			  MslTools::getOneLetterCode(sys.getResidue(i+1).getResidueName()) != "X"  && 
			  MslTools::getOneLetterCode(sys.getResidue(i-1).getResidueName()) != "X")){
		    Residue & nm1 = sys.getResidue(i-1);
		    Residue & np1 = sys.getResidue(i+1);
		    if (nm1.getChainId() == n.getChainId() && np1.getChainId() == n.getChainId()){
			    phi     = PhiPsiStatistics::getPhi(nm1,n);
			    psi     = PhiPsiStatistics::getPsi(n,np1);

			    if (opt.phiPsiTable != ""){
				    counts  = pps.getCounts(nm1,n,np1);
				    prob    = pps.getProbability(nm1,n,np1);
				    probAll = pps.getProbabilityAll(nm1,n,np1);
				    prop    = pps.getPropensity(nm1,n,np1);
			    }
		    }
		    
	    }



	    fprintf(stdout, "%8.3f %8.3f ",phi,psi);
	    phiAngles.push_back(phi);
	    psiAngles.push_back(psi);

	    if (opt.phiPsiTable != ""){
		    fprintf(stdout,"%5d  %5.2f  %5.2f  %5.2f", counts, prob, probAll, prop);
	    }

	    if (chi.getNumberChis(n) == -1) {
		    fprintf(stdout, "\n");
		    continue;
	    }
	    for (uint c = 0; c < chi.getNumberChis(n);c++){
		    if (!(chi.atomsExist(n,c+1))) {
			    fprintf(stdout, " ---- MISSING ATOMS ---- ");
			    break;
		    }

		    double angle = chi.getChi(n,c+1);
		    if (angle != MslTools::doubleMax){
			    fprintf(stdout,"%8.2f ",angle);
		    }
	    }
	    fprintf(stdout,"\n");
    }


    
#ifdef __R__
    // Start instance of R
    RInside R;

    // Pass stl vectors phiAngles,psiAngles to R as phi,psi
    R.assign(phiAngles, "phi");
    R.assign(psiAngles, "psi");

    string plotStr = "color=densCols(cbind(phi,psi));plot(phi,psi,col=color,pch=20,cex=1.5);";

    // Boundary for strict alpha-helix -90° ≤ phi ≤ -42°; -70° ≤ psi ≤ -15°; -125° ≤ phi+psi ≤ -77°
    plotStr       += "segments(x0=-90,y0=-15,x1=-90,y1=-35,col=\"red\",lty=2,lwd=2);";
    plotStr       += "segments(x0=-90,y0=-35,x1=-55,y1=-70,col=\"red\",lty=2,lwd=2);";
    plotStr       += "segments(x0=-55,y0=-70,x1=-42,y1=-70,col=\"red\",lty=2,lwd=2);";
    plotStr       += "segments(x0=-42,y0=-70,x1=-42,y1=-35,col=\"red\",lty=2,lwd=2);";
    plotStr       += "segments(x0=-42,y0=-35,x1=-62,y1=-15,col=\"red\",lty=2,lwd=2);";
    plotStr       += "segments(x0=-62,y0=-15,x1=-90,y1=-15,col=\"red\",lty=2,lwd=2);";

    // Boundary for loose alpha-helix -90° ≤ phi ≤ -35°; -70° ≤ psi ≤ 0°
    plotStr       += "segments(x0=-90,y0=-70,x1=-35,y1=-70,col=\"yellow\",lty=2,lwd=2);";
    plotStr       += "segments(x0=-35,y0=-70,x1=-35,y1=0,col=\"yellow\",lty=2,lwd=2);";
    plotStr       += "segments(x0=-35,y0=0,x1=-90,y1=0,col=\"yellow\",lty=2,lwd=2);";
    plotStr       += "segments(x0=-90,y0=0,x1=-90,y1=-70,col=\"yellow\",lty=2,lwd=2);";

    // Boundary for loose beta sheet -180 <= phi <= -30 ; 60 <= psi <=180 (used 175 to see boundary on plot)
    plotStr       += "segments(x0=-175,y0=60,x1=-30,y1=60,col=\"orange\",lty=2,lwd=2);";
    plotStr       += "segments(x0=-30,y0=60,x1=-30,y1=175,col=\"orange\",lty=2,lwd=2);";
    plotStr       += "segments(x0=-30,y0=175,x1=-175,y1=175,col=\"orange\",lty=2,lwd=2);";
    plotStr       += "segments(x0=-175,y0=175,x1=-175,y1=60,col=\"orange\",lty=2,lwd=2);";

    // Boundary for left-handed alpha helix 20 <= phi <= 125 ; 45 <= psi <= 90
    plotStr       += "segments(x0=20,y0=90,x1=125,y1=90,col=\"green\",lty=2,lwd=2);";
    plotStr       += "segments(x0=125,y0=90,x1=125,y1=45,col=\"green\",lty=2,lwd=2);";
    plotStr       += "segments(x0=125,y0=45,x1=20,y1=45,col=\"green\",lty=2,lwd=2);";
    plotStr       += "segments(x0=20,y0=45,x1=20,y1=90,col=\"green\",lty=2,lwd=2);";

    // Plot store using png
    string txt = "png(filename=\"rama.png\");"+plotStr+"dev.off()";
    R.parseEvalQ(txt);

    // Plot strore using svg
    txt = "library(RSvgDevice);svg(filename=\"rama.svg\");"+plotStr+"dev.off()";
    R.parseEvalQ(txt);


#endif

        
}


Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;

	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "getDihedrals --pdb PDB [ --phiPsiTable TABLE --debug ]\n";
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	// This is not implemented yet, but is a good idea (dwkulp 3/28/10)
	opt.selection = OP.getString("selection");
	if (OP.fail()){
		cerr << "WARNING 1111 selection not specified.\n";
	}


	opt.phiPsiTable = OP.getString("phiPsiTable");
	if (OP.fail()){
		//opt.phiPsiTable = "/home/dwkulp/software/mslib/trunk/tables/phiPsiCounts.txt";
		//cerr << "WARNING no phiPsiTable set, using: "<<opt.phiPsiTable<<endl;
		opt.phiPsiTable = "";
	}

	opt.debug = OP.getBool("debug");
	if (OP.fail()){
		opt.debug = false;
	}


	return opt;
}
