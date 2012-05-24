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


// MSL Includes
#include "System.h"
#include "MslTools.h"
#include "OptionParser.h"
#include "AtomSelection.h"
#include "Frame.h"
#include "PDBWriter.h"
#include "Line.h"
#include "Helanal.h"
#include "PyMolVisualization.h"
#include "HelixFit.h"
#include "helixWindowAnalysis.h"

// STL Includes
#include<iostream>
#include<fstream>
using namespace std;

int main(int argc, char *argv[]){
    
    // Option Parser
    Options opt = setupOptions(argc,argv);

    // Read PDB
    System sys;
    sys.readPdb(opt.pdb);
    sys.getAtoms().saveCoor("pre");


    if (sys.size() != 2){
	    cerr << "ERROR 2222 PDB does not have exactly 2 chains\n";
	    exit(2222);
    }


    // Make a selection object for selecting windows on each helix
    AtomSelection sel(sys.getAtoms());

    // Store chain references
    Chain &A = sys.getChain(0);
    Chain &B = sys.getChain(1);

    // Pymol visualization flag...
    PyMolVisualization py;

    double dist1 = A.getResidueByIndex(0).getAtom("CA").distance(B.getResidueByIndex(0).getAtom("CA"));
    double dist2 = A.getResidueByIndex(0).getAtom("CA").distance(B.getResidueByIndex(B.size()-1).getAtom("CA"));


    double dist3 = A.getResidueByIndex(A.size()-1).getAtom("CA").distance(B.getResidueByIndex(B.size()-1).getAtom("CA"));
    double dist4 = A.getResidueByIndex(A.size()-1).getAtom("CA").distance(B.getResidueByIndex(0).getAtom("CA"));

    bool parallelDimer = true;
    if (dist2 < dist1 || dist4 < dist3){
	    parallelDimer = false;
    }

    // Header..
    cout << "PDB CHAIN1 RESSTART1 RESEND1 CHAIN2 RESSTART2 RESEND2 ORI Z AVE CONTACT HELDIST HELCROSS HELPHASEA HELPHASEB MUTDIST WINDOWSZ MINCADIST MINRMSD MINRESA MINRESB HELIXAOFFSET HELIXAAMP HELIXAFREQ HELIXAPHASE HELIXARESIDUALS HELIXBOFFSET HELIXBAMP HELIXBFREQ HELIXBPHASE HELIXBRESIDUALS\n";




    // For each window
    HelixFit hfit;	
    if (opt.debug){
	    hfit.setVerbose();
    }

    for (uint i = 0; i < A.size()-opt.window+1;i++){
	    

	    // Select atoms on chain A
	    stringstream selStr1;
	    selStr1 << "((chain "<<A.getChainId()<<" and resi "<<A.getResidueByIndex(i).getResidueNumber()<<"-"<<A.getResidueByIndex(i).getResidueNumber()+opt.window-1<<") and name CA)";
	    AtomVector sel1 = sel.select(selStr1.str());





	    for (uint j = 0; j < B.size()-opt.window+1;j++){



		    // Select atoms on chain B
		    stringstream selStr2;
		    selStr2 << "((chain "<<B.getChainId()<<" and resi "<<B.getResidueByIndex(j).getResidueNumber()<<"-"<<B.getResidueByIndex(j).getResidueNumber()+opt.window-1<<") and name CA)";
		    AtomVector sel2 = sel.select(selStr2.str());

		    if (sel1.size() != sel2.size()){
			    cerr << "ERROR selection sizes do not match: "<<sel1.size()<<","<<sel2.size()<<endl;
			    cerr << "SEL1:"<<endl<<sel1.toString()<<endl;
			    cerr << "SEL2:"<<endl<<sel2.toString()<<endl;
			    continue;
		    }


		    if (opt.windowAtA != -1 && sel1(0).getResidueNumber() != opt.windowAtA) continue;
		    if (opt.windowAtB != -1 && sel2(0).getResidueNumber() != opt.windowAtB) continue;

		    vector<CartesianPoint *> projectionsOntoBundleAxis;
		    
		    Helanal h;
		    CartesianPoint helanalAxisA(0.0,0.0,0.0);
		    CartesianPoint helanalCenterA(0.0,0.0,0.0);

		    CartesianPoint helanalAxisB(0.0,0.0,0.0);
		    CartesianPoint helanalCenterB(0.0,0.0,0.0);


		    for (uint n = 0; n < sel1.size()-4;n++){

			    // Compute helanal 
			    h.update(sel1(n).getCoor(),sel1(n+1).getCoor(),sel1(n+2).getCoor(),sel1(n+3).getCoor());
			    helanalAxisA   += h.getAxis();
			    helanalCenterA += h.getCenter();


			    h.update(sel2(n).getCoor(),sel2(n+1).getCoor(),sel2(n+2).getCoor(),sel2(n+3).getCoor());
			    helanalAxisB   += h.getAxis();
			    helanalCenterB += h.getCenter();

		    }

		    // Compute Helanal Axis
		    CartesianPoint centerA = helanalCenterA /  (sel1.size()-4);
		    CartesianPoint centerB = helanalCenterB /  (sel1.size()-4);
		    CartesianPoint dir1 = helanalAxisA.getUnit();
		    CartesianPoint dir2 = helanalAxisB.getUnit();

		    Line helixA(centerA, dir1);
		    Line helixB(centerB, dir2);


		    // Retain a direction on helix B that is parallel to helix A.
		    CartesianPoint dir3 = dir2;
		    if (!parallelDimer){
			    dir3 *= -1;
		    }
		    Line helixBpar(centerB,dir3);

		    if (opt.pymol){
			    py.createArrow(centerA,dir1,"HelA");
			    py.createArrow(centerB,dir2,"HelB");
		    }
		    
		    // Bundle Axis creation
		    CartesianPoint bundleAxis     = (dir1 + dir3) / 2;
		    CartesianPoint bundleMidpoint = (centerA + centerB ) / 2;
		    
		    Line bundleZ(bundleMidpoint,bundleAxis);

		    if (opt.pymol){
			    py.createArrow(bundleMidpoint,bundleAxis,"BundleZ");
		    }

		    // Compute helanal paramters
		    vector<double> paramsHelanal = getGeometricParameters(helixA,helixB,bundleZ,sel1((int)sel1.size()/2).getCoor(),sel2((int)sel2.size()/2).getCoor());


		    // Get mutual perpendicular for helixA,helixB lines (helanal lines);
		    CartesianPoint mutA = helixA.pointOfMinDistanceToLine(helixBpar);
		    CartesianPoint mutB = helixBpar.pointOfMinDistanceToLine(helixA);
		    CartesianPoint mutMid = (mutA+mutB) / 2;

		    Line toMut(bundleMidpoint, mutMid - bundleMidpoint);


		    if (opt.pymol){
			    CartesianPoint c = toMut.getCenter();
			    CartesianPoint d = toMut.getDirection();
			    py.createArrow(c,d,"toMut");
		    }

		    // Project each point onto line, find max distance to bundleMidpoint
		    double maxDist = -MslTools::doubleMax;
		    double minMidpointDist =  MslTools::doubleMax;
		    Atom *closestCA = NULL;
		    char tmpWinA[80];
		    sprintf(tmpWinA,"/tmp/%s_windowA_%04d_%04d",MslTools::getFileName(opt.pdb).c_str(),sel1(0).getResidueNumber(),sel2(0).getResidueNumber());

		    ofstream foutA;
		    foutA.open(tmpWinA);

		    char tmpWinB[80];
		    sprintf(tmpWinB,"/tmp/%s_windowB_%04d_%04d",MslTools::getFileName(opt.pdb).c_str(),sel1(0).getResidueNumber(),sel2(0).getResidueNumber());
		    ofstream foutB;
		    foutB.open(tmpWinB);

		    
		    for (uint n = 0; n < sel1.size();n++){
			    int index2 = n;
			    if (!parallelDimer){
				    index2 = sel1.size() - (n+1);
			    }

			    CartesianPoint proj = toMut.projection(sel1(n).getCoor());
			    Line tmp(bundleMidpoint,(proj-bundleMidpoint));
			    double dist = proj.distance(bundleMidpoint);

			    if (opt.pymol){
				    CartesianPoint c = tmp.getCenter();
				    CartesianPoint d = tmp.getDirection();
				    char gg[12];
				    sprintf(gg,"projA_%03d",n);
				    py.createArrow(c,d,gg);
			    }


			    if (abs(tmp.getDirection().angle(toMut.getDirection())-180.0) > 1){
				    if (dist > maxDist){
					    maxDist = dist;
				    }

			    }
		    

			    proj = toMut.projection(sel2(index2).getCoor());
			    tmp.setCenter(bundleMidpoint);
			    tmp.setDirection(proj-bundleMidpoint);
			    dist = proj.distance(bundleMidpoint);

			    if (opt.pymol){
				    CartesianPoint c = tmp.getCenter();
				    CartesianPoint d = tmp.getDirection();
				    char gg[12];
				    sprintf(gg,"projB_%03d",n);
				    py.createArrow(c,d,gg);
			    }


			    if (abs(tmp.getDirection().angle(toMut.getDirection())-180.0) > 1){
				    if (dist > maxDist){
					    maxDist = dist;
				    }

			    } 


			    // Get distance to bundleAxis..
			    double d1 = bundleZ.projection(sel1(n).getCoor()).distance(sel1(n).getCoor());
			    double d2 = bundleZ.projection(sel2(index2).getCoor()).distance(sel2(index2).getCoor());

			    foutA << d1<<endl;
			    foutB << d2<<endl;


			    double midDist = sel1(n).getCoor().distance(bundleMidpoint);
			    if (midDist < minMidpointDist){
				    minMidpointDist = midDist;
				    closestCA = sel1[n];
			    }

			    midDist = sel2(index2).getCoor().distance(bundleMidpoint);
			    if (midDist < minMidpointDist){
				    minMidpointDist = midDist;
				    closestCA = sel2[index2];
			    }
			    
		    }
		    foutA.close();
		    foutB.close();
		    
		    stringstream cmd;
		    cmd << "Rscript "<<opt.rscript<<" --data "<<tmpWinA<<" --size "<<opt.window<<" > "<<tmpWinA<<".out";
		    int error = system(cmd.str().c_str());
		    vector<string> tmpLines;
		    MslTools::readTextFile(tmpLines,((string)tmpWinA)+".out");
		    string helixAsineParams = tmpLines[0];

		    cmd.str("");
		    cmd << "Rscript "<<opt.rscript<<" --data "<<tmpWinB<<" --size "<<opt.window<<" > "<<tmpWinB<<".out";
		    error = system(cmd.str().c_str());
		    tmpLines.clear();
		    MslTools::readTextFile(tmpLines,((string)tmpWinB)+".out");
		    string helixBsineParams = tmpLines[0];

		    if (!opt.debug){
			    cmd.str("");
			    
			    cmd <<"rm "<<tmpWinA<<" "<<tmpWinB<<" "<<tmpWinA<<".out "<<tmpWinB<<".out";
			    error = system(cmd.str().c_str());
		    } else {
			    exit(0);
		    }




		    // Now compute Z-offset, by transforming bundle and taking differences in Z-coordinate
		    Transforms trans;
		    /*
		    CartesianPoint z(0.0,0.0,1.0);
		    CartesianPoint pt  = (bundleMidpoint+bundleAxis);
		    trans.align(sys.getAtoms(),pt, z+bundleMidpoint,bundleMidpoint);		    
		    sys.getAtoms().translate(bundleMidpoint*-1);
		    */

		    if (closestCA == NULL){
			    cerr << "ERROR No closest CA? "<<i<<j<<endl;
			    exit(0);
		    }
		    CartesianPoint transInClosestCADir = closestCA->getCoor() - bundleZ.projection(closestCA->getCoor());
		    transInClosestCADir.getUnit();
		    bundleZ.setDirection(bundleZ.getDirection().getUnit());
		    Line x(bundleMidpoint,transInClosestCADir);
		    
		    Frame f;
		    f.computeFrameFrom2Lines(bundleZ,x);

		    ofstream frameOut;
		    frameOut.open("/tmp/frame.py");
		    frameOut << f.toString()<<endl;
		    frameOut.close();

		    f.transformToGlobalBasis(sys.getAtoms());

		    // Get average z-coordinate diference
		    double aveZDiff       = 0.0;
		    double aveDist        = 0.0;
		    int numCloseContacts  = 0;  // any 2 c-alphas are within 9 Angstroms

		    // Get any 2 C-alphas that are the closest
		    double minCalphaDist        = MslTools::doubleMax;
		    for (uint n = 0; n < sel1.size();n++){
			    int index2 = n;
			    if (!parallelDimer){
				    index2 = sel1.size() - (n+1);
			    }
			    aveZDiff += abs(sel1(n).getCoor().getZ() - sel2(index2).getCoor().getZ());
			    
			    
			    double dist = sel1(n).distance(sel2(index2));

			    aveDist  += dist;


			    for (uint n2 = 0; n2 < sel2.size();n2++){

				    double dist = sel1(n).distance(sel2(n2));

				    // Increment close contact counter
				    if ( dist < 9){
					    numCloseContacts++;
				    }

				    if (dist < minCalphaDist){
					    minCalphaDist = dist;
				    }
			    }

		    }

		    // Get averages
		    aveZDiff /= sel1.size();
		    aveDist  /= sel1.size();


		    // Fit Helical parameters 
		    hfit.setAtoms(sel1);
		    hfit.setStepSize(1.0);
		    hfit.setNumberSteps(1000);
 		    bool convergedFit = hfit.fit(HelixFit::NELDERMEAD1);
		    vector<double> paramsFitA = hfit.getFittedParameters();
		    if (!convergedFit){

			    hfit.setParameters(paramsFitA);
			    hfit.setNumberSteps(10000);
			    hfit.setStepSize(1.0);
			    hfit.fit(HelixFit::STEEPEST_DESCENT);
			    paramsFitA = hfit.getFittedParameters();
		    }

		    if (opt.debug){
			    AtomVector &test = hfit.fittedHelix(paramsFitA);
    
			    char fitA[80];
			    sprintf(fitA,"/tmp/fitA-%05d.pdb",i);
			    PDBWriter fitout;
			    fitout.open(fitA);
			    fitout.write(test);
			    fitout.close();
		    }

		    hfit.setAtoms(sel2);
		    hfit.setStepSize(1.0);
		    hfit.setNumberSteps(1000);
 		    convergedFit = hfit.fit(HelixFit::NELDERMEAD1);
		    
		    vector<double> paramsFitB = hfit.getFittedParameters();

		    if (!convergedFit){
			    hfit.setParameters(paramsFitB);
			    hfit.setNumberSteps(10000);
			    hfit.setStepSize(1.0);
			    hfit.fit(HelixFit::STEEPEST_DESCENT);
			    paramsFitB = hfit.getFittedParameters();
		    }

		    if (opt.debug){
			    AtomVector &test = hfit.fittedHelix(paramsFitB);
    
			    char fitB[80];
			    sprintf(fitB,"/tmp/fitB-%05d.pdb",i);
			    PDBWriter fitout;
			    fitout.open(fitB);
			    fitout.write(test);
			    fitout.close();
		    }

		    
		    // Get rmsd from A onto B. (best segment of 10 residues?)
		    if (!parallelDimer){
			    std::reverse(sel2.begin(),sel2.end());
		    }
		    double minRMSD = MslTools::doubleMax;
		    int minResidueA = 0;
		    int minResidueB = 0;
		    for (uint n = 0; n < sel1.size()-9;n++){


			    for (uint n2 = 0; n2 < sel2.size()-9;n2++){

				    AtomVector align1;
				    for (uint t = 0; t < 10;t++){
					    align1.push_back(new Atom(sel1(n+t)));
					    align1.push_back(new Atom(sel2(n2+t)));
				    }

				    AtomVector align2;
				    for (uint t = 0; t < 10;t++){
					    align2.push_back(new Atom(sel2(n2+t)));
					    align2.push_back(new Atom(sel1(n+t)));
				    }    


				    
				    trans.align(align1,align2);
				    double rmsd = align1.rmsd(align2);

				    if (rmsd < minRMSD){
					    minRMSD = rmsd;
					    minResidueA = sel1(n).getResidueNumber();
					    minResidueB = sel2(n2).getResidueNumber();

					    /*
					    PDBWriter pdbout;
					    pdbout.open("/tmp/min1.pdb");
					    pdbout.write(align1);
					    pdbout.close();
					    pdbout.open("/tmp/min2.pdb");
					    pdbout.write(align2);
					    pdbout.close();
					    sys.writePdb("/tmp/zalign.pdb");
					    */
					    
				    }
				    

				    align1.deletePointers();
				    align1.clear();

				    align2.deletePointers();
				    align2.clear();

			    }



		    }

		    if (!parallelDimer){
			    std::reverse(sel2.begin(),sel2.end());
		    }

		    string ori = "P";
		    if (!parallelDimer){
			    ori = "A";
		    }

		    // Finally print something out...
		    fprintf(stdout, "%-30s %1s %3d %3d %1s %3d %3d %1s %8.3f %8.3f %3d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %4d %4d %s %s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
			             MslTools::getFileName(opt.pdb).c_str(),
			             A.getChainId().c_str(),sel1(0).getResidueNumber(),sel1(sel1.size()-1).getResidueNumber(),
                                     B.getChainId().c_str(),sel2(0).getResidueNumber(),sel2(sel2.size()-1).getResidueNumber(),
			             ori.c_str(),aveZDiff,aveDist,numCloseContacts,paramsHelanal[0],paramsHelanal[1],paramsHelanal[2],paramsHelanal[3],
			    (bundleMidpoint.distance(mutMid)-maxDist),
			    bundleMidpoint.getZ(),
			    minCalphaDist,
			    minRMSD,
			    minResidueA,
			    minResidueB,
			    helixAsineParams.c_str(),
			    helixBsineParams.c_str(),

			    paramsFitA[0],
			    paramsFitA[1],
			    paramsFitA[2],
			    paramsFitA[3],
			    paramsFitA[4],
			    paramsFitA[5],
			    paramsFitA[6],
			    paramsFitA[7],
			    paramsFitA[8],
			    paramsFitA[9],

			    paramsFitB[0],
			    paramsFitB[1],
			    paramsFitB[2],
			    paramsFitB[3],
			    paramsFitB[4],
			    paramsFitB[5],
			    paramsFitB[6],
			    paramsFitB[7],
			    paramsFitB[8],
			    paramsFitB[9]
			    );

		    

		    if (opt.pymol){
			    py.createAtom(bundleMidpoint,"bundleMid");
			    py.createAtom(mutMid,"mutMid");
			    ofstream mout;
			    mout.open("tmp.py");
			    mout << py.toString();
			    mout.close();

		    }

		    if (opt.debug){
			    sys.writePdb("/tmp/foo.pdb");
			    exit(0);
		    }
		    
		    // Revert atoms back to original positionx
		    sys.getAtoms().applySavedCoor("pre");

	    }
    }

    if (opt.pymol){
	    ofstream fout;
	    fout.open("/tmp/helanal.py");

	    fout << py<<endl;
	    fout.close();
    }
}

vector<double> getGeometricParameters(Line &_axis1, Line &_axis2, Line &_axisCenter, CartesianPoint &_phasePoint1, CartesianPoint &_phasePoint2){
	
	vector<double> results;


	results.push_back(_axis1.getCenter().distance(_axis2.getCenter()));
	results.push_back(_axis1.segmentDihedral(_axis2));

	CartesianPoint projH = _axis1.projection(_phasePoint1);
	CartesianPoint projC = _axisCenter.projection(_phasePoint1);

	results.push_back(projC.angle(projH,_phasePoint1));

	projH = _axis2.projection(_phasePoint2);
	projC = _axisCenter.projection(_phasePoint2);

	results.push_back(projC.angle(projH,_phasePoint2));

	return results;

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
		cout << "helixPairAnalysis\n";
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.sel = OP.getMultiString("sel");
	if (OP.fail()){
	} else {

		if (opt.sel.size() != 2){
			cerr << "ERROR  1111 2 selections required, one per helix\n";
			exit(1111);
		}
	}

	opt.window = OP.getInt("window");
	if (OP.fail()){
		cerr << "WARNING window not specified, default = 10.\n";
		opt.window = 10;
	}

	opt.pymol = OP.getBool("pymol");
	if (OP.fail()){
		opt.pymol = false;
	}
	opt.debug = OP.getBool("debug");
	if (OP.fail()){
		opt.debug = false;
	}

	opt.rscript = OP.getString("sinFit");
	if (OP.fail()){
		opt.rscript = "/export/home/dwkulp/tmp/franken/opm101309/bin/sinFit.r";
		cerr << "WARNING RSCRIPT not set using: "<<opt.rscript<<endl;
	}

	opt.windowAtA = OP.getInt("windowAtA");
	if (OP.fail()){
		opt.windowAtA = -1;
	}

	opt.windowAtB = OP.getInt("windowAtB");
	if (OP.fail()){
		opt.windowAtB = -1;
	}
	return opt;
}
