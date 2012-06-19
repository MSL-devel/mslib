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



// MSL Includes
#include "System.h"
#include "ChiStatistics.h"
#include "MslTools.h"
#include "OptionParser.h"
#include "PhiPsiReader.h"
#include "PhiPsiWriter.h"
#include "PhiPsiStatistics.h"
#include "SurfaceAreaAndVolume.h"
#include "SasaCalculator.h"
#include "getDihedrals.h"

// STL Includes
#include<iostream>
#include<map>
#include<string>
#include<vector>

// R Includes
#ifdef __R__
   #include "RInside.h"
#endif


using namespace std;
using namespace MSL;

class BBDEPRotamer {
	public:
	BBDEPRotamer() {
	}


	void setChiVals(vector<double>& _chis) {
		chiVals.resize(_chis.size());
		for(int i = 0; i < _chis.size(); i++) {
			chiVals[i]= _chis[i];
		}
	}

	vector<double>& getChiVals() {
		return chiVals;
	}

	void setProb(double _prob) {
		prob = _prob;
	}

	double getProb() {
		return prob;
	}

	double msd(BBDEPRotamer* _other) {
		vector<double>& otherChis = _other->getChiVals();
		return msd(otherChis);
	}

	double msd(vector<double>& _otherChis) {
		double msd = 0.0;
		for(int i = 0; i < chiVals.size(); i++) {
			double diff = chiVals[i]-_otherChis[i];
			// we need to adjust 
			if(diff > 180) {
				diff = diff - 360;
			} else if (diff < -180) {
				diff = diff + 360;
			}
			msd += 	diff * diff;
		}
		return msd/double(chiVals.size());
	}

	private:

	double prob;
	vector<double> chiVals;
	/* Not needed for now
	unsigned count;
	vector<double> chiSigs;
	vector<double> rotBins;
	*/
};


class BBDEPResidue {
	public:
	BBDEPResidue() {
	}

	~BBDEPResidue() {
		for(map<int,map<int,vector<BBDEPRotamer*> > >::iterator it1 = rotamers.begin(); it1 != rotamers.end(); it1++) {
			for(map<int,vector<BBDEPRotamer*> >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
				for(vector<BBDEPRotamer*>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {
					delete *it3;
				}
				it2->second.clear();
			}
			it1->second.clear();
		}
	}

	void addRotamer(int _phiBin, int _psiBin, vector<double>& _chiVals,double _prob) {
		BBDEPRotamer* rot = new BBDEPRotamer;
		rotamers[_phiBin][_psiBin].push_back(rot);
		rot->setChiVals(_chiVals);
		rot->setProb(_prob);
	}



	double getProb(int _phiBin, int _psiBin, vector<double>& _chiVals) {
		if(rotamers.find(_phiBin) != rotamers.end()) {
			if(rotamers[_phiBin].find(_psiBin) != rotamers[_phiBin].end()) {
				vector<BBDEPRotamer*> rots =  rotamers[_phiBin][_psiBin];
				double bestMsd = MslTools::doubleMax;
				double prob = -1.0;
				for(int i = 0; i < rots.size(); i++) {
					double msd = rots[i]->msd(_chiVals);	
					if(msd < bestMsd) {
						bestMsd = msd;
						prob = rots[i]->getProb();
					}
				}
				return prob;
			} else {
				//cerr << "ERROR: PsiBin " << _psiBin << " not present in rotamer info file"<< endl;
			}
		} else {
			//cerr << "ERROR: PhiBin " << _phiBin << " not present in rotamer info file"<< endl;
		}
		return -1.0;
	}

	private:
	// phipsibins are from -18 to 17
	map<int,map<int,vector<BBDEPRotamer*> > > rotamers;
};

class BBDEPLibrary {

	public:
	BBDEPLibrary() {
	}

	BBDEPLibrary(string _rotinfofile) {
		read(_rotinfofile);  
	}

	~BBDEPLibrary() {
		for(map<string,BBDEPResidue*>::iterator it = residues.begin(); it != residues.end(); it++) {
			delete it->second;
		}
		residues.clear();
	}

	bool read(string _rotinfofile) {
		ifstream libFile;	
		libFile.open(_rotinfofile.c_str());
		if(!libFile.is_open()) {
			return false;
		}
		/* Format 
		 # T  Phi  Psi  Count    r1 r2 r3 r4 Probabil  chi1Val chi2Val chi3Val chi4Val   chi1Sig chi2Sig chi3Sig chi4Sig
		 ARG  -180 -180    10     1  2  2  1  0.249730    62.5   176.9   176.6    85.7       6.9    11.1    10.5     9.9
		*/
		string line;
		while(libFile.good()) {
			getline(libFile,line);
			line = MslTools::uncomment(line);
			vector<string> toks = MslTools::tokenizeAndTrim(line);
			if(toks.size() < 13) {
				continue;
			}
			// Reading in just the resName, phi psi prob and chiVals
			string resName = toks[0];
			double phi = MslTools::toDouble(toks[1]);
			double psi = MslTools::toDouble(toks[2]);
			double prob = MslTools::toDouble(toks[8]);
			//unsigned count = MslTools::toUnsignedInt(toks[3]);
			//cout << resName << " " << phi << " " << psi << " " << prob << " ";
			vector<double> chiVals;
			for(int i = 9; i < 13; i++) {
				chiVals.push_back(MslTools::toDouble(toks[i]));
				//cout << chiVals.back() << " ";
			}
			//cout << endl;
			if(residues.find(resName) == residues.end()) {
				BBDEPResidue* res = new BBDEPResidue;
				residues[resName] = res;
			}

			residues[resName]->addRotamer(getBackboneBin(phi), getBackboneBin(psi),chiVals,prob);
		}
		return true;
	}

	int getBackboneBin(double _angle) {
		// phipsibins are from -18 to 17
		int bin = 0;
		if(_angle > 0) {
			bin = (_angle + 5)/ 10; 
			if(bin == 18) {
				bin = -18;
			}
		} else {
			bin = (_angle - 5)/ 10; 
		}
		return bin;
	}

	double getProb(string _resName, double _phi, double _psi, vector<double>& _chiVals) {
		if(residues.find(_resName) != residues.end()) {
			// get the phi psi bins - 10*10 degrees
			return residues[_resName]->getProb(getBackboneBin(_phi), getBackboneBin(_psi), _chiVals);
		} else {
			//cerr << "ERROR Residue " << _resName << " not present in rotamer info file" << endl;
		}
		return -1.0;
	}

	private:

	map<string,BBDEPResidue*> residues;
};


int main(int argc, char *argv[]){
	
	// Option Parser
	Options opt = setupOptions(argc,argv);

	// Add pdb to pdblist
	if (opt.pdb != ""){
		opt.pdblist.push_back(opt.pdb);
	}
    
	// Phi-Psi statistics read in from file..
	PhiPsiStatistics pps;
	if (opt.phiPsiTable != ""){
		PhiPsiReader ppr(opt.phiPsiTable);
		ppr.open();
		ppr.read();
		ppr.close();
		pps = ppr.getPhiPsiStatistics();
	}

	// Create new phi-psi statistcs
	PhiPsiStatistics new_pps;
	if (opt.createNewPhiPsiTable != ""){
		new_pps.setGridSize(10);
	}

	// ChiStatistics
	ChiStatistics chi;
	if(opt.dofFile != "") {
		chi.read(opt.dofFile);
	}

	// Read the standard rotamer file if required
	BBDEPLibrary stdRotamers;

	if(opt.rotInfoFile != "") {
		if(!stdRotamers.read(opt.rotInfoFile)) {
			cerr << "Unable to read " << opt.rotInfoFile << endl;
			exit(0);
		}
	}


	for (uint s = 0; s < opt.pdblist.size();s++){

		// Read PDB
		System sys;
		sys.readPdb(opt.pdblist[s]);

		string filename = MslTools::getFileName(opt.pdb);
		// Each chain, PositionId -> dNSASA
		std::map<std::string, std::map<std::string,double > > deltaNormalizedSASA;
		if (opt.computeDeltaNormalizedSASA){
			SasaCalculator complex(sys.getAtomPointers(),1.4,2000);
			complex.calcSasa();
			cout << complex.getSasaTable()<<endl;
			cout << "DONE COMPLEX"<<endl;
		}

		// Compute Position-specific SASA in bound state
		for (uint i = 0; i < sys.chainSize();i++){

			for (uint p = 0; p < sys(i).positionSize();p++){
				double posComplexSASA = 0.0;

				for (uint a = 0; a < sys(i).getPosition(p).identitySize();a++){
					posComplexSASA += sys(i).getPosition(p).getAtom(a).getSasa();
				}
				deltaNormalizedSASA[sys(i).getChainId()][sys(i).getPosition(p).getPositionId()] = posComplexSASA; //(posChainSASA - posComplexSASA); // / refSasa[MslTools::getOneLetterCode(sys(i).getPosition(p).getResidueName())] ;

			}
		}

		// Compute Position-specific SASA in unbound state
		for (uint i = 0; i < sys.chainSize();i++){
			SasaCalculator chain(sys(i).getAtomPointers(),1.4,2000);
			chain.calcSasa();
			cout << chain.getSasaTable()<<endl;

			for (uint p = 0; p < sys(i).positionSize();p++){
				double posSasa = 0.0;
				for (uint a = 0; a < sys(i).getPosition(p).identitySize();a++){
					posSasa += sys(i).getPosition(p).getAtom(a).getSasa();
				}

				deltaNormalizedSASA[sys(i).getChainId()][sys(i).getPosition(p).getPositionId()] = posSasa - (deltaNormalizedSASA[sys(i).getChainId()][sys(i).getPosition(p).getPositionId()]);
	
			}

		}

		fprintf(stdout, "FILE CHAIN RESNUM RESICODE RESNAME PHI PSI ");
		if (opt.phiPsiTable != ""){
			fprintf(stdout,"PP-COUNTS PP-PROB PP-PROBALL PP-PROP ");
		}
		if (opt.rotInfoFile != ""){
			fprintf(stdout,"ROT-PROB ");
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


			double phi       = 0.0;
			double psi    	 = 0.0;
			int counts    	 = 0;
			double prob   	 = 0.0;   
			double probAll	 = 0.0;
			double prop      = 0.0; 
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

			if (opt.computeDeltaNormalizedSASA){
				fprintf(stdout, "%8.3f ",deltaNormalizedSASA[n.getChainId()][n.getPositionId()]);
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

			if (opt.createNewPhiPsiTable != ""){
				if (phi != 0.000 && psi != 0.000){
					cout << "Adding : "<<n.getPositionId()<<" "<<MslTools::stringf("%d",new_pps.getPhiPsiBin(phi))<<" "<<MslTools::stringf("%d",new_pps.getPhiPsiBin(psi))<<endl;
					new_pps.addStatisitics(n.getResidueName(),MslTools::stringf("%d",new_pps.getPhiPsiBin(phi)),MslTools::stringf("%d",new_pps.getPhiPsiBin(psi)),1);
				}
			}

			vector<double> chiVals;
			for (uint c = 0; c < chi.getNumberChis(n);c++){
				if (!(chi.atomsExist(n,c+1))) {
					fprintf(stdout, " ---- MISSING ATOMS ---- ");
					break;
				}

				double angle = chi.getChi(n,c+1);
				chiVals.push_back(angle);
			}
			if(opt.rotInfoFile != "") {
				double prob = stdRotamers.getProb(n.getResidueName(),phi,psi,chiVals);
				if(prob != -1.0) {
					fprintf(stdout,"%8.6f ",prob);
				} else {
					fprintf(stdout,"     UNK");
				}
			}
			for (uint c = 0; c < chiVals.size();c++){
				if (chiVals[c] != MslTools::doubleMax){
					fprintf(stdout,"%8.2f ",chiVals[c]);
				}
			}
			fprintf(stdout,"\n");
	
		}


		if (opt.createNewPhiPsiTable != ""){
			cout << "Writing: "<<opt.createNewPhiPsiTable<<endl;
			PhiPsiWriter phipsiout;
			phipsiout.open(opt.createNewPhiPsiTable);
			phipsiout.write(new_pps);
			phipsiout.close();
		}


#ifdef __R__
		// Start instance of R
		RInside R;

		// Pass stl vectors phiAngles,psiAngles to R as phi,psi
		R.assign(phiAngles, "phi");
		R.assign(psiAngles, "psi");


		string plotStr = "color=densCols(cbind(phi,psi));plot(phi,psi,col=color,pch=20,cex=1.5,main=\""+MslTools::getFileName(opt.pdb)+"\");";

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
		cout << "getDihedrals --pdb PDB [--doffile <DEGREE_OF_FREEDOM_FILE> --rotinfofile <ROTAMER_INFO_FILE> --phiPsiTable TABLE --debug ]\n";
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		opt.pdblist = OP.getStringVector("pdblist");
		if (OP.fail()){
			cerr << "ERROR 1111 pdb or pdblist not specified.\n";
			exit(1111);
		}
	} else {
		opt.pdblist = OP.getStringVector("pdblist");
	}
	opt.dofFile = OP.getString("doffile");
	if (OP.fail()){
		cerr << "WARNING 1111 doffile not specified.\n";
		opt.dofFile = "";
	}
	opt.rotInfoFile = OP.getString("rotinfofile");
	if (OP.fail()){
		cerr << "WARNING 1111 rotinfofile not specified.\n";
		opt.rotInfoFile = "";
	}

	// This is not implemented yet, but is a good idea (dwkulp 3/28/10)
	opt.selection = OP.getString("selection");
	if (OP.fail()){
		cerr << "WARNING 1111 selection not specified.\n";
	}


	opt.computeDeltaNormalizedSASA = OP.getBool("computeDeltaNormalizedSASA");
	if (OP.fail()){	
		opt.computeDeltaNormalizedSASA = false;
	} else {
		opt.charmmParameterFile = OP.getString("charmmPar");
		if (OP.fail()){
			cerr << "WARNING: you must specify a CHARMM parameter file (-charmmPar FILENAME) when computing SASA (for radii) and you did not."<<endl;
			exit(1111);
		}

		/*
		  SASA reference:
		  Protein Engineering vol.15 no.8 pp.659–667, 2002
		  Quantifying the accessible surface area of protein residues in their local environment
		  Uttamkumar Samanta Ranjit P.Bahadur and  Pinak Chakrabarti
		*/

		opt.refSasa["G"] = 83.91;
		opt.refSasa["A"] = 116.40;
		opt.refSasa["S"] = 125.68;
		opt.refSasa["C"] = 141.48;
		opt.refSasa["P"] = 144.80;
		opt.refSasa["T"] = 148.06;
		opt.refSasa["D"] = 155.37;
		opt.refSasa["V"] = 162.24;
		opt.refSasa["N"] = 168.87;
		opt.refSasa["E"] = 187.16;
		opt.refSasa["Q"] = 189.17;
		opt.refSasa["I"] = 189.95;
		opt.refSasa["L"] = 197.99;
		opt.refSasa["H"] = 198.51;
		opt.refSasa["K"] = 207.49;
		opt.refSasa["M"] = 210.55;
		opt.refSasa["F"] = 223.29;
		opt.refSasa["Y"] = 238.30;
		opt.refSasa["R"] = 249.26;
		opt.refSasa["W"] = 265.42;

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

	opt.createNewPhiPsiTable = OP.getString("newPhiPsiTable");
	if (OP.fail()){
		opt.createNewPhiPsiTable = "";
	}

	return opt;
}


/*
  CharmmParameterReader par;
  par.reset();
  par.open(opt.charmmParameterFile);
  par.read();
  par.close();

  SurfaceAreaAndVolume savComplex;
  savComplex.setProbeRadius(1.4); 
  savComplex.addAtomsAndCharmmRadii(sys.getAtomPointers(),par);
  savComplex.computeSurfaceAreaAndVolume();
  cout << "Complex SASA: "<<savComplex.getSurfaceArea()<<endl;

  for (uint i = 0; i < sys.size();i++){
  map<string,double> deltaNSASAforchain;

  SurfaceAreaAndVolume sav;
  sav.setProbeRadius(1.4); 
  sav.addAtomsAndCharmmRadii(sys(i).getAtomPointers(),par);
  sav.computeSurfaceAreaAndVolume();      

  cout << "Chain SASA: "<<sav.getSurfaceArea()<<endl;

  // Store norm-DSASA
  for (uint p = 0; p < sys(i).size();p++){

  double posComplexSASA = 0.0;
  double posChainSASA   = 0.0;
  for (uint a = 0; a < sys(i).getPosition(p).size();a++){
	
  Atom &at = sys(i).getPosition(p).getAtom(a);
  double complexSasa = savComplex.getRadiiSurfaceAreaAndVolume(&at)[1];
  double chainSasa   = sav.getRadiiSurfaceAreaAndVolume(&at)[1];

  posComplexSASA += complexSasa;
  posChainSASA   += chainSasa;
  }

  deltaNSASAforchain[sys(i).getPosition(p).getPositionId()] = posComplexSASA; //(posChainSASA - posComplexSASA); // / refSasa[MslTools::getOneLetterCode(sys(i).getPosition(p).getResidueName())] ;

  }

  deltaNormalizedSASA[sys(i).getChainId()] = deltaNSASAforchain;
  }
  }
    
*/
