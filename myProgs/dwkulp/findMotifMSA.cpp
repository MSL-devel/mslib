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
#include "MslTools.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "FastaReader.h"
#include "SasaCalculator.h"
#include "findMotifMSA.h"
#include "MslOut.h"

// STL Includes
#include<iostream>
#include<vector>
#include<map>

using namespace std;

using namespace MSL;


// MslOut 
static MslOut MSLOUT("findMotifMSA");


int main(int argc, char *argv[]) {	


	// Option Parser
	Options opt = setupOptions(argc,argv);
	
	System pdb;
	pdb.readPdb(opt.pdb);

	// Find exposed positions

     /*
	SASA reference:
	Protein Engineering vol.15 no.8 pp.659–667, 2002
	Quantifying the accessible surface area of protein residues in their local environment
	Uttamkumar Samanta Ranjit P.Bahadur and  Pinak Chakrabarti
      */
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

     // Identify Exposed positions
     SasaCalculator scalc(pdb.getAtomPointers());
     scalc.calcSasa();

     FastaReader fin(opt.msa);
     fin.open();
     if (!fin.read()){
       cerr << "ERROR reading "<<opt.msa<<endl;
       exit(143);
     }
     fin.close();

     string pdbSeq = fin.getSequence(MslTools::getFileName(opt.pdb));
     if (pdbSeq == ""){
       cerr << "ERROR couldn't find "<<MslTools::getFileName(opt.pdb)<<" as a name in : "<<opt.msa<<endl;
       exit(2342);
     }

     // Store sequence position = structural position
     map<int,int> pdbSeqToStructIndex;

     int nonDashPdbCounter = 0;
     for (uint i = 0; i < pdbSeq.size();i++){
         // Skip over '-' in pdbSeq
         if (pdbSeq[i] == '-')  continue;

	 pdbSeqToStructIndex[i] = nonDashPdbCounter;
	 nonDashPdbCounter++;
     }

     map<string,string> seqs = fin.getSequences();
     map<string,string>::iterator it;
     map<int,vector<pair<int,int> > > validatedMatches; // what is what?
     for (it = seqs.begin();it != seqs.end();it++){
       
       // Get non-blank sequences  (non-blank index, FASTA-original index)
       map<int,int> blanks;
       string removedBlanks ="";
       int noBlankIndex = 0;
       for (uint i = 0; i < it->second.size();i++){
	 if (it->second[i] != '-'){
	   removedBlanks += it->second[i];
	   blanks[noBlankIndex] = i;
	   noBlankIndex++;
	 }
       }

       // Search for regex
       vector<pair<int,int> > matches;
       MslTools::regex(removedBlanks,opt.regex,matches);
       
       // Correlate to structural position and store info
       for (uint i = 0; i < matches.size();i++){

	 bool validMatch = true;
	 for (uint j = matches[i].first; j <= matches[i].second;j++){
	   if (pdbSeq[ blanks[j] ] == '-'){
	     validMatch = false;
	     break;
	   }
	 }

	 if (!validMatch) continue;

	 // We have a valid matches[i] .. meaning no '-' in equivlanet PDB sequence
	 validatedMatches[blanks[matches[i].first]].push_back(pair<int,int>(blanks[matches[i].first],blanks[matches[i].second]));
	 
       }
     }


     // Tabulate answers
     map<int,vector<pair<int,int> > >::iterator vit;
     for (vit = validatedMatches.begin();vit != validatedMatches.end();vit++){
       int index1 = (vit->second)[0].first;
       int index2 = (vit->second)[0].second;
       string positionId1 = pdb.getPosition(pdbSeqToStructIndex[index1]).getPositionId();
       string positionId2 = pdb.getPosition(pdbSeqToStructIndex[index2]).getPositionId();

       double normSasa = scalc.getResidueSasa(positionId1) / refSasa[MslTools::getOneLetterCode(pdb.getPosition(positionId1).getResidueName())];
       if (normSasa > 1.0){
	 normSasa = 1.0;
       }

       
       double percent = ((double)vit->second.size()) / ( (double)seqs.size()) * 100;
       fprintf(stdout, "%-30s %-45s %20s %8u %8u %8.2f %8d %10s %10s %8.3f %8.3f\n", 
	       MslTools::getFileName(opt.pdb).c_str(),
	       MslTools::getFileName(opt.msa).c_str(),
	       opt.regex.c_str(),
	       vit->second.size(),
	       seqs.size(),
	       percent,
	       vit->first,
	       positionId1.c_str(),
	       positionId2.c_str(),
	       normSasa,
	       scalc.getResidueSasa(positionId1));
	       
     }

	
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
		cout << "findMotifMSA --pdb foo.pdb --msa foo.msa --regex N[^P][ST] \n";
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}


	opt.msa = OP.getString("msa");
	if (OP.fail()){
		cerr << "ERROR 1111 msa not specified.\n";
		exit(1111);
	}

	opt.regex = OP.getString("regex");
	if (OP.fail()){
		cerr << "ERROR 1111 regex not specified.\n";
		exit(1111);
	}

	return opt;
}



