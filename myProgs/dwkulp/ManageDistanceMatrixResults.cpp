#include "ManageDistanceMatrixResults.h"
#include "Transforms.h"
#include "PDBReader.h"
#include "PDBWriter.h"
#include "PolymerSequence.h"
#include "AtomSelection.h"

using namespace std;
using namespace MSL;

ManageDistanceMatrixResults::ManageDistanceMatrixResults(){
	alignPdbs = false;
	rmsdTol   = 2.0;
}

ManageDistanceMatrixResults::~ManageDistanceMatrixResults(){

}

void ManageDistanceMatrixResults::printResults(){

	cout << "Number of results to print: "<<allResults.size()<<endl;
 
	int numberAligned = 0;
	ofstream seqOut1;
	ofstream seqOut2;
	if (alignPdbs){
		seqOut1.open("seqAlignment1.txt", ios::out|ios::app);
		seqOut2.open("seqAlignment2.txt", ios::out|ios::app);
	}

	for (int i=0; i<allResults.size(); i++){


		//retrieve distance matrices
		DistanceMatrix dm1 = allResults[i][0].getDistanceMatrix1();
		DistanceMatrix dm2 = allResults[i][0].getDistanceMatrix2();
	
		//retrieve lists of matrix windows
		//print the name of PDBs we're comparing
		string PDBName1 = getFileName(dm1.getPDBid());
		string PDBName1Short = PDBName1.substr(0,17);

		string PDBName2 = getFileName(dm2.getPDBid());
		string PDBName2Short = PDBName2.substr(0,17);

		cout<<"Comparing PDBs: "<<PDBName1Short<<", "<<PDBName2Short<<endl;

	
		//for two fixed DM's (dm1, dm2), access and print the multiCompare results
		for(int j=0; j<allResults[i].size(); j++){

			//retrieve winning MWs and score from the Result
			DistanceMatrixResult & currentResult = allResults[i][j];
	    
			double minLikeness = currentResult.getLikeness(); 
			// compare function failed. Nothing to print--result has null windows. 
			if(minLikeness==1000000){
				cout<<"Badness, badness everywhere."<<endl;
				continue; //the compareFunction didn't find any matches
			} 

			MatrixWindow & minWindow1 = currentResult.getMatrixWindow1();
			MatrixWindow & minWindow2 = currentResult.getMatrixWindow2();
	    
			//print information
			int i1 = (minWindow1).getLeftR();
			int j1 = (minWindow1).getLeftC();
			int i2 = (minWindow2).getLeftR();
			int j2 = (minWindow2).getLeftC();
     
			string i1ID = dm1.getAtomVector()[i1]->getSegID().c_str();
			string j1ID = dm1.getAtomVector()[j1]->getSegID().c_str();
			string i2ID = dm2.getAtomVector()[i2]->getSegID().c_str();
			string j2ID = dm2.getAtomVector()[j2]->getSegID().c_str();
    
			if(i1ID=="" || j1ID=="" || i2ID=="" ||j2ID==""){
   
				i1ID = dm1.getAtomVector()[i1]->getChainId().c_str();
				j1ID = dm1.getAtomVector()[j1]->getChainId().c_str();
				i2ID = dm2.getAtomVector()[i2]->getChainId().c_str();
				j2ID = dm2.getAtomVector()[j2]->getChainId().c_str();
			}//end if

			int i1res = dm1.getAtomVector()[i1]->getResidueNumber();
			int j1res = dm1.getAtomVector()[j1]->getResidueNumber();
			int i2res = dm2.getAtomVector()[i2]->getResidueNumber();
			int j2res = dm2.getAtomVector()[j2]->getResidueNumber();


			fprintf(stdout, "\t\tWindow1 %3d,%3d (Residues: %1s %3d, %1s %3d)\tWindow2 %3d,%3d (Residues: %1s %3d, %1s %3d)\t%8.3f\n", i1, j1, i1ID.c_str(), i1res, j1ID.c_str(), j1res, i2, j2, i2ID.c_str(), i2res, j2ID.c_str(), j2res, minLikeness);


			if (alignPdbs){

			
				// Get AtomVectors of matching windows
				AtomPointerVector ca1 = minWindow1.getSmallAVec();
				AtomPointerVector ca2 = minWindow2.getSmallAVec();
//for(int i=0; i<ca2.size();i++) cout<<*ca2[i]<<endl;

				// Read PDB1 in 
				PDBReader r1;
				r1.open(dm1.getPDBid());
				r1.read();
				r1.close();
				AtomPointerVector &ref = r1.getAtomPointers();
				if (ref.size() == 0){
					cerr << "ERRROR 3453 in ManageDistanceMatrixResults::printResults while aligning pdbs, ref pdb not found?: "<<dm1.getPDBid()<<endl;
					exit(3453);
				}


				

				PDBReader r2;
				r2.open(dm2.getPDBid());
				r2.read();
				r2.close();
				AtomPointerVector &vec = r2.getAtomPointers();



				Transforms t;		    
//for(int i=0; i<vec.size();i++) cout<<*vec[i]<<endl;
//cout<<endl;
//				cout << ca2.size() << "\t" << ca1.size() << "\t" << vec.size() << endl;
				bool result = t.rmsdAlignment(ca2,ca1, vec);
// This align is to a get vec (full pdb coordinates)
//ca2,ca1 did not changed; vec changed
				if (!result){
					cerr << "Alignment has failed!"<<endl;
					exit(1211);
				}

				
				ca2.saveCoor("pre");
// here record old ca2 coordinates

				result = t.rmsdAlignment(ca2,ca1);
// this align is to get a new ca2 to calculate rmsd
				if (!result){
					cerr << "Alignment has failed!"<<endl;
					exit(1212);
				}
				
				double rmsd=ca1.rmsd(ca2);
				cout << "RMSD: "<<rmsd<<endl;
				if (rmsd <= rmsdTol){
// a smaller rmsd is to get optimized structures

					numberAligned++;
/*
					// Write out aligned PDB for dm2
					char a[80];
					sprintf(a, "%s.aligned.%5.3f.pdb",MslTools::getFileName(dm2.getPDBid()).c_str(),rmsd);

					PDBWriter w;
					w.open(a);
					w.write(vec);
					w.close();
*/

					// Setup a Polymer Sequence to print out the alignment...

					/*
					   Since chain A in dm1 can match chain A or chain B in dm2 , everything is named chain1 and chain2.
					   
					   Procedure:
					   1. Get chain1 and chain2 from all AtomVectors (ca1,ca2,ref) ; ref = ca1 , but includes all of the residues.
					   2. Create sequences for chain1s
					   3. Create a PolymerSequence for chain1s and print out.
					   .. chain2s

					 */
					AtomSelection selWin1(ca1);
					AtomPointerVector win1Ch1 = selWin1.select("chain "+ca1(0).getChainId()); 
					AtomPointerVector win1Ch2 = selWin1.select("chain "+ca1(ca1.size()/2).getChainId()); 

					AtomSelection selWin2(ca2);
					AtomPointerVector win2Ch1 = selWin2.select("chain "+ca2(0).getChainId()); 
					AtomPointerVector win2Ch2 = selWin2.select("chain "+ca2(ca2.size()/2).getChainId()); 

					AtomSelection selRef(ref);
					AtomPointerVector refCh1 = selRef.select("chain "+ca1(0).getChainId()); 
					AtomPointerVector refCh2 = selRef.select("chain "+ca1(ca2.size()/2).getChainId()); 

                                        // Write out aligned PDB for dm2
                                        char a[80];
                                        sprintf(a, "%s.aligned.%5.3f.pdb",MslTools::getFileName(dm2.getPDBid()).c_str(),rmsd);
//                                        sprintf(a, "SEQ.%s.SEQ.%s.%s.aligned.%5.3f.pdb",refCh1(0).getChainId().c_str(),refCh2(0).getChainId().c_str(),MslTools::getFileName(dm2.getPDBid()).c_str(),rmsd);

                                        PDBWriter w;
                                        w.open(a);
                                        w.write(vec);
//					w.write(ca2);
                                        w.close();

					// Create an input sequence for chain 1
					stringstream seqStr1;
					seqStr1 << win2Ch1(0).getChainId()<<" "<<win2Ch1(0).getResidueNumber()<<": ";
					seqStr1 << PolymerSequence::toThreeLetterCode(win2Ch1);
//seqStr1 << PolymerSequence::toOneLetterCode(win2Ch1); // it does not work as the whole one-letter code and "X" will be printed out as an unrecognized amino acid
//cout<<seqStr1.str()<<endl;
//for(int i=0; i<win2Ch1.size();i++) cout<<win2Ch1[i]->getResidueNumber()<<" "<<*win2Ch1[i]<<endl;
					// Create a polymer sequence for chain A
					PolymerSequence poly1;
					poly1.setName("SEQ"+refCh1(0).getChainId()+"-"+MslTools::getFileName(dm2.getPDBid()));
					poly1.setSequence(seqStr1.str());
poly1.setReferenceSequence(PolymerSequence::toOneLetterCode(refCh1).c_str(),"SEQ"+refCh1(0).getChainId()+"-REF-"+MslTools::getFileName(dm1.getPDBid()).c_str(),1,win1Ch1[0]->getResidueNumber(),win2Ch1[0]->getResidueNumber());
//cout <<"win1Ch1[0]->getResidueNumber(): "<<win1Ch1[0]->getResidueNumber()<<" win2Ch1([0]->getResidueNumber(): "<<win2Ch1[0]->getResidueNumber()<<endl;
//					poly1.setReferenceSequence(MslTools::getOneLetterCode(refCh1),"SEQ"+refCh1(0).getChainId().c_str()+"-REF",refCh1(0).getResidueNumber(),win1Ch1(0).getResidueNumber());


//cout<<poly1<<endl;
/*
vector<vector<vector<string> > > sequences = poly1.getSequence();
for(vector<vector<vector<string> > >::iterator it = sequences.begin(); it != sequences.end(); ++it) {
	for(vector<vector<string> >::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
		for(vector<string>::iterator it3 = it2->begin(); it3 != it2->end(); ++it3) {
			cout << *it3 <<" ";
		}
	}
}
cout <<endl;
*/
					stringstream seqStr2;
					seqStr2 << win2Ch2(0).getChainId()<<" "<<win2Ch2(0).getResidueNumber()<<": ";
					seqStr2 << PolymerSequence::toThreeLetterCode(win2Ch2);
//seqStr2 << PolymerSequence::toOneLetterCode(win2Ch2);

					// Create a polymer sequence for chain A
					PolymerSequence poly2;
					poly2.setName("SEQ"+refCh2(0).getChainId()+"-"+MslTools::getFileName(dm2.getPDBid()));
					poly2.setSequence(seqStr2.str());
//					poly2.setReferenceSequence(MslTools::getOneLetterCode(refCh2),"SEQ"+refCh2(0).getChainId()+"-REF", refCh2(0).getResidueNumber(),win1Ch2(0).getResidueNumber());
poly2.setReferenceSequence(PolymerSequence::toOneLetterCode(refCh2),"SEQ"+refCh2(0).getChainId()+"-REF-"+MslTools::getFileName(dm1.getPDBid()).c_str(),1,win1Ch2[0]->getResidueNumber(),win2Ch2[0]->getResidueNumber());

					if (numberAligned == 1){
						seqOut1 << poly1.getReferenceHeader();
						seqOut2 << poly2.getReferenceHeader();
					}
					

					seqOut1 << poly1;
					seqOut2 << poly2;


				}

				ca2.applySavedCoor("pre");
		    // ca2 is changed to its former value
			}

		}//end for on j


	}//end for on i

	seqOut1.close();
	seqOut2.close();

}


