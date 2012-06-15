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


#include <vector>
#include <map>
#include <string>
#include <stdlib.h>

#define CA_CA_DISTANCE          ((Real)3.78f)
#define CA_CA_DISTANCE_PROLINE  ((Real)3.2f)
#define CA_CA_TOLERANCE         ((Real)0.16f)

#include "BBQTableReader.h"
#include "BBQTable.h"
#include "PDBReader.h"
#include "PDBWriter.h"
#include "System.h"
#include "Chain.h"
#include "Position.h"
#include "Residue.h"

using namespace std;

using namespace MSL;


void processPdb(string &pdbFileName, string chainOfInterest, map<string, bool> &atomsOfInterest, BBQTable &bbqTable, fstream &outputFilestream, bool writePdbs, bool writeOutputFile);
void getAtomsOfInterest(map<string, bool> &atomsOfInterest);
void cloneResidueVector(vector<Residue *> &resVec, vector<Residue *> &resVecCopy);
void deleteResiduesInVector(vector<Residue *> &resVec);
void calcRMS(vector<Residue *> &resVec, vector<Residue *> &resVecCopy, fstream &outputFilestream);
void writePdb(vector<Residue *> &resVec, string pdbFilename);
bool resHasProperAtoms(Residue &res, map<string, bool> &atomsOfInterest);
void removeUninterestingAtoms(Residue &res, map<string, bool> &atomsOfInterest);
void printDihedralAngles(vector<Residue *> &resVec);
bool passDistanceTest(Residue *res1, Residue *res2);

/**
 * This program will take in a list of pdb files, loop over them, and
 * recreate the N, C, and O atoms for each quad.
 *
 * @param The first argument is the name of the text file with the list (including path) of pdbs.
 * @param The second argument is the name of the BBQ Table.
 */
int main(int argc, char **argv) {
    fstream inputFilestream, outputFilestream;
    BBQTable bbqTable;
    string pdbFileName;
    map<string, bool> atomsOfInterest;
    bool writePdbs = false;
    bool writeOutputFile = true;

    if(argc < 4) {
        cout << "Usage: testBBQ2 <pdb file> <BBQ Table> <Write Pdbs? T/F> <Optional Output File>\n";
        exit(1);
    }
   
    if(argc == 4)
        writeOutputFile = false;

    if( (argv[3][0] == 'T') || (argv[3][0] == 't') )
        writePdbs = true;

    // Parse input arguments and setup for the test.
    if(writeOutputFile)
        outputFilestream.open(argv[4], std::ios::out);
    
    BBQTableReader bbqTableReader(argv[2]);
    bbqTableReader.open();
    bbqTableReader.read( bbqTable );
    bbqTableReader.close();
    //getAtomsOfInterest(atomsOfInterest);

    // Loop over each of the pdbs listed in our input file.
    cout << "Now working on pdbID: " << argv[1] << ".\n";
    pdbFileName = argv[1];
    processPdb(pdbFileName, "", atomsOfInterest, bbqTable, outputFilestream, writePdbs, writeOutputFile);
};

/**
 * This function will open the given pdb file, read in all of the chains, and
 * add quadrilateral information for each of the chains.
 *
 * @param pdbFileName The filename, including path info, of the pdb.
 * @param chainOfInterest Which chain we want to look at for the given pdb.  Specify "" for all chains.
 * @param atomsOfInterest The atoms we are interested in collecting information about.  Usually C, O, and N.
 * @param bbqTable The BBQTable that we are adding our information to.
 * @param writePdbs A boolean indicating whether we should dump output pdbs or not.
 */
void processPdb(string &pdbFileName, string chainOfInterest, map<string, bool> &atomsOfInterest, BBQTable &bbqTable, fstream &outputFilestream, bool writePdbs, bool writeOutputFile) {
    PDBReader pdbReader(pdbFileName);
    System sys;
    vector<Chain *> allChains;

    bool success = pdbReader.open();
    success = success && pdbReader.read();
    if(!success) {
        cout << "Unable to read file: " << pdbFileName << ".\nSorry.  I'm quitting now.\n";
        exit(1);
    }

    sys.addAtoms( pdbReader.getAtomPointers() );
    pdbReader.close();
    if(chainOfInterest == "")
        allChains = sys.getChains();
    else
        allChains.push_back( &(sys.getChain(chainOfInterest)) );

    for (unsigned int chainIndex = 0; chainIndex < allChains.size(); ++chainIndex) {
        Chain *currChain = allChains[chainIndex];
        vector<Position *> posVec = currChain->getPositions();
        vector<Residue *> resVec, resVecCopy;

        for (unsigned int positionIndex = 0; positionIndex < posVec.size(); ++positionIndex) {
            Residue &currRes = posVec[positionIndex]->getCurrentIdentity();
            // We're only interested in testing residues that have CA, N, C, and O.
            if( resHasProperAtoms(currRes, atomsOfInterest) ) {
                removeUninterestingAtoms(currRes, atomsOfInterest);
                resVec.push_back(&currRes);
            }
        }

        cloneResidueVector(resVec, resVecCopy);
        bbqTable.fillInMissingBBAtoms(resVecCopy);

        printDihedralAngles(resVecCopy);

        if(writePdbs) {
            writePdb(resVec, pdbFileName + "_orig.pdb");
            writePdb(resVecCopy, pdbFileName + "_bbq.pdb");
        }

        if(writeOutputFile) {
            outputFilestream << pdbFileName << "\n";
            calcRMS(resVec, resVecCopy, outputFilestream);
        }
        
        deleteResiduesInVector(resVecCopy);
    }
}

/**
 * This function will create a map with a key for each atom type we are interested
 * in, and setting the bool value to true.  Usually we are interested in atoms
 * N, C, and O.
 */
void getAtomsOfInterest(map<string, bool> &atomsOfInterest) {
    atomsOfInterest["N"] = true;
    atomsOfInterest["C"] = true;
    atomsOfInterest["O"] = true;
}

void cloneResidueVector(vector<Residue *> &resVec, vector<Residue *> &resVecCopy) {
    for(vector<Residue *>::iterator currIter = resVec.begin(); currIter != resVec.end(); ++currIter) {
        Residue *currResidue = *currIter;
        Residue *newResidue = new Residue(*currResidue);
        newResidue->removeAllAtoms();

        if(currResidue->atomExists("CA"))
            newResidue->addAtom( currResidue->getAtom("CA") );
        resVecCopy.push_back(newResidue);
    }
}

void deleteResiduesInVector(vector<Residue *> &resVec) {
    for(vector<Residue *>::iterator currIter = resVec.begin(); currIter != resVec.end(); ++currIter) {
        delete *currIter;
    }
}

void calcRMS(vector<Residue *> &resVec, vector<Residue *> &resVecCopy, fstream &outputFilestream) {
    if(resVec.size() != resVecCopy.size()) {
        cout << "Why are these residue vectors different sizes? resVec = " << resVec.size() << " resVecCopy = " << resVecCopy.size() << "\n";
        return;
    }

    for(int index = 0; index < resVec.size(); ++index) {
        Real totalDistance = 0.0f;
        Residue *originalRes = resVec[index];
        Residue *bbqRes = resVecCopy[index];
        AtomPointerVector &bbqAV = bbqRes->getAtomPointers();

        // If we were unable to fill in the atoms for this residue,
        // skip calculating the RMS.  Reasons for not being able to
        // fill in the atoms include: there was a gap between residues,
        // the distance between neighboring C-alpha atoms in the quad
        // was greater than expected (often due to a Proline), this Residue
        // was not a natural amino acid (no C-alpha), etc.
        if(bbqAV.size() <= 1)
            continue;

        outputFilestream << bbqRes->getChainId() << "-" << bbqRes->getResidueNumber() << "_" << bbqRes->getResidueName() << ",";
        outputFilestream << originalRes->getChainId() << "-" << originalRes->getResidueNumber() << "_";
        outputFilestream  << originalRes->getResidueName() << ",";

        for(AtomPointerVector::iterator currIter = bbqAV.begin(); currIter != bbqAV.end(); ++currIter) {
            Atom *bbqAtom = *currIter;
            Atom *origAtom = &(originalRes->getAtom( bbqAtom->getName() ));

            CartesianPoint &bbqCP = bbqAtom->getCoor();
            CartesianPoint &origCP = origAtom->getCoor();

            Real distance = bbqCP.distance(origCP);
            totalDistance += distance;
            outputFilestream << origCP.getX() << "," << origCP.getY() << "," << origCP.getZ() << ",";
            outputFilestream << bbqCP.getX() << "," << bbqCP.getY() << "," << bbqCP.getZ() << ",";
   
            //cout << distance << ",";
            // outputFilestream << bbqCP << "," << origCP << ",";
        }
        int numAtoms = bbqAV.size();

        // We only want to calculate the RMS of the non-C-alpha atoms, so subtract
        // 1 from the total number of atoms.  Of course, if we only had a C-alpha
        // atom, then we don't want to divide by 0, so just leave it as 1.
        if(numAtoms > 1)
            --numAtoms;

        outputFilestream << totalDistance / (Real)numAtoms << "\n";
    }
}

void writePdb(vector<Residue *> &resVec, string pdbFileName) {
    PDBWriter pdbWriter;
    AtomPointerVector av;
    
    for(vector<Residue *>::iterator currIter = resVec.begin(); currIter != resVec.end(); ++currIter) {
        av += (*currIter)->getAtomPointers();
    }
    
    pdbWriter.open(pdbFileName, 2);
    pdbWriter.write(av);
    pdbWriter.close();
}

bool resHasProperAtoms(Residue &res, map<string, bool> &atomsOfInterest) {
    // All residues must have a C-alpha.
    bool hasProperAtoms = res.atomExists("CA");
    
    for(map<string, bool>::iterator currIter = atomsOfInterest.begin(); currIter != atomsOfInterest.end(); ++currIter) {
        hasProperAtoms = hasProperAtoms && res.atomExists(currIter->first);
    }
    
    return hasProperAtoms;
}

void removeUninterestingAtoms(Residue &res, map<string, bool> &atomsOfInterest) {
    Residue resCopy;
    
    resCopy.addAtom( res.getAtom("CA") );
    for(map<string, bool>::iterator currIter = atomsOfInterest.begin(); currIter != atomsOfInterest.end(); ++currIter) {
        if( res.atomExists(currIter->first) )
            resCopy.addAtom( res.getAtom(currIter->first) );
    }
    
    res.removeAllAtoms();
    res.addAtoms( resCopy.getAtomPointers() );
}

void printDihedralAngles(vector<Residue *> &resVec) {
    // We will only look at dihedral info if there are at least 4 residues.
    if(resVec.size() < 4)
        return;
    Residue *prevRes, *currRes, *nextRes;
    bool firstRes = true;


    for(vector<Residue *>::iterator currIter = (resVec.begin()+1); currIter != (resVec.end()-2); ++currIter ){
        prevRes = *(currIter - 1);
        currRes = *currIter;
        nextRes = *(currIter + 1);

        // Make sure that the prev residue and this residue are connected.
        // The first residue (res+1) only has CA, C, and O filled in, no N.
        if(!firstRes) {
            if( passDistanceTest(prevRes, currRes) ) {
                float phi = CartesianGeometry::dihedral(prevRes->getAtom("C").getCoor(),
                                                                    currRes->getAtom("N").getCoor(),
                                                                    currRes->getAtom("CA").getCoor(),
                                                                    currRes->getAtom("C").getCoor());

                cout << "Phi," << phi << "\n";
            }
        }

        // Make sure that the next residue and this residue are connected.
        if( passDistanceTest(currRes, nextRes) ) {
            // The first residue (res+1) only has CA, C, and O filled in, no N.
            if(!firstRes) {
                float psi = CartesianGeometry::dihedral(currRes->getAtom("N").getCoor(),
                                                                    currRes->getAtom("CA").getCoor(),
                                                                    currRes->getAtom("C").getCoor(),
                                                                    nextRes->getAtom("N").getCoor());

                cout << "Psi," << psi << "\n";
            }
            
            float omega = CartesianGeometry::dihedral(currRes->getAtom("CA").getCoor(),
                                                                currRes->getAtom("C").getCoor(),
                                                                nextRes->getAtom("N").getCoor(),
                                                                nextRes->getAtom("CA").getCoor());

            cout << "Omega," << omega << "\n";

        }

        firstRes = false;
    }
}

bool passDistanceTest(Residue *res1, Residue *res2) {
    float dist, diff, diffPro;

    dist = res1->getAtom("CA").distance( res2->getAtom("CA") );
    diff = fabs(dist - CA_CA_DISTANCE);
    diffPro = fabs(dist - CA_CA_DISTANCE_PROLINE);

    return( (diff <= CA_CA_TOLERANCE) || (diffPro <= CA_CA_TOLERANCE) );
}
