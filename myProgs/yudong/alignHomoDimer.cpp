#include <iostream>
#include <fstream>
#include <cstdlib>
#include "System.h"
#include "MslTools.h"
#include "Transforms.h"
#include "OptionParser.h"
#include "PDBReader.h"
#include "AtomSelection.h"
#include "Transforms.h"

using namespace std;
using namespace MSL;



string programName = "alignHomoDimer";
string programDescription = "This program do something?";
struct Solution
{
        double rmsd;
        int file_name_index;
        bool operator < (const Solution& sol) const
        {
                return(rmsd < sol.rmsd);
        }
};
time_t startTime, endTime;
double diffTime;
vector<Solution> solutions;
/********************************************
 *  Align a dimer with a perfect dimer
 ********************************************/

        

int main(int argc, char **argv){
        time(&startTime);
        // Get command line args
        OptionParser OP;
        OP.readArgv(argc,argv);
        string pdbFile = argv[2];
        string perfectHelicesFiles = argv[4];
        string outputFile = argv[6];
        string helixFile = OP.getString("helixFile");
        string outputDir = OP.getString("outputDir");

        // Read the imperfect helix 
        System ip_helix;
        if(!ip_helix.readPdb(pdbFile))
                cerr << "Fail to read the file: " << pdbFile << endl;
        AtomSelection ip_helix_sl(ip_helix.getAtomPointers());
        AtomPointerVector ip_helix_ca = ip_helix_sl.select("bb, name CA");
        
        //cout << ip_helix_ca;


        // Read in a perfect helix
        ifstream pf_helices_list(perfectHelicesFiles.c_str());
        vector<string> pf_helices;
        string line;
        while(getline(pf_helices_list,line)){
                        pf_helices.push_back(line);
        }
        cout << "Total Number of files that is going to process: " << pf_helices.size() << endl;
        ofstream file_list; file_list.open("file_list.dat");
        for(int i = 0;i< pf_helices.size(); i++){
                file_list << pf_helices[i] << endl;
        }
        for(int i = 0;i < pf_helices.size();i++){
                System pf_helix;
                if(!pf_helix.readPdb(pf_helices[i]))
                        cerr << "Fail to read file: " << pf_helices[i] << endl;
                AtomSelection pf_helix_sl(pf_helix.getAtomPointers());
                AtomPointerVector pf_helix_ca = pf_helix_sl.select("bb, name CA");
                // Align it with the imperfect helix
                Transforms trans;
                trans.rmsdAlignment(ip_helix_ca,pf_helix_ca);
                trans.resetHistory();
                //Check the rmsd of current alignet
                solutions.push_back(Solution {ip_helix_ca.rmsd(pf_helix_ca),i});
                 
        }
        ofstream log;
        log.open(outputFile.c_str());
        sort(solutions.begin(),solutions.end());
        for(int i = 0; i < solutions.size(); i++){
                log << "RMSD: " << solutions[i].rmsd << " File name: " << pf_helices[solutions[i].file_name_index] << endl;
        }
        time(&endTime);
        diffTime = difftime (endTime, startTime);
        log << "Total Time: " << diffTime << " seconds" << endl;
        return 0;
}
