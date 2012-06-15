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

#include "PhiPsiStatistics.h"

using namespace MSL;
using namespace std;



PhiPsiStatistics::PhiPsiStatistics(){
    gridSize = 0.0f;
}
PhiPsiStatistics::PhiPsiStatistics(const PhiPsiStatistics &_phiPsiStat){
	copy(_phiPsiStat);
}
PhiPsiStatistics::~PhiPsiStatistics(){
}

void PhiPsiStatistics::operator=(const PhiPsiStatistics &_phiPsiStat){
	copy(_phiPsiStat);
}
void PhiPsiStatistics::copy(const PhiPsiStatistics &_phiPsiStat){

    phiPsiTable = _phiPsiStat.getPhiPsiCounts();
    gridSize = _phiPsiStat.gridSize;
}

void PhiPsiStatistics::addStatisitics(string _residueType, string _phiBin, string _psiBin, int _count){
	stringstream ss;
	ss << _residueType<<":"<<_phiBin<<":"<<_psiBin;

    double psi = MslTools::toDouble(_psiBin);
  
    // So the first time through step will = 0.
    // We assume that the first psi will be a negative number.
    // Then the second time through, we'll calculate the diff
    // between the present value of psi and previous value
    // and that will be the step size.
    if(gridSize == 0.0f)
        gridSize = psi;
    else if(gridSize < 0.0f)
        gridSize = psi - gridSize;
	
     //cout << "Adding: "<<ss.str()<<" "<<_count<<endl;
     phiPsiTable[ss.str()]     += _count;
     phiPsiTable[_residueType] += _count;

    
}


int PhiPsiStatistics::operator()(string _key){

	map<string,int>::iterator it;
	it = phiPsiTable.find(_key);
	if (it == phiPsiTable.end()){
		return MslTools::intMax;
	}

	return it->second;
}

int PhiPsiStatistics::getCounts(string resName, double phi, double psi){
    stringstream ss;
	ss << resName <<":"<<getPhiPsiBin(phi)<<":"<<getPhiPsiBin(psi);

	//cout << "Key: "<<ss.str()<<"."<<endl;
	map<string, int>::iterator it;
	it = phiPsiTable.find(ss.str());
	if (it == phiPsiTable.end()){
		cout << "No Phi/Psi entry found for " << ss.str() << "." <<endl;
		return MslTools::intMax;
	}
	
	return it->second;
}

int PhiPsiStatistics::getCounts(const Residue &nMinus1, const Residue &n, const Residue &nPlus1){
	double phi = getPhi(nMinus1, n);
	double psi = getPsi(n, nPlus1);
	string resName = n.getResidueName();

    return getCounts(resName, phi, psi);
}

double PhiPsiStatistics::getFreqInBin(string resName, double phi, double psi){

        // Counts for residue type in phi,psi bin
	int AAxyz = getCounts(resName, phi, psi);
	if (AAxyz == MslTools::intMax) {
		return MslTools::doubleMax;
	}
	double AAxyzDouble = (double)AAxyz;

	stringstream allkey;
	allkey << "ALL" << ":" << getPhiPsiBin(phi) << ":" << getPhiPsiBin(psi);

	map<string,int>::iterator it;
	it = phiPsiTable.find(allkey.str());
	if (it == phiPsiTable.end()){
		return MslTools::doubleMax;
	}
	int AAallyz = it->second;
	double AAallyzDouble = (double)AAallyz;


	return (AAxyzDouble  / AAallyzDouble);

	
}
double PhiPsiStatistics::getProbability(string &resName, double phi, double psi){

        // Counts for residue type in phi,psi bin
	int AAxyz = getCounts(resName, phi, psi);
	if (AAxyz == MslTools::intMax) {
		return MslTools::doubleMax;
	}

	// Total number of counts for residue type "resName"
	map<string,int>::iterator it;
	it = phiPsiTable.find(resName);
	if (it == phiPsiTable.end()) {
		return MslTools::doubleMax;
	}
	int AAx   = it->second;
	
	double AAxyzDouble = (double)AAxyz;
	double AAxDouble   = (double)AAx;

	//double mlogprob = 0.0;
	if (AAxyzDouble > 0.00001 && AAxDouble > 0.00001){
		//mlogprob = -log(AAxyzDouble / AAxDouble);
                //cout << "Vals: "<<phi<<","<<psi<<" "<<AAxyzDouble<<" "<<AAxDouble<<endl;
		return AAxyzDouble/ AAxDouble;
	}

	//return mlogprob;
	return 0.0;
}

double PhiPsiStatistics::getProbability(const Residue &nMinus1, const Residue &n, const Residue &nPlus1){
    double phi = getPhi(nMinus1, n);
    double psi = getPsi(n, nPlus1);
    string resName = n.getResidueName();
	
    return getProbability(resName, phi, psi);
}

double PhiPsiStatistics::getProbabilityAll(double phi, double psi){
	// (#AA-Phi-Psi(all,y,z) / #AA(all)) 
	int AAallyz = 0;
	int AAall   = 0;

	stringstream allkey;
	allkey << "ALL" << ":" << getPhiPsiBin(phi) << ":" << getPhiPsiBin(psi);

	map<string,int>::iterator it;
	it = phiPsiTable.find(allkey.str());
	if (it == phiPsiTable.end()){
		return MslTools::doubleMax;
	}
	AAallyz = it->second;
	double AAallyzDouble = (double)AAallyz;

	it = phiPsiTable.find("ALL");
	if (it == phiPsiTable.end()){
		return MslTools::doubleMax;
	}
	AAall = it->second;

	double AAallDouble   = (double)AAall;

	return AAallyzDouble / AAallDouble;
}

double PhiPsiStatistics::getProbabilityAll(const Residue &nMinus1, const Residue &n, const Residue &nPlus1){
    double phi = getPhi(nMinus1, n);
	double psi = getPsi(n, nPlus1);
	
	return getProbabilityAll(phi, psi);
}

double PhiPsiStatistics::getPropensity(string &resName, double phi, double psi){
	//  propensity(x,y,z) = 
	//	  (#AA-Phi-Psi(x,y,z)   / #AA(x)) 
	//	  -------------------------------
	//	  (#AA-Phi-Psi(all,y,z) / #AA(all)) 
    double small = 0.0001;
    double bigProp = 30.0f;

	double probRes = getProbability(resName, phi, psi);
	double probAll = getProbabilityAll(phi, psi);

    if( (probRes == MslTools::doubleMax) || (probAll == MslTools::doubleMax))
        return MslTools::doubleMax;

	// If this Phi/Psi combination is rare in general and for this
    // AA in particular, return 1.  In other words, this Amino Acid
    // is no more or less likely than the average to have this Phi/Psi comb.
    // Also, cap the probability of this Phi/Psi combination for the average
    // AA to some small number.
    if( (probRes < small) && (probAll < small) )
        return 1.0;
    else if( (probAll < small) ) {
        // Note, we should never have a case where probAll == 0 but probRes
        // doesn't.  That wouldn't really make sense.  However, I suppose
        // due to some float imprecision, it is better to be safe here and
        // explicitly check for that.
        if(probAll == 0.0)
            return bigProp;
            
        double prop = probRes/probAll;
        if(prop > bigProp)
            return bigProp;
    }

	return probRes/probAll;
}

double PhiPsiStatistics::getPropensity(const Residue &nMinus1, const Residue &n, const Residue &nPlus1){
    double phi = getPhi(nMinus1, n);
	double psi = getPsi(n, nPlus1);
	string resName = n.getResidueName();
	
	return getPropensity(resName, phi, psi);
}


void PhiPsiStatistics::computeTotalCounts(){
	map<string,int>::iterator phiPsiIt;
	map<string,int> runningTotalByRes;
	int runningTotal = 0;
	int runningTotalbySubtotal = 0;
	for (phiPsiIt = phiPsiTable.begin(); phiPsiIt != phiPsiTable.end();phiPsiIt++){
		
		// Get key, get tokens. if token size = 3, then do counts
		string key = phiPsiIt->first;
		vector<string> toks = MslTools::tokenize(key,":");
		if (toks.size() == 3){
			if (toks[0] != "ALL"){
				stringstream newkey;
				newkey << "ALL"<<":"<<toks[1]<<":"<<toks[2];
				phiPsiTable[newkey.str()] += phiPsiIt->second;
				runningTotal += phiPsiIt->second;
				runningTotalByRes[toks[0]] += phiPsiIt->second;
			}

		} else {
			if (toks[0] != "ALL") {
				//cout << "Subtotal-by-read-in for "<<phiPsiIt->first<<" is "<< phiPsiIt->second<<endl;
				runningTotalbySubtotal += phiPsiIt->second;
			}
		}
	}

	//for (phiPsiIt = runningTotalByRes.begin(); phiPsiIt != runningTotalByRes.end();phiPsiIt++){
		//cout << "Subtotal-by-res for "<<phiPsiIt->first<<" is "<< phiPsiIt->second<<endl;
	//}
	if (runningTotalbySubtotal != runningTotal){
		cout << "ERROR 3333 subtotals don't match: "<<runningTotalbySubtotal<<" and "<<runningTotal<<endl;
	}

	phiPsiTable["ALL"] = runningTotal;
}

double PhiPsiStatistics::getPhiPsiBin(double _in){
	double out;
	int cint = int(_in*100);
	int gint = int(gridSize *100);
	double remainder = gint - (cint % gint);

	out = ( int((_in*100 - (gridSize*100 - remainder))/100));
	if (out < 0) {
		out -= (gridSize / 2);
	} else {
		out += (gridSize / 2);
	}
	
	return out;
}
double PhiPsiStatistics::getPhi(const Residue &nMinus1, const Residue &n){
    Residue &ncnMinus1 = const_cast<Residue&>(nMinus1);
    Residue &ncn = const_cast<Residue&>(n);

    if (!(ncnMinus1.atomExists("C") && ncn.atomExists("N") && ncn.atomExists("CA") && ncn.atomExists("C"))){
        return MslTools::doubleMax;
    }
    return CartesianGeometry::dihedral(ncnMinus1("C").getCoor(), ncn("N").getCoor(), ncn("CA").getCoor(), ncn("C").getCoor());
}

double PhiPsiStatistics::getPsi(const Residue &n, const Residue &nPlus1){
    Residue &ncn = const_cast<Residue&>(n);
    Residue &ncnPlus1 = const_cast<Residue&>(nPlus1);

    if (!(ncn.atomExists("N") && ncn.atomExists("CA") && ncn.atomExists("C") && ncnPlus1.atomExists("N"))){
        return MslTools::doubleMax;
    }
    return CartesianGeometry::dihedral(ncn("N").getCoor(),ncn("CA").getCoor(), ncn("C").getCoor(), ncnPlus1("N").getCoor());
}

double PhiPsiStatistics::getOmega(const Residue &n, const Residue &nPlus1){
    Residue &ncn = const_cast<Residue&>(n);
    Residue &ncnPlus1 = const_cast<Residue&>(nPlus1);

    if (!(ncn.atomExists("CA") && ncn.atomExists("C") && ncnPlus1.atomExists("N") && ncnPlus1.atomExists("CA"))){
        return MslTools::doubleMax;
    }
    return CartesianGeometry::dihedral(ncn("CA").getCoor(),ncn("C").getCoor(), ncnPlus1("N").getCoor(), ncnPlus1("CA").getCoor());
}



pair<double,double> PhiPsiStatistics::getRandomPhiPsi(std::string _resType){


    std::map<std::string,PhiPsiRNG *>::iterator itRand = phiPsiRandom.find(_resType);


    // If we do not have a random number generator, then build one for this residue type
    if (itRand == phiPsiRandom.end()){
	    std::map<std::string,int>::iterator itPhiPsi = phiPsiTable.find(_resType);

	    // If we have no data on this residue type then error.
	    if (itPhiPsi == phiPsiTable.end()){
		cerr << "ERROR 6495 PhiPsiStatistics::getRandomPhiPsi(string resType) could not find resType in dataset("<<_resType<<")\n";
	        exit(6495);
	    }

	    // We don't have a random number generator, but we do have data for this residueType.
	    PhiPsiRNG *ppRNG = new PhiPsiRNG();

	    // Iterate over keys in phiPsiTable, look for ones with data for this residue, store
	    for (itPhiPsi = phiPsiTable.begin();itPhiPsi != phiPsiTable.end();itPhiPsi++){

		std::string key = itPhiPsi->first;
		std::vector<std::string> tokens = MslTools::tokenize(key, ":");

		// Skip keys that are not in this format 'resName:phiBin:psiBin"
		if (tokens.size() != 3) continue;

		double phiMid = MslTools::toDouble(tokens[1]) + gridSize/2;
		double psiMid = MslTools::toDouble(tokens[2]) + gridSize/2;

		ppRNG->counts.push_back((double)itPhiPsi->second);
		ppRNG->phiPsiValues.push_back(std::pair<double,double>(phiMid,psiMid));

	    }	    


	/*  CONVERSION NO LONGER NEEDED
	    // Convert to c-style array for GSL
	    double arr[ppRNG->counts.size()];
	    for (uint i = 0; i  < ppRNG->counts.size();i++){
		arr[i] = ppRNG->counts[i];
	    }
	*/

	    //ppRNG->rng.setDiscreteProb(arr, ppRNG->counts.size());
	    ppRNG->rng.setDiscreteProb(ppRNG->counts);


	    phiPsiRandom[_resType] = ppRNG;
	    ppRNG = NULL;


	    // Now find it!
	    itRand = phiPsiRandom.find(_resType);
	    
    }


    PhiPsiRNG *ppRNG = itRand->second;    
    
    
    return ppRNG->getRandomAngles();
}
    




