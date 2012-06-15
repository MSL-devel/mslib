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
#include "PDBTopology.h"
#include "IcTable.h"
using namespace MSL;


#include "MslOut.h"
static MslOut MSLOUT("PDBTopology");

PDBTopology::PDBTopology(){
  setup();
}
PDBTopology::~PDBTopology(){
}

PDBTopology::PDBTopology(const PDBTopology & _top){
  copy(_top);
}

void PDBTopology::operator=(const PDBTopology & _top){
  copy(_top);
}

bool PDBTopology::readCharmmTopology(std::string _filename){
  reset();
  setup();
  charmmReader.setCharmmFileName(_filename);
  charmmReader.read();


  return true;
}

bool PDBTopology::readRotamerLibrary(std::string _filename){
  rotlib.readFile(_filename);

  return true;
}

bool PDBTopology::residueExists(std::string _name){
  std::map<std::string,std::vector<std::string> >::iterator it;
  it = atoms.find(_name);
  if (it == atoms.end()){
    return false;
  }
  return true;

}

AtomContainer PDBTopology::getGenericResidue(std::string _identityId, int _numRotamers){

  Atom N("N",-10.780,2.171,9.692,"N");
  Atom CA("CA",-10.799,0.924,10.492,"C");
  Atom C("C",-9.384,0.388,10.623,"C");

  AtomContainer ats;
  ats.addAtom(N);
  ats.addAtom(CA);
  ats.addAtom(C);

  AtomContainer result = getResidue(_identityId,ats.getAtomPointers(),_numRotamers);

  return result;

}
AtomContainer PDBTopology::getResidue(std::string _identityId) {
  AtomContainer foo;
  return foo;
}

AtomContainer PDBTopology::getResidue(std::string _identityId, AtomPointerVector &_seedAtoms, int _numRotamers){

  /*
    Build using a CharmmTopologyReader defintion 
  */	     

  // Parse the identityId into components...
  string chainId;
  string resName;
  int resNum;
  string icode;
  MslTools::parseIdentityId(_identityId,chainId,resNum,icode,resName);


   // Create storage for the new residue
   AtomContainer newResidue;

   // Add seed atoms to newResidue make 
   map<string,vector<string> >::iterator atomIt;
   atomIt = atoms.find(resName);
   if (atomIt == atoms.end()){
     cerr << "ERROR PDBTopology::getResidue(), can't find residue: '"<<resName<<"'"<<endl;
     exit(1234);
   }
   std::vector<std::string> atomsInResidue = atomIt->second;
   if (addAtomsFromRotLib){
     std::vector<std::string> mobileAtoms = rotlib.getMobileAtoms(rotlib.getDefaultLibrary(),resName);
     atomsInResidue.insert(atomsInResidue.end(),mobileAtoms.begin(),mobileAtoms.end());
     std::sort(atomsInResidue.begin(),atomsInResidue.end());
     atomsInResidue.erase(std::unique(atomsInResidue.begin(),atomsInResidue.end()),atomsInResidue.end());
   }
   for (unsigned int i = 0; i < atomsInResidue.size();i++){


     int foundAtom = -1;
     for (unsigned int a = 0; a < _seedAtoms.size();a++){
       if (atomsInResidue[i] == _seedAtoms(a).getName()){
	 foundAtom = a;
	 break;
       }
     }

     Atom *tmp = new Atom();
     if (foundAtom != -1){
       tmp->setCoor(_seedAtoms(foundAtom).getCoor());
       tmp->setHasCoordinates(true);
     } else {
       tmp->setHasCoordinates(false);
     }

     tmp->setChainId(chainId);
     tmp->setResidueNumber(resNum);
     tmp->setResidueIcode(icode);
     tmp->setResidueName(resName);
     tmp->setName(atomsInResidue[i]);
     newResidue.addAtom(*tmp);
     delete tmp;
   }

   MSLOUT.stream() << "Add: "<< endl<<newResidue.toString() <<endl;

   // First try to use the rotamer library approach
   bool builtResidue = false;

   buildRotamers(newResidue, resName, _numRotamers);
   builtResidue = true;
   MSLOUT.stream() << "Done: "<< endl<<newResidue.toString() <<endl;

   

   if (!builtResidue){
     cerr << "ERROR 9243 PDBTopology::getResidue() residue was not built!\n";
     exit(9243);
   }

  return newResidue;

}




void PDBTopology::reset() {
  atoms.clear();
  chis.clear();
}

void PDBTopology::setup(){

  addAtomsFromRotLib = false;

  // Hard coded atom names for each residue type
  //if (charmmReader.getCharmmFileName() == ""){
  	atoms["ALA"].push_back("N");
	atoms["ALA"].push_back("CA");
	atoms["ALA"].push_back("CB");
	atoms["ALA"].push_back("C");
	atoms["ALA"].push_back("O");

	atoms["CYS"].push_back("N");
	atoms["CYS"].push_back("CA");
	atoms["CYS"].push_back("CB");
	atoms["CYS"].push_back("SG");
	atoms["CYS"].push_back("C");
	atoms["CYS"].push_back("O");

 	atoms["ASP"].push_back("N");
	atoms["ASP"].push_back("CA");
	atoms["ASP"].push_back("CB");
	atoms["ASP"].push_back("CG");
	atoms["ASP"].push_back("OD1");
	atoms["ASP"].push_back("OD2");
	atoms["ASP"].push_back("C");
	atoms["ASP"].push_back("O");

 	atoms["D1Z"].push_back("N");
	atoms["D1Z"].push_back("CA");
	atoms["D1Z"].push_back("CB");
	atoms["D1Z"].push_back("CG");
	atoms["D1Z"].push_back("OD1");
	atoms["D1Z"].push_back("OD2");
	atoms["D1Z"].push_back("ZN");
	atoms["D1Z"].push_back("C");
	atoms["D1Z"].push_back("O");

	
	atoms["GLU"].push_back("N"); 
	atoms["GLU"].push_back("CA"); 
	atoms["GLU"].push_back("CB"); 
	atoms["GLU"].push_back("CG"); 
	atoms["GLU"].push_back("CD"); 
	atoms["GLU"].push_back("OE1"); 
	atoms["GLU"].push_back("OE2"); 
	atoms["GLU"].push_back("C"); 
	atoms["GLU"].push_back("O"); 
	
	atoms["PHE"].push_back("N");
	atoms["PHE"].push_back("CA");
	atoms["PHE"].push_back("CB");
	atoms["PHE"].push_back("CG");
	atoms["PHE"].push_back("CD1");
	atoms["PHE"].push_back("CD2");
	atoms["PHE"].push_back("CE1");
	atoms["PHE"].push_back("CE2");
	atoms["PHE"].push_back("CZ");
	atoms["PHE"].push_back("C");
	atoms["PHE"].push_back("O");
	
	atoms["GLY"].push_back("N"); 
	atoms["GLY"].push_back("CA"); 
	atoms["GLY"].push_back("C"); 
	atoms["GLY"].push_back("O"); 
 
	atoms["HIS"].push_back("N");	
	atoms["HIS"].push_back("CA");	
	atoms["HIS"].push_back("ND1");	
	atoms["HIS"].push_back("CG");	
	atoms["HIS"].push_back("CB");	
	atoms["HIS"].push_back("NE2");	
	atoms["HIS"].push_back("CD2");	
	atoms["HIS"].push_back("CE1");	
	atoms["HIS"].push_back("C");	
	atoms["HIS"].push_back("O");	


	atoms["HSE"].push_back("N");	
	atoms["HSE"].push_back("CA");	
	atoms["HSE"].push_back("ND1");	
	atoms["HSE"].push_back("CG");	
	atoms["HSE"].push_back("CB");	
	atoms["HSE"].push_back("NE2");	
	atoms["HSE"].push_back("CD2");	
	atoms["HSE"].push_back("CE1");	
	atoms["HSE"].push_back("C");	
	atoms["HSE"].push_back("O");	

	atoms["HSD"].push_back("N");	
	atoms["HSD"].push_back("CA");	
	atoms["HSD"].push_back("ND1");	
	atoms["HSD"].push_back("CG");	
	atoms["HSD"].push_back("CB");	
	atoms["HSD"].push_back("NE2");	
	atoms["HSD"].push_back("CD2");	
	atoms["HSD"].push_back("CE1");	
	atoms["HSD"].push_back("C");	
	atoms["HSD"].push_back("O");	

	atoms["HSP"].push_back("N");	
	atoms["HSP"].push_back("CA");	
	atoms["HSP"].push_back("ND1");	
	atoms["HSP"].push_back("CG");	
	atoms["HSP"].push_back("CB");	
	atoms["HSP"].push_back("NE2");	
	atoms["HSP"].push_back("CD2");	
	atoms["HSP"].push_back("CE1");	
	atoms["HSP"].push_back("C");	
	atoms["HSP"].push_back("O");	


	atoms["ILE"].push_back("N");
	atoms["ILE"].push_back("CA");
	atoms["ILE"].push_back("CB");
	atoms["ILE"].push_back("CG2");
	atoms["ILE"].push_back("CG1");
	atoms["ILE"].push_back("CD");
	atoms["ILE"].push_back("C");
	atoms["ILE"].push_back("O");

	atoms["LYS"].push_back("N");
	atoms["LYS"].push_back("CA");
	atoms["LYS"].push_back("CB");
	atoms["LYS"].push_back("CG");
	atoms["LYS"].push_back("CD");
	atoms["LYS"].push_back("CE");
	atoms["LYS"].push_back("NZ");
	atoms["LYS"].push_back("C");
	atoms["LYS"].push_back("O");
	

 	atoms["LEU"].push_back("N");
	atoms["LEU"].push_back("CA");
	atoms["LEU"].push_back("CB");
	atoms["LEU"].push_back("CG");
	atoms["LEU"].push_back("CD1");
	atoms["LEU"].push_back("CD2");
	atoms["LEU"].push_back("C");
 	atoms["LEU"].push_back("O");

	atoms["MET"].push_back("N");
	atoms["MET"].push_back("CA");
	atoms["MET"].push_back("CB");
	atoms["MET"].push_back("CG");
	atoms["MET"].push_back("SD");
	atoms["MET"].push_back("CE");
	atoms["MET"].push_back("C");
	atoms["MET"].push_back("O");


	atoms["ASN"].push_back("N");
	atoms["ASN"].push_back("CA");
	atoms["ASN"].push_back("CB");
	atoms["ASN"].push_back("CG");
	atoms["ASN"].push_back("OD1");
	atoms["ASN"].push_back("ND2");
	atoms["ASN"].push_back("C");
	atoms["ASN"].push_back("O");


	atoms["PRO"].push_back("N");
	atoms["PRO"].push_back("CA");
	atoms["PRO"].push_back("CD");
	atoms["PRO"].push_back("CB");
	atoms["PRO"].push_back("CG");
	atoms["PRO"].push_back("C");
	atoms["PRO"].push_back("O");
	

	atoms["GLN"].push_back("N");
	atoms["GLN"].push_back("CA");
	atoms["GLN"].push_back("CB");
	atoms["GLN"].push_back("CG");
	atoms["GLN"].push_back("CD");
	atoms["GLN"].push_back("OE1");
	atoms["GLN"].push_back("NE2");
	atoms["GLN"].push_back("C");
	atoms["GLN"].push_back("O");


	atoms["ARG"].push_back("N");
	atoms["ARG"].push_back("CA");
	atoms["ARG"].push_back("CB");
	atoms["ARG"].push_back("CG");
	atoms["ARG"].push_back("CD");
	atoms["ARG"].push_back("NE");
	atoms["ARG"].push_back("CZ");
	atoms["ARG"].push_back("NH1");
	atoms["ARG"].push_back("NH2");
	atoms["ARG"].push_back("C");
	atoms["ARG"].push_back("O");


	atoms["SER"].push_back("N");
	atoms["SER"].push_back("CA");
	atoms["SER"].push_back("CB");
	atoms["SER"].push_back("OG");
	atoms["SER"].push_back("C");
	atoms["SER"].push_back("O");
	
 	atoms["THR"].push_back("N");
	atoms["THR"].push_back("CA");
	atoms["THR"].push_back("CB");
	atoms["THR"].push_back("OG1");
	atoms["THR"].push_back("CG2");
	atoms["THR"].push_back("C");
	atoms["THR"].push_back("O");
	

	atoms["VAL"].push_back("N");
	atoms["VAL"].push_back("CA");
	atoms["VAL"].push_back("CB");
	atoms["VAL"].push_back("CG1");
	atoms["VAL"].push_back("CG2");
	atoms["VAL"].push_back("C");
	atoms["VAL"].push_back("O");

	atoms["TYR"].push_back("N");
	atoms["TYR"].push_back("CA");
	atoms["TYR"].push_back("CB");
	atoms["TYR"].push_back("CG");
	atoms["TYR"].push_back("CD1");
	atoms["TYR"].push_back("CD2");
	atoms["TYR"].push_back("CE1");
	atoms["TYR"].push_back("CE2");
	atoms["TYR"].push_back("CZ");
	atoms["TYR"].push_back("OH");
	atoms["TYR"].push_back("C");
	atoms["TYR"].push_back("O");
	
	atoms["TYS"].push_back("N");
	atoms["TYS"].push_back("CA");
	atoms["TYS"].push_back("CB");
	atoms["TYS"].push_back("CG");
	atoms["TYS"].push_back("CD1");
	atoms["TYS"].push_back("CD2");
	atoms["TYS"].push_back("CE1");
	atoms["TYS"].push_back("CE2");
	atoms["TYS"].push_back("CZ");
	atoms["TYS"].push_back("OH");
	atoms["TYS"].push_back("S");
	atoms["TYS"].push_back("O1");
	atoms["TYS"].push_back("O2");
	atoms["TYS"].push_back("O3");
	atoms["TYS"].push_back("C");
	atoms["TYS"].push_back("O");


	atoms["TRP"].push_back("N");
	atoms["TRP"].push_back("CA");
	atoms["TRP"].push_back("CB");
	atoms["TRP"].push_back("CG");
	atoms["TRP"].push_back("CD2");
	atoms["TRP"].push_back("CD1");
	atoms["TRP"].push_back("NE1");
	atoms["TRP"].push_back("CE2");
	atoms["TRP"].push_back("CE3");
	atoms["TRP"].push_back("CZ2");
	atoms["TRP"].push_back("CZ3");
	atoms["TRP"].push_back("CH2");
	atoms["TRP"].push_back("C");
	atoms["TRP"].push_back("O");
	


	atoms["HDZ"].push_back("N");	
	atoms["HDZ"].push_back("CA");	
	atoms["HDZ"].push_back("ND1");	
	atoms["HDZ"].push_back("CG");	
	atoms["HDZ"].push_back("CB");	
	atoms["HDZ"].push_back("NE2");	
	atoms["HDZ"].push_back("CD2");	
	atoms["HDZ"].push_back("CE1");	
	atoms["HDZ"].push_back("ZN");	
	atoms["HDZ"].push_back("C");	
	atoms["HDZ"].push_back("O");	

	atoms["HEZ"].push_back("N");	
	atoms["HEZ"].push_back("CA");	
	atoms["HEZ"].push_back("ND1");	
	atoms["HEZ"].push_back("CG");	
	atoms["HEZ"].push_back("CB");	
	atoms["HEZ"].push_back("NE2");	
	atoms["HEZ"].push_back("CD2");	
	atoms["HEZ"].push_back("CE1");	
	atoms["HEZ"].push_back("ZN");	
	atoms["HEZ"].push_back("C");	
	atoms["HEZ"].push_back("O");	




	std::vector<std::string> tmp;
	tmp.clear();
	
	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG"); 

	chis["ARG"].push_back(tmp);
	
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD");
	chis["ARG"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD");
	tmp.push_back("NE");
	chis["ARG"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CG");
	tmp.push_back("CD");
	tmp.push_back("NE");
	tmp.push_back("CZ");
	chis["ARG"].push_back(tmp);
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["ASN"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("OD1");
	chis["ASN"].push_back(tmp);

	      
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG"); 
	chis["ASP"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("OD1");
	chis["ASP"].push_back(tmp);

	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("SG");
	chis["CYS"].push_back(tmp);
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["GLN"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD");
	chis["GLN"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD");
	tmp.push_back("OE1");
	chis["GLN"].push_back(tmp);
	       
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["GLU"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD");
	chis["GLU"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD");
	tmp.push_back("OE1");
	chis["GLU"].push_back(tmp);
	       
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["HIS"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("ND1");
	chis["HIS"].push_back(tmp);
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["HSD"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("ND1");
	chis["HSD"].push_back(tmp);
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["HDF"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("ND1");
	chis["HDF"].push_back(tmp);
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["HEF"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("ND1");
	chis["HEF"].push_back(tmp);
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["HSE"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("ND1");
	chis["HSE"].push_back(tmp);
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG1");
	chis["ILE"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG1");
	tmp.push_back("CD1");
	chis["ILE"].push_back(tmp);
	       
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["LEU"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD1");
	chis["LEU"].push_back(tmp);
	       
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["LYS"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD");
	chis["LYS"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD");
	tmp.push_back("CE");
	chis["LYS"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CG");
	tmp.push_back("CD");
	tmp.push_back("CE");
	tmp.push_back("NZ");
	chis["LYS"].push_back(tmp);
	       
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["MET"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("SD");
	chis["MET"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("SD");
	tmp.push_back("CE");
	chis["MET"].push_back(tmp);
	       
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["PHE"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD1");
	chis["PHE"].push_back(tmp);
	       
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["PRO"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD");
	chis["PRO"].push_back(tmp);
	       
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("OG");
	chis["SER"].push_back(tmp);
	       
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("OG1");
	chis["THR"].push_back(tmp);
	       
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["TRP"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD1");
	chis["TRP"].push_back(tmp);
	       
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	chis["TYR"].push_back(tmp);
	tmp.clear();

	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG");
	tmp.push_back("CD1");
	chis["TYR"].push_back(tmp);
	       
	tmp.clear();

	tmp.push_back("N");
	tmp.push_back("CA");
	tmp.push_back("CB");
	tmp.push_back("CG1");
	chis["VAL"].push_back(tmp);


	backboneAtoms["N"]  = true;
	backboneAtoms["CA"] = true;
	backboneAtoms["C"]  = true;
	backboneAtoms["O"]  = true;

	//} else {
    
    // Process CharmmTopologyResidues into atoms, chis and backboneAtoms structure...
    
	//}
  
}

void PDBTopology::copy(const PDBTopology & _top){
  // Nothing to copy
}



void PDBTopology::buildRotamers(AtomContainer &_newResidue, std::string _resName, int _numRotamers){
  
  
   // For each initialized atom defined in rotamer library...
   std::vector<std::string> mobileAtoms = rotlib.getMobileAtoms(rotlib.getDefaultLibrary(),_resName);
  
   if (mobileAtoms.size() == 0){
     cerr << "ERROR 4533 PDBTopology::buildRotamers(..) no init atoms found for library "<<rotlib.getDefaultLibrary()<<" and residue "<<_resName<<endl;
     exit(4533);
   }

   // Setup two atom pointer vectors, one for seed atoms one for the rest
   AtomPointerVector seeds;
   AtomPointerVector rest;
   for (unsigned int i = 0; i < _newResidue.size();i++){
     if (_newResidue[i].hasCoor()){
       seeds.push_back(&_newResidue[i]);
     } else {
       rest.push_back(&_newResidue[i]);
     }
    }

   // Extract DEFI for this rotamer library/residue type
   vector<RotamerLibrary::InternalCoorDefi> defi = rotlib.getInternalCoorDefinition(rotlib.getLibraryNames()[0], _resName);

   // Extract ICVALUES for this rotamer library/residue type
   vector<vector<double> > icValues = rotlib.getInternalCoor(rotlib.getLibraryNames()[0], _resName);


    // Store a table of IcEntries so we can build this residue
     IcTable icTab;
   
   // Each rotamer
   for (uint r = 0; r < _numRotamers;r++){
     if (r > icValues.size()) { break; }


     // Add altConformation for seed atoms....
     if (r > 0){
       for (unsigned int s = 0; s < seeds.size();s++){
	 seeds[s]->addAltConformation();
         seeds[s]->setActiveConformation(seeds[s]->getNumberOfAltConformations()-1);
       }
     }


     stringstream atomId;
     atomId << _newResidue[0].getChainId()<<","<<_newResidue[0].getResidueNumber()<<",";

     // Each initAtom
     for (uint i = 0; i < mobileAtoms.size();i++){

       // Find atom in residue
       if (! _newResidue.atomExists(atomId.str()+mobileAtoms[i])){
	 cerr << "ERROR PDBTopology::buildRotamers() atom "<<atomId.str()+mobileAtoms[i]<<" does not exist in residue "<<_newResidue[0].getResidueName()<<endl;
	 cerr << "\tcontinue without building this atom\n";
	 continue;
       }

       // Pull out ic-values and 4 atom pointers
       Atom *p1 = NULL;
       Atom *p2 = NULL;
       Atom *p3 = NULL;
       Atom *p4 = NULL;
       double bond = 0.0;
       double angle = 0.0;
       double dihedral = 0.0;
       bool improper_flag = false;

       stringstream ss;
       for (unsigned int d = 0; d < defi.size();d++){

        if (defi[d].atomNames.size() == 4 && defi[d].atomNames[3] == mobileAtoms[i]){



	 // Get dihedral value and other 3 atoms
	 if (_newResidue.atomExists(atomId.str()+defi[d].atomNames[0])){
	   p1 = &_newResidue.getLastFoundAtom();
	 } else {
	   MSLOUT.stream() << "ERROR PDBTopology atom "<<atomId.str()+defi[d].atomNames[0]<< " does not exist in resdiue "<<endl;
	 }
	 if (p1 == NULL){
	   MSLOUT.stream() << "P1 is still  NULL"<<endl; 
	 }
	 if (_newResidue.atomExists(atomId.str()+defi[d].atomNames[1])){
	   p2 = &_newResidue.getLastFoundAtom();
	 }else {
	   MSLOUT.stream() << "ERROR PDBTopology atom "<<atomId.str()+defi[d].atomNames[0]<< "does not exist in residue "<<endl;
	 }
	 if (p2 == NULL){
	   MSLOUT.stream() << "P2 is still  NULL"<<endl; 
	 }
	 if (_newResidue.atomExists(atomId.str()+defi[d].atomNames[2])){
	   p3 = &_newResidue.getLastFoundAtom();
	 }else {
	   MSLOUT.stream() << "ERROR PDBTopology atom "<<atomId.str()+defi[d].atomNames[0]<< "does not exist in resdiue "<<endl;
	 }
	 if (p3 == NULL){
	   MSLOUT.stream() << "P3 is still  NULL"<<endl; 
	 }
	 if (_newResidue.atomExists(atomId.str()+defi[d].atomNames[3])){
	   p4 = &_newResidue.getLastFoundAtom();
	 }else {
	   MSLOUT.stream() << "ERROR PDBTopology atom "<<atomId.str()+defi[d].atomNames[0]<< "does not exist in resdiue "<<endl;
	 }
	 if (p4 == NULL){
	   MSLOUT.stream() << "P4 is still  NULL"<<endl; 
	 }
	 dihedral = icValues[r][d];
	 if (defi[d].type == 3){
	   improper_flag = true;
	 }
	 MSLOUT.stream() << "Dihedral found: "<<dihedral<<" "<<improper_flag<<" ats: "<<p4->toString()<<endl;
       }

       if (defi[d].atomNames.size() == 3 && defi[d].atomNames[2] == mobileAtoms[i]){

	 angle = icValues[r][d];
	 MSLOUT.stream() << "Angle found: "<< angle<<endl;
       }

       if (defi[d].atomNames.size() == 2 && defi[d].atomNames[1] == mobileAtoms[i]){

	 bond = icValues[r][d];
	 MSLOUT.stream() << "Bond found: "<<bond<<endl;
       }
       
     }

     // Create IcEntries into IcTable if we have all our atoms defined
     if (p1 != NULL && p2 != NULL && p3 != NULL && p4 != NULL){
       //icTab.push_back(new IcEntry(*p1,*p2,*p3,*p4,bond,angle,dihedral,bond,angle,improper_flag));


	 if (improper_flag){
	   if (r == 0){
	     //	     icTab.push_back(new IcEntry(*p4,*p2,*p3,*p1,bond,angle,dihedral,angle,bond,improper_flag));
	     icTab.push_back(new IcEntry(*p1,*p2,*p3,*p4,bond,angle,dihedral,angle,bond,improper_flag));
	     icTab.push_back(new IcEntry(*p4,*p2,*p3,*p1,bond,angle,dihedral,angle,bond,improper_flag));
	   } else {
	     icTab.editBond(p4,p3,bond);
	     icTab.editBond(p3,p4,bond);
	     icTab.editAngle(p4,p3,p2,angle);
	     icTab.editAngle(p2,p3,p4,angle);
	     icTab.editDihedral(p4,p3,p2,p1,dihedral);
	     icTab.editDihedral(p1,p3,p2,p4,-dihedral);
	   }
	 } else {
	   if (r == 0){
	     icTab.push_back(new IcEntry(*p4,*p3,*p2,*p1,bond,angle,dihedral,angle,bond,improper_flag));
	   } else {
	     icTab.editBond(p4,p3,bond);
	     icTab.editAngle(p4,p3,p2,angle);
	     icTab.editDihedral(p4,p3,p2,p1,dihedral);
	   }
	 }

     } else {
       cerr << "ERROR PDBToplogy an atom was not found in this residue\n";
       exit(1234);
     }
       
     } // FOR (i) mobileAtoms


     // Manually Add Carbonyl Oxygen Ic (Usually not in rotlib)
     if (_newResidue.atomExists(atomId.str()+"O") &&_newResidue.atomExists(atomId.str()+"C") && _newResidue.atomExists(atomId.str()+"CA") && _newResidue.atomExists(atomId.str()+"CB")){
       MSLOUT.stream() << "ADDING CARBONYL OXYGEN"<<endl;
       icTab.push_back(new IcEntry(_newResidue(atomId.str()+"O"), _newResidue(atomId.str()+"C"),_newResidue(atomId.str()+"CA"),_newResidue(atomId.str()+"CB"),1.2,120,-98,120,1.2,false));
     } else {
       MSLOUT.stream() << "NO IC FOR CARBONYL OXYGEN!\n";
     }

     // Add conformation to each non-seed atom in residue
     for (unsigned int j = 0; j < rest.size();j++){

        if (r > 0){
	 rest[j]->addAltConformation();
       }

       rest[j]->setActiveConformation(rest[j]->getNumberOfAltConformations()-1);
       rest[j]->wipeCoordinates();
     }

     // Build from IcTable
     int numAtomsNotBuilt = 0;
     for (unsigned int j = 0; j < rest.size();j++){
       if (!rest[j]->buildFromIc()) {
	 numAtomsNotBuilt++;
       }
     }
     
     MSLOUT.stream() << "NUMBER ATOMS NOT BUILT: "<<numAtomsNotBuilt<<endl;
//     PDBWriter pout;
//     char tmp[80];
//     sprintf(tmp,"/tmp/rotamer-%04d.pdb",r);
//     pout.open(tmp);
//     pout.write(_newResidue.getAtomPointers());
//     pout.close();

   } // FOR (r) NUM_ROTAMERS

   icTab.deletePointers();
   seeds.clear();
   rest.clear();

}




AtomPointerVector PDBTopology::getBackboneAtoms(Residue &_res){

    AtomPointerVector results;
    map<string,bool>::iterator it;
    for (it = backboneAtoms.begin();it != backboneAtoms.end();it++){
            if (it->second && _res.atomExists(it->first)){
	         results.push_back(&_res.getAtom(it->first));
	    } else {
	         if (it->second){
		     MSLOUT.stream() << " getBackboneAtoms(Residue) , residue "<<_res.toString()<<" does not have backbone atom '"<<it->second<<"'"<<endl;
	         }
	    }
    }

    return results;
}

Atom * PDBTopology::getPseudoCbeta(Residue &_glycine) {

  if (!_glycine.atomExists("N") || !_glycine.atomExists("CA") || !_glycine.atomExists("C")){
    cerr << "ERROR 34234 PDBTopology::getPseudoCbeta N,CA or C does not exist in : "<<_glycine.toString()<<endl;
    exit(34234);
  }
  
  using namespace MSL::CartesianGeometry;

  Atom *a = NULL;

  a = new Atom(_glycine.getAtom("CA"));

  a->setName("CB");
  a->setCoor(build(_glycine.getAtom("CA").getCoor(),_glycine.getAtom("N").getCoor(),_glycine.getAtom("C").getCoor(),1.521, 110.5, -122.5)); 

  return a;
 
}
