/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

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

#include "FormatConverter.h"

using namespace std;
using namespace MSL;

FormatConverter::FormatConverter() {
	orig = "CHARMM22";
	tgt = "PDB2.3";
}


FormatConverter::FormatConverter(string _orig, string _tgt) {
	setNamespaces(_orig,_tgt);
}

bool FormatConverter::setNamespaces(string _orig, string _tgt) {
	orig = _orig;
	tgt = _tgt;
	return isConversionSupported(orig,tgt);
}

FormatConverter::~FormatConverter() {
}


bool FormatConverter::isConversionSupported(string _orig, string _tgt) {
	if(_orig == "PDB2.3" && (_tgt == "CHARMM19" || _tgt == "CHARMM20" || _tgt == "CHARMM22" || _tgt == "CHARMM27" )) {
		return true;
	} else  if (_tgt == "PDB2.3" && (_orig == "CHARMM19" || _orig == "CHARMM20" || _orig == "CHARMM22" || _orig == "CHARMM27" )) {
		return true;
	} else {
		return false;
	}
}

string FormatConverter::getResidueName(string _resName, bool _protonatedOnD, bool _protonatedOnE ) {
	if(orig == "PDB2.3") {
		return getCharmmResName(_resName,tgt,_protonatedOnE,_protonatedOnE);
	} else if (orig == "CHARMM19" || orig == "CHARMM20" || orig == "CHARMM22" || orig == "CHARMM27") {
		return getPdbResName(_resName);
	} else {
		return _resName;
	}
}

string FormatConverter::getAtomName(string _atomName,string _resName, bool _NTerminal, bool _CTerminal) {
	if(orig == "PDB2.3") {
		//cout << "Converting pdb to charmm" << endl;
		return getCharmmAtomName(_atomName,_resName,tgt,_NTerminal,_CTerminal);
	} else if (orig == "CHARMM19" || orig == "CHARMM19" || orig == "CHARMM22" || orig == "CHARMM27") {
		//cout << "Converting charmm to pdb" << endl;
		return getPdbAtomName(_atomName,_resName,orig);
	} else {
		return _atomName;
	}
}

void FormatConverter::convert(Atom& _atom,bool _NTerminal,bool _CTerminal, bool _protonatedOnD , bool _protonatedOnE ) {
	//cout << "INP: " << _atom.getName() << " "  << _atom.getResidueName() << endl;
	string resName = getResidueName(_atom.getResidueName(),_protonatedOnE,_protonatedOnE);
	_atom.setResidueName(resName);
	_atom.setName(getAtomName(_atom.getName(),resName,_NTerminal,_CTerminal));
	//cout << "OUT: " << _atom.getName() << " "  << _atom.getResidueName() << endl;
}


void FormatConverter::convert(AtomPointerVector& _apV) {

	string oldChainId = "";
	string oldPosId = "";
	map<string,bool> nTer;
	map<string,bool> cTer;

	// collect the cTer and nTer posIds
	for(AtomPointerVector::iterator it = _apV.begin(); it != _apV.end(); it++) {
		string thisChainId = (*it)->getChainId();
		string thisPosId = (*it)->getPositionId();
		
		if(oldChainId != thisChainId ) {
			if(oldPosId != "") {
				cTer[oldPosId] = 1;
			}
			nTer[thisPosId] = 1;
		}
		oldPosId = thisPosId;
		oldChainId = thisChainId;
	}

	// convert all the atoms
	 
	for(AtomPointerVector::iterator it = _apV.begin(); it != _apV.end(); it++) {
		string thisPosId = (*it)->getPositionId();
		bool nTerm = false;
		bool cTerm = false;
		if(nTer.find(thisPosId) != nTer.end()) {
			nTerm = true;
		}
		if(cTer.find(thisPosId) != cTer.end()) {
			cTerm = true;
		}

		convert(**it,nTerm,cTerm);
	}

}


/**************************************************************************************************************
*PRIVATE FUNCTIONS
*
***************************************************************************************************************/


string FormatConverter::getCharmmResName(string _pdbResName, string _charmmVersion,bool _protonatedOnD, bool _protonatedOnE) {

	if(_pdbResName == "HIS") {
		if (_protonatedOnD && _protonatedOnE) {
			if (_charmmVersion == "CHARMM19" || _charmmVersion == "CHARMM20") {
				return "HSC";
			} else if (_charmmVersion == "CHARMM22" || _charmmVersion == "CHARMM27") {
				return "HSP";
			}
		} else if (_protonatedOnD) {
			if (_charmmVersion == "CHARMM19" || _charmmVersion == "CHARMM20") {
				return "HIS";
			} else if (_charmmVersion == "CHARMM22" || _charmmVersion == "CHARMM27") {
				return "HSD";
			}
		} else {
			if (_charmmVersion == "CHARMM19" || _charmmVersion == "CHARMM20") {
				return "HSD";
			} else if (_charmmVersion == "CHARMM22" || _charmmVersion == "CHARMM27") {
				return "HSE";
			}
		}
		if(_charmmVersion == "CHARMM22" || _charmmVersion == "CHARMM27") {
			return "HSE";
		} else if (_charmmVersion == "CHARMM19" || _charmmVersion == "CHARMM20") {
			return "HSD";
		}
	} else if (_pdbResName == "HOH") {
		return "TIP3";
	}

	return _pdbResName;
}

string FormatConverter::getCharmmAtomName(string _pdbName, string _resName, string _charmmVersion, bool _Nterminal, bool _Cterminal) {
	if (_charmmVersion == "CHARMM22" || _charmmVersion == "CHARMM27") {

		// WATER
		//HETATM 5765  O   HOH B 443       6.728 -16.138 -64.807  1.00 36.50      1HTB6081
		
		if (_resName == "HOH") {
			if (_pdbName == "O") {
				return("OH2");
				
			}
			if (_pdbName == "1H") {
				return("H1");
				
			}
			if (_pdbName == "2H") {
				return("H2");
				
			}
		}

		
		// CONVERT TERMINAL PATCHES
		if (_Nterminal) {
			if (_resName == "PRO") {
				/* STANDARD N-TERMINUS proline PROP */
				if (_pdbName == "1H") {
					return("HN1");
					
				}
				if (_pdbName == "2H") {
					return("HN2");
					
				}
			} else {
				/* STANDARD N-TERMINUS NTER and GLYP */
				if (_pdbName == "1H") {
					return("HT1");
					
				}
				if (_pdbName == "2H") {
					return("HT2");
					
				}
				if (_pdbName == "3H") {
					return("HT3");
					
				}
			}
			/* ACETYL GROUP ACE (or ACP for PRO) N-terminal patch
			 *
			 *  TO DO THINGS PROPERGLY WE'D HAVE TO SET IT AS A
			 *  PATCH AND NOT A SEPARATE RESIDUE IN CHARMM 22    */
			
			if (_resName == "ACE") {
				if (_pdbName == "CH3") {
					return("CAY");
					
				}
				if (_pdbName == "1H") {
					return("HY1");
					
				}
				if (_pdbName == "2H") {
					return("HY2");
					
				}
				if (_pdbName == "3H") {
					return("HY3");
					
				}
				if (_pdbName == "C") {
					return("CY");
					
				}
				if (_pdbName == "O") {
					return("OY");
					
				}
			}
		}
		if (_Cterminal) {
			/* Standard CTER */
			if (_pdbName == "O") {
				//return("OT1");
				return("OT2");
				
			}
			if (_pdbName == "OXT") {
				//return("OT2");
				return("OT1");
				
			}
		}

/*
HETATM    1  C   ACE A  22      -8.441  -0.600   8.995  1.00  0.00           C  
HETATM    2  O   ACE A  22      -7.625  -1.400   9.447  1.00  0.00           O  
HETATM    3  CH3 ACE A  22      -9.930  -0.806   9.243  1.00  0.00           C  
HETATM    4 1H   ACE A  22     -10.086  -1.705   9.839  1.00  0.00           H  
HETATM    5 2H   ACE A  22     -10.333   0.053   9.780  1.00  0.00           H  
HETATM    6 3H   ACE A  22     -10.448  -0.914   8.290  1.00  0.00           H  

PRES ACE          0.00 ! acetylated N-terminus
GROUP                  ! use in generate statement
ATOM CAY  CT3    -0.27 !
ATOM HY1  HA      0.09 ! HY1 HY2 HY3
ATOM HY2  HA      0.09 !    \ | /
ATOM HY3  HA      0.09 !     CAY
GROUP                  !      |
ATOM CY   C       0.51 !      CY=OY
ATOM OY   O      -0.51 !      |
!
BOND CY CAY CY N CAY HY1 CAY HY2 CAY HY3
DOUBLE OY CY  
IMPR CY CAY N OY    
IMPR N CY CA HN    
ACCEPTOR OY CY   
IC CY   N    CA   C     0.0000  0.0000  -60.0000  0.0000  0.0000
IC CY   CA   *N   HN    0.0000  0.0000  180.0000  0.0000  0.0000
IC CAY  CY   N    CA    0.0000  0.0000  180.0000  0.0000  0.0000
IC N    CAY  *CY  OY    0.0000  0.0000  180.0000  0.0000  0.0000
IC OY   CY   CAY  HY1   0.0000  0.0000  180.0000  0.0000  0.0000
IC OY   CY   CAY  HY2   0.0000  0.0000   60.0000  0.0000  0.0000
IC OY   CY   CAY  HY3   0.0000  0.0000  -60.0000  0.0000  0.0000



*/

		/************************************************
		 *  Make the necessary atom name changes for each residue
		 *  
		 *  NEED TO ADD NUCLEIC ACIDS!
		 ************************************************/
		if (_resName == "ALA") {
			/*
			   N    N   
			*  H    HN  
			   CA   CA  
			   HA   HA  
			   CB   CB  
			* 1HB   HB1 
			* 2HB   HB2 
			* 3HB   HB3 
			   C    C   
			   O    O   
			*/  
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			if (_pdbName == "3HB") {
				return("HB3");
				
			}
			if (_pdbName == "N" || _pdbName == "CA" || _pdbName == "HA" || _pdbName == "CB" || _pdbName == "C" || _pdbName == "O") {
				return(_pdbName);
				
			}
			return(_pdbName);
			
		}
		if (_resName == "CYS") {
			/*
			   N    N   
			*  H    HN  
			   CA   CA  
			   HA   HA  
			   CB   CB  
			* 1HB   HB1 
			* 2HB   HB2 
			   SG   SG  
			*  HG   HG1 
			   C    C   
			   O    O   
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			if (_pdbName == "HG") {
				return("HG1");
				
			}
			
		}
		if (_resName == "ASP") {
			/*
			   N     N   
			*  H     HN  
			   CA    CA  
			   HA    HA  
			   CB    CB  
			* 1HB    HB1 
			* 2HB    HB2 
			   CG    CG  
			   OD1   OD1 
			   OD2   OD2 
			   C     C   
			   O     O   
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			
		}
		if (_resName == "GLU") {
			/*
			   N     N   
			*  H     HN  
			   CA    CA  
			   HA    HA  
			   CB    CB  
			* 1HB    HB1 
			* 2HB    HB2 
			   CG    CG  
			* 1HG    HG1 
			* 2HG    HG2 
			   CD    CD  
			   OE1   OE1 
			   OE2   OE2 
			   C     C   
			   O     O   
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			if (_pdbName == "1HG") {
				return("HG1");
				
			}
			if (_pdbName == "2HG") {
				return("HG2");
				
			}
			
		}
		if (_resName == "PHE") {
			/*
			   N      N  
			*  H      HN 
			   CA     CA 
			   HA     HA 
			   CB     CB 
			* 1HB     HB1
			* 2HB     HB2
			   CG     CG 
			   CD1    CD1
			   HD1    HD1
			   CE1    CE1
			   HE1    HE1
			   CZ     CZ 
			   HZ     HZ 
			   CD2    CD2
			   HD2    HD2
			   CE2    CE2
			   HE2    HE2
			   C      C  
			   O      O  
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			
		}
		if (_resName == "GLY") {
			/*
			   N    N    
			*  H   HN   
			   CA   CA   
			* 1HA   HA1  
			* 2HA   HA2  
			   C    C    
			   O    O    
			*/       
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HA") {
				return("HA1");
				
			}
			if (_pdbName == "2HA") {
				return("HA2");
				
			}
			
		}
		if (_resName == "HSD" || _resName == "HSE" || _resName == "HSP" ) {
			/*
			   N      N  
			*  H      HN 
			   CA     CA 
			   HA     HA 
			   CB     CB 
			* 1HB     HB1
			* 2HB     HB2
			   ND1    ND1
			   HD1    HD1
			   CG     CG 
			   CE1    CE1
			   HE1    HE1
			   NE2    NE2
			   HE2    HE2
			   CD2    CD2
			   HD2    HD2
			   C      C  
			   O      O  
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			
		}
		if (_resName == "ILE") {
			/*
			   N      N    
			*  H      HN   
			   CA     CA   
			   HA     HA   
			   CB     CB   
			   HB     HB   
			   CG2    CG2  
			* 1HG2    HG21 
			* 2HG2    HG22 
			* 3HG2    HG23 
			   CG1    CG1  
			* 1HG1    HG11 
			* 2HG1    HG12 
			*  CD1    CD  
			* 1HD1    HD1  
			* 2HD1    HD2  
			* 3HD1    HD3  
			   C      C    
			   O      O    
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HG1") {
				return("HG11");
				
			}
			if (_pdbName == "2HG1") {
				return("HG12");
				
			}
			if (_pdbName == "1HG2") {
				return("HG21");
				
			}
			if (_pdbName == "2HG2") {
				return("HG22");
				
			}
			if (_pdbName == "3HG2") {
				return("HG23");
				
			}
			if (_pdbName == "1HD1") {
				return("HD1");
				
			}
			if (_pdbName == "2HD1") {
				return("HD2");
				
			}
			if (_pdbName == "3HD1") {
				return("HD3");
				
			}
			if (_pdbName == "CD1") {
				return("CD");
				
			}
			
		}
		if (_resName == "LYS") {
			/*
			    N    N   
			*   H    HN  
			    CA   CA  
			    HA   HA  
			    CB   CB  
			*  1HB   HB1 
			*  2HB   HB2 
			    CG   CG  
			*  1HG   HG1 
			*  2HG   HG2 
			    CD   CD  
			*  1HD   HD1 
			*  2HD   HD2 
			    CE   CE  
			*  1HE   HE1 
			*  2HE   HE2 
			    NZ   NZ  
			*  1HZ   HZ1 
			*  2HZ   HZ2 
			*  3HZ   HZ3 
			    C    C   
			    O    O   
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			if (_pdbName == "1HG") {
				return("HG1");
				
			}
			if (_pdbName == "2HG") {
				return("HG2");
				
			}
			if (_pdbName == "1HD") {
				return("HD1");
				
			}
			if (_pdbName == "2HD") {
				return("HD2");
				
			}
			if (_pdbName == "1HE") {
				return("HE1");
				
			}
			if (_pdbName == "2HE") {
				return("HE2");
				
			}
			if (_pdbName == "1HZ") {
				return("HZ1");
				
			}
			if (_pdbName == "2HZ") {
				return("HZ2");
				
			}
			if (_pdbName == "3HZ") {
				return("HZ3");
				
			}
			
		}
		if (_resName == "LEU") {
			/*
			   N      N    
			*  H      HN   
			   CA     CA   
			   HA     HA   
			   CB     CB   
			* 1HB     HB1  
			* 2HB     HB2  
			   CG     CG   
			   HG     HG   
			   CD1    CD1  
			* 1HD1    HD11 
			* 2HD1    HD12 
			* 3HD1    HD13 
			   CD2    CD2  
			* 1HD2    HD21 
			* 2HD2    HD22 
			* 3HD2    HD23 
			   C      C    
			   O      O    
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			if (_pdbName == "1HD1") {
				return("HD11");
				
			}
			if (_pdbName == "2HD1") {
				return("HD12");
				
			}
			if (_pdbName == "3HD1") {
				return("HD13");
				
			}
			if (_pdbName == "1HD2") {
				return("HD21");
				
			}
			if (_pdbName == "2HD2") {
				return("HD22");
				
			}
			if (_pdbName == "3HD2") {
				return("HD23");
				
			}
			
		}
		if (_resName == "MET") {
			/*
			   N      N   
			*  H      HN  
			   CA     CA  
			   HA     HA  
			   CB     CB  
			* 1HB     HB1 
			* 2HB     HB2 
			   CG     CG  
			* 1HG     HG1 
			* 2HG     HG2 
			   SD     SD  
			   CE     CE  
			* 1HE     HE1 
			* 2HE     HE2 
			* 3HE     HE3 
			   C      C   
			   O      O   
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			if (_pdbName == "1HG") {
				return("HG1");
				
			}
			if (_pdbName == "2HG") {
				return("HG2");
				
			}
			if (_pdbName == "1HE") {
				return("HE1");
				
			}
			if (_pdbName == "2HE") {
				return("HE2");
				
			}
			if (_pdbName == "3HE") {
				return("HE3");
				
			}
			
		}
		if (_resName == "ASN") {
			/*	         
			  NOTE: H numbering is inverted for terminal Hs

			   N     N   
			*  H     HN  
			   CA    CA  
			   HA    HA  
			   CB    CB  
			* 1HB    HB1 
			* 2HB    HB2 
			   CG    CG  
			   OD1   OD1 
			   ND2   ND2 
			* 1HD2   HD22
			* 2HD2   HD21
			   C     C   
			   O     O   
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			if (_pdbName == "2HD2") {
				return("HD21");
				
			}
			if (_pdbName == "1HD2") {
				return("HD22");
				
			}
			
		}
		if (_resName == "PRO") {
			/*
			  NOTE: H numbering is inverted

			  N     N  
			  CD    CD 
			 *1HD   HD2
			 *2HD   HD1
			  CA    CA 
			  HA    HA 
			  CB    CB 
			 *1HB   HB2
			 *2HB   HB1
			  CG    CG 
			 *1HG   HG2
			 *2HG   HG1
			  C     C  
			  O     O  
			*/
			if (_pdbName == "2HB") {
				return("HB1");
				
			}
			if (_pdbName == "1HB") {
				return("HB2");
				
			}
			if (_pdbName == "2HG") {
				return("HG1");
				
			}
			if (_pdbName == "1HG") {
				return("HG2");
				
			}
			if (_pdbName == "2HD") {
				return("HD1");
				
			}
			if (_pdbName == "1HD") {
				return("HD2");
				
			}
			
		}
		if (_resName == "GLN") {
			/*
			  NOTE: H numbering is inverted for terminal Hs

			   N    N    
			*  H    HN   
			   CA   CA   
			   HA   HA   
			   CB   CB   
			* 1HB   HB1  
			* 2HB   HB2  
			   CG   CG   
			* 1HG   HG1  
			* 2HG   HG2  
			   CD   CD   
			   OE1  OE1  
			   NE2  NE2  
			* 1HE2  HE22 
			* 2HE2  HE21 
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			if (_pdbName == "1HG") {
				return("HG1");
				
			}
			if (_pdbName == "2HG") {
				return("HG2");
				
			}
			if (_pdbName == "2HE2") {
				return("HE21");
				
			}
			if (_pdbName == "1HE2") {
				return("HE22");
				
			}
			
		}
		if (_resName == "ARG") {
			/*
			  NOTE: H numbering is inverted for terminal Hs

			   N     N   
			*  H     HN  
			   CA    CA  
			   HA    HA  
			   CB    CB  
			* 1HB    HB1 
			* 2HB    HB2 
			   CG    CG  
			* 1HG    HG1 
			* 2HG    HG2 
			   CD    CD  
			* 1HD    HD1 
			* 2HD    HD2 
			   NE    NE  
			   HE    HE  
			   CZ    CZ  
			   NH1   NH1 
			* 1HH1   HH12
			* 2HH1   HH11
			   NH2   NH2 
			* 1HH2   HH22
			* 2HH2   HH21
			   C     C   
			   O     O   
			*/  
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			if (_pdbName == "1HG") {
				return("HG1");
				
			}
			if (_pdbName == "2HG") {
				return("HG2");
				
			}
			if (_pdbName == "1HD") {
				return("HD1");
				
			}
			if (_pdbName == "2HD") {
				return("HD2");
				
			}
			if (_pdbName == "2HH1") {
				return("HH11");
				
			}
			if (_pdbName == "1HH1") {
				return("HH12");
				
			}
			if (_pdbName == "2HH2") {
				return("HH21");
				
			}
			if (_pdbName == "1HH2") {
				return("HH22");
				
			}
			
		}
		if (_resName == "SER") {
			/*
			   N    N   
			*  H    HN  
			   CA   CA  
			   HA   HA  
			   CB   CB  
			* 1HB   HB1 
			* 2HB   HB2 
			   OG   OG  
			*  HG   HG1 
			   C    C   
			   O    O   
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			if (_pdbName == "HG") {
				return("HG1");
				
			}
			
		}
		if (_resName == "THR") {
			/*
			   N      N    
			*  H      HN   
			   CA     CA   
			   HA     HA   
			   CB     CB   
			   HB     HB   
			   OG1    OG1  
			   HG1    HG1  
			   CG2    CG2  
			* 1HG2    HG21 
			* 2HG2    HG22 
			* 3HG2    HG23 
			   C      C    
			   O      O    
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HG2") {
				return("HG21");
				
			}
			if (_pdbName == "2HG2") {
				return("HG22");
				
			}
			if (_pdbName == "3HG2") {
				return("HG23");
				
			}
			
		}
		if (_resName == "VAL") {
			/*
			   N      N    
			*  H      HN   
			   CA     CA   
			   HA     HA   
			   CB     CB   
			   HB     HB   
			   CG1    CG1  
			* 1HG1    HG11 
			* 2HG1    HG12 
			* 3HG1    HG13 
			   CG2    CG2  
			* 1HG2    HG21 
			* 2HG2    HG22 
			* 3HG2    HG23 
			   C      C    
			   O      O    
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HG1") {
				return("HG11");
				
			}
			if (_pdbName == "2HG1") {
				return("HG12");
				
			}
			if (_pdbName == "3HG1") {
				return("HG13");
				
			}
			if (_pdbName == "1HG2") {
				return("HG21");
				
			}
			if (_pdbName == "2HG2") {
				return("HG22");
				
			}
			if (_pdbName == "3HG2") {
				return("HG23");
				
			}
			
		}
		if (_resName == "TRP") {
			/*
			   N      N  
			*  H      HN 
			   CA     CA 
			   HA     HA 
			   CB     CB 
			* 1HB     HB1
			* 2HB     HB2
			   CG     CG 
			   CD1    CD1
			   HD1    HD1
			   NE1    NE1
			   HE1    HE1
			   CE2    CE2
			   CD2    CD2
			   CE3    CE3
			   HE3    HE3
			   CZ3    CZ3
			   HZ3    HZ3
			   CZ2    CZ2
			   HZ2    HZ2
			   CH2    CH2
			   HH2    HH2
			   C      C  
			   O      O  
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			
		}
		if (_resName == "TYR") {
			/*
			   N     N   
			*  H    HN  
			   CA    CA  
			   HA    HA  
			   CB    CB  
			* 1HB    HB1 
			* 2HB    HB2 
			   CG    CG  
			   CD1   CD1 
			   HD1   HD1 
			   CE1   CE1 
			   HE1   HE1 
			   CZ    CZ  
			   OH    OH  
			   HH    HH  
			   CD2   CD2 
			   HD2   HD2 
			   CE2   CE2 
			   HE2   HE2 
			   C     C   
			   O     O   
			*/
			if (_pdbName == "H") {
				return("HN");
				
			}
			if (_pdbName == "1HB") {
				return("HB1");
				
			}
			if (_pdbName == "2HB") {
				return("HB2");
				
			}
			
		}

		/************************************************
		 *  Done making atom name changes
		 ************************************************/
	}

	if (_charmmVersion == "CHARMM19" || _charmmVersion == "CHARMM20" ) {

		// WATER
		if (_resName == "HOH") {
			if (_pdbName == "O") {
				return("OH2");
				
			}
			if (_pdbName == "1H") {
				return("H1");
				
			}
			if (_pdbName == "2H") {
				return("H2");
				
			}
		}

		// CONVERT TERMINAL PATCHES
		if (_Nterminal) {
			/* STANDARD N-TERMINUS NTER, GLYP */
			if (_resName == "PRO") {
				/* STANDARD N-TERMINUS proline PROP */
				if (_pdbName == "H1") {
					return("HN1");
					
				}
				if (_pdbName == "H2") {
					return("HN2");
					
				}
			} else {

				if (_pdbName == "H1") {
					return("HT1");
					
				}
				if (_pdbName == "H2") {
					return("HT2");
					
				}
				if (_pdbName == "H3") {
					return("HT3");
					
				}
			}
			/* ACETYL GROUP ACE (or ACP for PRO) N-terminal patch*/
			/* ACE IS A RESIDUE IN CHARMM 19, NOT A PATCH AS IN CHARMM 22 */
		}
		if (_Cterminal) {
			/* Standard CTER */
			if (_pdbName == "O") {
				//return("OT1");
				return("OT2");
				
			}
			if (_pdbName == "OXT") {
				//return("OT2");
				return("OT1");
				
			}
		}


		/************************************************
		 *  Make the necessary atom name changes for each residue
		 *  
		 *  NEED TO ADD NUCLEIC ACIDS!
		 ************************************************/
		
		if (_resName == "ILE") {
			/*
			   N      N    
			   H      H   
			   CA     CA   
			   CB     CB   
			   CG2    CG2  
			   CG1    CG1  
			   CD1    CD  
			   C      C    
			   O      O    
			*/
			if (_pdbName == "CD1") {
				return("CD");
				
			}
			
		}
		if (_resName == "LYS") {
			/*
			    N    N   
			    H    H  
			    CA   CA  
			    CB   CB  
			    CG   CG  
			    CD   CD  
			    CE   CE  
			    NZ   NZ  
			*  1HZ   HZ1 
			*  2HZ   HZ2 
			*  3HZ   HZ3 
			    C    C   
			    O    O   
			*/
			if (_pdbName == "1HZ") {
				return("HZ1");
				
			}
			if (_pdbName == "2HZ") {
				return("HZ2");
				
			}
			if (_pdbName == "3HZ") {
				return("HZ3");
				
			}
			
		}
		if (_resName == "ASN") {
			/*	         
			  NOTE: H numbering is inverted for terminal Hs

			   N     N   
			   H     H   
			   CA    CA  
			   CB    CB  
			   CG    CG  
			   OD1   OD1 
			   ND2   ND2 
			* 1HD2   HD22
			* 2HD2   HD21
			   C     C   
			   O     O   
			*/
			if (_pdbName == "1HD2") {
				return("HD22");
				
			}
			if (_pdbName == "2HD2") {
				return("HD21");
				
			}
			
		}
		
		if (_resName == "GLN") {
			/*
			  NOTE: H numbering is inverted for terminal Hs

			   N    N    
			   H    H   
			   CA   CA   
			   CB   CB   
			   CG   CG   
			   CD   CD   
			   OE1  OE1  
			   NE2  NE2  
			* 1HE2  HE22 
			* 2HE2  HE21 
			*/
			if (_pdbName == "2HE2") {
				return("HE21");
				
			}
			if (_pdbName == "1HE2") {
				return("HE22");
				
			}
			
		}
		if (_resName == "ARG") {
			/*
			  NOTE: H numbering is inverted for terminal Hs

			   N     N   
			   H     H   
			   CA    CA  
			   CB    CB  
			   CG    CG  
			   CD    CD  
			   NE    NE  
			   HE    HE  
			   CZ    CZ  
			   NH1   NH1 
			* 1HH1   HH12
			* 2HH1   HH11
			   NH2   NH2 
			* 1HH2   HH22
			* 2HH2   HH21
			   C     C   
			   O     O   
			*/  
			if (_pdbName == "2HH1") {
				return("HH11");
				
			}
			if (_pdbName == "1HH1") {
				return("HH12");
				
			}
			if (_pdbName == "2HH2") {
				return("HH21");
				
			}
			if (_pdbName == "1HH2") {
				return("HH22");
				
			}
			
		}
		if (_resName == "SER") {
			/*
			   N    N   
			   H    H   
			   CA   CA  
			   CB   CB  
			   OG   OG  
			*  HG   HG1 
			   C    C   
			   O    O   
			*/
			if (_pdbName == "HG") {
				return("HG1");
				
			}
			
		}
	}
	return _pdbName;
}




/**************************************************************************************************************/

string FormatConverter::getPdbResName(string _charmmResName) {
	if (_charmmResName == "HSD" || _charmmResName == "HSE" || _charmmResName == "HSP" || _charmmResName == "HSC") {
		return "HIS";
	} else if (_charmmResName == "TIP3") {
		return "HOH";
	} else {
		if (_charmmResName.size() > 3) {
				cerr << "WARNING 2912: charmm residue name " << _charmmResName << " too long for converting to pdb residue name: truncating it to the first 3 characters " << _charmmResName.substr(0, 3) << endl;
				_charmmResName = _charmmResName.substr(0,3);
		}
		return _charmmResName;
	}
}

string FormatConverter::getPdbAtomName(string _charmmName, string _resName, string _charmmVersion) {
	 
	if (_charmmVersion == "CHARMM22" || _charmmVersion == "CHARMM27") {

		if (_resName == "TIP3") {
			if (_charmmName == "OH2") {
				return("O");
				
			}
			if (_charmmName == "H1") {
				return("1H");
				
			}
			if (_charmmName == "H2") {
				return("2H");
				
			}
		}
		
		// CONVERT TERMINAL PATCHES
		/* STANDARD N-TERMINUS NTER and GLYP */
		if (_charmmName == "HT1") {
			return("1H");
			
		}
		if (_charmmName == "HT2") {
			return("2H");
			
		}
		if (_charmmName == "HT3") {
			return("3H");
			
		}
		/* STANDARD N-TERMINUS proline PROP */
		if (_charmmName == "HN1") {
			return("1H");
			
		}
		if (_charmmName == "HN2") {
			return("2H");
			
		}
		/* ACETYL GROUP ACE (or ACP for PRO) N-terminal patch*/
		if (_charmmName == "CAY") {
			return("CH3");
			/* What to do for the patch? 
			_rAtom.setAltPdbResName("ACE");
			_rAtom.setAltPdbResnum(_rAtom.getPdbResnum() - 1);
			_rAtom.setAltIcode("");
			_rAtom.setHetAtom(true);
			*/
			
		}
		if (_charmmName == "HY1") {
			return("1H");
			
		}
		if (_charmmName == "HY2") {
			return("2H");
			
		}
		if (_charmmName == "HY3") {
			return("3H");
			
		}
		if (_charmmName == "CY") {
			return("C");
			
		}
		if (_charmmName == "OY") {
			return("O");
			
		}
		/* Standard CTER */
		if (_charmmName == "OT1") {
			//return("O");
			return("OXT");
			
		}
		if (_charmmName == "OT2") {
			//return("OXT");
			return("O");
			
		}

/*
HETATM    1  C   ACE A  22      -8.441  -0.600   8.995  1.00  0.00           C  
HETATM    2  O   ACE A  22      -7.625  -1.400   9.447  1.00  0.00           O  
HETATM    3  CH3 ACE A  22      -9.930  -0.806   9.243  1.00  0.00           C  
HETATM    4 1H   ACE A  22     -10.086  -1.705   9.839  1.00  0.00           H  
HETATM    5 2H   ACE A  22     -10.333   0.053   9.780  1.00  0.00           H  
HETATM    6 3H   ACE A  22     -10.448  -0.914   8.290  1.00  0.00           H  

PRES ACE          0.00 ! acetylated N-terminus
GROUP                  ! use in generate statement
ATOM CAY  CT3    -0.27 !
ATOM HY1  HA      0.09 ! HY1 HY2 HY3
ATOM HY2  HA      0.09 !    \ | /
ATOM HY3  HA      0.09 !     CAY
GROUP                  !      |
ATOM CY   C       0.51 !      CY=OY
ATOM OY   O      -0.51 !      |
!
BOND CY CAY CY N CAY HY1 CAY HY2 CAY HY3
DOUBLE OY CY  
IMPR CY CAY N OY    
IMPR N CY CA HN    
ACCEPTOR OY CY   
IC CY   N    CA   C     0.0000  0.0000  -60.0000  0.0000  0.0000
IC CY   CA   *N   HN    0.0000  0.0000  180.0000  0.0000  0.0000
IC CAY  CY   N    CA    0.0000  0.0000  180.0000  0.0000  0.0000
IC N    CAY  *CY  OY    0.0000  0.0000  180.0000  0.0000  0.0000
IC OY   CY   CAY  HY1   0.0000  0.0000  180.0000  0.0000  0.0000
IC OY   CY   CAY  HY2   0.0000  0.0000   60.0000  0.0000  0.0000
IC OY   CY   CAY  HY3   0.0000  0.0000  -60.0000  0.0000  0.0000



*/

		/************************************************
		 *  Make the necessary atom name changes for each residue
		 *  
		 *  NEED TO ADD NUCLEIC ACIDS!
		 ************************************************/
		if (_resName == "ALA") {
			/*
			   N    N   
			*  H    HN  
			   CA   CA  
			   HA   HA  
			   CB   CB  
			* 1HB   HB1 
			* 2HB   HB2 
			* 3HB   HB3 
			   C    C   
			   O    O   
			*/  
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			if (_charmmName == "HB3") {
				return("3HB");
				
			}
			
		}
		if (_resName == "CYS") {
			/*
			   N    N   
			*  H    HN  
			   CA   CA  
			   HA   HA  
			   CB   CB  
			* 1HB   HB1 
			* 2HB   HB2 
			   SG   SG  
			*  HG   HG1 
			   C    C   
			   O    O   
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			if (_charmmName == "HG1") {
				return("HG");
				
			}
			
		}
		if (_resName == "ASP") {
			/*
			   N     N   
			*  H     HN  
			   CA    CA  
			   HA    HA  
			   CB    CB  
			* 1HB    HB1 
			* 2HB    HB2 
			   CG    CG  
			   OD1   OD1 
			   OD2   OD2 
			   C     C   
			   O     O   
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			
		}
		if (_resName == "GLU") {
			/*
			   N     N   
			*  H     HN  
			   CA    CA  
			   HA    HA  
			   CB    CB  
			* 1HB    HB1 
			* 2HB    HB2 
			   CG    CG  
			* 1HG    HG1 
			* 2HG    HG2 
			   CD    CD  
			   OE1   OE1 
			   OE2   OE2 
			   C     C   
			   O     O   
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			if (_charmmName == "HG1") {
				return("1HG");
				
			}
			if (_charmmName == "HG2") {
				return("2HG");
				
			}
			
		}
		if (_resName == "PHE") {
			/*
			   N      N  
			*  H      HN 
			   CA     CA 
			   HA     HA 
			   CB     CB 
			* 1HB     HB1
			* 2HB     HB2
			   CG     CG 
			   CD1    CD1
			   HD1    HD1
			   CE1    CE1
			   HE1    HE1
			   CZ     CZ 
			   HZ     HZ 
			   CD2    CD2
			   HD2    HD2
			   CE2    CE2
			   HE2    HE2
			   C      C  
			   O      O  
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			
		}
		if (_resName == "GLY") {
			/*
			   N    N    
			*  H   HN   
			   CA   CA   
			* 1HA   HA1  
			* 2HA   HA2  
			   C    C    
			   O    O    
			*/       
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HA1") {
				return("1HA");
				
			}
			if (_charmmName == "HA2") {
				return("2HA");
				
			}
			
		}
		if (_resName == "HIS") {
			/*
			   N      N  
			*  H      HN 
			   CA     CA 
			   HA     HA 
			   CB     CB 
			* 1HB     HB1
			* 2HB     HB2
			   ND1    ND1
			   HD1    HD1
			   CG     CG 
			   CE1    CE1
			   HE1    HE1
			   NE2    NE2
			   HE2    HE2
			   CD2    CD2
			   HD2    HD2
			   C      C  
			   O      O  
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			
		}
		if (_resName == "ILE") {
			/*
			   N      N    
			*  H      HN   
			   CA     CA   
			   HA     HA   
			   CB     CB   
			   HB     HB   
			   CG2    CG2  
			* 1HG2    HG21 
			* 2HG2    HG22 
			* 3HG2    HG23 
			   CG1    CG1  
			* 1HG1    HG11 
			* 2HG1    HG12 
			*  CD1    CD  
			* 1HD1    HD1  
			* 2HD1    HD2  
			* 3HD1    HD3  
			   C      C    
			   O      O    
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HG11") {
				return("1HG1");
				
			}
			if (_charmmName == "HG12") {
				return("2HG1");
				
			}
			if (_charmmName == "HG21") {
				return("1HG2");
				
			}
			if (_charmmName == "HG22") {
				return("2HG2");
				
			}
			if (_charmmName == "HG23") {
				return("3HG2");
				
			}
			if (_charmmName == "HD1") {
				return("1HD1");
				
			}
			if (_charmmName == "HD2") {
				return("2HD1");
				
			}
			if (_charmmName == "HD3") {
				return("3HD1");
				
			}
			if (_charmmName == "CD") {
				return("CD1");
				
			}
			
		}
		if (_resName == "LYS") {
			/*
			    N    N   
			*   H    HN  
			    CA   CA  
			    HA   HA  
			    CB   CB  
			*  1HB   HB1 
			*  2HB   HB2 
			    CG   CG  
			*  1HG   HG1 
			*  2HG   HG2 
			    CD   CD  
			*  1HD   HD1 
			*  2HD   HD2 
			    CE   CE  
			*  1HE   HE1 
			*  2HE   HE2 
			    NZ   NZ  
			*  1HZ   HZ1 
			*  2HZ   HZ2 
			*  3HZ   HZ3 
			    C    C   
			    O    O   
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			if (_charmmName == "HG1") {
				return("1HG");
				
			}
			if (_charmmName == "HG2") {
				return("2HG");
				
			}
			if (_charmmName == "HD1") {
				return("1HD");
				
			}
			if (_charmmName == "HD2") {
				return("2HD");
				
			}
			if (_charmmName == "HE1") {
				return("1HE");
				
			}
			if (_charmmName == "HE2") {
				return("2HE");
				
			}
			if (_charmmName == "HZ1") {
				return("1HZ");
				
			}
			if (_charmmName == "HZ2") {
				return("2HZ");
				
			}
			if (_charmmName == "HZ3") {
				return("3HZ");
				
			}
			
		}
		if (_resName == "LEU") {
			/*
			   N      N    
			*  H      HN   
			   CA     CA   
			   HA     HA   
			   CB     CB   
			* 1HB     HB1  
			* 2HB     HB2  
			   CG     CG   
			   HG     HG   
			   CD1    CD1  
			* 1HD1    HD11 
			* 2HD1    HD12 
			* 3HD1    HD13 
			   CD2    CD2  
			* 1HD2    HD21 
			* 2HD2    HD22 
			* 3HD2    HD23 
			   C      C    
			   O      O    
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			if (_charmmName == "HD11") {
				return("1HD1");
				
			}
			if (_charmmName == "HD12") {
				return("2HD1");
				
			}
			if (_charmmName == "HD13") {
				return("3HD1");
				
			}
			if (_charmmName == "HD21") {
				return("1HD2");
				
			}
			if (_charmmName == "HD22") {
				return("2HD2");
				
			}
			if (_charmmName == "HD23") {
				return("3HD2");
				
			}
			
		}
		if (_resName == "MET") {
			/*
			   N      N   
			*  H      HN  
			   CA     CA  
			   HA     HA  
			   CB     CB  
			* 1HB     HB1 
			* 2HB     HB2 
			   CG     CG  
			* 1HG     HG1 
			* 2HG     HG2 
			   SD     SD  
			   CE     CE  
			* 1HE     HE1 
			* 2HE     HE2 
			* 3HE     HE3 
			   C      C   
			   O      O   
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			if (_charmmName == "HG1") {
				return("1HG");
				
			}
			if (_charmmName == "HG2") {
				return("2HG");
				
			}
			if (_charmmName == "HE1") {
				return("1HE");
				
			}
			if (_charmmName == "HE2") {
				return("2HE");
				
			}
			if (_charmmName == "HE3") {
				return("3HE");
				
			}
			
		}
		if (_resName == "ASN") {
			/*	         
			  NOTE: H numbering is inverted for terminal Hs

			   N     N   
			*  H     HN  
			   CA    CA  
			   HA    HA  
			   CB    CB  
			* 1HB    HB1 
			* 2HB    HB2 
			   CG    CG  
			   OD1   OD1 
			   ND2   ND2 
			* 1HD2   HD22
			* 2HD2   HD21
			   C     C   
			   O     O   
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			if (_charmmName == "HD21") {
				return("2HD2");
				
			}
			if (_charmmName == "HD22") {
				return("1HD2");
				
			}
			
		}
		if (_resName == "PRO") {
			/*
			  NOTE: H numbering is inverted

			  N     N  
			  CD    CD 
			 *1HD   HD2
			 *2HD   HD1
			  CA    CA 
			  HA    HA 
			  CB    CB 
			 *1HB   HB2
			 *2HB   HB1
			  CG    CG 
			 *1HG   HG2
			 *2HG   HG1
			  C     C  
			  O     O  
			*/
			if (_charmmName == "HB1") {
				return("2HB");
				
			}
			if (_charmmName == "HB2") {
				return("1HB");
				
			}
			if (_charmmName == "HG1") {
				return("2HG");
				
			}
			if (_charmmName == "HG2") {
				return("1HG");
				
			}
			if (_charmmName == "HD1") {
				return("2HD");
				
			}
			if (_charmmName == "HD2") {
				return("1HD");
				
			}
			
		}
		if (_resName == "GLN") {
			/*
			  NOTE: H numbering is inverted for terminal Hs

			   N    N    
			*  H    HN   
			   CA   CA   
			   HA   HA   
			   CB   CB   
			* 1HB   HB1  
			* 2HB   HB2  
			   CG   CG   
			* 1HG   HG1  
			* 2HG   HG2  
			   CD   CD   
			   OE1  OE1  
			   NE2  NE2  
			* 1HE2  HE22 
			* 2HE2  HE21 
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			if (_charmmName == "HG1") {
				return("1HG");
				
			}
			if (_charmmName == "HG2") {
				return("2HG");
				
			}
			if (_charmmName == "HE21") {
				return("2HE2");
				
			}
			if (_charmmName == "HE22") {
				return("1HE2");
				
			}
			
		}
		if (_resName == "ARG") {
			/*
			  NOTE: H numbering is inverted for terminal Hs

			   N     N   
			*  H     HN  
			   CA    CA  
			   HA    HA  
			   CB    CB  
			* 1HB    HB1 
			* 2HB    HB2 
			   CG    CG  
			* 1HG    HG1 
			* 2HG    HG2 
			   CD    CD  
			* 1HD    HD1 
			* 2HD    HD2 
			   NE    NE  
			   HE    HE  
			   CZ    CZ  
			   NH1   NH1 
			* 1HH1   HH12
			* 2HH1   HH11
			   NH2   NH2 
			* 1HH2   HH22
			* 2HH2   HH21
			   C     C   
			   O     O   
			*/  
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			if (_charmmName == "HG1") {
				return("1HG");
				
			}
			if (_charmmName == "HG2") {
				return("2HG");
				
			}
			if (_charmmName == "HD1") {
				return("1HD");
				
			}
			if (_charmmName == "HD2") {
				return("2HD");
				
			}
			if (_charmmName == "HH11") {
				return("2HH1");
				
			}
			if (_charmmName == "HH12") {
				return("1HH1");
				
			}
			if (_charmmName == "HH21") {
				return("2HH2");
				
			}
			if (_charmmName == "HH22") {
				return("1HH2");
				
			}
			
		}
		if (_resName == "SER") {
			/*
			   N    N   
			*  H    HN  
			   CA   CA  
			   HA   HA  
			   CB   CB  
			* 1HB   HB1 
			* 2HB   HB2 
			   OG   OG  
			*  HG   HG1 
			   C    C   
			   O    O   
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			if (_charmmName == "HG1") {
				return("HG");
				
			}
			
		}
		if (_resName == "THR") {
			/*
			   N      N    
			*  H      HN   
			   CA     CA   
			   HA     HA   
			   CB     CB   
			   HB     HB   
			   OG1    OG1  
			   HG1    HG1  
			   CG2    CG2  
			* 1HG2    HG21 
			* 2HG2    HG22 
			* 3HG2    HG23 
			   C      C    
			   O      O    
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HG21") {
				return("1HG2");
				
			}
			if (_charmmName == "HG22") {
				return("2HG2");
				
			}
			if (_charmmName == "HG23") {
				return("3HG2");
				
			}
			
		}
		if (_resName == "VAL") {
			/*
			   N      N    
			*  H      HN   
			   CA     CA   
			   HA     HA   
			   CB     CB   
			   HB     HB   
			   CG1    CG1  
			* 1HG1    HG11 
			* 2HG1    HG12 
			* 3HG1    HG13 
			   CG2    CG2  
			* 1HG2    HG21 
			* 2HG2    HG22 
			* 3HG2    HG23 
			   C      C    
			   O      O    
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HG11") {
				return("1HG1");
				
			}
			if (_charmmName == "HG12") {
				return("2HG1");
				
			}
			if (_charmmName == "HG13") {
				return("3HG1");
				
			}
			if (_charmmName == "HG21") {
				return("1HG2");
				
			}
			if (_charmmName == "HG22") {
				return("2HG2");
				
			}
			if (_charmmName == "HG23") {
				return("3HG2");
				
			}
			
		}
		if (_resName == "TRP") {
			/*
			   N      N  
			*  H      HN 
			   CA     CA 
			   HA     HA 
			   CB     CB 
			* 1HB     HB1
			* 2HB     HB2
			   CG     CG 
			   CD1    CD1
			   HD1    HD1
			   NE1    NE1
			   HE1    HE1
			   CE2    CE2
			   CD2    CD2
			   CE3    CE3
			   HE3    HE3
			   CZ3    CZ3
			   HZ3    HZ3
			   CZ2    CZ2
			   HZ2    HZ2
			   CH2    CH2
			   HH2    HH2
			   C      C  
			   O      O  
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			
		}
		if (_resName == "TYR") {
			/*
			   N     N   
			*  H    HN  
			   CA    CA  
			   HA    HA  
			   CB    CB  
			* 1HB    HB1 
			* 2HB    HB2 
			   CG    CG  
			   CD1   CD1 
			   HD1   HD1 
			   CE1   CE1 
			   HE1   HE1 
			   CZ    CZ  
			   OH    OH  
			   HH    HH  
			   CD2   CD2 
			   HD2   HD2 
			   CE2   CE2 
			   HE2   HE2 
			   C     C   
			   O     O   
			*/
			if (_charmmName == "HN") {
				return("H");
				
			}
			if (_charmmName == "HB1") {
				return("1HB");
				
			}
			if (_charmmName == "HB2") {
				return("2HB");
				
			}
			
		}

		/************************************************
		 *  Done making atom name changes
		 ************************************************/
	}
	if (_charmmVersion == "CHARMM19" || _charmmVersion == "CHARMM20" ) {
		// CONVERT TERMINAL PATCHES
		/* STANDARD N-TERMINUS NTER, GLYP and PROP */
		if (_charmmName == "HT1") {
			return("H1");
			
		}
		if (_charmmName == "HT2") {
			return("H2");
			
		}
		if (_charmmName == "HT3") {
			return("H3");
			
		}
		/* STANDARD N-TERMINUS proline PROP */
		if (_charmmName == "HN1") {
			return("H1");
			
		}
		if (_charmmName == "HN2") {
			return("H2");
			
		}
		/* ACETYL GROUP ACE (or ACP for PRO) N-terminal patch*/
		/* ACE IS A RESIDUE IN CHARMM 19, NOT A PATCH AS IN CHARMM 22 */
		/*
		if (_resName == "ACE") {
			   CH3  CH3 
			   C    C   
			   O    O   
			if (_charmmName == "CH3" || _charmmName == "C" || _charmmName == "O") {
				return(_charmmName);
				
			}
			return(_charmmName);
			
		}
		*/  
		/* Standard CTER */
		if (_charmmName == "OT1") {
			//return("O");
			return("OXT");
			
		}
		if (_charmmName == "OT2") {
			//return("OXT");
			return("O");
			
		}


		/************************************************
		 *  Make the necessary atom name changes for each residue
		 *  
		 *  NEED TO ADD NUCLEIC ACIDS!
		 ************************************************/
		if (_resName == "ILE") {
			/*
			   N      N    
			   H      H   
			   CA     CA   
			   CB     CB   
			   CG2    CG2  
			   CG1    CG1  
			   CD1    CD  
			   C      C    
			   O      O    
			*/
			if (_charmmName == "CD") {
				return("CD1");
				
			}
			
		}
		if (_resName == "LYS") {
			/*
			    N    N   
			    H    H  
			    CA   CA  
			    CB   CB  
			    CG   CG  
			    CD   CD  
			    CE   CE  
			    NZ   NZ  
			*  1HZ   HZ1 
			*  2HZ   HZ2 
			*  3HZ   HZ3 
			    C    C   
			    O    O   
			*/
			if (_charmmName == "HZ1") {
				return("1HZ");
				
			}
			if (_charmmName == "HZ2") {
				return("2HZ");
				
			}
			if (_charmmName == "HZ3") {
				return("3HZ");
				
			}
			
		}
		
		if (_resName == "ASN") {
			/*	         
			  NOTE: H numbering is inverted for terminal Hs

			   N     N   
			   H     H   
			   CA    CA  
			   CB    CB  
			   CG    CG  
			   OD1   OD1 
			   ND2   ND2 
			* 1HD2   HD22
			* 2HD2   HD21
			   C     C   
			   O     O   
			*/
			if (_charmmName == "HD22") {
				return("1HD2");
				
			}
			if (_charmmName == "HD21") {
				return("2HD2");
				
			}
			
		}
		
		if (_resName == "GLN") {
			/*
			  NOTE: H numbering is inverted for terminal Hs

			   N    N    
			   H    H   
			   CA   CA   
			   CB   CB   
			   CG   CG   
			   CD   CD   
			   OE1  OE1  
			   NE2  NE2  
			* 1HE2  HE22 
			* 2HE2  HE21 
			*/
			if (_charmmName == "HE21") {
				return("2HE2");
				
			}
			if (_charmmName == "HE22") {
				return("1HE2");
				
			}
			
		}
		if (_resName == "ARG") {
			/*
			  NOTE: H numbering is inverted for terminal Hs

			   N     N   
			   H     H   
			   CA    CA  
			   CB    CB  
			   CG    CG  
			   CD    CD  
			   NE    NE  
			   HE    HE  
			   CZ    CZ  
			   NH1   NH1 
			* 1HH1   HH12
			* 2HH1   HH11
			   NH2   NH2 
			* 1HH2   HH22
			* 2HH2   HH21
			   C     C   
			   O     O   
			*/  
			if (_charmmName == "HH11") {
				return("2HH1");
				
			}
			if (_charmmName == "HH12") {
				return("1HH1");
				
			}
			if (_charmmName == "HH21") {
				return("2HH2");
				
			}
			if (_charmmName == "HH22") {
				return("1HH2");
				
			}
			
		}
	
	}
	return _charmmName;
}


