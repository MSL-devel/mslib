#include "FormatConverter.h"

FormatConverter::FormatConverter() {
}


FormatConverter::FormatConverter(System * _pSys) {
	pSys = _pSys;
}

FormatConverter::~FormatConverter() {
}


void FormatConverter::setCharmmFromPdb(string _charmmVersion) {
	if(pSys != NULL) {
		for(unsigned int i = 0; i < pSys->chainSize(); i++) {
			setCharmmFromPdb(pSys->getChain(i),_charmmVersion);
		}
	}
}

void FormatConverter::setCharmmFromPdb(Chain & _rChain, string _charmmVersion) {
	/*********************************
	 *  TRASNLATE THE RESIDUE NAMES
	 *********************************/
	bool nterm = false;
	bool cterm = false;
	for (int j=0; j<_rChain.positionSize(); j++) {
		Position& pos = _rChain.getPosition(j);
		if (j==0 && j==_rChain.positionSize() - 1) {
			// N-terminal & C-terminal
			nterm = cterm = true;
		} else if (j==0) {
			// N-terminal 
			nterm = true;
			cterm = false;
		} else if (j==_rChain.positionSize() - 1) {
			//C-terminal
			nterm = false;
			cterm = true;
		} else {
			nterm = cterm = false;
		}
		for (int i=0; i<pos.identitySize(); i++) {
			Residue & pRes = pos.getIdentity(i);
			setCharmmFromPdb(pRes, _charmmVersion, nterm, cterm);
		}
	}
}

void FormatConverter::setCharmmFromPdb(Residue & _rRes, string _charmmVersion, bool _Nterminal, bool _Cterminal) { 
	
	/*********************************
	 *  PDB RESNAME
	 *  
	 *  For the time being doing only the
	 *  protein residues is straightforward:
	 *  just change HIS to HIS/HSC/HSD/HSE/HSP 
	 *********************************/
	string pdbResName = _rRes.getResidueName();
	string resName = pdbResName;

	if (pdbResName == "HIS") {
		bool protonatedOnD = false;
		bool protonatedOnE = false;
		for (unsigned int i=0; i<_rRes.size(); i++) {
			if (_rRes.getAtom(i).getName() == "HD1") {
				protonatedOnD = true;
			} else if (_rRes.getAtom(i).getName() == "HE2") {
				protonatedOnE = true;
			}
		}
		if (protonatedOnD && protonatedOnE) {
			if (_charmmVersion == "19" || _charmmVersion == "20") {
				resName = "HSC";
			} else if (_charmmVersion == "22" || _charmmVersion == "27") {
				resName = "HSP";
			}
		} else if (protonatedOnD) {
			if (_charmmVersion == "19" || _charmmVersion == "20") {
				resName = "HIS";
			} else if (_charmmVersion == "22" || _charmmVersion == "27") {
				resName = "HSD";
			}
		} else {
			if (_charmmVersion == "19" || _charmmVersion == "20") {
				resName = "HSD";
			} else if (_charmmVersion == "22" || _charmmVersion == "27") {
				resName = "HSE";
			}
		}

	}

	if (pdbResName == "HOH") {
		resName = "TIP3";
	}

	// Here we set the CHARMM resName
	_rRes.setResidueName(resName);
	/*********************************
	 *  ATOM NAMES
	 *  
	 *  this ref has the XPLOR nomeclature,
	 *  http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl
	 *  need to check if the pro-S/R is consistent with
	 *  XPLOR or PDB
	 *
	 *  A bunch or hydrogen to change and
	 *  ILE CD to CD1
	 *********************************/
	for (int l=0; l<_rRes.size(); l++) {
		setCharmmFromPdb(_rRes.getAtom(l), _charmmVersion, _Nterminal, _Cterminal);
	}
}

/**************************************************************************************************************/

void FormatConverter::setCharmmFromPdb(Atom & _rAtom, string _charmmVersion, bool _Nterminal, bool _Cterminal) {
	string atomName = _rAtom.getName();

	string resName = _rAtom.getResidueName(); // This better be the CHARMM name
	
	if (_charmmVersion == "22" || _charmmVersion == "27") {

		// WATER
		//HETATM 5765  O   HOH B 443       6.728 -16.138 -64.807  1.00 36.50      1HTB6081
		
		if (resName == "HOH") {
			if (atomName == "O") {
				_rAtom.setName("OH2");
				return;
			}
			if (atomName == "1H") {
				_rAtom.setName("H1");
				return;
			}
			if (atomName == "2H") {
				_rAtom.setName("H2");
				return;
			}
		}

		
		// CONVERT TERMINAL PATCHES
		if (_Nterminal) {
			if (resName == "PRO") {
				/* STANDARD N-TERMINUS proline PROP */
				if (atomName == "1H") {
					_rAtom.setName("HN1");
					return;
				}
				if (atomName == "2H") {
					_rAtom.setName("HN2");
					return;
				}
			} else {
				/* STANDARD N-TERMINUS NTER and GLYP */
				if (atomName == "1H") {
					_rAtom.setName("HT1");
					return;
				}
				if (atomName == "2H") {
					_rAtom.setName("HT2");
					return;
				}
				if (atomName == "3H") {
					_rAtom.setName("HT3");
					return;
				}
			}
			/* ACETYL GROUP ACE (or ACP for PRO) N-terminal patch
			 *
			 *  TO DO THINGS PROPERGLY WE'D HAVE TO SET IT AS A
			 *  PATCH AND NOT A SEPARATE RESIDUE IN CHARMM 22    */
			
			if (resName == "ACE") {
				if (atomName == "CH3") {
					_rAtom.setName("CAY");
					return;
				}
				if (atomName == "1H") {
					_rAtom.setName("HY1");
					return;
				}
				if (atomName == "2H") {
					_rAtom.setName("HY2");
					return;
				}
				if (atomName == "3H") {
					_rAtom.setName("HY3");
					return;
				}
				if (atomName == "C") {
					_rAtom.setName("CY");
					return;
				}
				if (atomName == "O") {
					_rAtom.setName("OY");
					return;
				}
			}
		}
		if (_Cterminal) {
			/* Standard CTER */
			if (atomName == "O") {
				//_rAtom.setName("OT1");
				_rAtom.setName("OT2");
				return;
			}
			if (atomName == "OXT") {
				//_rAtom.setName("OT2");
				_rAtom.setName("OT1");
				return;
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
		if (resName == "ALA") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			if (atomName == "3HB") {
				_rAtom.setName("HB3");
				return;
			}
			if (atomName == "N" || atomName == "CA" || atomName == "HA" || atomName == "CB" || atomName == "C" || atomName == "O") {
				_rAtom.setName(atomName);
				return;
			}
			_rAtom.setName(atomName);
			return;
		}
		if (resName == "CYS") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			if (atomName == "HG") {
				_rAtom.setName("HG1");
				return;
			}
			return;
		}
		if (resName == "ASP") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			return;
		}
		if (resName == "GLU") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			if (atomName == "1HG") {
				_rAtom.setName("HG1");
				return;
			}
			if (atomName == "2HG") {
				_rAtom.setName("HG2");
				return;
			}
			return;
		}
		if (resName == "PHE") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			return;
		}
		if (resName == "GLY") {
			/*
			   N    N    
			*  H   HN   
			   CA   CA   
			* 1HA   HA1  
			* 2HA   HA2  
			   C    C    
			   O    O    
			*/       
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HA") {
				_rAtom.setName("HA1");
				return;
			}
			if (atomName == "2HA") {
				_rAtom.setName("HA2");
				return;
			}
			return;
		}
		if (resName == "HSD" || resName == "HSE" || resName == "HSP" ) {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			return;
		}
		if (resName == "ILE") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HG1") {
				_rAtom.setName("HG11");
				return;
			}
			if (atomName == "2HG1") {
				_rAtom.setName("HG12");
				return;
			}
			if (atomName == "1HG2") {
				_rAtom.setName("HG21");
				return;
			}
			if (atomName == "2HG2") {
				_rAtom.setName("HG22");
				return;
			}
			if (atomName == "3HG2") {
				_rAtom.setName("HG23");
				return;
			}
			if (atomName == "1HD1") {
				_rAtom.setName("HD1");
				return;
			}
			if (atomName == "2HD1") {
				_rAtom.setName("HD2");
				return;
			}
			if (atomName == "3HD1") {
				_rAtom.setName("HD3");
				return;
			}
			if (atomName == "CD1") {
				_rAtom.setName("CD");
				return;
			}
			return;
		}
		if (resName == "LYS") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			if (atomName == "1HG") {
				_rAtom.setName("HG1");
				return;
			}
			if (atomName == "2HG") {
				_rAtom.setName("HG2");
				return;
			}
			if (atomName == "1HD") {
				_rAtom.setName("HD1");
				return;
			}
			if (atomName == "2HD") {
				_rAtom.setName("HD2");
				return;
			}
			if (atomName == "1HE") {
				_rAtom.setName("HE1");
				return;
			}
			if (atomName == "2HE") {
				_rAtom.setName("HE2");
				return;
			}
			if (atomName == "1HZ") {
				_rAtom.setName("HZ1");
				return;
			}
			if (atomName == "2HZ") {
				_rAtom.setName("HZ2");
				return;
			}
			if (atomName == "3HZ") {
				_rAtom.setName("HZ3");
				return;
			}
			return;
		}
		if (resName == "LEU") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			if (atomName == "1HD1") {
				_rAtom.setName("HD11");
				return;
			}
			if (atomName == "2HD1") {
				_rAtom.setName("HD12");
				return;
			}
			if (atomName == "3HD1") {
				_rAtom.setName("HD13");
				return;
			}
			if (atomName == "1HD2") {
				_rAtom.setName("HD21");
				return;
			}
			if (atomName == "2HD2") {
				_rAtom.setName("HD22");
				return;
			}
			if (atomName == "3HD2") {
				_rAtom.setName("HD23");
				return;
			}
			return;
		}
		if (resName == "MET") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			if (atomName == "1HG") {
				_rAtom.setName("HG1");
				return;
			}
			if (atomName == "2HG") {
				_rAtom.setName("HG2");
				return;
			}
			if (atomName == "1HE") {
				_rAtom.setName("HE1");
				return;
			}
			if (atomName == "2HE") {
				_rAtom.setName("HE2");
				return;
			}
			if (atomName == "3HE") {
				_rAtom.setName("HE3");
				return;
			}
			return;
		}
		if (resName == "ASN") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			if (atomName == "2HD2") {
				_rAtom.setName("HD21");
				return;
			}
			if (atomName == "1HD2") {
				_rAtom.setName("HD22");
				return;
			}
			return;
		}
		if (resName == "PRO") {
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
			if (atomName == "2HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB2");
				return;
			}
			if (atomName == "2HG") {
				_rAtom.setName("HG1");
				return;
			}
			if (atomName == "1HG") {
				_rAtom.setName("HG2");
				return;
			}
			if (atomName == "2HD") {
				_rAtom.setName("HD1");
				return;
			}
			if (atomName == "1HD") {
				_rAtom.setName("HD2");
				return;
			}
			return;
		}
		if (resName == "GLN") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			if (atomName == "1HG") {
				_rAtom.setName("HG1");
				return;
			}
			if (atomName == "2HG") {
				_rAtom.setName("HG2");
				return;
			}
			if (atomName == "2HE2") {
				_rAtom.setName("HE21");
				return;
			}
			if (atomName == "1HE2") {
				_rAtom.setName("HE22");
				return;
			}
			return;
		}
		if (resName == "ARG") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			if (atomName == "1HG") {
				_rAtom.setName("HG1");
				return;
			}
			if (atomName == "2HG") {
				_rAtom.setName("HG2");
				return;
			}
			if (atomName == "1HD") {
				_rAtom.setName("HD1");
				return;
			}
			if (atomName == "2HD") {
				_rAtom.setName("HD2");
				return;
			}
			if (atomName == "2HH1") {
				_rAtom.setName("HH11");
				return;
			}
			if (atomName == "1HH1") {
				_rAtom.setName("HH12");
				return;
			}
			if (atomName == "2HH2") {
				_rAtom.setName("HH21");
				return;
			}
			if (atomName == "1HH2") {
				_rAtom.setName("HH22");
				return;
			}
			return;
		}
		if (resName == "SER") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			if (atomName == "HG") {
				_rAtom.setName("HG1");
				return;
			}
			return;
		}
		if (resName == "THR") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HG2") {
				_rAtom.setName("HG21");
				return;
			}
			if (atomName == "2HG2") {
				_rAtom.setName("HG22");
				return;
			}
			if (atomName == "3HG2") {
				_rAtom.setName("HG23");
				return;
			}
			return;
		}
		if (resName == "VAL") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HG1") {
				_rAtom.setName("HG11");
				return;
			}
			if (atomName == "2HG1") {
				_rAtom.setName("HG12");
				return;
			}
			if (atomName == "3HG1") {
				_rAtom.setName("HG13");
				return;
			}
			if (atomName == "1HG2") {
				_rAtom.setName("HG21");
				return;
			}
			if (atomName == "2HG2") {
				_rAtom.setName("HG22");
				return;
			}
			if (atomName == "3HG2") {
				_rAtom.setName("HG23");
				return;
			}
			return;
		}
		if (resName == "TRP") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			return;
		}
		if (resName == "TYR") {
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
			if (atomName == "H") {
				_rAtom.setName("HN");
				return;
			}
			if (atomName == "1HB") {
				_rAtom.setName("HB1");
				return;
			}
			if (atomName == "2HB") {
				_rAtom.setName("HB2");
				return;
			}
			return;
		}

		/************************************************
		 *  Done making atom name changes
		 ************************************************/
	}

	if (_charmmVersion == "19" || _charmmVersion == "20" ) {

		// WATER
		if (resName == "HOH") {
			if (atomName == "O") {
				_rAtom.setName("OH2");
				return;
			}
			if (atomName == "1H") {
				_rAtom.setName("H1");
				return;
			}
			if (atomName == "2H") {
				_rAtom.setName("H2");
				return;
			}
		}

		// CONVERT TERMINAL PATCHES
		if (_Nterminal) {
			/* STANDARD N-TERMINUS NTER, GLYP */
			if (resName == "PRO") {
				/* STANDARD N-TERMINUS proline PROP */
				if (atomName == "H1") {
					_rAtom.setName("HN1");
					return;
				}
				if (atomName == "H2") {
					_rAtom.setName("HN2");
					return;
				}
			} else {

				if (atomName == "H1") {
					_rAtom.setName("HT1");
					return;
				}
				if (atomName == "H2") {
					_rAtom.setName("HT2");
					return;
				}
				if (atomName == "H3") {
					_rAtom.setName("HT3");
					return;
				}
			}
			/* ACETYL GROUP ACE (or ACP for PRO) N-terminal patch*/
			/* ACE IS A RESIDUE IN CHARMM 19, NOT A PATCH AS IN CHARMM 22 */
		}
		if (_Cterminal) {
			/* Standard CTER */
			if (atomName == "O") {
				//_rAtom.setName("OT1");
				_rAtom.setName("OT2");
				return;
			}
			if (atomName == "OXT") {
				//_rAtom.setName("OT2");
				_rAtom.setName("OT1");
				return;
			}
		}


		/************************************************
		 *  Make the necessary atom name changes for each residue
		 *  
		 *  NEED TO ADD NUCLEIC ACIDS!
		 ************************************************/
		
		if (resName == "ILE") {
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
			if (atomName == "CD1") {
				_rAtom.setName("CD");
				return;
			}
			return;
		}
		if (resName == "LYS") {
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
			if (atomName == "1HZ") {
				_rAtom.setName("HZ1");
				return;
			}
			if (atomName == "2HZ") {
				_rAtom.setName("HZ2");
				return;
			}
			if (atomName == "3HZ") {
				_rAtom.setName("HZ3");
				return;
			}
			return;
		}
		if (resName == "ASN") {
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
			if (atomName == "1HD2") {
				_rAtom.setName("HD22");
				return;
			}
			if (atomName == "2HD2") {
				_rAtom.setName("HD21");
				return;
			}
			return;
		}
		
		if (resName == "GLN") {
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
			if (atomName == "2HE2") {
				_rAtom.setName("HE21");
				return;
			}
			if (atomName == "1HE2") {
				_rAtom.setName("HE22");
				return;
			}
			return;
		}
		if (resName == "ARG") {
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
			if (atomName == "2HH1") {
				_rAtom.setName("HH11");
				return;
			}
			if (atomName == "1HH1") {
				_rAtom.setName("HH12");
				return;
			}
			if (atomName == "2HH2") {
				_rAtom.setName("HH21");
				return;
			}
			if (atomName == "1HH2") {
				_rAtom.setName("HH22");
				return;
			}
			return;
		}
		if (resName == "SER") {
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
			if (atomName == "HG") {
				_rAtom.setName("HG1");
				return;
			}
			return;
		}
	}

}




/**************************************************************************************************************/
void FormatConverter::setPdbFromCharmm(string _charmmVersion) {
	if(pSys != NULL) {
		for (unsigned int i=0; i<pSys->chainSize(); i++) {
			setPdbFromCharmm(pSys->getChain(i), _charmmVersion);
		}
	}
}
	
void FormatConverter::setPdbFromCharmm(Chain & _rChain, string _charmmVersion) {
	/*********************************
	 *  TRASNLATE THE RESIDUE NUMBERS AND NAMES
	 *********************************/
	bool nterm = false;
	bool cterm = false;
	for (int j=0; j<_rChain.positionSize(); j++) {
		Position& pos = _rChain.getPosition(j);
		if (j==0 && j==_rChain.positionSize() - 1) {
			// N-terminal & C-terminal
			nterm = cterm = true;
		} else if (j==0) {
			// N-terminal 
			nterm = true;
			cterm = false;
		} else if (j==_rChain.positionSize() - 1) {
			//C-terminal
			nterm = false;
			cterm = true;
		} else {
			nterm = cterm = false;
		}
		for (int i=0; i<pos.identitySize(); i++) {
			Residue & pRes = pos.getIdentity(i);
			setPdbFromCharmm(pRes, _charmmVersion, nterm, cterm);
		}
	}

}

void FormatConverter::setPdbFromCharmm(Residue & _rRes, string _charmmVersion, bool _Nterminal, bool _Cterminal) { 
	/*********************************
	 *  PDB RESNAME
	 *  
	 *  For the time being doing only the
	 *  protein residues is straightforward:
	 *  just change HIS/HSC/HSD/HSE/HSP to HIS
	 *********************************/
	string charmmResName = _rRes.getResidueName();
	string resName = charmmResName;

	if (charmmResName == "HSD" || charmmResName == "HSE" || charmmResName == "HSP" || charmmResName == "HSC") {
		resName = "HIS";
	} else if (charmmResName == "TIP3") {
		resName = "HOH";
	}

	if (resName.size() > 3) {
		cerr << "WARNING 2912: charmm residue name " << resName << " too long for converting to pdb residue name: truncating it to the first 3 characters " << resName.substr(0, 3) << endl;
		resName = resName.substr(0,3);
	}

	_rRes.setResidueName(resName);
	/*********************************
	 *  ATOM NAMES
	 *  
	 *  this ref has the XPLOR nomeclature,
	 *  http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl
	 *  need to check if the pro-S/R is consistent with
	 *  XPLOR or PDB
	 *
	 *  A bunch or hydrogen to change and
	 *  ILE CD to CD1
	 *********************************/
	for (int l=0; l<_rRes.size(); l++) {
		setPdbFromCharmm(_rRes.getAtom(l), _charmmVersion, _Nterminal, _Cterminal);
	}
}

void FormatConverter::setPdbFromCharmm(Atom & _rAtom, string _charmmVersion, bool _Nterminal, bool _Cterminal) {
	string atomName = _rAtom.getName();

	// default: use the charmm atom name and set it
	
	string resName = _rAtom.getResidueName(); // This will be the PDB name
	
	if (_charmmVersion == "22" || _charmmVersion == "27") {

		if (resName == "TIP3") {
			if (atomName == "OH2") {
				_rAtom.setName("O");
				return;
			}
			if (atomName == "H1") {
				_rAtom.setName("1H");
				return;
			}
			if (atomName == "H2") {
				_rAtom.setName("2H");
				return;
			}
		}
		
		// CONVERT TERMINAL PATCHES
		if (_Nterminal) {
			/* STANDARD N-TERMINUS NTER and GLYP */
			if (atomName == "HT1") {
				_rAtom.setName("1H");
				return;
			}
			if (atomName == "HT2") {
				_rAtom.setName("2H");
				return;
			}
			if (atomName == "HT3") {
				_rAtom.setName("3H");
				return;
			}
			/* STANDARD N-TERMINUS proline PROP */
			if (atomName == "HN1") {
				_rAtom.setName("1H");
				return;
			}
			if (atomName == "HN2") {
				_rAtom.setName("2H");
				return;
			}
			/* ACETYL GROUP ACE (or ACP for PRO) N-terminal patch*/
			if (atomName == "CAY") {
				_rAtom.setName("CH3");
				/* What to do for the patch? 
				_rAtom.setAltPdbResName("ACE");
				_rAtom.setAltPdbResnum(_rAtom.getPdbResnum() - 1);
				_rAtom.setAltIcode("");
				_rAtom.setHetAtom(true);
				*/
				return;
			}
			if (atomName == "HY1") {
				_rAtom.setName("1H");
				return;
			}
			if (atomName == "HY2") {
				_rAtom.setName("2H");
				return;
			}
			if (atomName == "HY3") {
				_rAtom.setName("3H");
				return;
			}
			if (atomName == "CY") {
				_rAtom.setName("C");
				return;
			}
			if (atomName == "OY") {
				_rAtom.setName("O");
				return;
			}
		}
		if (_Cterminal) {
			/* Standard CTER */
			if (atomName == "OT1") {
				//_rAtom.setName("O");
				_rAtom.setName("OXT");
				return;
			}
			if (atomName == "OT2") {
				//_rAtom.setName("OXT");
				_rAtom.setName("O");
				return;
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
		if (resName == "ALA") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			if (atomName == "HB3") {
				_rAtom.setName("3HB");
				return;
			}
			return;
		}
		if (resName == "CYS") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			if (atomName == "HG1") {
				_rAtom.setName("HG");
				return;
			}
			return;
		}
		if (resName == "ASP") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			return;
		}
		if (resName == "GLU") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			if (atomName == "HG1") {
				_rAtom.setName("1HG");
				return;
			}
			if (atomName == "HG2") {
				_rAtom.setName("2HG");
				return;
			}
			return;
		}
		if (resName == "PHE") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			return;
		}
		if (resName == "GLY") {
			/*
			   N    N    
			*  H   HN   
			   CA   CA   
			* 1HA   HA1  
			* 2HA   HA2  
			   C    C    
			   O    O    
			*/       
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HA1") {
				_rAtom.setName("1HA");
				return;
			}
			if (atomName == "HA2") {
				_rAtom.setName("2HA");
				return;
			}
			return;
		}
		if (resName == "HIS") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			return;
		}
		if (resName == "ILE") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HG11") {
				_rAtom.setName("1HG1");
				return;
			}
			if (atomName == "HG12") {
				_rAtom.setName("2HG1");
				return;
			}
			if (atomName == "HG21") {
				_rAtom.setName("1HG2");
				return;
			}
			if (atomName == "HG22") {
				_rAtom.setName("2HG2");
				return;
			}
			if (atomName == "HG23") {
				_rAtom.setName("3HG2");
				return;
			}
			if (atomName == "HD1") {
				_rAtom.setName("1HD1");
				return;
			}
			if (atomName == "HD2") {
				_rAtom.setName("2HD1");
				return;
			}
			if (atomName == "HD3") {
				_rAtom.setName("3HD1");
				return;
			}
			if (atomName == "CD") {
				_rAtom.setName("CD1");
				return;
			}
			return;
		}
		if (resName == "LYS") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			if (atomName == "HG1") {
				_rAtom.setName("1HG");
				return;
			}
			if (atomName == "HG2") {
				_rAtom.setName("2HG");
				return;
			}
			if (atomName == "HD1") {
				_rAtom.setName("1HD");
				return;
			}
			if (atomName == "HD2") {
				_rAtom.setName("2HD");
				return;
			}
			if (atomName == "HE1") {
				_rAtom.setName("1HE");
				return;
			}
			if (atomName == "HE2") {
				_rAtom.setName("2HE");
				return;
			}
			if (atomName == "HZ1") {
				_rAtom.setName("1HZ");
				return;
			}
			if (atomName == "HZ2") {
				_rAtom.setName("2HZ");
				return;
			}
			if (atomName == "HZ3") {
				_rAtom.setName("3HZ");
				return;
			}
			return;
		}
		if (resName == "LEU") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			if (atomName == "HD11") {
				_rAtom.setName("1HD1");
				return;
			}
			if (atomName == "HD12") {
				_rAtom.setName("2HD1");
				return;
			}
			if (atomName == "HD13") {
				_rAtom.setName("3HD1");
				return;
			}
			if (atomName == "HD21") {
				_rAtom.setName("1HD2");
				return;
			}
			if (atomName == "HD22") {
				_rAtom.setName("2HD2");
				return;
			}
			if (atomName == "HD23") {
				_rAtom.setName("3HD2");
				return;
			}
			return;
		}
		if (resName == "MET") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			if (atomName == "HG1") {
				_rAtom.setName("1HG");
				return;
			}
			if (atomName == "HG2") {
				_rAtom.setName("2HG");
				return;
			}
			if (atomName == "HE1") {
				_rAtom.setName("1HE");
				return;
			}
			if (atomName == "HE2") {
				_rAtom.setName("2HE");
				return;
			}
			if (atomName == "HE3") {
				_rAtom.setName("3HE");
				return;
			}
			return;
		}
		if (resName == "ASN") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			if (atomName == "HD21") {
				_rAtom.setName("2HD2");
				return;
			}
			if (atomName == "HD22") {
				_rAtom.setName("1HD2");
				return;
			}
			return;
		}
		if (resName == "PRO") {
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
			if (atomName == "HB1") {
				_rAtom.setName("2HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HG1") {
				_rAtom.setName("2HG");
				return;
			}
			if (atomName == "HG2") {
				_rAtom.setName("1HG");
				return;
			}
			if (atomName == "HD1") {
				_rAtom.setName("2HD");
				return;
			}
			if (atomName == "HD2") {
				_rAtom.setName("1HD");
				return;
			}
			return;
		}
		if (resName == "GLN") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			if (atomName == "HG1") {
				_rAtom.setName("1HG");
				return;
			}
			if (atomName == "HG2") {
				_rAtom.setName("2HG");
				return;
			}
			if (atomName == "HE21") {
				_rAtom.setName("2HE2");
				return;
			}
			if (atomName == "HE22") {
				_rAtom.setName("1HE2");
				return;
			}
			return;
		}
		if (resName == "ARG") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			if (atomName == "HG1") {
				_rAtom.setName("1HG");
				return;
			}
			if (atomName == "HG2") {
				_rAtom.setName("2HG");
				return;
			}
			if (atomName == "HD1") {
				_rAtom.setName("1HD");
				return;
			}
			if (atomName == "HD2") {
				_rAtom.setName("2HD");
				return;
			}
			if (atomName == "HH11") {
				_rAtom.setName("2HH1");
				return;
			}
			if (atomName == "HH12") {
				_rAtom.setName("1HH1");
				return;
			}
			if (atomName == "HH21") {
				_rAtom.setName("2HH2");
				return;
			}
			if (atomName == "HH22") {
				_rAtom.setName("1HH2");
				return;
			}
			return;
		}
		if (resName == "SER") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			if (atomName == "HG1") {
				_rAtom.setName("HG");
				return;
			}
			return;
		}
		if (resName == "THR") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HG21") {
				_rAtom.setName("1HG2");
				return;
			}
			if (atomName == "HG22") {
				_rAtom.setName("2HG2");
				return;
			}
			if (atomName == "HG23") {
				_rAtom.setName("3HG2");
				return;
			}
			return;
		}
		if (resName == "VAL") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HG11") {
				_rAtom.setName("1HG1");
				return;
			}
			if (atomName == "HG12") {
				_rAtom.setName("2HG1");
				return;
			}
			if (atomName == "HG13") {
				_rAtom.setName("3HG1");
				return;
			}
			if (atomName == "HG21") {
				_rAtom.setName("1HG2");
				return;
			}
			if (atomName == "HG22") {
				_rAtom.setName("2HG2");
				return;
			}
			if (atomName == "HG23") {
				_rAtom.setName("3HG2");
				return;
			}
			return;
		}
		if (resName == "TRP") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			return;
		}
		if (resName == "TYR") {
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
			if (atomName == "HN") {
				_rAtom.setName("H");
				return;
			}
			if (atomName == "HB1") {
				_rAtom.setName("1HB");
				return;
			}
			if (atomName == "HB2") {
				_rAtom.setName("2HB");
				return;
			}
			return;
		}

		/************************************************
		 *  Done making atom name changes
		 ************************************************/
	}
	if (_charmmVersion == "19" || _charmmVersion == "20" ) {
		// CONVERT TERMINAL PATCHES
		if (_Nterminal) {
			/* STANDARD N-TERMINUS NTER, GLYP and PROP */
			if (atomName == "HT1") {
				_rAtom.setName("H1");
				return;
			}
			if (atomName == "HT2") {
				_rAtom.setName("H2");
				return;
			}
			if (atomName == "HT3") {
				_rAtom.setName("H3");
				return;
			}
			/* STANDARD N-TERMINUS proline PROP */
			if (atomName == "HN1") {
				_rAtom.setName("H1");
				return;
			}
			if (atomName == "HN2") {
				_rAtom.setName("H2");
				return;
			}
			/* ACETYL GROUP ACE (or ACP for PRO) N-terminal patch*/
			/* ACE IS A RESIDUE IN CHARMM 19, NOT A PATCH AS IN CHARMM 22 */
			/*
			if (resName == "ACE") {
				   CH3  CH3 
				   C    C   
				   O    O   
				if (atomName == "CH3" || atomName == "C" || atomName == "O") {
					_rAtom.setName(atomName);
					return;
				}
				_rAtom.setName(atomName);
				return;
			}
			*/  
		}
		if (_Cterminal) {
			/* Standard CTER */
			if (atomName == "OT1") {
				//_rAtom.setName("O");
				_rAtom.setName("OXT");
				return;
			}
			if (atomName == "OT2") {
				//_rAtom.setName("OXT");
				_rAtom.setName("O");
				return;
			}
		}


		/************************************************
		 *  Make the necessary atom name changes for each residue
		 *  
		 *  NEED TO ADD NUCLEIC ACIDS!
		 ************************************************/
		if (resName == "ILE") {
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
			if (atomName == "CD") {
				_rAtom.setName("CD1");
				return;
			}
			return;
		}
		if (resName == "LYS") {
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
			if (atomName == "HZ1") {
				_rAtom.setName("1HZ");
				return;
			}
			if (atomName == "HZ2") {
				_rAtom.setName("2HZ");
				return;
			}
			if (atomName == "HZ3") {
				_rAtom.setName("3HZ");
				return;
			}
			return;
		}
		
		if (resName == "ASN") {
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
			if (atomName == "HD22") {
				_rAtom.setName("1HD2");
				return;
			}
			if (atomName == "HD21") {
				_rAtom.setName("2HD2");
				return;
			}
			return;
		}
		
		if (resName == "GLN") {
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
			if (atomName == "HE21") {
				_rAtom.setName("2HE2");
				return;
			}
			if (atomName == "HE22") {
				_rAtom.setName("1HE2");
				return;
			}
			return;
		}
		if (resName == "ARG") {
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
			if (atomName == "HH11") {
				_rAtom.setName("2HH1");
				return;
			}
			if (atomName == "HH12") {
				_rAtom.setName("1HH1");
				return;
			}
			if (atomName == "HH21") {
				_rAtom.setName("2HH2");
				return;
			}
			if (atomName == "HH22") {
				_rAtom.setName("1HH2");
				return;
			}
			return;
		}
	
	}
}


