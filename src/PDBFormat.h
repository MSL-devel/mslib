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

#ifndef PDBFORMAT_H
#define PDBFORMAT_H
/*
  PDB Formating information goes here.  As well as a simple structure for containing a single line
  of a PDB.

  For each line type define start, stop, length, justification of fixed-width fields.
  The format should be:
  S_NAME_OF_FIELD  = NUM,  for starting positions
  E_NAME_OF_FIELD  = NUM,  for ending positions
  L_NAME_OF_FIELD  = NUM,  for length of field
  J_NAME_OF_FIELD  = NUM,  for justification ; 0 = right, 1 = left jusftified
  
  E.G.

  enum NAME {
   S_ATOM_NAME = 12,  
   E_ATOM_NAME = 15,  
   L_ATOM_NAME = 4,  
   J_ATOM_NAME = 1,  
  };


Example PDB lines (Put standard + any strange ones here)

          1         2         3         4         5         6         7                     
01234567890123456789012345678901234567890123456789012345678901234567890123456789
ATOM     15  CD  LYS A   1      17.517  50.534  29.567  1.00 29.97           C
HETATM 1101  O   HOH   117      30.011  51.798   1.270  1.00 43.26           O


Example of Remark 290 Cyrstal Symmmetry + Remark 350 Bio Unit

          1         2         3         4         5         6         7                     
01234567890123456789012345678901234567890123456789012345678901234567890123456789
REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000            
REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000            
REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000            
REMARK 290   SMTRY1   2 -0.500000 -0.866020  0.000000        35.0000            
REMARK 290   SMTRY2   2  0.866030 -0.500000  0.000000        10.0000            
REMARK 290   SMTRY3   2  0.000000  0.000000  1.000000        2.00000            


Example of SCALE line
          1         2         3         4         5         6         7         8
012345678901234567890123456789012345678901234567890123456789012345678901234567890
SCALE1      0.019231  0.000000  0.000000        0.00000               
SCALE2      0.000000  0.017065  0.000000        0.00000               
SCALE3      0.000000  0.000000  0.016155        0.00000               


Example of CRYST1
          1         2         3         4         5         6         7         8
012345678901234567890123456789012345678901234567890123456789012345678901234567890
CRYST1   52.000   58.600   61.900  90.00  90.00  90.00 P 21 21 21    8 

Example of MODEL and ENDMDL
          1         2         3         4         5         6         7         8
012345678901234567890123456789012345678901234567890123456789012345678901234567890
MODEL        1                                                                   
ENDMDL                                                                           
 */

//MSL Includes
#include "MslTools.h"
#include "Atom.h"

//STL Includes
#include <string>
#include <iostream>
#include <cstring>


namespace MSL { 
class PDBFormat {
	public:

		enum Cryst {
			S_CRYSTRECORD     = 0,
			S_CRYSTINDEX      = 5,
			S_CRYSTA          = 6,
			S_CRYSTB          = 15,
			S_CRYSTC          = 25,
			S_CRYSTALPHA      = 33,
			S_CRYSTBETA       = 40,
			S_CRYSTGAMMA      = 47,
			S_CRYSTSPACEGROUP = 55,
			S_CRYSTZ          = 66,


			E_CRYSTRECORD     = 4,
			E_CRYSTINDEX      = 5,
			E_CRYSTA          = 14,
			E_CRYSTB          = 23,
			E_CRYSTC          = 32,
			E_CRYSTALPHA      = 39,
			E_CRYSTBETA       = 46,
			E_CRYSTGAMMA      = 53,
			E_CRYSTSPACEGROUP = 65,
			E_CRYSTZ          = 70,

			L_CRYSTRECORD     = 5,
			L_CRYSTINDEX      = 1,
			L_CRYSTA          = 7,
			L_CRYSTB          = 7,
			L_CRYSTC          = 7,
			L_CRYSTALPHA      = 5,
			L_CRYSTBETA       = 5,
			L_CRYSTGAMMA      = 5,
			L_CRYSTSPACEGROUP = 9,
			L_CRYSTZ          = 3,

		};
		struct CrystData {
			char D_CRYSTRECORD[sizeof(char)*(L_CRYSTRECORD+1)];
			char D_CRYSTSPACEGROUP[sizeof(char)*(L_CRYSTSPACEGROUP+1)];

			int       D_CRYSTINDEX;
			int       D_CRYSTZ;			

			double    D_CRYSTA;			
			double    D_CRYSTB;			
			double    D_CRYSTC;			
			double    D_CRYSTALPHA;			
			double    D_CRYSTBETA;			
			double    D_CRYSTGAMMA;			

			CrystData() {
				clear();
			};

			CrystData(const CrystData &_cryst){

				strcpy(D_CRYSTRECORD,_cryst.D_CRYSTRECORD);
				strcpy(D_CRYSTSPACEGROUP,_cryst.D_CRYSTSPACEGROUP);

				D_CRYSTINDEX  = _cryst.D_CRYSTINDEX;
				D_CRYSTZ      = _cryst.D_CRYSTZ;
				D_CRYSTA      = _cryst.D_CRYSTA;
				D_CRYSTB      = _cryst.D_CRYSTB;
				D_CRYSTC      = _cryst.D_CRYSTC;
				D_CRYSTALPHA  = _cryst.D_CRYSTALPHA;
				D_CRYSTBETA   = _cryst.D_CRYSTBETA;
				D_CRYSTGAMMA  = _cryst.D_CRYSTGAMMA;
				

			};

			void clear() {

				// Create a blank std::string of size L_FIELD_NAME
				strncpy(D_CRYSTRECORD   ,"                           ",L_CRYSTRECORD);     
				strncpy(D_CRYSTSPACEGROUP   ,"                           ",L_CRYSTSPACEGROUP);     


				// Terminate c-style std::strings properly..
				D_CRYSTRECORD[L_CRYSTRECORD]              = '\0';
				D_CRYSTSPACEGROUP[L_CRYSTSPACEGROUP]       = '\0';

				// Initialize numeric variables
				D_CRYSTINDEX = 0;
				D_CRYSTZ     = 0;
				D_CRYSTA     = 0.0;
				D_CRYSTB     = 0.0;
				D_CRYSTC     = 0.0;
				D_CRYSTALPHA = 0.0;
				D_CRYSTBETA  = 0.0;
				D_CRYSTGAMMA = 0.0;
			};
			
		};


		enum Scale {
			S_SCALERECORD   = 0,
			S_SCALELINE     = 5,
			S_SCALEX        = 10,
			S_SCALEY        = 20,
			S_SCALEZ        = 30,
			S_SCALETRANS    = 45,
			
			E_SCALERECORD   = 4,
			E_SCALELINE     = 5,
			E_SCALEX        = 19,
			E_SCALEY        = 29,
			E_SCALEZ        = 39,
			E_SCALETRANS    = 54,

			L_SCALERECORD   = 5,
			L_SCALELINE     = 1,
			L_SCALEX        = 10,
			L_SCALEY        = 10,
			L_SCALEZ        = 10,
			L_SCALETRANS    = 10,
		};

		struct ScaleData {
			char D_SCALERECORD[sizeof(char)*(L_SCALERECORD+1)];
			int       D_SCALELINE;
			double    D_SCALEX;			
			double    D_SCALEY;			
			double    D_SCALEZ;			
			double    D_SCALETRANS;			

			ScaleData() {
				clear();
			};

			ScaleData(const ScaleData &_scale){

				strcpy(D_SCALERECORD,_scale.D_SCALERECORD);

				D_SCALELINE   = _scale.D_SCALELINE;
				D_SCALEX      = _scale.D_SCALEX;
				D_SCALEY      = _scale.D_SCALEY;
				D_SCALEZ      = _scale.D_SCALEZ;
				D_SCALETRANS  = _scale.D_SCALETRANS;

			};

			void clear() {

				// Create a blank std::string of size L_FIELD_NAME
				strncpy(D_SCALERECORD   ,"                           ",L_SCALERECORD);     


				// Terminate c-style std::strings properly..
				D_SCALERECORD[L_SCALERECORD]       = '\0';

				// Initialize numeric variables
				D_SCALELINE  = 0;
				D_SCALEX     = 0.0;
				D_SCALEY     = 0.0;
				D_SCALEZ     = 0.0;
				D_SCALETRANS = 0.0;
			};
			
		};

		enum Remark290 {
			S_SYMMRECORD    =  13,
			S_SYMMLINE      =  18,
			S_SYMMINDEX     =  20,
			S_SYMMX         =  24,
			S_SYMMY         =  34,
			S_SYMMZ         =  44,
			S_SYMTRANS      =  55,

			E_SYMMRECORD    =  17,
			E_SYMMLINE      =  18,
			E_SYMMINDEX     =  22,
			E_SYMMX         =  32,
			E_SYMMY         =  42,
			E_SYMMZ         =  52,
			E_SYMTRANS      =  67,

			L_SYMMRECORD    =   5,
			L_SYMMLINE      =   1,
			L_SYMMINDEX     =   3,
			L_SYMMX         =   9,
			L_SYMMY         =   9,
			L_SYMMZ         =   9,
			L_SYMTRANS      =  13
		};

		struct SymData {
			char   D_SYMMRECORD[sizeof(char)*(L_SYMMRECORD+1)];
			int    D_SYMMLINE;
			int    D_SYMMINDEX;
			double D_SYMMX;
			double D_SYMMY;
			double D_SYMMZ;
			double D_SYMTRANS;

			SymData() {
				clear();
			};

			SymData(const SymData &_sym){

				strcpy(D_SYMMRECORD,_sym.D_SYMMRECORD);

				D_SYMMLINE  = _sym.D_SYMMLINE;
				D_SYMMINDEX = _sym.D_SYMMINDEX;
				D_SYMMX     = _sym.D_SYMMX;
				D_SYMMY     = _sym.D_SYMMY;
				D_SYMMZ     = _sym.D_SYMMZ;
				D_SYMTRANS  = _sym.D_SYMTRANS;

			}
			void clear() {

				// Create a blank std::string of size L_FIELD_NAME
				strncpy(D_SYMMRECORD   ,"                           ",L_SYMMRECORD);     


				// Terminate c-style std::strings properly..
				D_SYMMRECORD[L_SYMMRECORD]       = '\0';

				// Initialize numeric variables
				D_SYMMLINE     = 0;
				D_SYMMINDEX    = 0;
				D_SYMMX        = 0.0;
				D_SYMMY        = 0.0;
				D_SYMMZ        = 0.0;
				D_SYMTRANS     = 0.0;
			};
		};

		enum Remark350 {
			S_BIOURECORD    =  13,
			S_BIOULINE      =  18,
			S_BIOUINDEX     =  20,
			S_BIOUX         =  24,
			S_BIOUY         =  34,
			S_BIOUZ         =  44,
			S_BIOUTRANS     =  55,

			E_BIOURECORD    =  17,
			E_BIOULINE      =  18,
			E_BIOUINDEX     =  22,
			E_BIOUX         =  32,
			E_BIOUY         =  42,
			E_BIOUZ         =  52,
			E_BIOUTRANS     =  67,

			L_BIOURECORD    =   5,
			L_BIOULINE      =   1,
			L_BIOUINDEX     =   3,
			L_BIOUX         =   9,
			L_BIOUY         =   9,
			L_BIOUZ         =   9,
			L_BIOUTRANS     =  13
		};

		struct BioUData {
			char   D_BIOURECORD[sizeof(char)*(L_BIOURECORD+1)];
			int    D_BIOULINE;
			int    D_BIOUINDEX;
			double D_BIOUX;
			double D_BIOUY;
			double D_BIOUZ;
			double D_BIOUTRANS;

			BioUData() {
				clear();
			};

			BioUData(const BioUData &_bioU){

				strcpy(D_BIOURECORD,_bioU.D_BIOURECORD);

				D_BIOULINE  = _bioU.D_BIOULINE;
				D_BIOUINDEX = _bioU.D_BIOUINDEX;
				D_BIOUX     = _bioU.D_BIOUX;
				D_BIOUY     = _bioU.D_BIOUY;
				D_BIOUZ     = _bioU.D_BIOUZ;
				D_BIOUTRANS = _bioU.D_BIOUTRANS;

			}
			void clear() {

				// Create a blank std::string of size L_FIELD_NAME
				strncpy(D_BIOURECORD   ,"                           ",L_BIOURECORD);     


				// Terminate c-style std::strings properly..
				D_BIOURECORD[L_BIOURECORD]       = '\0';

				// Initialize numeric variables
				D_BIOULINE     = 0;
				D_BIOUINDEX    = 0;
				D_BIOUX        = 0.0;
				D_BIOUY        = 0.0;
				D_BIOUZ        = 0.0;
				D_BIOUTRANS    = 0.0;
			};
		};

		enum AtomField { 

			// Start Record Position
			S_RECORD_NAME    =  0,
			S_SERIAL         =  6,
			S_ATOM_NAME      = 12,
			S_ALT_LOC        = 16,
			S_RES_NAME       = 17,
			S_CHAIN_ID       = 21,
			S_RES_SEQ        = 22,
			S_I_CODE         = 26,
			S_X              = 30,
			S_Y              = 38,
			S_Z              = 46,
			S_OCCUP          = 54,
			S_TEMP_FACT      = 60,
			S_SEG_ID         = 72,
			S_ELEMENT_SYMBOL = 76,
			S_CHARGE         = 78,

			// End Record Position
			E_RECORD_NAME    =  5,
			E_SERIAL         = 10,
			E_ATOM_NAME      = 15,
			E_ALT_LOC        = 16,
			E_RES_NAME       = 19,
			E_CHAIN_ID       = 21,
			E_RES_SEQ        = 25,
			E_I_CODE         = 26,
			E_X              = 37,
			E_Y              = 45,
			E_Z              = 53,
			E_OCCUP          = 59,
			E_TEMP_FACT      = 65,
			E_SEG_ID         = 75,
			E_ELEMENT_SYMBOL = 77,
			E_CHARGE         = 79,


			// Length of Record
			L_RECORD_NAME    =  6,
			L_SERIAL         =  5,
			L_ATOM_NAME      =  4,
			L_ALT_LOC        =  1,
			L_RES_NAME       =  3,
			L_CHAIN_ID       =  1,
			L_RES_SEQ        =  4,
			L_I_CODE         =  1,
			L_X              =  8,
			L_Y              =  8,
			L_Z              =  8,
			L_OCCUP          =  6,
			L_TEMP_FACT      =  6,
			L_SEG_ID         =  4,
			L_ELEMENT_SYMBOL =  2,
			L_CHARGE         =  2,


			// Justification within record 0 = right-justified, 1 = left-justified.
			// TODO: Make this correct!
			J_RECORD_NAME    = 1,
			J_SERIAL         = 1,
			J_ATOM_NAME      = 1,
			J_ALT_LOC        = 1,
			J_RES_NAME       = 1,
			J_CHAIN_ID       = 1,
			J_RES_SEQ        = 1,
			J_I_CODE         = 1,
			J_X              = 1,
			J_Y              = 1,
			J_Z              = 1,
			J_OCCUP          = 1,
			J_TEMP_FACT      = 1,
			J_SEG_ID         = 1,
			J_ELEMENT_SYMBOL = 1,
			J_CHARGE         = 1

			
		};

		struct AtomData {
			char   D_RECORD_NAME[sizeof(char)*(L_RECORD_NAME+1)];
			int    D_SERIAL;
			char   D_ATOM_NAME[sizeof(char)*(L_ATOM_NAME+1)];
			char   D_ALT_LOC[sizeof(char)*(L_ALT_LOC+1)];
			char   D_RES_NAME[sizeof(char)*(L_RES_NAME+1)];
			char   D_CHAIN_ID[sizeof(char)*(L_CHAIN_ID+1)];
			int    D_RES_SEQ;
			char   D_I_CODE[sizeof(char)*(L_I_CODE+1)];
			double D_X;
			double D_Y;
			double D_Z;
			double D_OCCUP;
			double D_TEMP_FACT;
			char   D_SEG_ID[sizeof(char)*(L_SEG_ID+1)];
			char   D_ELEMENT_SYMBOL[sizeof(char)*(L_ELEMENT_SYMBOL+1)];
			double D_CHARGE;

			AtomData() {
				clear();
			};

			AtomData(const AtomData &_at){

				strcpy(D_RECORD_NAME,_at.D_RECORD_NAME);
				strcpy(D_ATOM_NAME,_at.D_ATOM_NAME);
				strcpy(D_ALT_LOC,_at.D_ALT_LOC);
				strcpy(D_RES_NAME,_at.D_RES_NAME);
				strcpy(D_CHAIN_ID,_at.D_CHAIN_ID);
				strcpy(D_I_CODE,_at.D_I_CODE);
				strcpy(D_SEG_ID,_at.D_SEG_ID);
				strcpy(D_ELEMENT_SYMBOL,_at.D_ELEMENT_SYMBOL);

				D_SERIAL = _at.D_SERIAL ;       
				D_RES_SEQ = _at.D_RES_SEQ ;      
				D_X = _at.D_X ;	       
				D_Y = _at.D_Y ;	       
				D_Z = _at.D_Z ;	       
				D_OCCUP = _at.D_OCCUP ;	       
				D_TEMP_FACT = _at.D_TEMP_FACT ;    


			}
			void clear() {

				// Create a blank std::string of size L_FIELD_NAME
				strncpy(D_RECORD_NAME   ,"                           ",L_RECORD_NAME);     
				strncpy(D_ATOM_NAME     ,"                           ",L_ATOM_NAME);      
				strncpy(D_ALT_LOC       ,"                           ",L_ALT_LOC);
				strncpy(D_RES_NAME      ,"                           ",L_RES_NAME);
				strncpy(D_CHAIN_ID      ,"                           ",L_CHAIN_ID);
				strncpy(D_I_CODE        ,"                           ",L_I_CODE);
				strncpy(D_SEG_ID        ,"                           ",L_SEG_ID);
				strncpy(D_ELEMENT_SYMBOL,"                           ",L_ELEMENT_SYMBOL); 	

				// Terminate c-style std::strings properly..
				D_RECORD_NAME[L_RECORD_NAME]       = '\0';
				D_ATOM_NAME[L_ATOM_NAME]           = '\0';
				D_ALT_LOC[L_ALT_LOC]               = '\0';       
				D_RES_NAME[L_RES_NAME]             = '\0';      
				D_CHAIN_ID[L_CHAIN_ID]             = '\0';      
				D_I_CODE[L_I_CODE]                 = '\0';        
				D_SEG_ID[L_SEG_ID]                 = '\0';        
				D_ELEMENT_SYMBOL[L_ELEMENT_SYMBOL] = '\0';

				// Initialize numeric variables
				D_SERIAL    = 0;
				D_RES_SEQ   = 0;        
				D_X         = 0.0;       
				D_Y         = 0.0;       
				D_Z         = 0.0;       
				D_OCCUP     = 0.0;       
				D_TEMP_FACT = 0.0;     
				D_CHARGE    = 0.0;     
			};
		};

		enum Model {
			S_MODEL_RECORD  =   0,
			S_MODEL_NUMBER  =  10,

			E_MODEL_RECORD  =   5,
			E_MODEL_NUMBER  =  13,

			L_MODEL_RECORD  =   6,
			L_MODEL_NUMBER  =   4
		};

		struct ModelData {
			unsigned int   D_MODEL_NUMBER;
			bool           D_ENDMODEL_FLAG;

			ModelData() {
				clear();
			};

			ModelData(const ModelData &_modelD){
				D_MODEL_NUMBER  = _modelD.D_MODEL_NUMBER ;
				D_ENDMODEL_FLAG = _modelD.D_ENDMODEL_FLAG;
			}
			void clear() {
				// Initialize numeric variables
				D_MODEL_NUMBER     = 0;
				D_ENDMODEL_FLAG    = false;
			};
		};

		static CrystData parseCrystLine(const std::string &_crystLine);
		static ScaleData parseScaleLine(const std::string &_scaleLine);
		static SymData  parseSymLine(const std::string &_symLine);
		static BioUData parseBioULine(const std::string &_bioULine);
		static AtomData parseAtomLine(const std::string &_pdbAtomLine);
		static ModelData parseModelLine(const std::string &_pdbModelLine);
		static AtomData createAtomData(const Atom &_at);
		static AtomData createAtomData(std::string _resName, Real &_x, Real &_y, Real &_z, std::string _element);
		static std::string createAtomLine(const AtomData &ad);
		static std::string createTerLine(const AtomData &ad);


	protected:
	private:


};
}

#endif
