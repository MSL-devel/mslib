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

#ifndef FORMATCONVERTER_H
#define FORMATCONVERTER_H

#include "AtomPointerVector.h"
#include "Residue.h"
#include "PDBFormat.h"


/*! \brief  This Object converts charmm variables from and to pdb variables
 *  (atom and residue names )
 * supports PDB3, PDB2 <==> CHARMM19, CHARMM20, CHARMM22, CHARMM27
 *   
 */

namespace MSL {
	class FormatConverter {
		public:

			FormatConverter();
			FormatConverter(std::string _orig, std::string _tgt);

			~FormatConverter();

			bool setNamespaces (std::string _orig,std::string _tgt);
			std::string getResidueName(std::string _resName, bool _protonatedOnD = false, bool _protonatedOnE = false);
			// _resName for getAtomName should be obtained using getResidueName above
			std::string getAtomName(std::string _atomName, std::string _resName, bool _NTerminal=false,bool _CTerminal=false);

			void convert(Atom& _atom,bool _NTerminal=false,bool _CTerminal=false, bool _protonatedOnD = false, bool _protonatedOnE = false);
			void convert(AtomPointerVector& _apV);
			void convert(PDBFormat::AtomData& _atomLine);
			bool isConversionSupported(std::string _orig, std::string _tgt);

		private:

			std::string getCharmmResName(std::string _pdbResName,std::string _charmmVersion, bool _protonatedOnD = false, bool _protonatedOnE = false);
			std::string getCharmmAtomName(std::string _pdbName, std::string _resName, std::string _charmmVersion, bool _Nterminal, bool _Cterminal, bool _PDB2_flag=false);

			std::string getPdbResName(std::string _charmmResName);
			std::string getPdbAtomName(std::string _charmmName, std::string _resName, std::string _charmmVersion, bool _PDB2_flag=false);
			
			std::string orig;
			std::string tgt;

	};
}
#endif

