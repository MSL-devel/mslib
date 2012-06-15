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

#ifndef HELIXFUSION_H
#define HELIXFUSION_H

#include "Chain.h"

namespace MSL { 
class HelixFusion {
	
	public:
		HelixFusion();
		~HelixFusion();

		// Code assumes nterm is fixed or reference
		void setChains(Chain &_nterm, Chain &_cterm);
		/*
		  fusionByAtomicAlignment
		     take the last 2 turns (7 residues)  from c-term and n-term chains
		     aligning sets of 4 of the C-alphas, minimizing RMSD (could be over all bb atoms?). Will return best fusion
                     under rmsdTolerance, or will return false. 
		 */
		bool fusionByAtomicAlignment(double _rmsdTolerance, std::string _newChainId="A");
		bool fusionByHelicalFrames();
		bool fusionByHydrogenBonding();

		int getNumberFusions();
		Chain& getFusedChain(int index); 

		std::vector<int> getCtermResiduesFromLastFusion();
		std::vector<int> getNtermResiduesFromLastFusion();

	private:
		// remove fused chains
		void deleteFusedChains();

		Chain *nterm; // Nterminal is defined as residue 0
		Chain *cterm; // Cterminal is defined as chain.size()-1

		std::vector<Chain *> fused;

		std::vector<int>    nIndex;
		std::vector<int>    cIndex;
};
inline HelixFusion::HelixFusion(){ nterm = NULL; cterm = NULL; }
inline HelixFusion::~HelixFusion(){ deleteFusedChains(); }
inline void HelixFusion::setChains(Chain &_nterm, Chain &_cterm) { nterm = &_nterm; cterm = &_cterm;}
inline Chain& HelixFusion::getFusedChain(int _index) { return *fused[_index];} 
inline int HelixFusion::getNumberFusions() { return fused.size(); }
inline void HelixFusion::deleteFusedChains() { for (uint f = 0; f < fused.size();f++) { delete(fused[f]); } fused.clear();}
inline std::vector<int> HelixFusion::getCtermResiduesFromLastFusion(){ return cIndex; }
inline std::vector<int> HelixFusion::getNtermResiduesFromLastFusion(){ return nIndex; }
}

#endif
