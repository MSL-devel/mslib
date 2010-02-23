/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

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
