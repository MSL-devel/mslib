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

#include "MslOut.h"
#include "Transforms.h"
#include "AtomPointerVector.h"
#include "AtomContainer.h"
#include "OptionParser.h"
#include "Timer.h"
#include <stdio.h>
using namespace std;
using namespace MSL;

static MslOut MSLOUT("testMSLOut2");
int main(int argc, char *argv[]){

	OptionParser OP;
	OP.readArgv(argc, argv);

	Transforms t;
	AtomPointerVector av;
	MSLOUT.stream(MslOut::GENERAL) << "Its a general statement, this will always output."<<std::endl;
	MSLOUT.stream(MslOut::SPECIFIC) << "Its a specific statement, this will only output when this object output is on"<<std::endl;
	MSLOUT.stream(MslOut::WARNING) << "Its a WARNING that is specific to this object"<<std::endl;
	MSLOUT.stream(MslOut::ERROR) << "Its an ERROR that is specific to this object"<<std::endl;
	MSLOUT.debug() << "A debug statement like this is only seen when __MSL_MSLOUT_DEBUG_OFF__ is unset!"<<std::endl;

	// Performance test for MSLOUT;
	Timer tme;
	double start = tme.getWallTime();
	AtomContainer a;
	for (uint i = 0; i < 100000;i++){
		Atom *b = new Atom();
		a.addAtom(*b);
	}
	double end = tme.getWallTime();
	fprintf(stdout, "Time: %8.3f\n",(end-start));
	
}

