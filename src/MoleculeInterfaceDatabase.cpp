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

#include "MoleculeInterfaceDatabase.h"

using namespace MSL;


/**
 * The deconstructor of our class.  The deconstructor
 * will iterate through all entries in our database
 * and delete each one.  Thus, any InterfaceResidueDescriptor pointer
 * that was obtained from our database will be unusable
 * after the database has been destroyed.
 */
MoleculeInterfaceDatabase::~MoleculeInterfaceDatabase(){
	for (uint i = 0; i < residueDescriptors.size();i++){
		delete(residueDescriptors[i]);
	}

	residueDescriptors.clear();
	databaseFile = "";
}

/**
 * This method will clone a second MoleculeInterfaceDatabase.
 *
 * @param _mid      The database to be cloned.
 */
void MoleculeInterfaceDatabase::copy(MoleculeInterfaceDatabase &_mid){

	residueDescriptors.clear();

	
	for (uint i = 0; i < _mid.size();i++){
		InterfaceResidueDescriptor *ird = _mid[i];
		residueDescriptors.push_back(new InterfaceResidueDescriptor(*ird));
	}

	databaseFile = _mid.getDatabaseFile();

	archiveType = _mid.getArchiveType();
}

/**
 * This method allows a user to add a new entry
 * to our database.
 *
 * @param _ird        The new InterfaceResidueDescriptor to add to our database.
 * @see InterfaceResidueDescriptor
 */
void MoleculeInterfaceDatabase::addInterafaceResidueDescriptor(InterfaceResidueDescriptor *_ird){
	InterfaceResidueDescriptor &ird = *_ird;
	residueDescriptors.push_back(new InterfaceResidueDescriptor(ird));
}
