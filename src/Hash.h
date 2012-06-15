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

#ifndef HASH_H
#define HASH_H


/*
  Allow for multiple implementations hashing...
  
  Usage:

  Hash<std::string,std::string>::Table hashTable; // hashes std::strings to std::strings

  
*/

template <class Key, class Value> struct Hash {


#if defined(__GOOGLE__)
	#include <google/dense_hash_map> 
	typedef google::dense_hash_map<Key,Value> Table;
#else
// Include better std::maps (unordered_map) if GCC 4.1.1 or better
#if defined(__GNUC__) && (__GNUC__ == 4) && (__GNUC_MINOR__ == 1) && (__GNUC_PATCHLEVEL == 1)

	#include <tr1/unordered_map>
	/*
	  In order to be compatibily with google dense hash we need
	  to create an object that inherits unordered_map, then 
	  implmenents a "set_empty_key" function.
	 */

	typedef std::tr1::unordered_map<Key,Value> Table;
#else

	#include <map>

	/*
	  In order to be compatibily with google dense hash we need
	  to create an object that inherits std::map, then 
	  implmenents a "set_empty_key" function.
	 */

	typedef std::map<Key,Value> Table;

#endif
#endif


};

#endif
