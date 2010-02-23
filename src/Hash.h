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
