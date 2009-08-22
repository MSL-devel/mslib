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

#ifndef SELECTABLE_H
#define SELECTABLE_H

// STL Includes
#include <string>
#include <map>

// MSL Includes
#include "Real.h"
#include "Hash.h"
#include "MslTools.h"

// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif

// Namespaces
using namespace std;

//Forward Declarations
class AtomVector;


// Abstract base class
class KeyLookup {
	public:
		KeyLookup() {};
		virtual ~KeyLookup() {};
		virtual string getString(string _key)=0;
		virtual Real   getReal(string _key)  =0;
		virtual int    getInt(string _key)   =0;
		virtual bool   getBool(string _key)  =0;

		virtual bool   getSelectionFlag(string _key) =0;
		virtual string isValidKeyword(string _key)=0;

};


template <class T> class Selectable : public KeyLookup {
	public:
		Selectable() { ptTObj = NULL; }
		Selectable(T* _tObj) { ptTObj = _tObj;}
		virtual ~Selectable() {}

		void setObjectInstance(T *_tObj) { ptTObj = _tObj;}


		virtual void addSelectableFunctions() {};

		
		void addStringFunction(string _key, string(T::*_fpt)() const){
			keyValuePairStrings[_key] = _fpt;
			validKeywords[_key] = "string";	
			//cout << "Added: "<<_key<<" to strings"<<endl;
			
		}

		/*
		  // This gives compilation errors...
		string (T::*)() const getStringFunction(string _key){
			return keyValuePairStrings[_key];
		}
		*/
		string getString(string _key){
			return ((*ptTObj).*(keyValuePairStrings[_key]))();
		}


		void addRealFunction(string _key, Real(T::*_fpt)() const){
			keyValuePairReals[_key] = _fpt;
			validKeywords[_key] = "real";
		}

		Real getReal(string _key){
			return ((*ptTObj).*(keyValuePairReals[_key]))();
		}

		void addIntFunction(string _key, int(T::*_fpt)() const){
			keyValuePairInts[_key] = _fpt;
			validKeywords[_key] = "int";
		}

		int getInt(string _key){
			return ((*ptTObj).*(keyValuePairInts[_key]))();
		}


		void addBoolFunction(string _key, bool(T::*_fpt)() const){
			keyValuePairBools[_key] = _fpt;
			validKeywords[_key] = "bool";
		}

		bool getBool(string _key){
			return ((*ptTObj).*(keyValuePairBools[_key]))();
		}


		string isValidKeyword(string _key){
			Hash<string,string>::Table::iterator it;
			it = validKeywords.find(_key);
			if (it == validKeywords.end()) {
				return "";
			} 
			return it->second;
		}


		void setSelectionFlag(string _key, bool _flag) { _key = MslTools::toUpper(_key); selectionFlags[_key] = _flag; }
		bool getSelectionFlag(string _key) { 
			_key = MslTools::toUpper(_key);
			Hash<string,bool>::Table::iterator it = selectionFlags.find(_key); 

			if (it != selectionFlags.end()){
				return selectionFlags[_key];  
			}


			return false;
		}


		void printAllFlags(){
			Hash<string,bool>::Table::iterator it;
			for (it = selectionFlags.begin(); it != selectionFlags.end(); it++){
				if (it->second){
					if (it != selectionFlags.begin()){
						cout << " , ";
					}
					cout <<it->first;
				}
			}
		}
	

	private:
		T *ptTObj;
		map<string,string (T::*)() const> keyValuePairStrings;
		map<string,Real (T::*)() const>   keyValuePairReals;
		map<string,int (T::*)() const>    keyValuePairInts;
		map<string,bool (T::*)() const>   keyValuePairBools;

		Hash<string,string>::Table validKeywords;
		Hash<string,bool>::Table selectionFlags;


#ifdef __BOOST__		
		friend class boost::serialization::access;		


		template<class Archive> void serialize(Archive & ar, const unsigned int version){
			ar & ptTObj;
			ar & keyValuePairStrings;
			ar & keyValuePairReals;
			ar & keyValuePairInts; 
			ar & keyValuePairBools;

			ar & validKeywords;
			ar & selectionFlags;
		}
#endif
		
};
class Atom;
class Residue;
class Truck;
template class Selectable<Atom>;
template class Selectable<Residue>;
template class Selectable<Truck>;
#endif

