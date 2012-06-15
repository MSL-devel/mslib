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

//Forward Declarations
namespace MSL { 
class AtomPointerVector;


// Abstract base class
class KeyLookup {
	public:
		KeyLookup() {};
		virtual ~KeyLookup() {};
		virtual std::string getString(std::string _key)=0;
		virtual Real   getReal(std::string _key)  =0;
		virtual int    getInt(std::string _key)   =0;
		virtual bool   getBool(std::string _key)  =0;
		virtual bool   getQueryBool(std::string _key, std::string _arg) =0;

		virtual bool   getSelectionFlag(std::string _key) =0;
		virtual std::string isValidKeyword(std::string _key)=0;

};


template <class T> class Selectable : public KeyLookup {
	public:
                inline Selectable()  { ptTObj = NULL; }
		inline Selectable(T* _tObj) { ptTObj = _tObj;}
                virtual inline ~Selectable() {}

		inline void setObjectInstance(T *_tObj) { ptTObj = _tObj;}


		virtual inline void addSelectableFunctions() {};

		
		inline void addStringFunction(std::string _key, std::string(T::*_fpt)() const){
			keyValuePairStrings[_key] = _fpt;
			validKeywords[_key] = "string";	
			//std::cout << "Added: "<<_key<<" to std::strings"<<std::endl;
			
		}

		/*
		  // This gives compilation errors...
		std::string (T::*)() const getStringFunction(std::string _key){
			return keyValuePairStrings[_key];
		}
		*/
		inline std::string getString(std::string _key){
			return ((*ptTObj).*(keyValuePairStrings[_key]))();
		}


		inline void addRealFunction(std::string _key, Real(T::*_fpt)() const){
			keyValuePairReals[_key] = _fpt;
			validKeywords[_key] = "real";
		}

		inline Real getReal(std::string _key){
			return ((*ptTObj).*(keyValuePairReals[_key]))();
		}

		inline void addIntFunction(std::string _key, int(T::*_fpt)() const){
			keyValuePairInts[_key] = _fpt;
			validKeywords[_key] = "int";
		}

		inline int getInt(std::string _key){
			return ((*ptTObj).*(keyValuePairInts[_key]))();
		}


		inline void addBoolFunction(std::string _key, bool(T::*_fpt)() const){
			keyValuePairBools[_key] = _fpt;
			validKeywords[_key] = "bool";
		}

		inline bool getBool(std::string _key){
			return ((*ptTObj).*(keyValuePairBools[_key]))();
		}

		inline void addQueryBoolFunction(std::string _key, bool(T::*_fpt)(std::string _arg)){
			keyValuePairQueryBools[_key] = _fpt;
			validKeywords[_key] = "queryBool";
		}

		inline bool getQueryBool(std::string _key, std::string _arg){
		  return ((*ptTObj).*(keyValuePairQueryBools[_key]))(_arg);
		}


		inline std::string isValidKeyword(std::string _key){
			Hash<std::string,std::string>::Table::iterator it;
			it = validKeywords.find(_key);
			if (it == validKeywords.end()) {
				return "";
			} 
			return it->second;
		}


		inline void setSelectionFlag(std::string _key, bool _flag) { _key = MslTools::toUpper(_key); selectionFlags[_key] = _flag; }
		inline bool getSelectionFlag(std::string _key) { 
			_key = MslTools::toUpper(_key);
			Hash<std::string,bool>::Table::iterator it = selectionFlags.find(_key); 

			if (it != selectionFlags.end()){
				return selectionFlags[_key];  
			}


			return false;
		}


		inline void printAllFlags(){
			Hash<std::string,bool>::Table::iterator it;
			for (it = selectionFlags.begin(); it != selectionFlags.end(); it++){
				if (it->second){
					if (it != selectionFlags.begin()){
						std::cout << " , ";
					}
					std::cout <<it->first;
				}
			}
		}
	
		inline void clearFlag(std::string _key){
			_key = MslTools::toUpper(_key);
			Hash<std::string,bool>::Table::iterator it = selectionFlags.find(_key); 

			if (it != selectionFlags.end()){
				selectionFlags.erase(it);  
			}
		}

		inline void clearAllFlags(){
			selectionFlags.clear();
		}
		

	private:
		T *ptTObj;
		std::map<std::string,std::string (T::*)() const> keyValuePairStrings;
		std::map<std::string,Real (T::*)() const>   keyValuePairReals;
		std::map<std::string,int (T::*)() const>    keyValuePairInts;
		std::map<std::string,bool (T::*)() const>   keyValuePairBools;

		std::map<std::string,bool (T::*)(std::string)>   keyValuePairQueryBools;

		Hash<std::string,std::string>::Table validKeywords;
		Hash<std::string,bool>::Table selectionFlags;


#ifdef __BOOST__		
		friend class boost::serialization::access;		


		template<class Archive> inline void serialize(Archive & ar, const unsigned int version){
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

#ifndef __MACOS__
  class Atom;
  class Residue;
  class Truck;
  template class Selectable<Atom>;
  template class Selectable<Residue>;
  template class Selectable<Truck>;
#endif


}

#endif

