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

#ifndef MSLOUT_H
#define MSLOUT_H

#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
using namespace std;

namespace MSL {

  class MslOut : public std::ostream {

	public:
	    MslOut();
	    MslOut(std::string _name);

	    enum MSG_TYPE { GENERAL=0,SPECIFIC=1,WARNING=2,ERROR=3,DEBUG=4};

	    void turnAllOff();
	    void turnAllOn();
	    void turnOff(std::string _name);
	    void turnOn(std::string _name);
	
	    void fprintf(FILE *_stream, const char * _format, ...);
	    
	    void printOutputOnFlags();
	    void printAllFlags();
	
	    std::ostream& stream(MSG_TYPE _type=SPECIFIC);
	    std::ostream& debug();

	    static std::map<std::string,double >& getStaticLookup();

	private:
	    // Helper functions to make the vector like a map.  (OSX doesn't like static maps, but likes static vectors?)
	    bool nameInVector(std::vector<std::string> &_vec,std::string _key);	     
	    void uniqueAdd(std::vector<std::string> &_vec,std::string _key);	     
	    void uniqueRemove(std::vector<std::string> &_vec,std::string _key);	     

	    //static std::map<std::string, bool> outputFlag; // problems with OSX and STL maps
	    static std::vector<std::string>& getAllFlags();
	    static std::vector<std::string>& getOutputOnFlags();


	    std::vector<std::string>::iterator it;
	    std::string name;
	    
};



inline bool MslOut::nameInVector(std::vector<string> &_vec,std::string _key){
	it = find(_vec.begin(),_vec.end(),_key);
	return (it != _vec.end());
}

inline void MslOut::uniqueAdd(std::vector<std::string> &_vec, string _key){

        if (!nameInVector(_vec,_key)){
	  _vec.push_back(_key);
	}	
}

inline void MslOut::uniqueRemove(std::vector<std::string> &_vec, string _key){

        if (nameInVector(_vec,_key)){
	  _vec.erase(it);
	}	
}

inline void MslOut::turnOff(std::string _name){
	static std::vector<std::string> &outputOnFlags = MslOut::getOutputOnFlags();

	uniqueRemove(outputOnFlags,_name);
}
inline void MslOut::turnOn(std::string _name){
	static std::vector<std::string> &allFlags = MslOut::getAllFlags();
	static std::vector<std::string> &outputOnFlags = MslOut::getOutputOnFlags();

	uniqueAdd(allFlags,_name);
        uniqueAdd(outputOnFlags,_name);
}


inline void MslOut::turnAllOff(){
	static std::vector<std::string> &outputOnFlags = MslOut::getOutputOnFlags();
	outputOnFlags.clear();
}

inline void MslOut::turnAllOn(){
	static std::vector<std::string> &allFlags = MslOut::getAllFlags();
	static std::vector<std::string> &outputOnFlags = MslOut::getOutputOnFlags();

	outputOnFlags = allFlags;

}

inline void MslOut::fprintf(FILE *_stream, const char * _format, ...){

	static std::vector<std::string> &outputOnFlags = MslOut::getOutputOnFlags();

	if (nameInVector(outputOnFlags,name)){
	    std::stringstream ss;
	    ss << name << ":\t"<< _format;

	    va_list args;
	    va_start (args, _format);
	    vfprintf (_stream, ss.str().c_str(), args);
	    va_end (args);
	}
}


inline std::ostream& MslOut::debug(){

#ifndef __MSL_MSLOUT_DEBUG_OFF__
  static std::vector<std::string> &outputOnFlags = MslOut::getOutputOnFlags();
  if (nameInVector(outputOnFlags,name)){
	  std::cout << name << " DEBUG:\t";
	  return std::cout;
  }
  return *this;
#else
  return *this;
#endif
}

inline std::ostream& MslOut::stream(MSG_TYPE _type){	
  static std::vector<std::string> &outputOnFlags = MslOut::getOutputOnFlags();

  // When MSG_TYPE is GENERAL, print regardless if outputflag is on
  if (_type == GENERAL){
    std::cout << "MSLGENERAL:\t";
    return std::cout;
  }

  // When MSG_TYPE is SPECIFIC/DEBUG, print only if flag is set.
  if (nameInVector(outputOnFlags,name)){
	  if (_type == SPECIFIC){
		  std::cout << name << ":\t";
	  }
	  if (_type == WARNING){
		  std::cout << name <<  " WARNING:\t";
	  }
	  if (_type == ERROR){
		  std::cout << name <<  " ERROR:\t";
	  }

	  return std::cout;
  }

  return *this;

}



inline void MslOut::printOutputOnFlags(){

  static std::vector<std::string> &outputOnFlags = MslOut::getOutputOnFlags();

  for (it = outputOnFlags.begin(); it != outputOnFlags.end();it++){
	    cout << "outputOnFlags -> "<<*it<<endl;
  }

}
inline void MslOut::printAllFlags(){
  static std::vector<std::string> &allFlags = MslOut::getAllFlags();

  for (it = allFlags.begin(); it != allFlags.end();it++){
	    cout << "allFlags -> "<<*it<<endl;
  }
}


}


#endif
