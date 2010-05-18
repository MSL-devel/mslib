#include "MslOut.h"
using namespace MSL;


// STL maps do not seem to work as static variables on OSX..
//std::map<std::string, bool> MslOut::outputFlag = std::map<std::string,bool>();

//std::vector<std::string> MslOut::outputFlagOn = std::vector<std::string>();
//std::vector<std::string> MslOut::allFlags  = std::vector<std::string>();

//std::vector<std::string> MslOut::outputFlagOn;
//std::vector<std::string> MslOut::allFlags;



std::vector<std::string>& MslOut::getAllFlags(){
  static std::vector<std::string>* allFlagsVar = new std::vector<std::string>();
  return *allFlagsVar;
}

std::vector<std::string>& MslOut::getOutputOnFlags(){
  static std::vector<std::string>* outputOnFlagsVar = new std::vector<std::string>();
  return *outputOnFlagsVar;
}
