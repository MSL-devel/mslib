#include "MslOut.h"
#include "Transforms.h"
#include "AtomPointerVector.h"
#include <stdio.h>
using namespace std;
using namespace MSL;

static MslOut MSLOUT("testMSLOut");
int main(){

  //MSLOUT.printOutputOnFlags();   
  //MSLOUT.printAllFlags();   
  //MSLOUT.turnAllOff();
  //MSLOUT.turnOff("AtomPointerVector");
  //MSLOUT.turnOn("AtomPointerVector");

  MSLOUT.turnOff("testMSLOut");
  MSLOUT.fprintf(stdout,"Hi there\n");

  Transforms t;
  AtomPointerVector av;


  MSLOUT.stream(MslOut::GENERAL) << "Hey general!"<<std::endl;
  MSLOUT.stream(MslOut::SPECIFIC) << "Hey specific!"<<std::endl;
}
