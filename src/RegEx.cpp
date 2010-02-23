#include "RegEx.h"
#include "PolymerSequence.h"

using namespace MSL;
using namespace std;


RegEx::RegEx(){
}

RegEx::RegEx(const RegEx &_regex){
	copy(_regex);
}

RegEx::~RegEx(){
}

void RegEx::operator=(const RegEx & _regex) {
    copy(_regex);
}

void RegEx::copy(const RegEx &_regex){
}


vector<pair<int,int> > RegEx::getResidueRanges(Chain &_ch, string _regex){

	vector<pair<int,int> > results;

	// Chain to 1-AA string...
	string aas = PolymerSequence::toOneLetterCode(_ch.getAtoms());
	
	cout << "CHAIN: "<<aas<<endl;
	// Iterative search storing indices.
	boost::regex expression(_regex);

	//boost::sregex_token_iterator r1(aas.begin(),aas.end(),expression,-1);
	boost::sregex_iterator r1(aas.begin(),aas.end(),expression);
	boost::sregex_iterator r2;
	/*
	cout << "r1: "<<(*r1).size()<<endl;
	for (uint i =0 ;i < (*r1).size();i++){
		cout << "M: "<<(*r1).position()<<endl;
		*r1++;
	}
	*/
	while (r1 != r2){
		cout << "MATCH: "<<*r1<<" "<<(*r1).position()<<" "<<(*r1).length()<<endl;

		pair<int,int> a;
		a.first = (*r1).position();
		a.second = (*r1).position()+(*r1).length()-1;
		results.push_back(a);

		*r1++;

		
	}

	return results;

}
