#ifndef DISTANCEMATRIX_H
#define DISTANCEMATRIX_H

// STL Includes
#include <iostream>
#include <string>
#include <map>
#include <math.h>

// MSL Includes
#include "AtomPointerVector.h"
#include "Atom.h"
#include "MslTools.h"


// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#endif




using namespace std;


namespace MSL{

using namespace MslTools;
class DistanceMatrixResult;
class MatrixWindow;


class DistanceMatrix : public Matrix{

	public:
		DistanceMatrix();
		~DistanceMatrix();
		DistanceMatrix(const DistanceMatrix &_dm);

		// Inline accessor functions
		AtomPointerVector& getAtomVector() {return atomVec;}

		int getDimension() const{return dimension;}

		void setGeneralWinSize(int _winSize){generalWinSize = _winSize;}
		int getGeneralWinSize() const {return generalWinSize;}

		vector<MatrixWindow*> getMatrixWindows() const {return listMW;}  //should I be returning a pointer instead?? (would be pointer to a vec of pointers)

		void setIntraChain(bool _inChain){intraChain = _inChain;}
		bool getIntraChain(){return intraChain;}

		void setPDBid(string _id){PDBid = _id;}
		string getPDBid(){return PDBid;}

		void setWeight(vector<double> _weight){ weight = _weight;}
		vector<double> getWeight() { return weight; }
		

		void setDebug(bool _flag) { debug = _flag; }
		bool getDebug() { return debug; }

		void setDiagnolMatrixWindowsOnlyFlag(bool _flag) { diagnolMatrixWindowsOnly = _flag;}
		bool getDiagnolMatrixWindowsOnlyFlag() { return diagnolMatrixWindowsOnly;}

		void addAtom(const Atom &_atm);

		void createDistanceMatrix();


		MatrixWindow* operator()(string _key);
		vector<MatrixWindow*> & getDiagnolMWs();

		// Comparison Functions
		pair<vector<int>, double> compareAllWindows(DistanceMatrix &_distMat, int choice);	// Returns the coordinates of our best matching windows in the list of matrix windows (listMW)
		pair<int, double> compareAllWindows(vector<double> minDistVec1,vector<double> minDistVec2, vector<string> allowedAminoAcids1,vector<string> allowedAminoAcids2, int choice);	// Using weighted minimal distances 
		//  					    				       		 Returns the index of the best match
		//									  	   	 in listMW and the likeness double. 
		
		vector<DistanceMatrixResult> multiCompareAllWindows(DistanceMatrix &_dstMat, int choice, int _numCompare);
		
		//print the compare info
		void printCompareInfo(vector<double> _inputVec1,vector<double> _inputVec2, pair<int, double> _result, int choice);//Prints the result from compareAllMinRowWeighted
		void printCompareInfo(DistanceMatrix &_distMat, pair<vector<int>, double> _result, int choice);//Prints the result from compareAllWindows
		
		void  createMatrixWindows();

		enum CompareType {standard, diag, doubleDiag, minDist, minDistRow, minDistCol};
		enum CompareTypeVec {minRowWeighted, minRowAndColWeighted, correlation, correlationRowCol, correlationHeteroRowCol};


		// Boost related, is ok if no BOOST libraries are being used, just a string.
		void setArchiveType(string _type) { archiveType = _type; }
		string getArchiveType() { return archiveType; }

		
	protected:
		void copy(const DistanceMatrix &_dm);
		bool filterWindow(int _i, int _j, int _dim);
		
		AtomPointerVector atomVec;
		int dimension;
		int generalWinSize;
		bool intraChain;
		string archiveType;

		vector<MatrixWindow*> diagMW;		
		vector<MatrixWindow*> listMW;		
		map<string,MatrixWindow*> mapMW;

		string PDBid;

		vector<double> weight;
		bool debug;
		bool diagnolMatrixWindowsOnly;

		// BOOST-RELATED FUNCTIONS , keep them away from main class def.
#ifdef __BOOST__
	public:

		void save_checkpoint(string filename) const{

			if (archiveType == "binary"){
				std::ofstream fout(filename.c_str(),std::ios::binary);
				boost::archive::binary_oarchive oa(fout);
				oa << (*this);
			} else if (archiveType == "xml"){
				std::ofstream fout(filename.c_str());
				boost::archive::xml_oarchive oa(fout);
				oa << boost::serialization::make_nvp("DistanceMatrix",*this);
			} else {
				std::ofstream fout(filename.c_str());
				boost::archive::text_oarchive oa(fout);
				oa << (*this);
			}

		}

		void load_checkpoint(string filename){

			if (archiveType == "binary"){
				std::ifstream fin(filename.c_str(), std::ios::binary);
				boost::archive::binary_iarchive ia(fin);
				ia >> (*this);
			} else if (archiveType == "xml"){
				std::ifstream fin(filename.c_str());
				boost::archive::xml_iarchive ia(fin);
				ia >> boost::serialization::make_nvp("DistanceMatrix",*this);
			} else {
				std::ifstream fin(filename.c_str());
				boost::archive::text_iarchive ia(fin);
				ia >> (*this);
			}
		}
	private:

		friend class boost::serialization::access;		


		template<class Archive> void serialize(Archive & ar, const unsigned int version){

			using boost::serialization::make_nvp;
			//ar & boost::serialization::base_object<vector<Atom *> >(*this);
			//ar & geometricCenter;
			ar & make_nvp("matrix",boost::serialization::base_object<Matrix>(*this));
			ar & make_nvp("atomVec",atomVec);
			ar & make_nvp("dimension",dimension);
			ar & make_nvp("generalWinSize",generalWinSize);
			ar & make_nvp("intraChain",intraChain);
			ar & make_nvp("archiveType",archiveType);
			ar & make_nvp("listMW", listMW);
			ar & make_nvp("diagMW", diagMW);
			ar & make_nvp("mapMW", mapMW);
			ar & make_nvp("PDBid",PDBid);
			ar & make_nvp("weight",weight);
			ar & make_nvp("debug",debug);
		}
#else
	public:
		void save_checkpoint(string filename) const{
			cout << "NO IMPLEMENTATION OF SAVE_CHECKPOINT WITHOUT BOOST LIBRARIES INSTALLED.\n";
		}
		void load_checkpoint(string filename) const{
			cout << "NO IMPLEMENTATION OF LOAD_CHECKPOINT WITHOUT BOOST LIBRARIES INSTALLED.\n";
		}
#endif
		
};
}
#endif
