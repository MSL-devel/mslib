#ifndef MATRIXWINDOW_H
#define MATRIXWINDOW_H

#include <iostream>

#include "DistanceMatrix.h"
#include "MslTools.h"
//#include "Matrix.h"

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
namespace MSL {

class MatrixWindow{
	public:
		MatrixWindow();
		MatrixWindow(int initLeftR, int initLeftC, DistanceMatrix &m, int _winSize=10);
		MatrixWindow(const MatrixWindow &_mw);
		~MatrixWindow();

		

		//inline accessor functions
		void setLeftR (int leftRow){upLeftR=leftRow;}
		int getLeftR() const {return upLeftR;}
		void setLeftC(int leftCol){upLeftC = leftCol;}
		int getLeftC() const {return upLeftC;}

		void setWinSize(int windowSize) {winSize=windowSize;}
		int getWinSize() const {return winSize;}

		void setMatrix(DistanceMatrix & newMatrix) {bigMatrix = &newMatrix;}
		DistanceMatrix & getMatrix() const {return *bigMatrix;}

		AtomPointerVector getSmallAVec() const {return smallAtomVec;}

		vector<int> getMinDistanceRow() const {return minDistanceRow;}
		vector<int> getMinDistanceCol() const {return minDistanceCol;}
		vector<double> getMinRowValues() const {return minRowValues;}
		vector<double> getMinColValues() const {return minColValues;}

		//set the small atomVector (used by constructor and initialize)
		void setSmallAVec();

		//set the minDistance vector (used by constructor and initialize)
		void setMinDistanceRow();					//sets both minDistanceRow and the value vector minRowValues
		void setMinDistanceCol();					//sets minDistanceCol and minColValues

		//compare functions
		double compare(MatrixWindow & mw);
		double compareDiagonal(MatrixWindow & _mw);
		double compareDoubleDiagonal(MatrixWindow & _mw);
		double compareMinDist(MatrixWindow & _mw);		//compares minimum values along each row AND col
		double compareMinRow(MatrixWindow &_mw);		//compares minimum values along each row
		double compareMinCol(MatrixWindow &_mw);		//compares minimum values along each col
		double compareMinRowWeighted(vector<double> _minDistVec);//compares weighted minimum values along each row (takes a vector)
		double compareMinRCWeighted(vector<double> _minDistVec); //compares weighted min values along each row AND col (takes a vector)
		double compareCorrelation(vector<double> _minDistVec);	 //compares using correlation
		double compareCorrelationRowCol(vector<double> _minDistVec);// compares correlation using both row and column minDist vecs
		double compareCorrelationRowColAverage(vector<double> _minDistVec);// compares correlation using both row and column minDist vecs, by averaging
		double compareCorrelationHeteroRowCol(vector<double> _minDistVec1, vector<double> _minDistVec2);// compares correlation using both row and column minDist vecs

		//initialize
		void initialize(int _leftR, int _leftC, DistanceMatrix &m, int _winSize=10); 

		string toString();

		vector<int> getUpLeftResidueNumbers(); // Gives Row,Col residue numbers vector[0] = row, vector[1] = col;

		bool doesSequenceMatchRow(vector<string> _allowed);
		bool doesSequenceMatchCol(vector<string> _allowed);

		// Boost related, is ok if no BOOST libraries are being used, just a string.
		void setArchiveType(string _type) { archiveType = _type; }
		string getArchiveType() { return archiveType; }

	private:
		void copy(const MatrixWindow &_mw);

		int upLeftR;
		int upLeftC;
		int winSize;
		DistanceMatrix *bigMatrix;
		AtomPointerVector smallAtomVec;
		vector<int> minDistanceRow; //minDistanceRow[i] gives column coordinate of the minimum distance along the ith ROW
		vector<int> minDistanceCol; //minDistanceCol[i] gives row coordinates of the mimum distance along the ith COLUMN
		vector<double> minRowValues;
		vector<double> minColValues;
		string archiveType;

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
				oa << boost::serialization::make_nvp("MatrixWindow",*this);
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
				ia >> boost::serialization::make_nvp("MatrixWindow",*this);
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

			ar & make_nvp("upLeftR",upLeftR);
			ar & make_nvp("upLeftC",upLeftC);
			ar & make_nvp("winSize",winSize);
			ar & make_nvp("DistanceMatrix",bigMatrix);
			ar & make_nvp("smallAtomVec",smallAtomVec);
			ar & make_nvp("minDistanceRow",minDistanceRow);
			ar & make_nvp("minDistanceCol",minDistanceCol);
			ar & make_nvp("minRowValues",minRowValues);
			ar & make_nvp("minColValues",minColValues);
			ar & make_nvp("archiveType",archiveType);			

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
