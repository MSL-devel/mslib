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

#ifndef INTERFACERESIDUEDESCRIPTOR_H
#define INTERFACERESIDUEDESCRIPTOR_H

// MSL Includes
#include "MslTools.h"

// STL Includes
#include <map>
#include <string>
#include <fstream>

// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#endif



/// This enum will define the various types of chains that can be found in complexes.
typedef enum{ Other = 0,    ///< A catch-all default value.
              Protein = 1,  ///< This residue is part of a protein.
              Peptide = 2,  ///< This residue is part of a peptide.
}ChainType_t; 

/**
 * This class is essentially a container object which holds
 * information about individual residues. This information includes
 * solvent accesibility, secondary structure, number of backbone
 * to backbone contacts etc.  These objects are used to hold
 * the entries in a MoleculeInterfaceDatabase.
 *
 * @see MoleculeInterfaceDatabase
 */
namespace MSL { 
class InterfaceResidueDescriptor {

    public:



        InterfaceResidueDescriptor();
        InterfaceResidueDescriptor(InterfaceResidueDescriptor & _ird);
        ~InterfaceResidueDescriptor();


	void operator=(InterfaceResidueDescriptor & _ird); // assignment

        void setPdbId(std::string _pdbid);
        std::string getPdbId();


        void setChainId(std::string _chainid);
        std::string getChainId();

        void setResidueNumber(int _residueNumber);
        int getResidueNumber();


	void setResidueName(std::string _residueName);
	std::string getResidueName();

	void setSecondaryStructure(std::string _secondaryStructure);
	std::string getSecondaryStructure();
		
	void setSingleChainDeltaSolventAccessibility(double _deltaSolventAcc);
	double getSingleChainDeltaSolventAccessibility();
        void setOtherChainDeltaSolventAccessibility(std::string _peptideChain, double _deltaSolventAcc);
	double getOtherChainDeltaSolventAccessibility(std::string _peptideChain);

        void setAAProb(std::string _aaName, double _prob);
        double getAAProb(std::string _aaName);

        void setConsScore(double _score);
        double  getConsScore();
        void setConsScoreConfInterval(double _low, double _high);
        std::pair<double, double> getConsScoreConfInterval();

        void setResolution(std::string _resolution);
        void setResolution(double _resolution);
        double getResolution();
        bool isResolutionValid();


	void setResidueIcode(std::string _icode);
	std::string getResidueIcode();

        void setChainType(ChainType_t _type);
        ChainType_t getChainType();


        void setNumberBackboneContacts(int _backboneContacts);
        int getNumberBackboneContacts();


        void setNumberSidechainContacts(int _sidechainContacts);
        int getNumberSidechainContacts();


        void setNumberMixedContacts(int _mixedContacts);
        int getNumberMixedContacts();

        void addToDescriptionTable(std::string & _key, std::string & _value);
        std::string getFromDescriptionTable(std::string &_key);
        void removefromDescriptionTable(std::string &_key);
        void clearDescriptionTable();


        bool isStringInResidueDescription(std::string & _query);


        std::map<std::string,std::string> & getDescriptionTable();

	std::string toString();

	// MOVED HERE FROM BOOST PRIVATE TO ALLOW COMPILATION WITHOUT BOOST
	std::string getArchiveType() { return archiveType; }
	void setArchiveType(std::string _type) { archiveType = _type; }

    private:



	void copy(InterfaceResidueDescriptor &_ird);


        std::string pdbId;
        std::string chainId;
        std::string residueIcode;

        double resolution;
        bool resolutionIsValid;
        std::string residueName;

	int residueNumber;
	std::string secondaryStructure;
	double singleChainDeltaSolventAccessibility;
        std::map<std::string, double> otherChainDeltaSolventAccessibility;
        std::map<std::string, double> aaProbs;
        double consScore;
        std::pair<double, double> consScoreConfInterval;

	// Chain Type (Protein,Peptide,etc..)
	ChainType_t chainType;
		
	int numberBackboneContacts;  // BB-to-BB
	int numberSidechainContacts; // SC-to-SC
	int numberMixedContacts;     // SC-to-BB or BB-to-SC
		
	// Could we use this for mol1desc, mol2desc, function, etc?
	std::map<std::string, std::string> descriptionTable;

	// MOVED HERE FROM BOOST PRIVATE TO ALLOW COMPILATION WITHOUT BOOST
	std::string archiveType;


        // BOOST-RELATED FUNCTIONS , keep them away from main class def.

#ifdef __BOOST__

    public:


	void save_checkpoint(std::string filename) const{

		if (archiveType == "binary"){
			std::ofstream fout(filename.c_str(),std::ios::binary);
			boost::archive::binary_oarchive oa(fout);
			oa << (*this);
		} else if (archiveType == "xml"){
			std::ofstream fout(filename.c_str());
			boost::archive::xml_oarchive oa(fout);
			oa << boost::serialization::make_nvp("InterfaceResidueDescriptor",*this);

		} else {
			std::ofstream fout(filename.c_str());
			boost::archive::text_oarchive oa(fout);
			oa << (*this);
		}

	}

	void load_checkpoint(std::string filename){

		if (archiveType == "binary"){
			std::ifstream fin(filename.c_str(), std::ios::binary);
			boost::archive::binary_iarchive ia(fin);
			ia >> (*this);
		} else if (archiveType == "xml"){
			std::ifstream fin(filename.c_str());
			boost::archive::xml_iarchive ia(fin);
			ia >> boost::serialization::make_nvp("InterfaceResidueDescriptor",*this);
			//ia >> (*this);
		} else {
			std::ifstream fin(filename.c_str());
			boost::archive::text_iarchive ia(fin);
			ia >> (*this);
		}


	}





    private:
        friend class boost::serialization::access;        


        template<class Archive> void serialize(Archive & ar, const unsigned int version){

		//		if (archiveType == "xml"){
			using boost::serialization::make_nvp;
			ar & make_nvp("pdbId",pdbId);
			ar & make_nvp("chainId",chainId);
			ar & make_nvp("resolution",resolution);
			ar & make_nvp("resolutionIsValid",resolutionIsValid);
			ar & make_nvp("residueName",residueName);
			ar & make_nvp("residueNumber",residueNumber);
			ar & make_nvp("residueIcode",residueIcode);
			ar & make_nvp("secondaryStructure",secondaryStructure);
			ar & make_nvp("singleChainDeltaSolventAccessibility",singleChainDeltaSolventAccessibility);
                        ar & make_nvp("otherChainDeltaSolventAccessibility",otherChainDeltaSolventAccessibility);
                        ar & make_nvp("aaProbs",aaProbs);
                        ar & make_nvp("consScore", consScore);
                        ar & make_nvp("consScoreConfInterval", consScoreConfInterval);
			ar & make_nvp("chainType",chainType);
			ar & make_nvp("numberBackboneContacts",numberBackboneContacts);
			ar & make_nvp("numberSidechainContacts",numberSidechainContacts);
			ar & make_nvp("numberMixedContacts",numberMixedContacts);
			ar & make_nvp("descriptionTable",descriptionTable);
/* 		} else { */
/* 			ar & pdbId; */
/* 			ar & chainId; */
/* 			ar & resolution; */
/* 			ar & resolutionIsValid; */
/* 			ar & residueName; */
/* 			ar & residueNumber; */
/* 			ar & residueIcode; */
/* 			ar & secondaryStructure; */
/* 			ar & deltaSolventAccessibility; */
/* 			ar & chainType; */
/* 			ar & numberBackboneContacts; */
/* 			ar & numberSidechainContacts; */
/* 			ar & numberMixedContacts; */
/* 			ar & descriptionTable; */
//		}
        }
#endif        
        
};

//Inlines go HERE
/**
 * The basic constructor.
 */
inline InterfaceResidueDescriptor::InterfaceResidueDescriptor() {
    resolutionIsValid = false;
    singleChainDeltaSolventAccessibility =0.0;
    consScore = 0.0;
    consScoreConfInterval.first = 0.0f;
    consScoreConfInterval.second = 0.0f;
    archiveType="text";
}
/**
 * The basic deconstructor.
 */
inline InterfaceResidueDescriptor::~InterfaceResidueDescriptor() {}
/**
 * A copy constructor which will clone a given InterfaceResidueDescriptor.
 *
 * @param _ird  An InterfaceResidueDescriptor to clone.
 */
inline InterfaceResidueDescriptor::InterfaceResidueDescriptor(InterfaceResidueDescriptor & _ird) { copy(_ird); }

/**
 * Overloading the = operator so that we can easily clone
 * a given InterfaceResidueDescriptor.
 */
inline void InterfaceResidueDescriptor::operator=(InterfaceResidueDescriptor & _ird) { copy(_ird); }

/**
 * This method will set the PDB ID that this residue came from.
 *
 * @param _pdbId  The PDB ID.  This is typically a 4 digit alpha-numeric code.
 */
	inline void   InterfaceResidueDescriptor::setPdbId(std::string _pdbId) { pdbId = _pdbId; }
/**
 * This method will get the PDB ID for this residue.
 *
 * @return The PDB ID.  This is typically a 4 digit alpha-numeric code.
 */
inline std::string InterfaceResidueDescriptor::getPdbId() { return pdbId; }

/**
 * This method will set the chain ID for this residue.  This information 
 * comes from the PDB.
 *
 * @param _chainId The Chain ID of this residue. Typically this is a single letter.
 */
inline void   InterfaceResidueDescriptor::setChainId(std::string _chainId) { chainId = _chainId; }
/**
 * This method will get the chain ID for this residue.
 *
 * @return The chain ID.  Typically this is a single letter.
 */
inline std::string InterfaceResidueDescriptor::getChainId() { return chainId; }

/**
 * This method will set the resolution of this PDB structure.  It will attempt
 * to translate the resolution value in the std::string into a double.  If it
 * is unable to do so (the resolution is "NOT DEFINED" for instance, then
 * it will set a flag indicating that the resolution is not applicable.
 * 
 * @param _resolution  A std::string giving the resolution of this PDB.  
 */
inline void InterfaceResidueDescriptor::setResolution(std::string _resolution) {
	try{
		resolution = MslTools::toDouble(_resolution, "Converting resolution std::string");
		resolutionIsValid = true;
	} catch(...) {
		resolutionIsValid = false;
	}
}

/**
 * This method will set the resolution of this PDB structure.
 *
 * @param _resolution A double giving the resolution of this PDB.
 */
inline void InterfaceResidueDescriptor::setResolution(double _resolution) {
	resolution = _resolution;
	resolutionIsValid = true;
}

/**
 * This method will get the resolution of the PDB structure that this
 * residue came from.  If no resolution was given, then this 
 * method will return a negative value.  Note that users should probably
 * call isResolutionValid before calling getResolution.
 *
 * @return Returns the resolution of the PDB structure if a resolution was given.
 * If no resolution was given, or resolution is not applicable, this method
 * returns a negative number.
 */
inline double InterfaceResidueDescriptor::getResolution() {
	if( resolutionIsValid ) {
		return resolution;
	} else {
		return -1.0f;
	}
}

/**
 * This method will tell the user whether the resolution
 * given for this residue is valid or not.  PDB's created
 * from NMR data for instance don't give resolution values,
 * so this method would return false.
 *
 * @return Returns true if the resolution information is valid
 * and false otherwise.
 */
inline bool InterfaceResidueDescriptor::isResolutionValid() { return resolutionIsValid; }

/**
 * This method sets the name of the Amino Acid for this residue.
 *
 * @param _residueName The name of the amino acid for this residue.
 */
inline void InterfaceResidueDescriptor::setResidueName(std::string _residueName) { residueName = _residueName; }

/**
 * This method gets the name of the Amino Acid for this residue.
 *
 * @return A std::string giving the one letter code for the Amino Acid of this residue.
 */
inline std::string InterfaceResidueDescriptor::getResidueName() { return residueName; }

/**
 * This method sets the residue number of this residue, as found in the PDB.
 *
 * @param _residueNumber  The residue number of this residue as found in the PDB.
 */
inline void InterfaceResidueDescriptor::setResidueNumber(int _residueNumber) { residueNumber = _residueNumber; }
/**
 * This method gets the residue number for this residue.  The residue number is the
 * number assigned to this residue in the PDB file.
 *
 * @return The residue number as found in the PDB.
 */
inline int  InterfaceResidueDescriptor::getResidueNumber() { return residueNumber; }

/**
 * This method sets the icode for this residue, as found in the PDB.
 * Some pdbs have multiple residues with the same residue number and
 * chain id.  The icode is used to further differentiate residues.
 *
 * @param _icode A std::string of the icode used in the pdb for this residue.
 *  icodes are usually single letters like a, b, or c.
 */
inline void InterfaceResidueDescriptor::setResidueIcode(std::string _icode) { residueIcode = _icode; }

/**
 * This method returns the icode for this residue.  Some pdbs
 * have multiple residues with the same residue number and chain id.
 * The icode is used by the pdb to further differentiate residues.
 *
 * @return A std::string giving the icode for this residue.  Usually the
 *  icode is a single letter like a, b, or c.
 */
inline std::string InterfaceResidueDescriptor::getResidueIcode() { return residueIcode; }

/**
 * This method sets the secondary structure of this residue. This information comes 
 * from running the DSSP program on the particular PDB file. Each residue is assigned 
 * a code which describes its secondary structure.
 *
 * @param _secondaryStructure A std::string indicating the secondary structure of this residue. Values include: 
 * H = alpha helix, 
 * B = residue in isolated beta-bridge, 
 * E = extended strand, participates in beta ladder, 
 * G = 3-helix, 
 * I = 5 helix (pi helix), 
 * T= hydrogen bonded turn, 
 * S = bend. 
 */

inline void   InterfaceResidueDescriptor::setSecondaryStructure(std::string _secondaryStructure) { secondaryStructure = _secondaryStructure; }
/**
 * This method gets the secondary structure of this residue.  This information
 * came from running the DSSP program on the particular PDB file.
 *
 * @return A std::string indicating the secondary structure of this residue.
 * @see getResidueNumber for a description of the codes used to 
 * indicate secondary structure.
 */
inline std::string InterfaceResidueDescriptor::getSecondaryStructure() { return secondaryStructure; }

/**
 * This method sets the delta solvent accessibility of this particular residue.
 * To calculate this value, DSSP was run on the PDB file and each residue's solvent
 * accessibility was noted.  Then, a separate pdb was created with only residues
 * from the given chain.  The solvent accessibility was again calculated for each residue
 * using DSSP.  The difference between the solvent accessibility of each residue
 * by itself and when its part of the full complex was noted.  This difference is stored
 * in the singleChainDeltaSolventAccesibility.
 *
 * @param _deltaSolventAcc  The delta solvent accessibility for this residue.
 */
inline void    InterfaceResidueDescriptor::setSingleChainDeltaSolventAccessibility(double _deltaSolventAcc) { singleChainDeltaSolventAccessibility = _deltaSolventAcc; }
/**
 * This method gets the single chain delta solvent accessibility for this residue.
 * 
 * @see setPeptideDeltaSolventAccessibility for a description on how this value is calculated.
 */
inline double  InterfaceResidueDescriptor::getSingleChainDeltaSolventAccessibility() { return singleChainDeltaSolventAccessibility; }

/**
 * This method sets the delta solvent accessibility of this particular residue,
 * when a given chain was withheld.  To calculate this value, DSSP was run on the
 * PDB file and each residue's solvent accessibility was noted.  Then separate
 * PDB's were made - each time holding out one of the chains.  The solvent
 * accessibility was again calculated for each residue using DSSP.  The difference
 * between the solvent accessibility of each residue
 * in this second constructed pdb and when its part of the full complex was noted.
 * This difference is stored in the otherChainDeltaSolventAccesibility.  Note, that
 * this value can obviously change depending upon which chain is being withheld.
 * For instance, let's say we have 3 chains, A, B, and C.  B is a globular
 * protein, while A and C are short peptides, but on opposite faces of B.  Then, when
 * A is the withheld chain, the residues of B which complex with A will have some neg. value
 * of delta solvent accessibility, while the residues of B which complex with C will have
 * a delta solvent accessibility of 0.  Likewise, all residues of C will have a delta
 * solvent accessibility of 0 when treating A as the withheld chain.
 *
 * @param _withheldChain Which chain is being treated as the peptide chain
 * @param _deltaSolventAcc  The delta solvent accessibility for this residue.
 */
inline void    InterfaceResidueDescriptor::setOtherChainDeltaSolventAccessibility(std::string _withheldChain, double _deltaSolventAcc) { otherChainDeltaSolventAccessibility[_withheldChain] = _deltaSolventAcc; }
/**
 * This method gets the delta solvent accessibility for this residue.
 *
 * @param _withheldChain Which chain is being treated as the peptide chain.
 * @see setOtherChainDeltaSolventAccessibility for a description on how this value is calculated.
 */
inline double  InterfaceResidueDescriptor::getOtherChainDeltaSolventAccessibility(std::string _withheldChain) { return otherChainDeltaSolventAccessibility[_withheldChain]; }

/**
 * Using the ConSurf database, we collected information about observed amino acid identities at each given position.
 * We store a std::map which gives the fraction of times a given amino acid is seen at each position.
 *
 * @param _aaName  The name of the amino acid.
 * @param _prob    The fraction of times that amino acid is seen at this position.
 */
inline void    InterfaceResidueDescriptor::setAAProb(std::string _aaName, double _prob) { aaProbs[_aaName] = _prob; }

/**
 * Using the ConSurf database, we collected information about observed amino acid identities at each given position.
 * We store a std::map which gives the fraction of times a given amino acid is seen at each position.
 *
 * @param _aaName  The name of the amino acid.
 */
inline double  InterfaceResidueDescriptor::getAAProb(std::string _aaName) { return aaProbs[_aaName]; }

/**
 * Using the ConSurf database, we collected information about how conserved each
 * residue is.  We'll store that info.
 *
 * @param _score    The ConSurf score for this residue.  The score is normalized so
 *                  that the average score is 0 with a std. dev. of 1.  A negative
 *                  score indicates conserved position while a positive score
 *                  indicates a variable position.
 */
inline void InterfaceResidueDescriptor::setConsScore(double _score) { consScore = _score; }

/**
 * Using the ConSurf database, we collected information about how conserved each
 * residue is.  
 *
 * @see setConsScore
 */
inline double  InterfaceResidueDescriptor::getConsScore() { return consScore; }

/**
 * Using the ConSurf database, we collected information about how conserved each
 * residue is.  ConSurf also calculates a confidence interval for the conservation score.
 * We'll store that info here.
 *
 * @param _low    The low range of the confidence interval for our conservation score.
 * @param _high   The high range of the confidence interval for our conservation score.
 */
inline void InterfaceResidueDescriptor::setConsScoreConfInterval(double _low, double _high) { consScoreConfInterval.first = _low; consScoreConfInterval.first = _high; }

/**
 * Using the ConSurf database, we collected information about how conserved each
 * residue is.  ConSurf also calculates a confidence interval for the conservation score.
 * We'll store that info here.
 *
 * @param _low    The low range of the confidence interval for our conservation score.
 * @param _high   The high range of the confidence interval for our conservation score.
 */
inline std::pair<double, double> InterfaceResidueDescriptor::getConsScoreConfInterval() { return consScoreConfInterval; }

/**
 * This method sets the descriptor of the chain type of which this residue is a part.
 *
 * @param _type  The chain type of this residue.
 * @see ChainType_t for a description of the possible chain types.
 */
inline void InterfaceResidueDescriptor::setChainType(ChainType_t _type) { chainType = _type; }
/**
 * This method gets the descriptor fo the chain type for this residue.
 *
 * @return The chain type of this residue.
 * @see ChainType_t for a description of the possible chain types.
 */
inline ChainType_t InterfaceResidueDescriptor::getChainType() { return chainType; }

/**
 * This method allows the user to set the number of backbone to 
 * backbone contacts this residue makes.  Note: Currently for residues
 * where this information is unavailable, the default value is left to 0.
 *
 * @param _numberBackboneContacts  The number of backbone-backbone contacts this residue makes.
 */
inline void InterfaceResidueDescriptor::setNumberBackboneContacts(int _numberBackboneContacts) { numberBackboneContacts = _numberBackboneContacts; }
/**
 * This method gets the number of backbone to backbone contacts
 * made by this residue.  Note:  Currently for residues where this information
 * is unavailable, the default value is left to 0.
 *
 * @return The number of backbone-backbone contacts.
 */
inline int  InterfaceResidueDescriptor::getNumberBackboneContacts() { return numberBackboneContacts; }

/**
 * This method sets the number of side-chain to side-chain contacts
 * made by this residue.  Note:  Currently for residues where this
 * information is unavailable, the default value is left to 0.
 *
 * @param _numberSidechainContacts  The number of side-chain to side-chain contacts this residue makes.
 */
inline void InterfaceResidueDescriptor::setNumberSidechainContacts(int _numberSidechainContacts) { numberSidechainContacts = _numberSidechainContacts; }
/**
 * This method gets the number of side-chain to side-chain contacts
 * made by this residue.  Note:  Currently for residues where this
 * information is unavailable, the default value is left to 0.
 *
 * @return The number of side-chain to side-chain contacts this residue makes.
 */
inline int  InterfaceResidueDescriptor::getNumberSidechainContacts() { return numberSidechainContacts; }

/**
 * This method sets the number of mixed contacts (backbone to side-chain
 * or side-chain to backbone) contacts made by this residue.  Note:
 * Currently for residues where this information is unavailable, the 
 * default value is left to 0.
 *
 * @param _numberMixedContacts  The number of mixed contacts this residue makes.
 */
inline void InterfaceResidueDescriptor::setNumberMixedContacts(int _numberMixedContacts) { numberMixedContacts = _numberMixedContacts; }
/**
 * This method gets the number of mixed contacts (backbone to side-chain
 * or side-chain to backbone) contacts made by this residue.  Note:
 * Currently for residues where this information is unavailable, the 
 * default value is left to 0.
 *
 * @return The number of mixed contacts this residue makes.
 */ 
inline int  InterfaceResidueDescriptor::getNumberMixedContacts() { return numberMixedContacts; }

/**
 * This method will add a misc. description of this residue to a table of
 * descriptions.  Each description has a name, and then a std::string giving
 * the description.
 *
 * @param _key  The name of this description.  For example "function"
 * may describe the function of this protein complex.
 * @param _value A std::string that holds the description for this _key.
 */
inline void InterfaceResidueDescriptor::addToDescriptionTable(std::string & _key, std::string & _value) { descriptionTable[_key] = _value; }
/**
 * This method will get the description for a given key from the description table
 * for this residue.
 * 
 * @param _key The key to use to look up our description.
 * @return The description stored for the given key.
 * @see addToDescriptionTable for more information on what information
 * the description table is expected to hold.
 */
inline std::string InterfaceResidueDescriptor::getFromDescriptionTable(std::string &_key) { return descriptionTable[_key]; }
/**
 * This method will remove the description with a given key from the description table.  Note
 * That there will still be an entry for the given key, but the description will now be an empty std::string.
 *
 * @param _key The name of the description to delete.
 */
inline void InterfaceResidueDescriptor::removefromDescriptionTable(std::string &_key) { descriptionTable[_key] = "";}
/**
 * This method will clear the description table.  All entries and keys will be deleted.
 */
inline void InterfaceResidueDescriptor::clearDescriptionTable() { descriptionTable.clear(); }
/**
 * This method will return the std::map object that holds the description table.  
 *
 * @return the std::map object which holds the description table.  The key is a std::string
 * which holds the name of the descriptions, and the value is also a std::string which
 * holds the actual description.
 */
inline std::map<std::string,std::string> & InterfaceResidueDescriptor::getDescriptionTable() { return descriptionTable; }

/**
 * Create a std::string representation of object
 *
 */
inline std::string InterfaceResidueDescriptor::toString() {

	char c[100];
	if (residueIcode != ""){
		sprintf(c,"%4s %1s %-4d %1s %1s %1s %8.3f",pdbId.c_str(), chainId.c_str(), residueNumber, residueName.c_str(), residueIcode.c_str(), secondaryStructure.c_str(), singleChainDeltaSolventAccessibility);
	} else {
		sprintf(c,"%4s %1s %-4d %1s  %1s %8.3f",pdbId.c_str(), chainId.c_str(), residueNumber, residueName.c_str(), secondaryStructure.c_str(), singleChainDeltaSolventAccessibility);
	}
	
	return std::string(c);
}
}

#endif
