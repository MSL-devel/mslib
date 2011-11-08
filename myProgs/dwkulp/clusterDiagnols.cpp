//STL Includes
#include<fstream>
#include<string>
#include<vector>
#include<iostream>

//MSL Includes

#include "System.h"
#include "MatrixWindow.h"
#include "DistanceMatrix.h"
#include "DistanceMatrixDatabase.h"
#include "OptionParser.h"
#include "PyMolVisualization.h"
#include "MslTools.h"
#include "clusterDiagnols.h"

using namespace std;
using namespace MSL;
using namespace MslTools;

// R Includes
#ifdef __R__
   #include "Rcpp.h"
   #include "RInside.h"
#endif

#ifdef __R__
   RInside R;
#endif


int main(int argc, char *argv[]){

    // Option Parser
    Options opt = setupOptions(argc,argv);

    DistanceMatrixDatabase dmd;
    dmd.load_checkpoint(opt.pdb_binary_database);

    vector<DistanceMatrix *> &dms = dmd.getDistanceMatrixList();

    vector<MatrixWindow *> allDiagMWs(0,NULL);
    map<int,bool> dmAlreadyUsed;
    RandomNumberGenerator rng;
    rng.setSeed(52341);
    for (uint i = 0;i < dms.size();i++){
    //for (uint i = 0;i < opt.numDMs;i++){
            //int randomIndex = rng.getRandomInt(dms.size()-1);
	    //while (!dmAlreadyUsed[randomIndex]){
	    //  randomIndex = rng.getRandomInt(dms.size()-1);
	    //}
           int randomIndex = i;
	    dmAlreadyUsed[randomIndex] = true;

	    cout << "DM["<<randomIndex<<"]: "<<dms[randomIndex]->getPDBid()<<endl;

	    vector<MatrixWindow*> &diags = dms[randomIndex]->getDiagnolMWs();
	    
	    allDiagMWs.insert(allDiagMWs.end(),diags.begin(),diags.end());

	    cout << "Diag size: "<<allDiagMWs.size()<<endl;
    }

    cout << "Create Distance Matrix:"<<endl;
    vector<vector<double> > distMatrix;

    for (uint i = 0; i < allDiagMWs.size();i++){
      vector<double> rowMatrix;

      for (uint j = 0; j < allDiagMWs.size();j++){
	double distSq = allDiagMWs[i]->compare(*allDiagMWs[j]);
	rowMatrix.push_back(distSq);
      }
      distMatrix.push_back(rowMatrix);
    }

    cout << "K-means clustering..."<<endl;
    vector<int> clusterInfo;
    clusterInfo = getClusters(distMatrix);

    cout << "Output.."<<endl;
    PyMolVisualization pymol;
    vector<vector<double> > colors;

    if (opt.pymol){


#ifdef __R__
    
      R["v"] = clusterInfo;
    
      string evalStr = "rgb1=col2rgb(densCols(v,nbin=10,colramp = colorRampPalette(c(\"red\",\"orange\",\"yellow\",\"blue\"))))[1,]/255;rgb2=col2rgb(densCols(v,nbin=10,colramp = colorRampPalette(c(\"red\",\"orange\",\"yellow\",\"blue\"))))[2,]/255;rgb3=col2rgb(densCols(v,nbin=10,colramp = colorRampPalette(c(\"red\",\"orange\",\"yellow\",\"blue\"))))[3,]/255;m=cbind(rgb1,rgb2);m = cbind(m,rgb3);mat=as.matrix(m);print(mat);mat";
    
      SEXP ans = R.parseEval(evalStr);

      Rcpp::NumericMatrix M(ans);

      for (uint i = 0; i < clusterInfo.size();i++){
	vector<double> color;
	for (uint j = 0; j < 3;j++){
	  color.push_back(M(i,j));
	}
	colors.push_back(color);
      }

#endif

    }

    int j = 0;
    for (uint i = 0; i < clusterInfo.size();i++){

      DistanceMatrix &dm    = allDiagMWs[i]->getMatrix();
      int residue1index     = allDiagMWs[i]->getLeftR();
      int residue2index     = allDiagMWs[i]->getLeftR()+dm.getGeneralWinSize()-1;
      AtomPointerVector &av = dm.getAtomVector();

      Atom *ca1 = av[residue1index];
      Atom *ca2 = av[residue2index];

      // Obtain a SSE string for this window?
      stringstream sse;
      for (uint k = residue1index; k <= residue2index;k++){
	if (av[k]->getSegID().length() >= 1){
	  sse << av[k]->getSegID().substr(0,1);
	}
      }

      fprintf(stdout,"%4s %1s %4d-%-4d %4d-%-4d %-4d %s\n",
	      MslTools::getFileName(dm.getPDBid()).c_str(),
	      ca1->getChainId().c_str(),
	      ca1->getResidueNumber(),
	      ca2->getResidueNumber(),
	      allDiagMWs[i]->getLeftR(),
	      allDiagMWs[i]->getLeftR()+dm.getGeneralWinSize()-1,
	      clusterInfo[i],
	      sse.str().c_str());
	      

      if (opt.pymol){
	  if (i > 0 && dm.getPDBid() == allDiagMWs[i-1]->getMatrix().getPDBid()){
	    j += 1;
	    char c[80];
	    sprintf(c,"Segment%03d",j);

	    pymol.createCylinder(ca1->getCoor(),ca2->getCoor(),(string)c,0.5,colors[i][0],colors[i][1],colors[i][2]);
	  }

      }
      
    }

    if (opt.pymol){
      ofstream fout;
      fout.open("clusterPymol.py");
      fout << pymol;
      fout <<endl;
      fout.close();
    }
    cout << "Done"<<endl;
}

vector<int> getClusters(vector<vector<double> > &_distMatrix) {
  vector<int> clusters;

  int minClusters = 5;
  int maxClusters = (int) _distMatrix.size()/3;
  int stepSize    = 1;

#ifdef __R__
    // Start instance of R
    //RInside R;
  cout << "Convert to NumericMatrix"<<endl;
    Rcpp::NumericMatrix M(_distMatrix.size(),_distMatrix.size());
    for (uint i = 0; i<_distMatrix.size();i++){
      for (uint j=0; j<_distMatrix.size();j++){
	M(i,j) = _distMatrix[i][j];
      }
    }
    R["M"] = M;
      cout << "Convert to create R-string"<<endl;
    stringstream evalStr;
    evalStr << "library(cluster);";

    //evalStr << "nClusters=c(seq(2,10,1),seq(20,500,10));";

    evalStr << "nClusters=c(seq("<<minClusters<<","<<maxClusters<<","<<stepSize<<"));";
    evalStr << "clust <- {};";
    evalStr << "k.width <- {};";
    evalStr << "for (k in 1:length(nClusters)) { ";
    evalStr <<     "clust <- rbind(clust,pam(M,nClusters[k]));";
    evalStr <<     "k.width <- rbind(k.width,clust[k,]$silinfo$avg.width);";
    evalStr << "};";
    evalStr << "k.best <- which.max(k.width);";
    evalStr << "print(paste(\"Best sized cluster:\",nClusters[k.best],sep=\" \"));";
    evalStr << "png(filename=\"clusterSilhouettes.png\");";
    evalStr << "plot(nClusters,k.width,type=\"h\",main=\"pam() clustering assesment\",xlab=\"k (#clusters)\",ylab=\"average silhouette width\");";
    evalStr << "dev.off();";
    //evalStr << "print(clust[k.best,]$clustering);";
    evalStr << "clust[k.best,]$clustering;";

    
      cout << "Run clustering..."<<endl;
    // Run clustering
    SEXP ans = R.parseEval(evalStr.str());

      cout << "Convert answer"<<endl;
    // Convert answer to proper format
    Rcpp::NumericVector v(ans);
    vector<int> result;
    for (uint i = 0; i < v.size();i++){
      result.push_back(v[i]);
    }
    return result;

#endif
    return clusters;
}

Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;

	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "clusterDiagnols --pdb_db BINARY_PDB_DB.bin\n";
		exit(0);
	}

	opt.pdb_binary_database = OP.getString("pdb_db");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb_db not specified.\n";
		exit(1111);
	}

	opt.numDMs = OP.getInt("numPdbs");
	if (OP.fail()){
	  opt.numDMs = 10;
	}

	opt.pymol = OP.getBool("pymol");
	if (OP.fail()){
	  opt.pymol = false;
	}
	return opt;
}


