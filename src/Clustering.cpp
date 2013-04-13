
/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2013 The MSL Developer Group (see README.TXT)
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


/*
  Clustering Object:
  K-medoid
 */

#include "Clustering.h"

using namespace MSL;

#include "MslOut.h"
static MslOut MSLOUT("Clustering");


#ifdef __GSL__
  #include <gsl/gsl_statistics.h>
#endif

Clustering::Clustering(){
  elementLabels = NULL;
  distMatrix    = NULL;
  markMaximumCutoff = 0;
  markVariarianceCutoff = 0;
  markFlag = false;
  nonMatchingValue = 0;
  noMatchFlag = false;
  RNG.setTimeBasedSeed();
  RNG.setRNGType("knuth2");
  cout << "\tRNG Type: "<<RNG.getRNGType()<<" Seed: "<<RNG.getSeed()<<endl;
}
Clustering::Clustering(vector<vector< double > > *_distMatrix,vector<string> *_labels){
   distMatrix    = _distMatrix;
   elementLabels = _labels;
   markMaximumCutoff = 0;
   markVariarianceCutoff = 0;
   markFlag = false;
   nonMatchingValue = 0;
   noMatchFlag = false;

   RNG.setTimeBasedSeed();	
   RNG.setRNGType("knuth2");
   cout << "\tRNG Type: "<<RNG.getRNGType()<<" Seed: "<<RNG.getSeed()<<endl;
}

Clustering::Clustering(vector<vector< double > > *_distMatrix){
   elementLabels = NULL;
   distMatrix = _distMatrix;
   markMaximumCutoff = 0;
   markVariarianceCutoff = 0;
   markFlag = false;
   nonMatchingValue = 0;
   markFlag = false;

   RNG.setTimeBasedSeed();
   RNG.setRNGType("knuth2");	
   cout << "\tRNG Type: "<<RNG.getRNGType()<<" Seed: "<<RNG.getSeed()<<endl;
}
Clustering::~Clustering(){
}


void Clustering::Kmedoids(int _nIterations, int _nClusters){

  cout << endl;
  nClusters = _nClusters;

  vector<int> tmpCentroids;
  tmpCentroids.resize(nClusters);
  vector<int> tmpClusterAssignments;
  double globalSumDistToCentroids = 100000;

  // Initially set the tmpClusterAssignments from clusterAssignments
  if (clusterAssignments.size() != distMatrix->size()) clusterAssignments.resize(distMatrix->size());
  tmpClusterAssignments.resize(distMatrix->size());

  int iteration = 0;
  do {

    double sumDistToCentroids = 100000;


    // Generate random cluster assignments
    //cout << "Generating Random Clusters"<<endl;
    generateRandomClusters(&tmpClusterAssignments);

    /*
    for (int iCl = 0; iCl < nClusters; iCl++){
      cout << "\t\tCluster "<<iCl<<endl;
      for (int iEl = 0; iEl < distMatrix->size(); iEl++){
      if (tmpClusterAssignments[iEl] != iCl) continue;
      cout << "\t\t\tMember: "<<iEl<<" "<<tmpClusterAssignments[iEl]<<endl;
      }
      cout << endl;
    }
    */
  

    while(true){

      // Store previous result and clear sumDistToCentroids.
      double previousDistToCentroids = sumDistToCentroids;
      sumDistToCentroids = 0;

      // Find centroids
       getClusterCentroids(&tmpClusterAssignments);

      // For each element
      for (uint iEl = 0; iEl < distMatrix->size();iEl++){

      if (elementFlag[iEl] == "frozen") continue;

	double distToCentroid = 100000;

	// Skip if element frozen 

	bool noMatch = true; 

	// Find closest cluster
	for (uint iCl = 0; iCl < nClusters; iCl++){

	  // No adding to frozen clusters
 	  if (clusterFlag[iCl] == "frozen") continue;

	  // Skip cluster assignment if we are a centroid
	  if (iEl == centroids[iCl]){
	    distToCentroid = 0;
	    tmpClusterAssignments[iEl] = iCl;
	    break;
	  }

	  
	  // Since we have a symmetric matrix, that is lower-triangular
	  double distToThisCentroid = getMatrixValue(iEl,centroids[iCl]);
	  

	  if (noMatchFlag && fabs(distToThisCentroid - nonMatchingValue) < 0.01){
		  continue;
	  }


	  // If closer centroid, update distToCentroid and update clusterAssignment table
	  if (distToThisCentroid < distToCentroid){

	    distToCentroid = distToThisCentroid;
	    tmpClusterAssignments[iEl] = iCl;
	    noMatch = false;
	  }
	}

	// Check for a non-matching entry
	// Don't add to sumDist
	// Add cluster Assignment unmatched cluster
	if (!noMatch){
		sumDistToCentroids += distToCentroid;
	} else {
		tmpClusterAssignments[iEl] = -1;
		sumDistToCentroids += nonMatchingValue;
	}

      } // END for(elements)

      

      if (sumDistToCentroids >= previousDistToCentroids) break;
    } // WHILE (TRUE)
        



    if (sumDistToCentroids < globalSumDistToCentroids) {
      cout << "  SumDistToCentroids";
      cout <<"*";
      fprintf(stdout," %4d % 6.10f\n",iteration,sumDistToCentroids);

      globalSumDistToCentroids = sumDistToCentroids;
      for (uint iEl = 0; iEl < distMatrix->size();iEl++) {
	clusterAssignments[iEl] = tmpClusterAssignments[iEl];
      }
      for (uint i = 0; i < nClusters;i++){	
	      tmpCentroids[i] = centroids[i];
      }
    }
    iteration++;
  }   while(iteration < _nIterations);
  
  for (uint i = 0; i < nClusters;i++){
	  centroids[i] = tmpCentroids[i];
  }

  // Set centroids for the 'winning' clusters  ... This can re-define the centroids. Which is bad because if an entity was placed in given cluster due to the centroid, then the centroid changes it may not match to the new centroid at aall..... hmmmm....
  // getClusterCentroids(&clusterAssignments);

}

void Clustering::setElementsFrozen(){

     if (elementFlag.size() != distMatrix->size()){
		elementFlag.resize(distMatrix->size());
     }

     for (uint iEl = 0; iEl < distMatrix->size();iEl++){
	     if (clusterAssignments[iEl] < 0) continue;
	      if (clusterFlag[clusterAssignments[iEl]] == "frozen"){
		      elementFlag[iEl] = "frozen";
	      }
     }

}
void Clustering::setClustersFrozen(double _maxValue, double _sizePercent){

	cout << "\n***********\nCLUSTER REPORT\n";
	getStatistics();
	if (clusterFlag.size() != nClusters){
		clusterFlag.resize(nClusters);
	}
	for (int iCl = 0; iCl < nClusters; iCl++){
		if (clusterStats[iCl].max <= _maxValue && (clusterStats[iCl].prob *100) >= _sizePercent){

			// Freeze this cluster (centroid).
			clusterFlag[iCl] = "frozen";

			fprintf(stdout, "\tFROZEN: Cluster %5d (%5d) max %8.3f prob %8.3f\n",iCl,centroids[iCl],clusterStats[iCl].max, (clusterStats[iCl].prob *100));
		} else {
			fprintf(stdout, "\tRETRY : Cluster %5d (%5d) max %8.3f prob %8.3f\n",iCl,centroids[iCl],clusterStats[iCl].max, (clusterStats[iCl].prob *100));
			clusterFlag[iCl] = "";

		}
	}
	cout << endl;	
}
void Clustering::setClustersFlag(string _flag){


	if (clusterFlag.size() != nClusters){
		clusterFlag.resize(nClusters);
	}
	for (int iCl = 0; iCl < nClusters; iCl++){
		clusterFlag[iCl] = _flag;
	}
}
double Clustering::getStatistics(){
  
  maxStats.variance  = 0;
  maxStats.max       = 0;
  maxStats.mean      = 0;
  clusterStats.clear();
  clusterStats.resize(nClusters);
  for (int iCl = 0; iCl < nClusters; iCl++){
  
    vector<double> vdata;
    for (int iEl = 0; iEl < distMatrix->size(); iEl++){
      if (clusterAssignments[iEl] != iCl) continue;
      vdata.push_back(getMatrixValue(iEl, centroids[iCl]));
    }


    
    int size = vdata.size();
    double data[size];
    for (uint i = 0 ; i < size;i++){
      data[i] = vdata[i];
    }
#ifdef __GSL__

    stats s;
    s.mean = gsl_stats_mean(data,1,size);


    s.variance = gsl_stats_variance(data,1,size);


    s.std = gsl_stats_sd(data,1,size);


    gsl_stats_minmax(&s.min,&s.max,data,1,size);

    
    s.prob = (double)vdata.size() / (double)distMatrix->size();


    clusterStats[iCl] = s;

    if (s.variance > maxStats.variance) maxStats.variance = s.variance;
    if (s.max > maxStats.max)           maxStats.max      = s.max;
    if (s.mean > maxStats.mean)		maxStats.mean     = s.mean;

#endif

  }

  return maxStats.variance;
}


void Clustering::printClusters(double numSTDs){

  cout << endl;
  for (int iCl = 0; iCl < nClusters; iCl++){
    cout << "\tCluster "<<iCl<<" centroid is ("<<centroids[iCl]<<")";
    if (elementLabels !=NULL) cout << " "<<(*elementLabels)[centroids[iCl]];
    cout <<endl;

    for (int iEl = 0; iEl < distMatrix->size(); iEl++){
    if (clusterAssignments[iEl] != iCl) continue;
      double distToCentroid = getMatrixValue(iEl, centroids[iCl]);
      if (iEl == centroids[iCl]) distToCentroid = 0;

      if ((numSTDs == 0) || 
	  (clusterStats.size() > iCl &&
	   distToCentroid > clusterStats[iCl].mean - numSTDs*clusterStats[iCl].std && 
	   distToCentroid < clusterStats[iCl].mean + numSTDs*clusterStats[iCl].std
	  )){
	cout << "\t\tMember: ("<<iEl<<")";
	if (elementLabels != NULL) cout << " "<<(*elementLabels)[iEl];
	cout<<" measurement to centroid: "<<distToCentroid<<endl;
      }
      
    }

    if (clusterStats.size() > iCl){
      cout << "\t\tStatistics "<< clusterStats[iCl].toString();
    }
    cout << endl<<endl;
  } // END nClusters

  // Now print out non-matching
  cout << "\tUnmatched: "<<endl;
  
  for (int iEl = 0; iEl < distMatrix->size(); iEl++){
	  if (clusterAssignments[iEl] != -1) continue;
	  cout << "\t\tMember: ("<<iEl<<")";
	  if (elementLabels != NULL) cout << " "<<(*elementLabels)[iEl];
	  cout <<endl;
  }
	
  


}
void Clustering::printClusters(){
  printClusters(0);
}

/* Algorithms yet to be implemented, take a shot if you want.. */
void Clustering::Kmeans(){
}
void Clustering::SingleLinkage(){
   // Written by Jason Donald, based on code of Cinque Soto
   uint numElements = distMatrix->size();
   vector<vector <double> > clusteringMatrix(numElements);

   MSLOUT.stream() << "Setup clustering matrix"<<endl;
   for (uint i = 0; i < numElements; i++)
   {
      clusteringMatrix[i].resize(i);
      for (uint j = 0; j < i; j++)
      {
         clusteringMatrix[i][j] = (*distMatrix)[i][j];
      }
   }

   MSLOUT.stream() << "Setup cluster id and markers"<<endl;   
   vector<int> clusterid(numElements, 0);
   for (uint i = 0; i < clusterid.size(); i++) { clusterid[i] = i; }
   vector<int> clusterMarkers(numElements, 0);
   for (uint i = 0; i < clusterMarkers.size(); i++) { clusterMarkers[i] = i; }

   vector<double> linkdist(numElements, 0.);

   MSLOUT.stream() << "Linking distance loop"<<endl;
   for (uint numNodes = numElements; numNodes > 1; numNodes--)
   {
      // Get the smallest distance left
      int isaved = 1;
      int jsaved = 0;
      double distance = clusteringMatrix[isaved][jsaved];
      for (uint i = 0; i < numNodes; i++)
      {
         for (uint j = 0; j < i; j++)
         {
            if (clusteringMatrix[i][j] < distance)
            {
               isaved = i;
               jsaved = j;
               distance = clusteringMatrix[i][j];
            }
         }
      }
      linkdist[(numElements-numNodes)] = distance;
      
      // Fix the distances by merging lines
      // Replace smaller cluster position row with statistics of the new cluster
      for (uint k = 0; k < jsaved; k++)
      {
         clusteringMatrix[jsaved][k] = min(clusteringMatrix[isaved][k], clusteringMatrix[jsaved][k]);
      }
      // Replace smaller cluster position column with statistics of the new cluster
      for (uint k = jsaved+1; k < isaved; k++)
      {
         clusteringMatrix[k][jsaved] = min(clusteringMatrix[isaved][k],clusteringMatrix[k][jsaved]);
      }
      // Replace smaller cluster position column with statistics of the new cluster, below new larger position
      for (uint k = isaved+1; k < numNodes; k++)
      {
        clusteringMatrix[k][jsaved] = min(clusteringMatrix[k][isaved],clusteringMatrix[k][jsaved]);
      }

      // Get rid of isaved row with last row
      for (uint k = 0; k < isaved; k++)
      {
         clusteringMatrix[isaved][k] = clusteringMatrix[numNodes-1][k];
      }
      // Get rid of isaved column with last column
      for (uint k = isaved+1; k < numNodes-1; k++)
      {
         clusteringMatrix[k][isaved] = clusteringMatrix[numNodes-1][k];
      }

      /* Update clusterids */
      // For each merger, move any clusters associated with i to number j
      for (uint k = 0; k < numElements; k++)
      {
         if (clusterMarkers[k] == isaved) { clusterMarkers[k] = jsaved; }
         else if (clusterMarkers[k] == (numNodes-1)) { clusterMarkers[k] = isaved; }
      }
      uint numClusters = numNodes - 1;

      MSLOUT.stream() << "Printing clusters"<<endl;
      // Arbitrary printing parameter
      //if (numClusters <= 15) {
         // Now print the merger
         cout << "\nNumber of clusters: " << numClusters << "\tAt cutoff: " << distance << endl;
         for (uint k = 0; k < numClusters; k++)
         {
            cout << "Cluster number " << k << endl;
            for (uint m = 0; m < numElements; m++)
            {
               if (clusterMarkers[m] == k) {

		       cout << "\t" << m;
		       if (elementLabels != NULL) cout << " "<<(*elementLabels)[m];
		       cout <<endl;

               }
            }
         }
      //}
      clusterid[jsaved] = numNodes-numElements-1;
      clusterid[isaved] = clusterid[numNodes-1];
   }
}
void Clustering::AverageLinkage(){
   // Written by Jason Donald, based on code of Cinque Soto
   uint numElements = distMatrix->size();
   vector<vector <double> > clusteringMatrix(numElements);
   for (uint i = 0; i < numElements; i++)
   {
      clusteringMatrix[i].resize(i);
      for (uint j = 0; j < i; j++)
      {
         clusteringMatrix[i][j] = (*distMatrix)[i][j];
      }
   }
   vector<int> clusterid(numElements, 0);
   for (uint i = 0; i < clusterid.size(); i++) { clusterid[i] = i; }
   vector<int> clusterMarkers(numElements, 0);
   for (uint i = 0; i < clusterMarkers.size(); i++) { clusterMarkers[i] = i; }
   vector<int> numMembers(numElements, 1);

   vector<double> linkdist(numElements, 0.);

   for (uint numNodes = numElements; numNodes > 1; numNodes--)
   {
      // Get the smallest distance left
      int isaved = 1;
      int jsaved = 0;
      double distance = clusteringMatrix[isaved][jsaved];
      for (uint i = 0; i < numNodes; i++)
      {
         for (uint j = 0; j < i; j++)
         {
            if (clusteringMatrix[i][j] < distance)
            {
               isaved = i;
               jsaved = j;
               distance = clusteringMatrix[i][j];
            }
         }
      }
      linkdist[(numElements-numNodes)] = distance;
      
      // Fix the distances by merging lines
      // Replace smaller cluster position row with statistics of the new cluster
      for (uint k = 0; k < jsaved; k++)
      {
         clusteringMatrix[jsaved][k] = avg(clusteringMatrix[isaved][k], clusteringMatrix[jsaved][k], (numMembers[k]*numMembers[isaved]), (numMembers[k]*numMembers[jsaved]));
      }
      // Replace smaller cluster position column with statistics of the new cluster
      for (uint k = jsaved+1; k < isaved; k++)
      {
         clusteringMatrix[k][jsaved] = avg(clusteringMatrix[isaved][k],clusteringMatrix[k][jsaved], (numMembers[k]*numMembers[isaved]), (numMembers[k]*numMembers[jsaved]));
      }
      // Replace smaller cluster position column with statistics of the new cluster, below new larger position
      for (uint k = isaved+1; k < numNodes; k++)
      {
        clusteringMatrix[k][jsaved] = avg(clusteringMatrix[k][isaved],clusteringMatrix[k][jsaved], (numMembers[k]*numMembers[isaved]), (numMembers[k]*numMembers[jsaved]));
      }

      // Get rid of isaved row with last row
      for (uint k = 0; k < isaved; k++)
      {
         clusteringMatrix[isaved][k] = clusteringMatrix[numNodes-1][k];
      }
      // Get rid of isaved column with last column
      for (uint k = isaved+1; k < numNodes-1; k++)
      {
         clusteringMatrix[k][isaved] = clusteringMatrix[numNodes-1][k];
      }

      /* Update clusterids */
      // For each merger, move any clusters associated with i to number j
      for (uint k = 0; k < numElements; k++)
      {
         if (clusterMarkers[k] == isaved) { clusterMarkers[k] = jsaved; }
         else if (clusterMarkers[k] == (numNodes-1)) { clusterMarkers[k] = isaved; }
      }
      uint numClusters = numNodes - 1;

      // Arbitrary printing parameter
      //if (numClusters <= 15) {
         // Now print the merger
         cout << "\nNumber of clusters: " << numClusters << "\tAt cutoff: " << distance << endl;
         for (uint k = 0; k < numClusters; k++)
         {
            cout << "Cluster number " << k << endl;
            for (uint m = 0; m < numElements; m++)
            {
               if (clusterMarkers[m] == k) {
		       cout << "\t" << m;
		       if (elementLabels != NULL) cout << " "<<(*elementLabels)[m];
		       cout <<endl;
               }
            }
         }
      //}
      clusterid[jsaved] = numNodes-numElements-1;
      clusterid[isaved] = clusterid[numNodes-1];
      numMembers[jsaved] += numMembers[isaved];
      numMembers[isaved] = numMembers[(numNodes-1)];
   }
}
void Clustering::CenterLinkage(){
}
void Clustering::CompleteLinkage(){
   // Written by Jason Donald, based on code of Cinque Soto
   uint numElements = distMatrix->size();
   vector<vector <double> > clusteringMatrix(numElements);
   for (uint i = 0; i < numElements; i++)
   {
      clusteringMatrix[i].resize(i);
      for (uint j = 0; j < i; j++)
      {
         clusteringMatrix[i][j] = (*distMatrix)[i][j];
      }
   }
   vector<int> clusterid(numElements, 0);
   for (uint i = 0; i < clusterid.size(); i++) { clusterid[i] = i; }
   vector<int> clusterMarkers(numElements, 0);
   for (uint i = 0; i < clusterMarkers.size(); i++) { clusterMarkers[i] = i; }

   vector<double> linkdist(numElements, 0.);

   for (uint numNodes = numElements; numNodes > 1; numNodes--)
   {
      // Get the smallest distance left
      int isaved = 1;
      int jsaved = 0;
      double distance = clusteringMatrix[isaved][jsaved];
      for (uint i = 0; i < numNodes; i++)
      {
         for (uint j = 0; j < i; j++)
         {
            if (clusteringMatrix[i][j] < distance)
            {
               isaved = i;
               jsaved = j;
               distance = clusteringMatrix[i][j];
            }
         }
      }
      linkdist[(numElements-numNodes)] = distance;
      
      // Fix the distances by merging lines
      // Replace smaller cluster position row with statistics of the new cluster
      for (uint k = 0; k < jsaved; k++)
      {
         clusteringMatrix[jsaved][k] = max(clusteringMatrix[isaved][k], clusteringMatrix[jsaved][k]);
      }
      // Replace smaller cluster position column with statistics of the new cluster
      for (uint k = jsaved+1; k < isaved; k++)
      {
         clusteringMatrix[k][jsaved] = max(clusteringMatrix[isaved][k],clusteringMatrix[k][jsaved]);
      }
      // Replace smaller cluster position column with statistics of the new cluster, below new larger position
      for (uint k = isaved+1; k < numNodes; k++)
      {
        clusteringMatrix[k][jsaved] = max(clusteringMatrix[k][isaved],clusteringMatrix[k][jsaved]);
      }

      // Get rid of isaved row with last row
      for (uint k = 0; k < isaved; k++)
      {
         clusteringMatrix[isaved][k] = clusteringMatrix[numNodes-1][k];
      }
      // Get rid of isaved column with last column
      for (uint k = isaved+1; k < numNodes-1; k++)
      {
         clusteringMatrix[k][isaved] = clusteringMatrix[numNodes-1][k];
      }

      /* Update clusterids */
      // For each merger, move any clusters associated with i to number j
      for (uint k = 0; k < numElements; k++)
      {
         if (clusterMarkers[k] == isaved) { clusterMarkers[k] = jsaved; }
         else if (clusterMarkers[k] == (numNodes-1)) { clusterMarkers[k] = isaved; }
      }
      uint numClusters = numNodes - 1;

      // Now print the merger
      cout << "\nNumber of clusters: " << numClusters << "\tAt cutoff: " << distance << endl;
      for (uint k = 0; k < numClusters; k++){

	  cout << "Cluster number " << k+1 << endl;
	  for (uint m = 0; m < numElements; m++)
            {
               if (clusterMarkers[m] == k) {
		       cout << "\t" << m+1;
		       if (elementLabels != NULL) cout << " "<<(*elementLabels)[m];
		       cout <<endl;
               }
            }
	}
      clusterid[jsaved] = numNodes-numElements-1;
      clusterid[isaved] = clusterid[numNodes-1];
   }
}

void Clustering::generateRandomClusters(vector<int> *assignments){

  if (clusterFlag.size() != nClusters){
	  clusterFlag.resize(nClusters);
  }

  if (elementFlag.size() != distMatrix->size()){
		elementFlag.resize(distMatrix->size());
  }

  if (assignments->size() != distMatrix->size()){
    assignments->resize(distMatrix->size());
  } 

  assignments->clear();


  bool continueLoop;
  do{
    continueLoop = false;
    vector<bool> clusterNonEmptyFlags(nClusters,bool(false));

    // Insure all elements have been assigned
    for (uint iEl = 0; iEl < distMatrix->size();iEl++){

       // Default assignment
      (*assignments)[iEl] = clusterAssignments[iEl];

      // Get random cluster
      int randomCluster;
      do {
	      randomCluster = RNG.getRandomInt()%nClusters;    
      } while(clusterFlag[randomCluster] == "frozen");


      clusterNonEmptyFlags[randomCluster] = true;

      // Don't assign if element or cluster is frozen
      if (elementFlag[iEl] == "frozen") continue;


      // Assign random cluster
      (*assignments)[iEl] = randomCluster;
    }

    // Insure atleast one entry per cluster
    for (int iCl = 0; iCl < nClusters; iCl++){
      if (!(clusterNonEmptyFlags[iCl]) && clusterFlag[iCl] != "frozen") {
	continueLoop = true;
	break;
      }
    }
    
  }while(continueLoop);
    
}



void Clustering::getClusterCentroids(vector<int> *assignments){
  
  if (centroids.size() != nClusters){
    centroids.resize(nClusters);
  }

  if (clusterFlag.size() != nClusters){
	  clusterFlag.resize(nClusters);
  }
  vector<double> minSumOfDist(nClusters,1000000);
  for (int iCl = 0; iCl < nClusters; iCl++){

    // Skip clusters that are marked as frozen
    if (clusterFlag[iCl] == "frozen") continue;


    // Try each element as centroid
    for (uint i = 0; i < distMatrix->size();i++){

      // Only elements within cluster iCl
      if ((*assignments)[i] != iCl) continue;

      // Sum dist to this temporary centroid
      double sumOfDist = 0;
      for (uint j = 0; j < distMatrix->size();j++){

	// Only elments in cluster iCL, not element i
	if (i == j || (*assignments)[j] != iCl ) continue;
	
	sumOfDist += getMatrixValue(i,j);

	if (sumOfDist > minSumOfDist[iCl]) break;
      }

      // Store if this temporary centroid has smaller sumOfDist the current centroid
      if (sumOfDist < minSumOfDist[iCl]){
	minSumOfDist[iCl] = sumOfDist;
	centroids[iCl]    = i;
      }
    }
  }
  
}

double Clustering::getMatrixValue(int row, int col){


  if (row > col)  return (*distMatrix)[row][col];
  
  return (*distMatrix)[col][row];
}

/*
void Clustering::outputClusters(vector<AtomVector<> *> & _structureList) {
  if (_structureList.size() != clusterAssignments.size()) { cerr << "ERROR 1857: _structureList is not of the correct size!" << endl; exit(1857); } 
  vector<uint> counters(clusterAssignments.size(),0);
  for (uint i = 0; i < clusterAssignments.size(); i++) {
    stringstream ss;
    ss << "cluster" << i << ".pdb";
    if (counters[clusterAssignments[i]] != 0) {
      ofstream fout(ss.str().c_str(), ios::app);
      fout << "MODEL\t" << counters[clusterAssignments[i]] << endl;
      fout << *(_structureList[i]);
      fout << "ENDMDL" << endl;
      fout.close();
    }
    else {
      ofstream fout(ss.str().c_str());
      fout << "MODEL\t" << counters[clusterAssignments[i]] << endl;
      fout << *(_structureList[i]);
      fout << "ENDMDL" << endl;
      fout.close();
    }
    counters[clusterAssignments[i]]++;
  }
}

void Clustering::outputClusters(vector<vector<BaseAtom *> > & _structureList) {
  if (_structureList.size() != clusterAssignments.size()) { cerr << "ERROR 1857: _structureList is not of the correct size!" << endl; exit(1857); } 
  vector<uint> counters(clusterAssignments.size(),0);
  for (uint i = 0; i < clusterAssignments.size(); i++) {
    stringstream ss;
    ss << "cluster" << clusterAssignments[i] << ".pdb";
    if (counters[clusterAssignments[i]] != 0) {
      ofstream fout(ss.str().c_str(), ios::app);
      fout << "MODEL\t" << counters[clusterAssignments[i]] << endl;
      for (uint j = 0; j < _structureList[i].size(); j++) {
        fout << _structureList[i][j] << endl;
      }
      fout << "ENDMDL" << endl;
      fout.close();
    }
    else {
      ofstream fout(ss.str().c_str());
      fout << "MODEL\t" << counters[clusterAssignments[i]] << endl;
      for (uint j = 0; j < _structureList[i].size(); j++) {
        fout << _structureList[i][j] << endl;
      }
      fout << "ENDMDL" << endl;
      fout.close();
    }
    counters[clusterAssignments[i]]++;
  }
}
*/
double Clustering::max(double i, double j) {
   if (i > j) { return i; }
   else { return j; }
}
double Clustering::min(double i, double j) {
   if (i < j) { return i; }
   else { return j; }
}
double Clustering::avg(double i, double j, uint k, uint m) {
   double thisAvg = i*k + j*m;
   thisAvg /= (k+m);
   return thisAvg;
}
