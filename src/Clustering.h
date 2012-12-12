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

/*
  Clustering Object:
  K-medoid
  Single-Linkage
  Complete-Linkage
  Average-Linkage
 */


#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <signal.h>
#include <math.h>
#include <map>
#include <vector>
#include "RandomNumberGenerator.h"
using namespace std;

namespace MSL { 
class Clustering {
  
 public:
  Clustering();
  Clustering(vector<vector<double> > *_distMatrix);
  Clustering(vector<vector<double> > *_distMatrix,vector<string> *_labels);
  ~Clustering();

  void   Kmedoids(int _nIterations,int _nClusters);
  void   Kmeans();

  void   SingleLinkage();
  void   AverageLinkage();
  void   CenterLinkage();
  void   CompleteLinkage();

  void   printClusters();
  void   printClusters(double numSTDs);
  double getStatistics(); // return max variance across clusters

  double getMaxVariance() { return maxStats.variance; }
  double getMaxValue()    { return maxStats.max;      }
  double getMaxMean()     { return maxStats.mean;     }

  void setMarkClustersWithVar(double _variance) { markVariarianceCutoff = _variance; markFlag = true;}
  void setMarkClustersWithMax(double _max)      { markMaximumCutoff = _max;  markFlag = true;}
  void setNonMatchingValue(double _noMatch)     { nonMatchingValue = _noMatch; noMatchFlag = true;}
  void setClustersFrozen(double _maxValue, double _sizePercent);   
  void setElementsFrozen();   
  void setClustersFlag(string _flag);
 private:

  void   generateRandomClusters(vector<int> *clusterAssignments);
  void   getClusterCentroids(vector<int> *clusterAssignments);
  double getMatrixValue(int i, int j);
  double max(double i, double j);
  double min(double i, double j);
  double avg(double i, double j, uint k, uint m);

  vector<vector< double > > *distMatrix;
  vector<string>            *elementLabels;

  // Variables for Kmedoids, Kmeans algorithms
  int                        nClusters;
  vector<int>                centroids;
  vector<int>                clusterAssignments;

  RandomNumberGenerator RNG;

  struct stats {
    stats() {
      mean = variance = std = min = prob = 999;
      max  = -999; 
    }
    double mean;
    double variance;
    double std;
    double min;
    double max;
    double prob; 
    string toString(){
      stringstream ss;
      //      ss << "Mean: "<<mean<<" Variance: "<<variance<<" STD: "<<std<<" Max: "<<max<<" Prob: "<<prob;
      char tmp[100];
      sprintf(tmp,"Mean: %6.4f Variance: %6.4f STDEV: %6.4f Max: %6.4f Prob: %6.4f", mean,variance, std,max,prob);
      ss << tmp;
      return ss.str();
    }
  };
  vector<stats>              clusterStats;
  stats maxStats;

  bool markFlag;
  double markMaximumCutoff;
  double markVariarianceCutoff;

  bool noMatchFlag;
  double nonMatchingValue;


  vector<string> clusterFlag;
  vector<string> elementFlag;

};
}
#endif
