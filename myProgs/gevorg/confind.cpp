/*
----------------------------------------------------------------------------
This file is part of MaDCaT.
Copyright (C) 2012 Gevorg Grigoryan (see README)
MaDCaT: http://www.grigoryanlab.org/madcat/

MaDCaT is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

MaDCaT is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with MaDCaT.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------------
*/

#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <typeinfo>
#include <stdio.h>
#include <getopt.h>
#include <map>
#include <iomanip>

#include "AtomContainer.h"
#include "System.h"
#include "PDBFormat.h"
#include "Frame.h"
#include "OptimalRMSDCalculator.h"
#include "AtomSelection.h"
#include "PrincipleComponentAnalysis.h"
#include "CrystalLattice.h"
#include "SystemRotamerLoader.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "PDBWriter.h"

using namespace std;
using namespace MSL;

// ---- Utility Functions
void openFile (fstream& fs, string lsofn, ios_base::openmode mode);
void file2array(string _filename, vector<string>& lines);
string fileBase(string fn);
bool isNameLegal(Position& p, vector<string>& legalNames);
void proteinOnly(System& S, AtomPointerVector& A, vector<string>& legalNames, bool renumber = false);
Atom* getCA(Position& p);
bool isBackbone(string an);
bool isBackbone(Atom& a);
int findAtomInRes(Residue& r, string& an);
void assert(bool cond, string mes);
template <class T> string toString (T val);
void error(string mes);
template <class T> bool sortAsendingHelper (const T& i, const T& j);
char* trim(char* str);
void display(string mes, int w = 80, string tab = "");
string contactString(System& S, int i, int j, double d, bool perm = false);

// ---- Data Structure
// nearest-neighbor searches
class nnclass {
  public:
    nnclass(double _xlo, double _ylo, double _zlo, double _xhi, double _yhi, double _zhi, int _N = 20) {
      xlo = _xlo; ylo = _ylo; zlo = _zlo;
      xhi = _xhi; yhi = _yhi; zhi = _zhi;
      reinitBuckets(_N);
    }

    nnclass(AtomPointerVector& _atoms, int _N, bool _addAtoms = true, vector<int>* tags = NULL) {
      calculateExtent(_atoms);
      reinitBuckets(_N);
      if (_addAtoms) {
        for (int i = 0; i < _atoms.size(); i++) {
          addPoint(_atoms[i]->getCoor(), (tags == NULL) ? i : (*tags)[i]);
        }
      }
    }

    nnclass(AtomPointerVector& _atoms, double _characteristicDistance, bool _addAtoms = true, vector<int>* tags = NULL) {
      calculateExtent(_atoms);
      if (xlo == xhi) { xlo -= _characteristicDistance/2; xhi += _characteristicDistance/2; }
      if (ylo == yhi) { ylo -= _characteristicDistance/2; yhi += _characteristicDistance/2; }
      if (zlo == zhi) { zlo -= _characteristicDistance/2; zhi += _characteristicDistance/2; }
      int _N = int(ceil(max(max((xhi - xlo), (yhi - ylo)), (zhi - zlo))/_characteristicDistance));
      reinitBuckets(_N);
      if (_addAtoms) {
        for (int i = 0; i < _atoms.size(); i++) {
          addPoint(_atoms[i]->getCoor(), (tags == NULL) ? i : (*tags)[i]);
        }
      }
    }

    void calculateExtent(AtomPointerVector& _atoms) {
      if (_atoms.size() == 0) { cout << "Error in nnclass::calculateExtent() -- empty atom vector passed!\n"; exit(-1); }
      xlo = xhi = _atoms[0]->getX();
      ylo = yhi = _atoms[0]->getY();
      zlo = zhi = _atoms[0]->getZ();
      for (int i = 0; i < _atoms.size(); i++) {
        if (xlo > _atoms[i]->getX()) xlo = _atoms[i]->getX();
        if (ylo > _atoms[i]->getY()) ylo = _atoms[i]->getY();
        if (zlo > _atoms[i]->getZ()) zlo = _atoms[i]->getZ();
        if (xhi < _atoms[i]->getX()) xhi = _atoms[i]->getX();
        if (yhi < _atoms[i]->getY()) yhi = _atoms[i]->getY();
        if (zhi < _atoms[i]->getZ()) zhi = _atoms[i]->getZ();
      }
    }

    double getXLow() { return xlo; }
    double getYLow() { return ylo; }
    double getZLow() { return zlo; }
    double getXHigh() { return xhi; }
    double getYHigh() { return yhi; }
    double getZHigh() { return zhi; }
    int pointSize() { return pointList.size(); }
    CartesianPoint& getPoint(int i) { return pointList[i]; }

    void reinitBuckets(int _N) {
      N = _N;
      buckets.resize(N);
      for (int i = 0; i < N; i++) {
        buckets[i].resize(N);
        for (int j = 0; j < N; j++) {
          buckets[i][j].resize(N);
          for (int k = 0; k < N; k++) { buckets[i][j][k].resize(0); }
        }
      }
      pointList.resize(0);
    }

    void addPoint(CartesianPoint& p, int tag) {
      int i, j, k;
      pointBucket(&p, &i, &j, &k);
      if ((i < 0) || (j < 0) || (k < 0) || (i > N-1) || (j > N-1) || (k > N-1)) { cout << "Error: point " << p << " out of range for nnclass object!\n"; exit(-1); }
      buckets[i][j][k].push_back(pair<CartesianPoint*, int>(&p, tag));
      pointList.push_back(p);
    }
    
    void pointBucket(CartesianPoint* p, int* i, int* j, int* k) {
      *i = min((int) floor(N*(p->getX() - xlo)/(xhi - xlo)), N-1); // xhi should technically map to N, but we will map it to N-1
      *j = min((int)floor(N*(p->getY() - ylo)/(yhi - ylo)), N-1);
      *k = min((int)floor(N*(p->getZ() - zlo)/(zhi - zlo)), N-1);
    }
    void pointBucket(CartesianPoint p, int* i, int* j, int* k) { pointBucket(&p, i, j, k); }
    void limitIndex(int *ind) {
      if (*ind < 0) *ind = 0;
      if (*ind > N-1) *ind = N-1;
    }

    double gridSpacingX() { return (xhi - xlo)/N; }
    double gridSpacingY() { return (yhi - ylo)/N; }
    double gridSpacingZ() { return (zhi - zlo)/N; }

    void pointsWithin(CartesianPoint& c, double dmin, double dmax, vector<int>& list) {
      double d2, dmin2, dmax2;
      int ci, cj, ck;
      int imax1, jmax1, kmax1, imax2, jmax2, kmax2; // external box (no point in looking beyond it, points there are too far)
      int imin1, jmin1, kmin1, imin2, jmin2, kmin2; // internal box (no point in looking within it, points there are too close)
      pointBucket(c, &ci, &cj, &ck);
      pointBucket(c - CartesianPoint(dmax, dmax, dmax), &imax1, &jmax1, &kmax1);
      pointBucket(c + CartesianPoint(dmax, dmax, dmax), &imax2, &jmax2, &kmax2);
      pointBucket(c - CartesianPoint(dmin, dmin, dmin)/sqrt(3), &imin1, &jmin1, &kmin1);
      pointBucket(c + CartesianPoint(dmin, dmin, dmin)/sqrt(3), &imin2, &jmin2, &kmin2);
      // need to trim the internal box to make sure it is fully contained within the sphere of radius dmin from the central point
      if (imin1 != ci) imin1++;
      if (jmin1 != cj) jmin1++;
      if (kmin1 != ck) kmin1++;
      if (imin2 != ci) imin2--;
      if (jmin2 != cj) jmin2--;
      if (kmin2 != ck) kmin2--;
      limitIndex(&imin1); limitIndex(&imin2); limitIndex(&jmin1); limitIndex(&jmin2); limitIndex(&kmin1); limitIndex(&kmin2);
      limitIndex(&imax1); limitIndex(&imax2); limitIndex(&jmax1); limitIndex(&jmax2); limitIndex(&kmax1); limitIndex(&kmax2);

      // search only within the boxes where points of interest can be, in principle
      dmin2 = dmin*dmin; dmax2 = dmax*dmax;
      for (int i = imax1; i <= imax2; i++) {
        bool insi = (i >= imin1) && (i <= imin2);
        for (int j = jmax1; j <= jmax2; j++) {
          bool ins = insi && (j >= jmin1) && (j <= jmin2);
          for (int k = kmax1; k <= kmax2; k++) {
            // check all points in bucket i, j, k
            for (int ii = 0; ii < buckets[i][j][k].size(); ii++) {
              d2 = c.distance2(*(buckets[i][j][k][ii].first));
              if ((d2 >= dmin2) && (d2 <= dmax2)) {
                list.push_back(buckets[i][j][k][ii].second);
              }
            }
            // skip the range from kmin1 to kmin2 (too close)
            if (ins && (k == kmin1) && (kmin1 != kmin2)) k = kmin2 - 1;
          }
        }
      }
    }

  private:
    int N; // dimension of bucket list is N x N x N
    double xlo, ylo, zlo, xhi, yhi, zhi;
    vector<vector<vector<vector<pair<CartesianPoint*, int> > > > > buckets;
    vector<CartesianPoint> pointList;
};

class rotamer {
  public:
    rotamer() { atoms = NULL; rP = aaP = 1.0; aaN = "XXX"; rID = -1; }
    rotamer(nnclass* _atoms, double _aaProp, double _rotProb, string _name, int _rotID) { atoms = _atoms; rP = _rotProb; aaP = _aaProp; aaN = _name; rID = _rotID; }
    nnclass* grid() { return atoms; }
    double aaProp() { return aaP; }
    double rotProb() { return rP; }
    string aaName() { return aaN; }
    int rotID() { return rID; }

  private:
    nnclass* atoms;
    double rP, aaP;
    string aaN;
    int rID;
};

class contact {
  public:
    int resi;
    int resj;
    double degree;
    string info;
    contact(int _resi, int _resj, double _degree, string _info) {
      resi = _resi;
      resj = _resj;
      degree = _degree;
      info = _info;
    }
    contact() { resi = 0; resj = 0; degree = 0; info = ""; }
};
bool contactOrder (const contact& i, const contact& j) { return (j.resi == i.resi) ? (j.resj > i.resj) : (j.resi > i.resi); }

class options {
  public:
    options() {
      dcut = 25.0;
      clashDist = 2.0; contDist = 3.0; rotLibFile = ""; rotOutFile = ""; rotLevel = "";
      verbose = false; renumPDB = false;
      aaProp["ALA"] = 7.73; aaProp["CYS"] = 1.84; aaProp["ASP"] = 5.82; aaProp["GLU"] = 6.61; aaProp["PHE"] = 4.05;
      aaProp["GLY"] = 7.11; aaProp["HIS"] = 2.35; aaProp["HSD"] = 2.35; aaProp["ILE"] = 5.66; aaProp["LYS"] = 6.27;
      aaProp["LEU"] = 8.83; aaProp["MET"] = 2.08; aaProp["ASN"] = 4.50; aaProp["PRO"] = 4.52; aaProp["GLN"] = 3.94;
      aaProp["ARG"] = 5.03; aaProp["SER"] = 6.13; aaProp["THR"] = 5.53; aaProp["VAL"] = 6.91; aaProp["TRP"] = 1.51; aaProp["TYR"] = 3.54;
    }
    vector<string> pdbfs, omapfs, opdbfs;
    bool verbose, renumPDB;
    string selection, rotLibFile, rotOutFile, rotLevel;
    double dcut, clashDist, contDist;
    map<string, double> aaProp; // amino-acid propensities (in percent)
    map<string, vector<double> > rotProbs; // rotamer probabilities, mapped by aa name and then indexed by rotamer index
};

string option(string opt, string mes, int w, int p1, int p2) {
  // first print the name of the option
  string text(p1, ' ');
  text += opt;
  if (p2 > text.size()) text += string(p2 - text.size(), ' ');

  // next print the description text
  int i = 0, k, L = text.size(), n;
  while (i < mes.size()) {
    k = mes.find_first_of(" ", i);
    if (k == string::npos) k = mes.size();
    n = k - i;
    if ((L + n >= w) && (L > 0)) { text += "\n" + string(p2, ' '); L = p2; }
    text += mes.substr(i, n) + " ";
    L += n + 1;
    i = k+1;
  }
  return text;
}

void usage() {
  int w = 80, p1 = 3, p2 = p1+8;
  cout << endl << option("", "Creates distance map(s) for use with MaDCaT (as either database or query maps) from input PDB file(s). Options:", w, 0, 0) << endl;
  cout << option("--p", "input PDB file.", w, p1, p2) << endl;
  cout << option("--pL", "a file with a list of PDB files. Either --p or --pL must be given.", w, p1, p2) << endl;
  cout << option("--o", "output file name for writing contacts. If not given, will write to standard output.", w, p1, p2) << endl;
  cout << option("--oL", "a file with a list of contact file names (one per input PDB structure).", w, p1, p2) << endl;
  cout << option("--opdb", "optional: output post-processed PDB file (useful for keeping track of how all PDB weirdnesses got parsed).", w, p1, p2) << endl;
  cout << option("--opdbL", "optional: a file with a list of file names for post-processed PDBs, one per input structure.", w, p1, p2) << endl;
  cout << option("--rLib", "a path to an MSL-formatter rotamer library. Needed only if --type is 2 (contact probability map).", w, p1, p2) << endl;
  cout << option("--rout", "name of a file into which to place PDB-formated coordinates of rotamers that ended up surviving at each considered position.", w, p1, p2) << endl;
  cout << option("--sel", "optional: selection string to apply before doing anything (only the selected part of structure will be considered). Will select a residue if its CA atom is included in the given selection.", w, p1, p2) << endl;
  cout << option("--verb", "optional: generate lots of detailed output (i.e., for --type 2, it will print the details of which rotamer pairs are in contact).", w, p1, p2) << endl << endl;
  cout << option("--ren", "optional: if flag specified, will renumber the structure before output. Useful for keeping track of residues in the output list of contacts if the input PDB file is strangely numbered.", w, p1, p2) << endl;
}

void parseCommandLine(int argc, char** argv, options& iopts) {
  map<string, bool> spec;

  while (1) {
    int oind = 0;
    static struct option opts[] = {
      {"p", 1, 0, 1},
      {"o", 1, 0, 2},
      {"opdb", 1, 0, 7},
      {"pL", 1, 0, 10},
      {"oL", 1, 0, 11},
      {"opdbL", 1, 0, 12},
      {"sel", 1, 0, 18},
      {"rLib", 1, 0, 19},
      {"verb", 0, 0, 20},
      {"rout", 1, 0, 21},
      {"ren", 0, 0, 22},
      {0, 0, 0, 0}
    };

    int c = getopt_long (argc, argv, "", opts, &oind);
    if (c == -1) break;

    switch (c) {
      case 1:
        if (iopts.pdbfs.size() > 0) { usage(); exit(-1); }
        iopts.pdbfs.push_back(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 2:
        if (iopts.omapfs.size() > 0) { usage(); exit(-1); }
        iopts.omapfs.push_back(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 7:
        if (iopts.opdbfs.size() > 0) { usage(); exit(-1); }
        iopts.opdbfs.push_back(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 10:
        if (iopts.pdbfs.size() > 0) { usage(); exit(-1); }
        file2array(optarg, iopts.pdbfs);
        spec[string(opts[oind].name)] = true;
        break;

      case 11:
        if (iopts.omapfs.size() > 0) { usage(); exit(-1); }
        file2array(optarg, iopts.omapfs);
        spec[string(opts[oind].name)] = true;
        break;

      case 12:
        if (iopts.opdbfs.size() > 0) { usage(); exit(-1); }
        file2array(optarg, iopts.opdbfs);
        spec[string(opts[oind].name)] = true;
        break;

      case 18:
        iopts.selection = string(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 19:
        iopts.rotLibFile = string(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 20:
        iopts.verbose = true;
        spec[string(opts[oind].name)] = true;
        break;

      case 21:
        iopts.rotOutFile = string(optarg);
        spec[string(opts[oind].name)] = true;
        break;

      case 22:
        iopts.renumPDB = true;
        spec[string(opts[oind].name)] = true;
        break;

      case '?':
        break;

      default:
        printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }

  // make sure all required options have been specified
  if (!(((spec.find(string("p")) != spec.end()) || (spec.find(string("pL")) != spec.end())) && (spec.find(string("rLib")) != spec.end()))) {
    cout << "Not all required options specified!\n"; usage(); exit(-1);
  }

  // error checking
  // make sure lists are of the proper size
  if (iopts.omapfs.size() > 1) assert(iopts.omapfs.size() == iopts.pdbfs.size(), "the number of input PDB files and output map files does not agree");
  if (iopts.opdbfs.size() > 1) assert(iopts.opdbfs.size() == iopts.pdbfs.size(), "the number of input PDB files and output PDB files does not agree");
  if (iopts.pdbfs.size() > 1) {
    if (iopts.omapfs.size() == 1) {
      iopts.omapfs.resize(iopts.pdbfs.size());
      string base = fileBase(iopts.omapfs[0]);
      for (int i = 0; i < iopts.pdbfs.size(); i++) { iopts.omapfs[i] = base + ".f" + toString(i+1) + ".cont"; }
    }
    if (iopts.opdbfs.size() == 1) {
      iopts.opdbfs.resize(iopts.pdbfs.size());
      string base = fileBase(iopts.opdbfs[0]);
      for (int i = 0; i < iopts.pdbfs.size(); i++) { iopts.opdbfs[i] = base + ".f" + toString(i+1) + ".pdb"; }
    }
  }
}

// ---- Number Crunchers
CartesianPoint getCentroid(vector<Atom*> atoms) {
  CartesianPoint centroid(0.0, 0.0, 0.0);
  if (atoms.size() == 0.0) error("cannot calculate the centroid of an empty atom set!");

  for (int i = 0; i < atoms.size(); i++) {
    centroid += atoms[i]->getCoor();
  }
  centroid /= (double) atoms.size();

  return centroid;
}

void readRotamerProbabilities(options& _opts) {
  bool verbose = false;
  string currAA = "XXX";
  ifstream fs;
  fs.open(_opts.rotLibFile.c_str());
  if (fs.fail()) error("could not open rotamer library file '" + _opts.rotLibFile + "' to read probabilities");

  while (!fs.fail()) {
    string line;
    getline(fs,line);
    vector<string> tokens = MslTools::tokenize(line);
    if (tokens.size() == 0) continue;
    if (tokens[0] == "RESI") {
      if (tokens.size() < 1) error("could not parse residue name from rotamer library " + _opts.rotLibFile + ", line '" + line + "'");
      currAA = tokens[1];
    } else if (tokens[0] == "PROBS") {
      if (tokens.size() < 1) error("could not parse rotamer probabilities from rotamer library " + _opts.rotLibFile + ", line '" + line + "'");
      vector<double> empty;
      _opts.rotProbs[currAA] = empty;
      for (int i = 1; i < tokens.size(); i++) {
        _opts.rotProbs[currAA].push_back(MslTools::toDouble(tokens[i]));
        if (verbose) cout << currAA << ", rotamer " << i << "/" << tokens.size()-1 << ", probability = " << tokens[i] << endl;
      }
      if (verbose) cout << "... read " << _opts.rotProbs[currAA].size() << " rotamer probabilities for " << currAA << endl;
    }
  }
  fs.close();
}

void filterRotamers(System& _sys, options& _opt, vector<int>& resIndex, vector<vector<rotamer> >& rotamers, vector<set<int> >& permanentContacts, vector<double>& fractionPruned) {
  fstream rof;
  PDBWriter pdbw;
  stringstream ss;
  int param = 19;
  rotamers.resize(resIndex.size());
  permanentContacts.resize(resIndex.size());
  fractionPruned.resize(resIndex.size(), 0);

  // all amino acid names to add rotamers for (all except Gly and Pro)
  vector<string> aaNames = MslTools::tokenize("ARG ASN ASP CYS GLN GLU HIS ILE LEU LYS MET PHE SER THR TRP TYR VAL ALA", " ");

  /* Instead of having to both with topologies, we will manually create all identities by putting in dummy side-chain atoms,
   * and simply copy the backbone. Thus we have to know which atoms are present in all amino acids. */
  map<string, vector<string> > aaAtomNames;
  if (param == 19) {
    aaAtomNames["ALA"] = MslTools::tokenize("N H CA CB C O", " ");
    aaAtomNames["ARG"] = MslTools::tokenize("N H CA CB CG CD NE HE CZ NH1 HH11 HH12 NH2 HH21 HH22 C O", " ");
    aaAtomNames["ASN"] = MslTools::tokenize("N H CA CB CG OD1 ND2 HD21 HD22 C O", " ");
    aaAtomNames["ASP"] = MslTools::tokenize("N H CA CB CG OD1 OD2 C O", " ");
    aaAtomNames["CYS"] = MslTools::tokenize("N H CA CB SG C O", " ");
    aaAtomNames["GLN"] = MslTools::tokenize("N H CA CB CG CD OE1 NE2 HE21 HE22 C O", " ");
    aaAtomNames["GLU"] = MslTools::tokenize("N H CA CB CG CD OE1 OE2 C O", " ");
    aaAtomNames["GLY"] = MslTools::tokenize("N H CA C O", " ");
    aaAtomNames["HIS"] = MslTools::tokenize("N H CA CB CG ND1 HD1 CD2 NE2 CE1 C O", " ");
    aaAtomNames["HSC"] = MslTools::tokenize("N H CA CB CG CD2 ND1 HD1 CE1 NE2 HE2 C O", " ");
    aaAtomNames["HSD"] = MslTools::tokenize("N H CA CB CG ND1 CE1 CD2 NE2 HE2 C O", " ");
    aaAtomNames["ILE"] = MslTools::tokenize("N H CA CB CG2 CG1 CD C O", " ");
    aaAtomNames["LEU"] = MslTools::tokenize("N H CA CB CG CD1 CD2 C O", " ");
    aaAtomNames["LYS"] = MslTools::tokenize("N H CA CB CG CD CE NZ HZ1 HZ2 HZ3 C O", " ");
    aaAtomNames["MET"] = MslTools::tokenize("N H CA CB CG SD CE C O", " ");
    aaAtomNames["MSE"] = MslTools::tokenize("N H CA CB CG SE CE C O", " ");
    aaAtomNames["PHE"] = MslTools::tokenize("N H CA CB CG CD1 CD2 CE1 CE2 CZ C O", " ");
    aaAtomNames["PRO"] = MslTools::tokenize("N CD CA CB CG C O", " ");
    aaAtomNames["SER"] = MslTools::tokenize("N H CA CB OG HG C O", " ");
    aaAtomNames["THR"] = MslTools::tokenize("N H CA CB OG1 HG1 CG2 C O", " ");
    aaAtomNames["TRP"] = MslTools::tokenize("N H CA CB CG CD2 CE2 CE3 CD1 NE1 HE1 CZ2 CZ3 CH2 C O", " ");
    aaAtomNames["TYR"] = MslTools::tokenize("N H CA CB CG CD1 CE1 CD2 CE2 CZ OH HH C O", " ");
    aaAtomNames["VAL"] = MslTools::tokenize("N H CA CB CG1 CG2 C O", " ");
  } else {
    aaAtomNames["ALA"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 HB3 C O", " ");
    aaAtomNames["ARG"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CG HG1 HG2 CD HD1 HD2 NE HE CZ NH1 HH11 HH12 NH2 HH21 HH22 C O", " ");
    aaAtomNames["ASN"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CG OD1 ND2 HD21 HD22 C O", " ");
    aaAtomNames["ASP"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CG OD1 OD2 C O", " ");
    aaAtomNames["CYS"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 SG HG1 C O", " ");
    aaAtomNames["GLN"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CG HG1 HG2 CD OE1 NE2 HE21 HE22 C O", " ");
    aaAtomNames["GLU"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CG HG1 HG2 CD OE1 OE2 C O", " ");
    aaAtomNames["GLY"] = MslTools::tokenize("N HN CA HA1 HA2 C O", " ");
    aaAtomNames["HIS"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 ND1 HD1 CG CE1 HE1 NE2 CD2 HD2 C O", " "); // arbitrarily make "HIS" HSD (as in CHARMM 19)
    aaAtomNames["HSD"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 ND1 HD1 CG CE1 HE1 NE2 CD2 HD2 C O", " ");
    aaAtomNames["HSE"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 ND1 CG CE1 HE1 NE2 HE2 CD2 HD2 C O", " ");
    aaAtomNames["HSP"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CD2 HD2 CG NE2 HE2 ND1 HD1 CE1 HE1 C O", " ");
    aaAtomNames["ILE"] = MslTools::tokenize("N HN CA HA CB HB CG2 HG21 HG22 HG23 CG1 HG11 HG12 CD HD1 HD2 HD3 C O", " ");
    aaAtomNames["LEU"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CG HG CD1 HD11 HD12 HD13 CD2 HD21 HD22 HD23 C O", " ");
    aaAtomNames["LYS"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CG HG1 HG2 CD HD1 HD2 CE HE1 HE2 NZ HZ1 HZ2 HZ3 C O", " ");
    aaAtomNames["MET"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CG HG1 HG2 SD CE HE1 HE2 HE3 C O", " ");
    aaAtomNames["MSE"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CG HG1 HG2 SE CE HE1 HE2 HE3 C O", " ");
    aaAtomNames["PHE"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CG CD1 HD1 CE1 HE1 CZ HZ CD2 HD2 CE2 HE2 C O", " ");
    aaAtomNames["PRO"] = MslTools::tokenize("N CD HD1 HD2 CA HA CB HB1 HB2 CG HG1 HG2 C O", " ");
    aaAtomNames["SER"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 OG HG1 C O", " ");
    aaAtomNames["THR"] = MslTools::tokenize("N HN CA HA CB HB OG1 HG1 CG2 HG21 HG22 HG23 C O", " ");
    aaAtomNames["TRP"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CG CD1 HD1 NE1 HE1 CE2 CD2 CE3 HE3 CZ3 HZ3 CZ2 HZ2 CH2 HH2 C O", " ");
    aaAtomNames["TYR"] = MslTools::tokenize("N HN CA HA CB HB1 HB2 CG CD1 HD1 CE1 HE1 CZ OH HH CD2 HD2 CE2 HE2 C O", " ");
    aaAtomNames["VAL"] = MslTools::tokenize("N HN CA HA CB HB CG1 HG11 HG12 HG13 CG2 HG21 HG22 HG23 C O", " ");
  }

  // read rotamer library
  SystemRotamerLoader sysRot;
  sysRot.setSystem(_sys);
  sysRot.defineRotamerSamplingLevels();
  if(!sysRot.readRotamerLibraryFile(_opt.rotLibFile)) { error("unable to load rotamer library from " + _opt.rotLibFile); }
  readRotamerProbabilities(_opt);

  // get backbone atoms, mark which residue each belongs to, and make up a grid for searchig them
  vector<Position*> & positions = _sys.getPositions();
  map<string, int> posIdMap;
  for (int i = 0; i < resIndex.size(); i++) posIdMap[positions[resIndex[i]]->getPositionId()] = i;
  AtomPointerVector backbone;
  vector<int> tags; // these tags will indicate with position index each backbone atom in the NNclass is from
  for (int i = 0; i < _sys.atomSize(); i++) {
    if (isBackbone(_sys[i]) && _sys[i].getName().compare(0, 1, "H")) {
      backbone.push_back(&(_sys[i]));
      string resId = _sys[i].getParentResidue()->getPositionId();
      if (posIdMap.find(resId) != posIdMap.end()) {
        tags.push_back(posIdMap[resId]);
      } else {
        tags.push_back(-1);
      }
    }
  }
  nnclass bbNN(backbone, _opt.clashDist, true, &tags);

  // load rotamers
  if (_opt.verbose) {
    printf("rotamer filtering...\n");
  }
  if (!_opt.rotOutFile.empty()) {
    openFile(rof, _opt.rotOutFile, fstream::out);
  }
  for (int i = 0; i < resIndex.size(); i++) {
    Position &posi = _sys.getPosition(resIndex[i]);
    if (_opt.verbose) {
      printf("position %s\n", posi.getPositionId().c_str());
    }
    // instead of erasing the native residue, rename it so that it is not used in counting rotamer probabilities
    // and so that it does not interfere with building novel identities
    // note: the native identity always stays as the first one
    for (int ii = 0; ii < posi.getNumberOfIdentities(); ii++) { // sometimes there is more than one native identity (in cases of multiple possible identities in crystal structures, for example)
      posi.getIdentity(ii).setResidueName("_NAT_" + posi.getIdentity(ii).getResidueName());
    }
    Residue& nativeRes = posi.getCurrentIdentity();

    int numRemRotsInPosition = 0; int totNumRotsInPosition = 0;
    for (int j = 0; j < aaNames.size(); j++) {
      if (_opt.aaProp.find(aaNames[j]) == _opt.aaProp.end()) error("no propensity defined for amino acid " + aaNames[j]);
      double aaProp = _opt.aaProp[aaNames[j]];
      if (_opt.verbose) {
        printf("%s %.3f: ", aaNames[j].c_str(), aaProp);
      }
      // build in the identity, copying the backbone from native
      AtomPointerVector ratoms;
      for (int k = 0; k < aaAtomNames[aaNames[j]].size(); k++) {
        string an = aaAtomNames[aaNames[j]][k];
        string atomId = posi.getPositionId() + "," + aaNames[j] + "," + an;
        int nidx = findAtomInRes(nativeRes, an);
        if (isBackbone(an) && (nidx >= 0) && nativeRes[nidx].hasCoor()) ratoms.push_back(new Atom(atomId, nativeRes[nidx].getCoor()));
        else { ratoms.push_back(new Atom(atomId)); ratoms.back()->setHasCoordinates(false); }
      }
      Residue tmpres(ratoms, aaNames[j], posi.getResidueNumber(), posi.getResidueIcode());
      posi.addIdentity(tmpres);
      posi.setActiveIdentity(aaNames[j]);
      if (!sysRot.loadRotamers(&posi, aaNames[j], sysRot.getRotamerLibrary()->size("", aaNames[j]), "", false)) {
        error("Cannot load rotamers for of " + aaNames[j] + " in " + posi.getPositionId());
      }
      // check to see if the rotamers were successfully built (i.e., make sure all side-chain atoms were correctly placed)
      bool succ = true;
      Residue& addedIdentity = posi.getCurrentIdentity();
      for (int k = 0; k < addedIdentity.atomSize(); k++) {
        if (isBackbone(addedIdentity.getAtom(k))) continue;
        if (!addedIdentity.getAtom(k).hasCoor()) {
          cout << "Warning: could not build rotamers for " << aaNames[j] << " at position " << posi.getPositionId() << ", skipping (this will affect the accuracy of the contact probability map)...\n";
          for (int ii = 0; ii < addedIdentity.atomSize(); ii++) {
            cout << addedIdentity[ii] << endl;
          }
          posi.removeIdentity(aaNames[j]);
          succ = false;
          break;
        }
      }
      if (!succ) continue;

      // if rotamers were successfully built, see which need to be pruned (clash with the backbone)
      int numRemRots = 0;
      for (int r = 0; r < addedIdentity.getNumberOfAltConformations(); r++) {
        bool prune = false;
        addedIdentity.setActiveConformation(r);
        for (int k = 0; k < addedIdentity.atomSize(); k++) {
          if (isBackbone(addedIdentity[k])) continue;
          vector<int> closeOnes;
          bbNN.pointsWithin(addedIdentity[k].getCoor(), 0.0, _opt.clashDist, closeOnes);
          for (int ci = 0; ci < closeOnes.size(); ci++) {
            // backbone atoms of the same residue do not count as clashing (the rotamer library should not allow true clashes with own backbone)
            if (closeOnes[ci] != i) {
              prune = true;
              // clashes with ALA have a special meaning (permanent "unavoidable" contacts; need to find all of them, though unlikely to have more than one)
              if (!aaNames[j].compare("ALA")) permanentContacts[i].insert(closeOnes[ci]);
              else break;
            }
          }
          if (prune) break;
        }
        // if rotamer not pruned, prepare a datastructure corresponding to it for use in NN searching
        if (!prune) {
          if (!_opt.rotOutFile.empty()) {
            ss.str("");
            pdbw.open(ss);
            pdbw.write(addedIdentity.getAtomPointers());
            string pdbs = ss.str();
            rof << "REM " << addedIdentity.getIdentityId() << ", rotamer " << r+1 << endl << pdbs << "END\n";
            pdbw.close();
          }
          AtomPointerVector heavySC;
          for (int ai = 0; ai < addedIdentity.atomSize(); ai++) {
            if ((!isBackbone(addedIdentity[ai])) && addedIdentity[ai].getName().compare(0, 1, "H")) heavySC.push_back(&(addedIdentity[ai]));
          }
          double rotProb = 1;
          if (_opt.rotProbs.find(aaNames[j]) == _opt.rotProbs.end()) rotProb = 1;
          else {
            if (_opt.rotProbs[aaNames[j]].size() <= r) error("no probability defined for rotamer " + MslTools::doubleToString(r) + " of amino acid " + aaNames[j]);
            rotProb = _opt.rotProbs[aaNames[j]][r];
          }
          rotamers[i].push_back(rotamer(new nnclass(heavySC, _opt.contDist), aaProp, rotProb, aaNames[j], r));
          if (_opt.verbose) {
            printf("%d %.3f; ", r, rotProb);
          }
          numRemRots++;
          numRemRotsInPosition++;
        }
      }
      totNumRotsInPosition += addedIdentity.getNumberOfAltConformations();
      if (_opt.verbose) cout << numRemRots << "/" << addedIdentity.getNumberOfAltConformations() << " remaining at position " << addedIdentity.toString() << endl;
    }
    fractionPruned[i] = (totNumRotsInPosition - numRemRotsInPosition)*1.0/totNumRotsInPosition;
    posi.setActiveIdentity(nativeRes.getResidueName()); // return the native as the active identity
    if (nativeRes.getResidueName().compare(0, 5, "_NAT_") == 0) { // rename native identity back to its original
      nativeRes.setResidueName(nativeRes.getResidueName().substr(5));
    }
  }
  if (_opt.verbose) {
    printf("end of rotamer filtering...\n");
  }
  if (!_opt.rotOutFile.empty()) {
    rof.close();
  }
}

double contactProbability (vector<rotamer>& posi, vector<rotamer>& posj, options& iopts, string* contactInfo = NULL) {
  double n = 0; double c = 0;
  double contDist = iopts.contDist;
  if ((posi.size() == 0) || (posj.size() == 0)) return 0.0;

  for (int i = 0; i < posi.size(); i++) {
    double ixlo = posi[i].grid()->getXLow();
    double iylo = posi[i].grid()->getYLow();
    double izlo = posi[i].grid()->getZLow();
    double ixhi = posi[i].grid()->getXHigh();
    double iyhi = posi[i].grid()->getYHigh();
    double izhi = posi[i].grid()->getZHigh();
    double p1 = posi[i].aaProp() * posi[i].rotProb();

    for (int j = 0; j < posj.size(); j++) {
      double p1p2 = p1 * posj[j].aaProp() * posj[j].rotProb();
      n += p1p2;
      double jxlo = posj[j].grid()->getXLow();
      double jylo = posj[j].grid()->getYLow();
      double jzlo = posj[j].grid()->getZLow();
      double jxhi = posj[j].grid()->getXHigh();
      double jyhi = posj[j].grid()->getYHigh();
      double jzhi = posj[j].grid()->getZHigh();

      // skip right away if boxes are farther appart than interaction distance
      if ((jxlo > ixhi + contDist) || (ixlo > jxhi + contDist) || (jylo > iyhi + contDist) || (iylo > jyhi + contDist) || (jzlo > izhi + contDist) || (izlo > jzhi + contDist)) {
        continue;
      }

      // otherwise, investigage point-by-point
      bool cont = false;
      for (int ai = 0; ai < posj[j].grid()->pointSize(); ai++) {
        vector<int> closeOnes;
        posi[i].grid()->pointsWithin(posj[j].grid()->getPoint(ai), 0.0, iopts.contDist, closeOnes);
        if (closeOnes.size() > 0) {
          cont = true;
          break;
        }
      }

      // count contacts
      if (cont) {
        if (contactInfo != NULL) {
          *contactInfo += posi[i].aaName() + " " + MslTools::intToString(posi[i].rotID()) + " -- " + posj[j].aaName() + " " + MslTools::intToString(posj[j].rotID()) + " : " + MslTools::doubleToString(p1p2) + "\n";
        }
        c += p1p2;
      }
    }
  }
  return c*1.0/n;
}

void computeContactDegrees(System& S, options& iopts, vector<int>& resIndex, vector<vector<rotamer> >& rotamers, vector<contact>& conts) {
  double d; int ii, jj;
  // extra CA atoms only
  AtomPointerVector R;
  for (int i = 0; i < resIndex.size(); i++) {
    Position &p = S.getPosition(resIndex[i]);
    R.push_back(getCA(p));
  }

  // aling CA atoms to have the principal component along X
  Frame O, F;
  CoordAxes xyz(CartesianPoint(1, 0, 0), CartesianPoint(0, 1, 0), CartesianPoint(0, 0, 1));
  O.computeFrameFromAxes(xyz);
  F.computeFrameFromPCA(R);
  AtomContainer Rt(R);
  Frame::transformAtoms(Rt.getAtomPointers(), F, O);

  // sort atoms in ascending order of the coordinate along the 1st principal component
  vector<triple<double, Atom*, int> > atoms;
  for (int i = 0; i < Rt.size(); i++) {
    atoms.push_back(triple<double, Atom*, int> (Rt[i].getZ(), &(Rt[i]), i));
  }
  sort(atoms.begin(), atoms.end(), sortAsendingHelper<triple<double, Atom*, int> >);

  int hi = 0;
  // get all contacts
  conts.resize(0);
  for (int i = 0; i < atoms.size(); i++) {
    // get window around the i-th atom of all atoms within dcut along 1st principal axis
    // (only the upper bound of this window is needed since we don't want to repeat contacts)
    while ((hi < atoms.size()-1) && (atoms[hi].second->getZ() < atoms[i].second->getZ() + iopts.dcut)) hi = hi + 1;
    // iterate over that window to find all atoms that are within dcut of the i-th
    for (int j = i+1; j <= hi; j++) {
      if (j == i) continue;
      ii = min(atoms[i].third, atoms[j].third);
      jj = max(atoms[i].third, atoms[j].third);
      string info;
      d = contactProbability(rotamers[ii], rotamers[jj], iopts, iopts.verbose ? &info : NULL);
      if (d > 0) {
        conts.push_back(contact(ii, jj, d, info));
      }
    }
  }

  // sort contacts in a nice manner (by residue indices)
  sort(conts.begin(), conts.end(), contactOrder);
}

// ---- Main Program
int main(int argc, char *argv[]) {
  int ii, jj; double d;

  // process input arguments
  options iopts;
  parseCommandLine(argc, argv, iopts);
  
  // legal residue names that are considered "protein" here
  vector<string> legalNames;
  legalNames.push_back("ALA"); legalNames.push_back("CYS"); legalNames.push_back("ASP"); legalNames.push_back("GLU"); legalNames.push_back("PHE"); legalNames.push_back("GLY");
  legalNames.push_back("HIS"); legalNames.push_back("ILE"); legalNames.push_back("LYS"); legalNames.push_back("LEU"); legalNames.push_back("MET"); legalNames.push_back("ASN");
  legalNames.push_back("PRO"); legalNames.push_back("GLN"); legalNames.push_back("ARG"); legalNames.push_back("SER"); legalNames.push_back("THR"); legalNames.push_back("VAL");
  legalNames.push_back("TRP"); legalNames.push_back("TYR"); legalNames.push_back("HSD"); legalNames.push_back("HSE"); legalNames.push_back("HSC"); legalNames.push_back("HSP");
  legalNames.push_back("MSE");
  legalNames.push_back("CSO"); legalNames.push_back("HIP"); legalNames.push_back("SEC"); legalNames.push_back("SEP"); legalNames.push_back("TPO"); legalNames.push_back("PTR");

  for (int si = 0; si < iopts.pdbfs.size(); si++) {
    AtomContainer C;                                     // original input PDB structure
    System S;                                            // just the region of the original structure corresponding to the map
    assert (C.readPdb(iopts.pdbfs[si]), (string) "Could not read PDB file " + iopts.pdbfs[si]);
    if (!iopts.selection.empty()) {
      AtomSelection selector(C.getAtomPointers());
      proteinOnly(S, selector.select(iopts.selection, true), legalNames, iopts.renumPDB);
    } else {
      proteinOnly(S, C.getAtomPointers(), legalNames, iopts.renumPDB);
    }

    // open output file and write header
    fstream of;
    streambuf * buf;
    if (!iopts.omapfs.empty()) {
      openFile(of, iopts.omapfs[si], fstream::out);
      buf = of.rdbuf();
    } else {
      cout << iopts.pdbfs[si] << endl;
      buf = cout.rdbuf();
    }
    ostream out(buf);
    
    // --- load and rotamers, create a data structure for each rotamer at each position for fast neighbor searching
    vector<vector<rotamer> > rotamers;
    vector<set<int> > permanentContacts;
    vector<double> fractionPruned;
    vector<int> resIndex; for (int i = 0; i < S.positionSize(); i++) resIndex.push_back(i); // can optionally limit to a subset of residues within the system
    filterRotamers(S, iopts, resIndex, rotamers, permanentContacts, fractionPruned);

    // --- compute contact degrees
    vector<contact> conts; computeContactDegrees(S, iopts, resIndex, rotamers, conts);

    // --- write contact degree information
    for (int i = 0; i < conts.size(); i++) {
      ii = conts[i].resi;
      jj = conts[i].resj;
      d = conts[i].degree;
      out << contactString(S, resIndex[ii], resIndex[jj], d) << endl;
      if (iopts.verbose) { printf("%s", conts[i].info.c_str()); }
    }

    // -- write out permanent contacts
    for (int i = 0; i < permanentContacts.size(); i++) {
      for (set<int>::iterator it = permanentContacts[i].begin(); it != permanentContacts[i].end(); ++it) {
        out << contactString(S, resIndex[i], resIndex[*it], -1, true) << endl;
      }
    }

    // -- write out the crowdedness parameter
    for (int i = 0; i < fractionPruned.size(); i++) {
      out << "crowd  \t" << S.getPosition(resIndex[i]).getPositionId() << "\t" << std::setprecision(6) << std::fixed << fractionPruned[i] << endl;
    }

    // --- free near-neighbor structures if was dealing with contact probability maps
    for (int i = 0; i < rotamers.size(); i++) {
      for (int j = 0; j < rotamers[i].size(); j++) delete(rotamers[i][j].grid());
    }

    // write sequence information
    out << "SEQUENCE:";
    for (int i = 0; i < resIndex.size(); i++) {
      Position &p = S.getPosition(resIndex[i]);
      out << " " << p.getResidueName();
    }
    out << endl;

    // close output file
    if (!iopts.omapfs.empty()) of.close();

    // write out the parsed region of interest
    if (!iopts.opdbfs.empty()) S.writePdb(iopts.opdbfs[si]);

  }

}


// ---- utility functions definitions
void openFile (fstream& fs, string lsofn, ios_base::openmode mode) {
  fs.open(lsofn.c_str(), mode);
  if (fs.fail()) error("could not open file " + lsofn);
}

void file2array(string _filename, vector<string>& lines) {
  FILE* ifp;
  int maxline = 1000;
  char *line, *tline;
  line = (char*) malloc(sizeof(char)*maxline);

  ifp = fopen(_filename.c_str(), "r");
  if (ifp == NULL) { error("unable to open file " + _filename); }

  while (fgets(line, maxline, ifp) != NULL) {
    if (line[strlen(line)-1] != '\n') { error("lines in file " + _filename + " are over " + toString(maxline) + "  long - increase max line limit and recompile."); }
    tline = trim(line);
    if (strlen(tline) > 0) { lines.push_back(line); }
    if (feof(ifp)) { break; }
  }
  if (ferror(ifp)) { error("an error occurred while reading file " + _filename); }
  free(line);
}

string fileBase(string fn) {
  if (fn.find_last_of(".") == string::npos) return fn;
  else return fn.substr(0, fn.find_last_of("."));
}

bool isNameLegal(Position& p, vector<string>& legalNames) {
  const char* name = p.getResidueName().c_str();
  for (int i = 0; i < legalNames.size(); i++) {
    if (strcasecmp(name, legalNames[i].c_str()) == 0) return true;
  }
  return false;
}

void proteinOnly(System& S, AtomPointerVector& A, vector<string>& legalNames, bool renumber) {
  AtomPointerVector v;
  for (int i = 0; i < A.size(); i++) {
    bool isNameLegal = false;
    for (int j = 0; j < legalNames.size(); j++) {
      if (A[i]->getResidueName().compare(legalNames[j]) == 0) { isNameLegal = true; break; }
    }
    if (isNameLegal) v.push_back(A[i]);
  }
  S.addAtoms(v);
  if (renumber) {
    for (int i = 0; i < S.chainSize(); i++) {
      S(i).renumberChain(1);
    }
  }
}

Atom* getCA(Position& p) {
  if (!p.atomExists("CA")) {
    Atom ca("A,1,XXX,CA");
    // attempt to find the centroid of the backbone
    vector<string> bban;
    bban.push_back("CA");
    bban.push_back("C");
    bban.push_back("N");
    bban.push_back("O");
    vector<Atom*> bba;
    for (int i = 0; i < bban.size(); i++) {
      if (p.atomExists(bban[i])) bba.push_back(&(p.getAtom(bban[i])));
    }
    if (bba.size() > 0) {
      ca.setCoor(getCentroid(bba));
    } else {
      // otherwise, get the overall centroid
      ca.setCoor(p(0).getCentroid());
    }
    p(0).addAtom(ca);
  }
  return &(p.getAtom("CA"));
}

bool isBackbone(string an) {
  return (!an.compare("N") || !an.compare("C") || !an.compare("CA") || !an.compare("H") || !an.compare("O") || !an.compare("NT") || !an.compare("HA") || !an.compare("HN"));
}

bool isBackbone(Atom& a) {
  return (isBackbone(a.getName()));
}

int findAtomInRes(Residue& r, string& an) {
  for (int i = 0; i < r.atomSize(); i++) {
    if (an.compare(r[i].getName()) == 0) return i;
  }
  return -1;
}

void assert(bool cond, string mes) {
  if (!cond) error(mes);
}

template <class T>
string toString (T val) {
  return static_cast<ostringstream*>( &(ostringstream() << val) )->str();
}

void error(string mes) {
  display("Error: " + mes + "\n");
  exit(-1);
}

template <class T>
bool sortAsendingHelper (const T& i, const T& j) { return (j.first > i.first); }

char* trim(char* str) {
  char* nptr; int i;
  for (i = 0; i < strlen(str); i++) {
    if ((str[i] != '\n') && (str[i] != '\t') && (str[i] != ' ')) { break; }
  }
  nptr = str + i;
  if (i == strlen(str)) { return nptr; }

  for (i = strlen(str)-1; i >= 0; i--) {
    if ((str[i] != '\n') && (str[i] != '\t') && (str[i] != ' ')) {
      str[i+1] = '\0';
      break;
    }
  }
  return nptr;
}

void display(string mes, int w, string tab) {
  int i = 0, k, L = 0, n;
  while (i < mes.size()) {
    k = mes.find_first_of(" ", i);
    if (k == string::npos) k = mes.size();
    n = k - i;
    if ((L + n > w) && (L > 0)) { cout << endl << tab; L = tab.size(); }
    cout << mes.substr(i, n) << " ";
    L += n + 1;
    i = k+1;
  }
}

string contactString(System& S, int i, int j, double d, bool perm) {
  stringstream ss;
  ss << (perm ? "percont" : "contact") << "\t" << S.getPosition(i).getPositionId() << "\t" << S.getPosition(j).getPositionId() << "\t" << std::setprecision(6) << std::fixed << d << "\t" << S.getPosition(i).getCurrentIdentity().getResidueName() << "\t" << S.getPosition(j).getCurrentIdentity().getResidueName();
  return ss.str();
}


// ---- OLD CODE

// double contactProbability (Position& posi, Position& posj, options& iopts) {
//   int n = 0; int c = 0;
//   for (int ri = 0; ri < posi.getTotalNumberOfRotamers(); ri++) {
//     posi.setActiveRotamer(ri);
//     for (int rj = 0; rj < posj.getTotalNumberOfRotamers(); rj++) {
//       posj.setActiveRotamer(rj);
// 
//       // are the two rotamers in contact?
//       bool cont = false;
//       for (int ai = 0; ai < posi.atomSize(); ai++) {
//         Atom& atomi = posi.getAtom(ai);
//         if (atomi.getName().compare(0, 1, "H") == 0) continue; // skip hydrogens
//         if (isBackbone(atomi)) continue;                       // skip backbone atoms
//         for (int aj = 0; aj < posj.atomSize(); aj++) {
//           Atom& atomj = posj.getAtom(aj);
//           if (atomj.getName().compare(0, 1, "H") == 0) continue;
//           if (isBackbone(atomj)) continue;
//           if (atomi.distance(atomj) < iopts.contDist) {
//             cont = true;
//             break;
//           }
//         }
//         if (cont) break;
//       }
// 
//       n++;
//       if (cont) {
//         c++;
//       }
//     }
//   }
//   return c*1.0/n;
// }
