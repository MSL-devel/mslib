#include "Atom3DGrid.h"

using namespace MSL;
using namespace std;


Atom3DGrid::Atom3DGrid() {
	AtomPointerVector a;
	setup(a, 1.0);
}

Atom3DGrid::Atom3DGrid(AtomPointerVector & _atoms, double _gridSize) {
	setup(_atoms, _gridSize);
	buildGrid();
}

Atom3DGrid::~Atom3DGrid() {
}

void Atom3DGrid::setup(AtomPointerVector & _atoms, double _gridSize) {
	gridSize = _gridSize;
	atoms = _atoms;
	xMin = 0.0;
	xMax = 0.0;
	yMin = 0.0;
	yMax = 0.0;
	zMin = 0.0;
	zMax = 0.0;
	xSize = 0;
	ySize = 0;
	zSize = 0;
}


void Atom3DGrid::buildGrid() {

	xMin = 0.0;
	xMax = 0.0;
	yMin = 0.0;
	yMax = 0.0;
	zMin = 0.0;
	zMax = 0.0;
	// calculate the max dimensions of the box
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		if (*k != NULL) {
			double x = (*k)->getX();
			double y = (*k)->getY();
			double z = (*k)->getZ();
			if (x < xMin) {
				xMin = x;
			}
			if (x > xMax) {
				xMax = x;
			}
			if (y < yMin) {
				yMin = y;
			}
			if (y > yMax) {
				yMax = y;
			}
			if (z < zMin) {
				zMin = z;
			}
			if (z > zMax) {
				zMax = z;
			}
		}
	}
	double xLen = xMax - xMin;
	double yLen = yMax - yMin;
	double zLen = zMax - zMin;
	double xCenter = (xMax + xMin)/2;
	double yCenter = (yMax + yMin)/2;
	double zCenter = (zMax + zMin)/2;
	xSize = (int)(xLen/gridSize) + 2;  // add another bin to catch overflow due to loss of precision
	ySize = (int)(yLen/gridSize) + 2;  // for eg. int bin = int(6.9999999999999999999999500) turns out to be 7 
	zSize = (int)(zLen/gridSize) + 2;
	double xStart = xCenter - (double)xSize/2 * gridSize;
	double yStart = yCenter - (double)ySize/2 * gridSize;
	double zStart = zCenter - (double)zSize/2 * gridSize;
	//cout << xMax << " " << xMin << " " << xLen << " " << xSize << " " << xStart << " " << xStart + gridSize * xSize << endl;
	//cout << yMax << " " << yMin << " " << yLen << " " << ySize << " " << yStart << " " << yStart + gridSize * ySize << endl;
	//cout << zMax << " " << zMin << " " << zLen << " " << zSize << " " << zStart << " " << zStart + gridSize * zSize << endl;

	grid = vector<vector<vector<AtomPointerVector> > >(xSize, vector<vector<AtomPointerVector> >(ySize, vector<AtomPointerVector>(zSize, AtomPointerVector(0))));

	atomIndeces = vector<vector<unsigned int> >(atoms.size(), vector<unsigned int>(3, 0));
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		if (*k != NULL) {
			double x = (*k)->getX();
			double y = (*k)->getY();
			double z = (*k)->getZ();
			
			unsigned int xBin = (int)((x - xStart)/gridSize);
			unsigned int yBin = (int)((y - yStart)/gridSize);
			unsigned int zBin = (int)((z - zStart)/gridSize);
			//cout << **k << " " << xBin << " " << yBin << " " << zBin << endl;
			grid[xBin][yBin][zBin].push_back(*k);
			atomIndeces[k-atoms.begin()][0] = xBin;
			atomIndeces[k-atoms.begin()][1] = yBin;
			atomIndeces[k-atoms.begin()][2] = zBin;
		}
	}

}

AtomPointerVector Atom3DGrid::getNeighbors(unsigned int _atomIndex) {
	AtomPointerVector out;
	unsigned int iMin = atomIndeces[_atomIndex][0];	
	unsigned int iMax = iMin;
	unsigned int jMin = atomIndeces[_atomIndex][1];	
	unsigned int jMax = jMin;
	unsigned int kMin = atomIndeces[_atomIndex][2];	
	unsigned int kMax = kMin;

	if (iMin > 0) {
		iMin--;
	}
	if (iMax < xSize - 1) {
		iMax++;
	}
	if (jMin > 0) {
		jMin--;
	}
	if (jMax < ySize - 1) {
		jMax++;
	}
	if (kMin > 0) {
		kMin--;
	}
	if (kMax < zSize - 1) {
		kMax++;
	}

	for (unsigned int i=iMin; i<=iMax; i++) {
		for (unsigned int j=jMin; j<=jMax; j++) {
			for (unsigned int k=kMin; k<=kMax; k++) {
				out.insert(out.end(), grid[i][j][k].begin(), grid[i][j][k].end());
			}
		}
	}
	for (AtomPointerVector::iterator k=out.begin(); k!=out.end(); k++) {
		if (atoms[_atomIndex] == *k) {
			out.erase(k);
			k--;
		}
	}
	return out;
}

