#include"System.h"
#include"CharmmSystemBuilder.h"
#include"AtomSelection.h"

using namespace std;
using namespace MSL;


ofstream dot;
class GraphNode {
	public:
		GraphNode(Atom* _a) {
			a = _a;
			visited = articulationPoint = false;
			dfn = l = -1; // largest possible value
			parent = NULL;
		}
		
		Atom* getAtom() {
			return a;
		}

		void addNeighbor(GraphNode* g) {
			neighbors.push_back(g);
		}

		vector<GraphNode*>& getNeighbors() {
			return neighbors;
		}

		void setVisited(bool _value) {
			visited = _value;
		}

		bool getVisited() {
			return visited;
		}

		int numNeighbors() {
			return neighbors.size();
		}
	
		unsigned int dfn,l; // for articulation point determination
		GraphNode* parent;
		bool articulationPoint;
	private:
		Atom* a;
		vector<GraphNode*> neighbors;
		bool visited;
};


void setArticulationPoints(vector<GraphNode*>& _nodes) {
	if(_nodes.size() < 3) {
		// No articulation points
		return;
	}

	// reset all visted to false
	for(int i = 0; i < _nodes.size(); i++) {
		_nodes[i]->setVisited(false);
	}

	// do a DFS and find articulation points
	// Algorithm here http://www.seas.gwu.edu/~ayoussef/cs212/graphsearch.html
	vector<GraphNode*> stack; // contains  graphnodes
	
	// visit the first node
	int num = 0;

	_nodes[0]->setVisited(true);
	stack.push_back(_nodes[0]);
	_nodes[0]->dfn = num;
	_nodes[0]->l = num;
	num++; 

	while(stack.size() > 0) {
		GraphNode * top = stack.back();
		vector<GraphNode*>& neighbors = top->getNeighbors();
		GraphNode* unvisited = NULL;
		for(int i = 0; i < neighbors.size(); i++) {
			if(neighbors[i]->getVisited()) {
			} else {
				unvisited = neighbors[i];
				break;
			}
		}

		if(unvisited != NULL) {
			// found an unvisited neighbor
			unvisited->setVisited(true);
			stack.push_back(unvisited);
			unvisited->dfn = num;
			unvisited->l = num;
			num++;
			unvisited->parent = top;
		} else {
			stack.erase(stack.end()-1);
			// no unvisited neighbor
			for(int i = 0; i < neighbors.size(); i++) {
				if(neighbors[i]->parent != top && neighbors[i]->dfn < top->dfn) {
					top->l = ((top->l > neighbors[i]->l) ? neighbors[i]->l : top->l); // top->l = min(top->l, neighbors[i]->l);
				} else {
					if (neighbors[i]->parent == top){
						top->l = ((top->l > neighbors[i]->l) ? neighbors[i]->l : top->l); // top->l = min(top->l, neighbors[i]->l);
						if(neighbors[i]->l >= top->dfn && top != stack[0]) { // top != root
							top->articulationPoint = true;
						}
					}
				}
			}

		}
	}
	
	vector<GraphNode*>& neighbors = _nodes[0]->getNeighbors();
	// if root is the parent of more than one node then it is an articulation point
	int children = 0;
	for(int i = 0; i < neighbors.size();i++) {
		if(neighbors[i]->parent == _nodes[0]) {
			children++;
			if(children > 1) {
				break;
			}
		}
	}

	if(children > 1) {
		_nodes[0]->articulationPoint = true;
	}
}
void printArticulationPoints(int compNum, char* pdbFile, vector<GraphNode*>& _nodes) {
	// print all the articulation points

	// print in the format pdb componentNumber componentSize <AP/N> <ResidueId> <neighborSize for APs> 
	for(int i = 0; i < _nodes.size(); i++) {
		if(_nodes[i]->articulationPoint) {
			dot << "INFO " << pdbFile << " " << compNum << " " << _nodes.size() << " AP " <<  (_nodes[i]->getAtom())->getIdentityId() << " " << _nodes[i]->numNeighbors() << endl;
		} else {
			dot << "INFO " << pdbFile << " " << compNum << " " << _nodes.size() << " N " <<  (_nodes[i]->getAtom())->getIdentityId() << endl;
		}
	}
	
}

// generates a graph using interCA distances as cutoffs and decomposes it into biconnected components
int main(int argc, char* argv[]) {

	if( argc != 4 ) {
		cout << "Usage: getInteractionGraph.cpp <pdbFile> <distCutoff> <graphFile>" << endl;
		cout << "Creates a graph in the graphFile based on interCA distance that can be plotted using graphviz. Also writes to stdout info regading each node in the graph in the format  <pdb> <componentNumber> <componentSize> <AP/N> <ResidueId> <neighborSize for APs>" << endl;
		exit(0);
	}

	double distCutoff = 6.0;
	distCutoff = MslTools::toDouble(argv[2]);

	System sys;
	if(!sys.readPdb(string(argv[1]))) {
		cout << "Unable to read " << argv[1] << endl;
		exit(0);
	}

//	cout << argv[1] << " " << distCutoff << " " << tertiary << " " << quaternary << endl;
	AtomSelection sel(sys.getAtomPointers());
	AtomPointerVector caAtoms = sel.select("CAs,name CA");
	//vector<Position*> positions = sys.getPositions();
	//cout << "Number of positions " << positions.size() << endl;

	//cout << "Number of CA Atoms " << caAtoms.size() << endl;

	vector<GraphNode*> graph(caAtoms.size());

	for(int i = 0; i < caAtoms.size(); i++) {
		graph[i] = new GraphNode(caAtoms[i]);
	}

	dot.open(argv[3]);
	if(!dot.is_open()) {
		cerr << "unable to write to " << argv[3] << endl;
		exit(0);
	}
	dot << "graph {" << endl;
	dot << "overlap=scale;" << endl;
	
	for(int i = 0; i < caAtoms.size(); i++) {
		for(int j = 0; j < i; j++) {
			double dist = caAtoms[i]->distance(*caAtoms[j]);
			string chain1 = caAtoms[i]->getChainId();
			string chain2 = caAtoms[j]->getChainId();
			int resNum1 = caAtoms[i]->getResidueNumber();
			int resNum2 = caAtoms[j]->getResidueNumber();
			int diff2 = (resNum1 - resNum2) * (resNum1 - resNum2);
			//cout << dist << endl;
			if(dist <= distCutoff && (diff2 > 16 || chain1 != chain2)) {
				//cout << chain1 << "," << chain2 << endl;
				dot << "\"" << caAtoms[i]->getAtomId() << "\" -- \"" << caAtoms[j]->getAtomId() << "\";" << endl;
				graph[i]->addNeighbor(graph[j]);
				graph[j]->addNeighbor(graph[i]);
				//cout << graph[i]->numNeighbors() << endl;
			}
		}
	}
	dot << "}" << endl;
/*
	vector<int> degrees;
	// get the degrees of all the nodes
	for(int i = 0; i < graph.size(); i++) {
		degrees.push_back(graph[i]->numNeighbors());
	}
	// sort degrees 
	sort(degrees.begin(),degrees.end());
	for(int i = 0; i < degrees.size(); i++) {
		cout << degrees[i] << endl;
	}

*/

	// and get all the connected components
	
	vector<vector<GraphNode*> > connectedComponents;
	while(graph.size() > 0) {
		//cout << "Remaining Nodes in Graph " << graph.size() << endl;
		connectedComponents.push_back(vector<GraphNode*>());
		vector<GraphNode*> queue;
		queue.push_back(graph[0]);
		queue[0]->setVisited(true);

		// Do BFS 
		while(queue.size() > 0) {
			//cout << "Remaining Nodes in queue" << queue.size() << endl;
			// visit node at head of queue and add it to current cluster
			GraphNode* head = queue[0];
			connectedComponents.back().push_back(head);
			// erase head node from queue
			queue.erase(queue.begin());

			vector<GraphNode*>& neighbors = head->getNeighbors();
			// add neighbors to queue if not visited
			//cout << "No of Neighbors " <<  neighbors.size() << endl;
			for(int i = 0; i < neighbors.size(); i++) {
				if(neighbors[i]->getVisited()) {
				} else {
					//cout << "Pushing" << endl;
					queue.push_back(neighbors[i]);
					(queue.back())->setVisited(true);
				}
			}

		}

		// erase all visited nodes from graph
		for(int i = 0; i < graph.size(); i++) {
			if(graph[i]->getVisited()) {
				graph.erase(graph.begin() + i);
				i--;
			}
		}

	}
	//dot << "//Number of connected components " << connectedComponents.size() << endl;

	// get the articulation points for each connected subgraph

	// print in the format INFO pdb <ResidueId> <componentNumber> <componentSize> <neighborSize>

	dot << "/*" << endl; // write other information within comments
	
	for(int i = 0; i < connectedComponents.size(); i++) {
		if(connectedComponents[i].size() > 2) {
			//dot << "//Component Size " << connectedComponents[i].size() << endl;
			setArticulationPoints(connectedComponents[i]);
		} else {
			// all are normal nodes - no articulation points
		}
		printArticulationPoints(i,argv[1],connectedComponents[i]);
	}
	dot << "*/" << endl; // write other information within comments

	cout << "Done" << endl;
} 
