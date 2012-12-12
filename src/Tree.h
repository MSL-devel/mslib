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

#ifndef TREE_H
#define TREE_H

#include "Predicate.h"

#include <vector>
#include <map>
#include <queue>
#include <string>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <sys/types.h>
#include <iostream>
using namespace std;

namespace MSL { 
template<class T>
class Tree {
	public:
		Tree() { parent = NULL; visitFlag = false; depth= 0; subTree.clear();}
		Tree(Tree *_parent) { parent = _parent; visitFlag = false; depth = _parent->depth-1;subTree.clear();}
		~Tree() {
			// NOTE: this is deleting a pointer that was passed from the outside
			//       not an ideal fix
			for(typename std::vector<Tree<T> *>::iterator it = subTree.begin(); it != subTree.end(); it++) {
				delete *it;
			}
		}

		void operator=(Tree<T> & _tree);
		//void iterateDepthFirst(function pointer takes type T* as arguement);
		//void iterateBreadthFirst(function pointer takes type T* as arguement);
		
		inline T* getData() { return &data; }
		inline void setData(T &_data) { data = _data;}


		inline Tree<T>* getParent()   { return parent; }
		inline void setParent(Tree<T> *_parent) { parent = _parent; }

		inline void addSubTree(Tree<T> *_child) { subTree.push_back(_child); }
		
		Tree<T>* operator[](size_t n);

		inline void setVisitFlag(bool _flag) { visitFlag = _flag; }
		inline bool getVisitFlag()    	     { return visitFlag;  }

		inline int getNumSubtrees()          { return subTree.size();}

		int maxDepth(Tree<T>* _t);
		void printTree(int middle);
		void printTreeDFS(Tree<T>* _t);

		void printNode(Tree<T> *cur, int numTabs);
		void printTreeRec(Tree<T> *cur, int level);
		void printTree();
	private:

		T data;
		Tree<T> *parent;
		std::vector<Tree<T> *> subTree;
		std::vector<int> subTreeSize; // storage for size of each subtree 

		bool visitFlag;
		int depth;
	
};
template<class T>
void Tree<T>::printNode(Tree<T> *cur, int numTabs){
	  for (uint i = 0; i < numTabs; i++){
	    cout << "\t";
	  }
	  cout << (*cur->getData())<<endl;
}
template<class T>
void Tree<T>::printTreeRec(Tree<T> *cur, int level){

  if (cur->getNumSubtrees() == 0){
    printNode(cur, level);
    return;
  }
  if (cur->getNumSubtrees() == 1){
    printNode(cur, level);
    level++;
    printTreeRec((*cur)[0],level);
  }
  if (cur->getNumSubtrees() == 2){
    printNode(cur,level);
    level++;
    printTreeRec((*cur)[0],level);
    printTreeRec((*cur)[1],level);
  }


}
template<class T>
void Tree<T>::printTree(){
	printTreeRec(this,0);
}
template<class T>
void Tree<T>::printTree(int middle){

	Tree<T> *cur  = this;
	std::queue<Tree<T> *> Q;
	std::map<Tree<T> *, int> visited;
	std::map<Tree<T> *, int> numChildrenVisited;

	std::vector< std::vector<Tree<T> *> > treesEachDepth;
	treesEachDepth.resize(cur->maxDepth(cur));

	Q.push(cur);
	visited[cur] = middle;
	numChildrenVisited[cur] = 0;
	cur->depth = 0;
	//std::cout << "START"<<std::endl;
	while (!Q.empty()){
		
		cur = Q.front(); Q.pop();
		if (cur == NULL) { Q.pop(); continue;} 

		//std::cout << "Trying :"<<(*cur->getData())<<std::endl;
		if (cur->getVisitFlag()) continue;
		//std::cout << "NEW!: "<<cur->getNumSubtrees()<<std::endl;
		
		// Visit Node
		cur->setVisitFlag(true);
		Tree<T> *parent = cur->getParent();

		if (parent != NULL){

			// Keep track of depth..
			cur->depth = parent->depth+1;
			
			// Keep track of number of child nodes visited
			numChildrenVisited[parent]++;
			numChildrenVisited[cur] = 0;

			// Figure out proper spacing
			visited[cur] = visited[parent];
			if (numChildrenVisited[parent] == 1){
				visited[cur] -= (middle/2)/(cur->depth+1);
			}
			if (numChildrenVisited[parent] == 2){
				visited[cur] += (middle/2)/(cur->depth+1);
			}

		} else {
			cur->depth = 0;
		}


		for (uint i = 0; i < cur->getNumSubtrees();i++){
			//std::cout << "ADDING: "<<*(cur->getData())<<" "<<(*(*cur)[i]->getData())<<std::endl;
			Q.push((*cur)[i]);
		}

		treesEachDepth[cur->depth].push_back(cur);
		

	}


	for (uint i = 0; i < treesEachDepth.size();i++){
		std::string line = "                                                                                                                                                                                                                                        ";
		std::string branchLine = "                                                                                                                                                                                                                                        ";

		for (uint j = 0; j < treesEachDepth[i].size();j++){

			Tree<T> *tmp = treesEachDepth[i][j];

			std::stringstream ss;
			ss << (*tmp->getData());
			line.insert(visited[tmp],ss.str());
			if (tmp->parent == NULL) continue;

			int pos = visited[tmp];
			int offset = middle/(5+tmp->depth+1);

			if (visited[tmp] < visited[tmp->parent]) {
				pos += offset;
				branchLine.insert(pos,"/");

			} else {
				pos -= offset;
				branchLine.insert(pos, "\\");
			}


			
			
		}

		std::cout << std::endl<<branchLine<<std::endl<<std::endl<<line<<std::endl;
		
	}

	
}

template<class T>
int Tree<T>::maxDepth(Tree<T> *_t){

   if(_t==NULL)return(0);

   
   if (_t->getNumSubtrees() == 0) return(1);

   std::vector<int> heightEachChild(_t->getNumSubtrees(),0);
   for (uint i = 0 ; i < _t->getNumSubtrees();i++){
	   heightEachChild[i] = maxDepth((*_t)[i]);
   }
   
   std::sort(heightEachChild.begin(), heightEachChild.end(), std::greater<int>());
   return(heightEachChild[0]+1);
}


template<class T>
void Tree<T>::printTreeDFS(Tree<T>* _t){
	if (_t == NULL) return;

	if (_t->getNumSubtrees() == 0) {
		std::cout << "subtree0" << *_t->getData()<<" "<<std::endl;
		return;
	}


	if (_t->getNumSubtrees() == 1){
		printTreeDFS((*_t)[0]);
		return; 
	}

	printTreeDFS((*_t)[0]);
	std::cout << data<<" "<<_t->getNumSubtrees()<<std::endl;
	printTreeDFS((*_t)[1]);
}

template<class T>
void Tree<T>::operator=(Tree<T> &_tree){ 
	 data = *_tree.getData();
	 parent = _tree.getParent(); 
	 visitFlag = _tree.getVisitFlag();
	 depth   = _tree.depth;
	 for (uint i =0; i < _tree.getNumSubtrees();i++){
		 subTree.push_back(new Tree<T>(_tree[i]));
	 }
 }
template<class T>
Tree<T>* Tree<T>::operator[](size_t n) { 
	if (subTree.size() <= n) return NULL; 
	return subTree[n]; 
}

#endif

#ifndef __MACOS__
  template class Tree<double>;
  template class Tree<Predicate>;
#endif

}
