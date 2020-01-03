/*********************************************
Language:	C++
Project:	CPX class for Branch and Bound List
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Sept 2017
*********************************************/
#pragma once
#include "vector"

using namespace std;

typedef DCP NODEPROB;

struct NODE {
	double LB; // lower bound value for node problem
	Matrix X; // current solution to make branch
	NODEPROB *NodeProb; // pointer to node problem
};

typedef NODE NODE;

class BBList
{
public:
	vector<NODE> tree; // using vector for list data

public:
	BBList();
	~BBList();
	// get size of BBList
	size_t size();
	// get a node at the end of list
	void push_back(NODE node);
	// check wheather the list is empty or not
	bool isempty();
	// pop and free the last node
	void pop_back();
	// update list via current best upper bound
	void updatelist(double UB);
	// add a new node at the end of list (no initialization)
	void emplace_back();
	// get last node of list
	NODE & back();
	// erase the node at given position (index from 0)
	void erase(size_t Pos);
	// get first node of the list
	NODE & front();
	// get the node at position Pos
	NODE & at(size_t Pos);
	// get index of the Node with minimal lower bound
	int getNode_MinLB();
	// empty BBlist
	void empty();
};

static NODE CopyNode(NODE& node) {
	NODE newnode;
	newnode.LB = node.LB;
	newnode.X = node.X;
	newnode.NodeProb = new NODEPROB(node.NodeProb->G, node.NodeProb->H, node.NodeProb->DH, node.NodeProb->DGC, node.NodeProb->C);
	return newnode;
}