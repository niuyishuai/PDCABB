/*********************************************
Language:	C++
Project:	CPX class for Branch and Bound List
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Sept 2017
*********************************************/
#include "stdafx.h"
#include "BBList.h"

BBList::BBList()
{
}

BBList::~BBList()
{
}


// get size of BBList
size_t BBList::size()
{
	return tree.size();
}

// get a node at the end of list
void BBList::push_back(NODE node)
{
	tree.push_back(node);
}

// check wheather the list is empty or not
bool BBList::isempty()
{
	if (tree.size() == 0)
		return true;
	else
		return false;
}

// pop and free the last node
void BBList::pop_back()
{
	delete tree.back().NodeProb;
	tree.pop_back();
}


// update list via current best upper bound
void BBList::updatelist(double UB)
{
	if (isempty() == false) {
		for (auto i = size(); i > 0; i--) {
			auto idx = i - 1;
			if (tree[idx].LB >= UB) {
				erase(idx);
			}
		}
	}
}

// add a new node at the end of list (no initialization)
void BBList::emplace_back()
{
	tree.emplace_back();
	tree.back().LB = -INFINITY;
	tree.back().NodeProb = new NODEPROB;
}

// get last node of list
NODE & BBList::back()
{
	return tree.back();
}


// erase the node at given position (index from 0)
void BBList::erase(size_t Pos)
{
	delete tree[Pos].NodeProb;
	tree.erase(tree.begin() + Pos);
}


// get first node of the list
NODE & BBList::front()
{
	return tree.front();
}


// get the node at position Pos
NODE & BBList::at(size_t Pos)
{
	return tree.at(Pos);
}

// get index of the Node with minimal lower bound (return -1 if exmpty list)
int BBList::getNode_MinLB()
{
	if (isempty() == false) {
		auto minlb = back().LB;
		auto idxminlb = size() - 1;
		for (auto i = size(); i > 0; i--) {
			auto idx = i - 1;
			if (at(idx).LB < minlb) {
				minlb = at(idx).LB;
				idxminlb = idx;
			}
		}
		return int(idxminlb);
	}
	return -1; // iff empty list
}


// empty BBlist
void BBList::empty()
{
	auto treesize = size();
	for (size_t i = 0; i < treesize; i++)
	{
		pop_back();
	}
}