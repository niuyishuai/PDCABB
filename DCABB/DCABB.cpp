/*********************************************
Language:	C++
Project:	DCABB class
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Sept 2017
*********************************************/

#include "stdafx.h"

// compute initial point
const Matrix DCABB::ComputeInitialPoint(const DCP& dcp) {
	int nvar = dcp.C.lc.A.cols();
	double _lb, _ub;
	_lb = min_m(dcp.C.lc.lb);
	_ub = max_m(dcp.C.lc.ub);
	if (InitMethod == 1) {
		return random(RANDOM_REAL, nvar, 1, _lb, _ub);
	}
	if (InitMethod == 2) {
		DCA dca;
		dca.SetIntParam(DCA::DCA_PARAM_VERBOSE, 0);
		double bestfopt = INF;
		Matrix bestsol;
		for (int i = 0; i < RandomRounds; i++) {
			auto X0 = random(RANDOM_REAL, nvar, 1, 0, 2);
			auto output = dca.Solve(dcp, X0);
			if (output.fopt < bestfopt) {
				bestfopt = output.fopt;
				bestsol = output.xopt;
			}
		}
		return bestsol;
	}
}

DCABB::DCABB()
{
	t = 1e+6;
	bblist = new BBList();
}


DCABB::~DCABB()
{
	delete dcp;
	delete bblist;
	dcp = NULL;
	bblist = NULL;
}


// Read a model from model file (support lp, mps)
bool DCABB::ReadModel(const char* modelname)
{
	// Read model
	CPX cpx;
	if (cpx.ReadModel((char*)modelname) == false) {
		throw exception("DCABB::ReadModel::Can't read model file");
		return false;
	}
	// Extract informations
	obj = cpx.ExtractObj();
	orgobjsense = obj.objsense;

	// Convert to minimization problem
	if (obj.objsense == CPX_MAX) {
		obj.Q = -1.0*obj.Q;
		obj.c = -1.0*obj.c;
	}
	obj.objsense = CPX_MIN;
	// Set linear part of objective function and linear constraints.
	f = obj.c;
	C.lc = cpx.ExtractLC();

	isQP = false;
	if (cpx.isQP() || cpx.isQCP()) {
		isQP = true;
		spdecomp(obj.Q, P, Q);
	}

	nvars = cpx.GetNumVars();
	auto ctypes = new char[nvars];
	nintvars = cpx.GetNumInt();
	idxintvars = new int[nintvars]; // index for integer variables
	cpx.ExtractVarType(ctypes);
	int idx = 0;
	for (auto i = 0; i < nvars; i++) {
		if (ctypes[i] == 'I' || ctypes[i] == 'B') {
			idxintvars[idx] = i;
			idx++;
		}
	}

	// Get variable names
	for (auto i = 0; i < nvars; i++) {
		varnames.push_back(cpx.GetVarName(i));
	}

	return true;
}

// Create a DCP from model
bool DCABB::CreateDCP(PTR_CVX_FUNC g, PTR_CVX_FUNC h, PTR_SUBDIFF_FUNC dh, PTR_SUBDIFF_FUNC dgc) {
	dcp = new DCP(g, h, dh, dgc, C);
	return true;
}

// Check feasibility of given point
bool DCABB::CheckFeas(Matrix & x) {
	for (auto i = 0; i < nintvars; i++) {
		if (fabs(x[idxintvars[i]][0] - round(x[idxintvars[i]][0])) > MATRIX_ZERO) {
			return false;
		}
	}
	return true;
}

// DCABB for global optimization
DCABB::OUTPUT DCABB::Optimize()
{
	clock_t time_start = clock();
	clock_t time_end;
	DCABB::OUTPUT output;

	lb = -INFINITY; // best lower bound
	ub = INFINITY; // best upper bound
	int iteration = 0; // iterations o
	int nodescounter = 0; // number of nodes
	int n = nintvars; // number of integer variables
	int m = nvars - nintvars; // number of continuous variables

	PrintCopyright();
	if (verbose == 1) {
		PrintProbInfo();
		cout << "Starting BBDCA ..." << endl;
		cout << "Branching method : " << BranchNodeMethod << ", Restart DCA Gap : " << TolGapRestart << endl;
	}

	// Step 0. Create root node
	root.LB = -INF;
	root.NodeProb = dcp;

	// Solve relaxtion at root
	RELAX_OUTPUT relax_output = SolveRelaxation(root);
	if (relax_output.stat != CPX_STAT_OPTIMAL) {
		output.fopt = INFINITY;
		output.nbiter = 0;
		output.nbnodes = 1; // only root node
		output.nbstartdca = 0;
		output.stat = -1; // infeasible or unbounded
		time_end = clock();
		output.time = TimeGap(time_start, time_end);
		return output;
	}

	// update lower bound
	if (lb < relax_output.fopt) {
		lb = relax_output.fopt;
	}

	// Check feasibility of computed solution at root
	if (CheckFeas(relax_output.xopt)) {
		// if feasible we get a global optimal solution
		output.xopt = relax_output.xopt;
		output.fopt = relax_output.fopt;
		if (orgobjsense==CPX_MAX){output.fopt = -relax_output.fopt;}
		output.nbiter = 0;
		output.nbnodes = 1; // nonly root
		output.nbstartdca = 0;
		output.stat = 0; // optimal
		time_end = clock();
		output.time = TimeGap(time_start, time_end);
		return output;
	}

	// Step 1. Start DCA
	// Compute initial point
	//auto X0 = ComputeInitialPoint(*dcp);
	auto X0 = relax_output.xopt;

	// Run DCA
	DCA dca;
	dca.SetIntParam(DCA::DCA_PARAM_VERBOSE, 0);
	DCA::OUTPUT dcaoutput = dca.Solve(*dcp, X0);
	nstartdca++;
	if (CheckFeas(dcaoutput.xopt)) {
		// if feasible we get an upper bound solution
		ub = dcp->EvalDC(dcaoutput.xopt);
		ubsol = dcaoutput.xopt;
	}

	// Starting B&B
	if (UseParal == false) {
		auto CurNode = CopyNode(root); // set current branch node
		CurNode.X = relax_output.xopt; // branch from lb solution at root
		NODE LeftNode, RightNode;

		while (iteration <= MaxBB || nodescounter <= MaxNodes) {
			// choose branch variable
			int chosenidx = BranchVar(CurNode.X);
			// nodes operations
			LeftNode = CopyNode(CurNode);
			RightNode = CopyNode(CurNode);

			// set bounds of left and right nodes
			LeftNode.NodeProb->C.lc.ub[chosenidx][0] = floor(CurNode.X[chosenidx][0]);
			RightNode.NodeProb->C.lc.lb[chosenidx][0] = ceil(CurNode.X[chosenidx][0]);

			// Left and Right Nodes Operations
			if (NodeOperation(LeftNode) == 0) {
				bblist->push_back(LeftNode);
			}
			if (NodeOperation(RightNode) == 0) {
				bblist->push_back(RightNode);
			}

			iteration++;
			nodescounter += 2;
			// update BBlist with current best upper bound (eliminate all nodes whose lower bound is greater than ub)
			bblist->updatelist(ub);

			//check if BBList is empty or not
			if (bblist->isempty()) {
				goto TERMINATE;
			}

			// update lb
			lb = bblist->at(bblist->getNode_MinLB()).LB;

			// compute gap
			double gap = ComputeGap(ub);
			if (gap <= TolGap) {
				goto TERMINATE;
			}
			// find branch node
			int pos = FindBranchNode();
			CurNode = CopyNode(bblist->at(pos));
			bblist->erase(pos);

			// print information
			if (verbose == 1) {
				if (iteration % 100 == 0)
					printf_s("  Iter = %5d, LB = %10.3f, UB = %10.3f, gap= %6.3f%%, nbnodes=%5d, nbdca=%5d, cputime=%6.3f\n", iteration, lb, ub, ComputeGap(ub) * 100, (int)bblist->size(), nstartdca, TimeGap(time_start, clock()));
			}
		}
	}
	// Parallel mode
	else if (UseParal == true) {
		auto CurNode = CopyNode(root); // set current branch node
		CurNode.X = relax_output.xopt; // branch from lb solution at root
		vector<NODE> SelectedNodes;
		SelectedNodes.push_back(CurNode);
		static int lastprintiter = 0;


		while (iteration <= MaxBB || nodescounter <= MaxNodes) {
			omp_set_num_threads(NbParalNodes);
#pragma omp parallel for
			for (auto i = 0; i < SelectedNodes.size(); i++) {
				// choose branch variable for each selected node
				int chosenidx = BranchVar(SelectedNodes[i].X);
				// nodes operations
				auto LeftNode = CopyNode(SelectedNodes[i]);
				auto RightNode = CopyNode(SelectedNodes[i]);

				// set bounds of left and right nodes
				LeftNode.NodeProb->C.lc.ub[chosenidx][0] = floor(SelectedNodes[i].X[chosenidx][0]);
				RightNode.NodeProb->C.lc.lb[chosenidx][0] = ceil(SelectedNodes[i].X[chosenidx][0]);

				// Left and Right Nodes Operations
				if (NodeOperation(LeftNode) == 0) {
#pragma omp critical
					{
						bblist->push_back(LeftNode);
					}
				}
				if (NodeOperation(RightNode) == 0) {
#pragma omp critical
					{
						bblist->push_back(RightNode);
					}
				}
#pragma omp critical
				{ // update shared variables
					iteration++;
					nodescounter += 2;
				}
			}

			// update BBlist with current best upper bound (eliminate all nodes whose lower bound is greater than ub)
			bblist->updatelist(ub);

			//check if BBList is empty or not
			if (bblist->isempty()) {
				goto TERMINATE;
			}

			// update lb
			lb = bblist->at(bblist->getNode_MinLB()).LB;

			// compute gap
			double gap = ComputeGap(ub);
			if (gap <= TolGap) {
				goto TERMINATE;
			}

			SelectedNodes.clear();
			// find branch nodes
			for (auto i = 0; i < NbParalNodes; i++) {
				int pos = FindBranchNode();
				if (pos < 0) break;
				SelectedNodes.push_back(CopyNode(bblist->at(pos)));
				bblist->erase(pos);
			}

			// print information
			if (verbose == 1) {
				if (iteration - lastprintiter > 100)
					lastprintiter = iteration;
					printf_s("  Iter = %5d, LB = %10.3f, UB = %10.3f, gap= %6.3f%%, nbnodes=%5d, nbdca=%5d, cputime=%6.3f\n", iteration, lb, ub, gap * 100, (int)bblist->size(), nstartdca, TimeGap(time_start, clock()));
			}
		}
	}
TERMINATE:
	if (verbose == 1) {
		printf_s("  Iter = %5d, LB = %10.3f, UB = %10.3f, gap= %6.3f%%, nbnodes=%5d, nbdca=%5d, cputime=%6.3f\n", iteration, lb, ub, ComputeGap(ub) * 100, (int)bblist->size(), nstartdca, TimeGap(time_start, clock()));
	}
	output.fopt = ub;
	if (orgobjsense == CPX_MAX) { output.fopt = -ub; }
	output.xopt = ubsol;
	output.nbiter = iteration;
	output.nbnodes = nodescounter;
	output.nbstartdca = nstartdca;
	output.stat = 0; // optimized
	time_end = clock();
	output.time = TimeGap(time_start, time_end);
	return output;
}

// Solve relaxed problem
DCABB::RELAX_OUTPUT DCABB::SolveRelaxation(NODE& node)
{
	// construct and solve convex relaxation problem
	CPX cpx;
	cpx.SetOptimizationTarget(CPX_OPTIMALITYTARGET_OPTIMALGLOBAL);
	int status;
	if (isQP) {
		status = cpx.SolveQP(obj, node.NodeProb->C.lc, node.NodeProb->C.qc);
		if (status == CPXMIP_OPTIMAL) {
			status = CPX_STAT_OPTIMAL;
		}
	}
	else {
		status = cpx.SolveLP(obj, node.NodeProb->C.lc);
	}
	RELAX_OUTPUT output;
	output.fopt = cpx.GetObjVal();
	output.xopt = cpx.GetOptx();
	output.stat = status;
	return output;
}


// Choose a branching variable
int DCABB::BranchVar(Matrix & x)
{
	int chosenidx = idxintvars[0];
	if (BranchVarMethod == 1) {
		for (auto i = 1; i < nintvars; i++) {
			if (fabs(x[chosenidx][0] - 0.5) > fabs(x[idxintvars[i]][0] - 0.5)) {
				chosenidx = idxintvars[i];
			}
		}
		return chosenidx;
	}
	return 0;
}


// Node operation (solve convex relaxation and restart dca if necessary)
int DCABB::NodeOperation(NODE & node)
{
	// solve relaxed problem
	auto node_output = SolveRelaxation(node);
	node.LB = node_output.fopt;
	node.X = node_output.xopt;
	if (node_output.stat != CPX_STAT_OPTIMAL) {
		// infeasible or unbounded, no need to add node
		return -1;
	}

	if (CheckFeas(node_output.xopt)) {
		// the lower bound solution of node is integer, then update upper bound, no need to add node
#pragma omp critical 
		{
			auto loc_ub = ub;
			ub = loc_ub;
		}
		if (ub > node_output.fopt) {
			ub = node_output.fopt;
			ubsol = node_output.xopt;
			return 1; // ub is updated, no need to add node into list
		}
		return 2; // ub is not updated, node problem is completly solved
	}

	// if the lower bound solution is not integer, try dca
	// if (ComputeGap(node_output.fopt) >= TolGapRestart) {
#pragma omp critical
	{
		auto loc_ub = ub;
		ub = loc_ub;
	}
	auto gap = ComputeGap(ub);
	if (gap >= TolGapRestart) {
		DCA dca;
		dca.SetIntParam(DCA::DCA_PARAM_VERBOSE, 0);
		//C = node.NodeProb->C;
		auto dca_output = dca.Solve(*node.NodeProb, node_output.xopt);
#pragma omp critical
		nstartdca++;
		if (CheckFeas(dca_output.xopt)) {
			// if dca found integer solution, then update upper bound
			if (ub > dca_output.fopt) {
				ub = dca_output.fopt;
				ubsol = dca_output.xopt; // ub is updated, but still need deeper search
			}
		}
	}
	return 0;
}

double DCABB::ComputeGap(double cur_ub) {
	double D = fmax(fabs(cur_ub), fabs(lb));
	if (D >= INF) { return INF; }
	return fabs(cur_ub - lb) / (1 + fmax(fabs(cur_ub), fabs(lb)));
}

// Find a branch node in BBList
int DCABB::FindBranchNode()
{
	if (BranchNodeMethod == 1) {
		return bblist->getNode_MinLB();
	}

	return 0;
}


// Show results
void DCABB::PrintResults(DCABB::OUTPUT & output) {
	cout << "The optimization is terminated" << endl;
	cout << "Optimal value : " << output.fopt << ", CPU time : " << output.time << "secs." << endl;
	cout << "Number of nodes : " << output.nbnodes << ", number of restarts DCA : " << output.nbstartdca << endl;
	if (output.stat == 0) {
		cout << "Solution status : Optimized." << endl;
	}
	else if (output.stat == -1) {
		cout << "Solution status : Infeasible or Unbounded." << endl;
	}
}

void DCABB::PrintProbInfo() {
	printf_s("**************************************************************************\n");
	printf_s("* The DC optimization model consists of \n");
	printf_s("* Variables: %d (Continuous variables: %d, Integer variables: %d)\n",nvars,nvars-nintvars,nintvars);
	printf_s("* Constraints: %d linear (%d equalities, %d inequalities), %d quadratic\n",
		(int)dcp->C.lc.b.rows() + dcp->C.lc.beq.rows(), (int)dcp->C.lc.beq.rows(), (int)dcp->C.lc.b.rows(), (int)dcp->C.qc.r.size());
	printf_s("* Objective: %s \n", (isQP==true)?"Quadratic":"Linear");
	printf_s("**************************************************************************\n");
	printf_s("* The machine has %d processor(s)\n", omp_get_num_procs());
	printf_s("* The parallel mode is: %s, with max number of threads %d\n", UseParal ? "ON" : "OFF", UseParal ? NbParalNodes : 0);
	printf_s("**************************************************************************\n");
}

void DCABB::PrintCopyright() {
	printf_s("**************************************************************************\n");
	cout << "* A C++ class of DCABB for solving general dc programming problem ver 1.0" << endl;
	cout << "* Author: 牛一帅 (Yi-Shuai Niu)" << endl;
	cout << "* Copyright since Sept 2017 by Yi-Shuai Niu, All right reserved." << endl;
	printf_s("**************************************************************************\n");
}

// Set int parameters
void DCABB::SetIntParam(enum DCABB::DCABB_PARAM_INT param, int val)
{
	switch (param)
	{
	case DCABB::MAXBB:
		MaxBB = val;
		break;
	case DCABB::MAXDCA:
		MaxDCA = val;
		break;
	case DCABB::MAXNODES:
		MaxNodes = val;
		break;
	case DCABB::BRANCHNODEMETHOD:
		BranchNodeMethod = val;
		break;
	case DCABB::BRANCHVARMETHOD:
		BranchVarMethod = val;
		break;
	case DCABB::INITMETHOD:
		InitMethod = val;
		break;
	case DCABB::NBPARALNODES:
		NbParalNodes = val;
		break;
	default:
		break;
	}
}


// Get int parameters
int DCABB::GetIntParam(enum DCABB::DCABB_PARAM_INT param)
{
	switch (param)
	{
	case DCABB::MAXBB:
		return MaxBB;
		break;
	case DCABB::MAXDCA:
		return MaxDCA;
		break;
	case DCABB::MAXNODES:
		return MaxNodes;
		break;
	case DCABB::BRANCHNODEMETHOD:
		return BranchNodeMethod;
		break;
	case DCABB::BRANCHVARMETHOD:
		return BranchVarMethod;
		break;
	case DCABB::INITMETHOD:
		return InitMethod;
		break;
	case DCABB::NBPARALNODES:
		return NbParalNodes;
		break;
	default:
		break;
	}
	return -1;
}



// Set double parameters
void DCABB::SetDblParam(enum DCABB::DCABB_PARAM_DBL param, double val)
{
	switch (param)
	{
	case DCABB::TOLF:
		TolF = val;
		break;
	case DCABB::TOLXOPT:
		TolXopt = val;
		break;
	case DCABB::TOLGAP:
		TolGap = val;
		break;
	case DCABB::TOLGAPRESTARTDCA:
		TolGapRestart = val;
		break;
	case DCABB::PENALTY:
		t = val;
		break;
	default:
		break;
	}
}

// Get double parameters
double DCABB::GetDblParam(enum DCABB::DCABB_PARAM_DBL param)
{
	switch (param)
	{
	case DCABB::TOLF:
		return TolF;
		break;
	case DCABB::TOLXOPT:
		return TolXopt;
		break;
	case DCABB::TOLGAP:
		return TolGap;
		break;
	case DCABB::TOLGAPRESTARTDCA:
		return TolGapRestart;
		break;
	case DCABB::PENALTY:
		return t;
		break;
	default:
		break;
	}
}


// Set bool parameters
void DCABB::SetBoolParam(enum DCABB::DCABB_PARAM_BOOL param, bool val)
{
	switch (param)
	{
	case DCABB::VERBOSE:
		verbose = val;
		break;
	case DCABB::USEPARAL:
		UseParal = val;
		break;
	default:
		break;
	}
}

// Get bool parameters
bool DCABB::GetBoolParam(enum DCABB::DCABB_PARAM_BOOL param)
{
	switch (param)
	{
	case DCABB::VERBOSE:
		return verbose;
		break;
	case DCABB::USEPARAL:
		return UseParal;
		break;
	default:
		break;
	}
	return false;
}
