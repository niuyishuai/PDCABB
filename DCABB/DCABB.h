/*********************************************
Language:	C++
Project:	DCABB class
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Sept 2017
*********************************************/
#pragma once

#ifndef __DCABB_H__
#define __DCABB_H__

#define INF INFINITY

class DCABB
{
public:
	// data for dc programming
	Matrix P; // quadratic coefficients for G
	Matrix Q; // quadratic coefficients for H
	Matrix f; // linear coefficients for G
	CVXC C; // convex constraints

public:
	// output for dcabb
	struct OUTPUT
	{
		Matrix xopt; // optimal solution
		double fopt; // optimal objective value
		int nbiter; // iteration of BB
		int stat; // soluation status: 0 optimized, -1 infeasible or unbounded
		double time; // computing time
		int nbnodes; // number of nodes
		int nbstartdca;
	};
	// output for relaxation problem
	struct RELAX_OUTPUT
	{
		Matrix xopt; // optimal solution
		double fopt; // optimal objective value
		int stat; // solution status CPX_STAT
	};
private:
	DCP * dcp = NULL;
	BBList *bblist = NULL;
	OBJ obj; // initial objective function (only for linear and quadratic objective)
	int orgobjsense;
	bool isQP; // true if QP false if LP
	NODE root; // root node

	Matrix ubsol; // upper bound solution
	double ub;
	double lb;
	int nstartdca = 0;

public:
	int *idxintvars = NULL; // index of integer variables
	int nvars; // number of variables
	int nintvars; // number of integers
	vector<char*> varnames; // variable names

				  // parameters
	bool verbose = true; // true: with verbose, false: silence mode
	double t = 1e+5; // penalty parameter
	int BranchNodeMethod = 1; // 1: node with smallest lower bound 
	int BranchVarMethod = 1; // 1: branch at variable whose value is closest to 1/2. 2: branch at fractional component with maximal absolute value of objective coef
	int MaxBB = (int)1e+5; // max iterations for BB
	int MaxDCA = (int)1e+5; // max iterations for DCA
	int MaxNodes = 2 * MaxBB; // max number of nodes (each iteration create two nodes)
	double TolF = 1e-4; // tolerence for stopping DCA via error of objective
	double TolXopt = 1e-4; // tolerance for stopping DCA via error of X
	double TolGap = 1e-4; // tolerance of Gap(for terminating Branch-and-Bound. 0.1% by default)
	double TolGapRestart = 0.01;// 2 * TolGap; // tolerance of Gap for restarting DCA. 
								// Restarting DCA when the current Gap is greater than TolGapRestart. set INF if don't use TolGapRestart rule. 
	int InitMethod = 1; // 1: random inital point 2: random many initial point for dca
	int RandomRounds = 10; //
						   // The following parameters are used for parallel mode
	bool UseParal = true; // activate parallel mode
	int NbParalNodes = 20; // number of nodes to be solved parallely

	enum DCABB_PARAM_INT {
		MAXBB,
		MAXDCA,
		MAXNODES,
		BRANCHNODEMETHOD,
		BRANCHVARMETHOD,
		INITMETHOD,
		NBPARALNODES,
	};
	enum DCABB_PARAM_DBL {
		TOLF,
		TOLXOPT,
		TOLGAP,
		TOLGAPRESTARTDCA,
		PENALTY
	};
	enum DCABB_PARAM_BOOL {
		VERBOSE,
		USEPARAL
	};

public:
	DCABB();
	~DCABB();
	// Read a model from model file (support lp, mps)
	bool ReadModel(const char* modelname);
	// Create a DCP from model
	bool CreateDCP(PTR_CVX_FUNC g, PTR_CVX_FUNC h, PTR_SUBDIFF_FUNC dh, PTR_SUBDIFF_FUNC dgc);
	// DCABB for global optimization
	DCABB::OUTPUT Optimize();
	// Check feasibility of given point
	bool CheckFeas(Matrix & x);
	void SetPenaltyParam(double param) { t = param; }
	// Solve relaxed problem
	DCABB::RELAX_OUTPUT SolveRelaxation(NODE& node);
	// Choose a branching variable
	int BranchVar(Matrix & x);
	// Node operation (solve convex relaxation and restart dca if necessary)
	int NodeOperation(NODE & node);
	// Compute Gap at current solution with respect to the best lower bound
	double ComputeGap(double cur_ub);
	// Find a branch node in BBList
	int FindBranchNode();
	// Show results
	void PrintResults(DCABB::OUTPUT & output);
	const Matrix ComputeInitialPoint(const DCP& dcp);
	void PrintProbInfo();
	void PrintCopyright();
	// Set int parameters
	void SetIntParam(enum DCABB::DCABB_PARAM_INT param, int val);
	// Get int parameters
	int GetIntParam(enum DCABB::DCABB_PARAM_INT param);
	// Set double parameters
	void SetDblParam(enum DCABB::DCABB_PARAM_DBL param, double val);
	// Get double parameters
	double GetDblParam(enum DCABB::DCABB_PARAM_DBL param);
	// Set bool parameters
	void SetBoolParam(enum DCABB::DCABB_PARAM_BOOL param, bool val);
	// Get bool parameters
	bool GetBoolParam(enum DCABB::DCABB_PARAM_BOOL param);
};

// Compute CPU time gap
static double TimeGap(clock_t start, clock_t end) { return (double)(end - start) / CLOCKS_PER_SEC; };

#endif // !__DCABB_H__
