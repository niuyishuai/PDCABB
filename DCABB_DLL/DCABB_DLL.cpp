/*********************************************
Language:	C++
Project:	DCABB_PARAL_DLL export functions
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Jan 2018
*********************************************/

#include "stdafx.h"

#ifndef  PDCABBLIB
#define PDCABBLIB extern "C" _declspec(dllexport)
#endif // ! PDCABBLIB

static DCABB* dcabb;
static DCABB::OUTPUT dcabb_output;

// declarations for dc decompositions
const double G(const Matrix X);
const double H(const Matrix X);
const Matrix DH(const Matrix X);
const Matrix DGC(const Matrix X);

// Optimization with DCABB algorithm
PDCABBLIB void Optimize() {
	dcabb_output = dcabb->Optimize();
}

// Set double parameter
PDCABBLIB void SetDblParam(int param, double value) {
	dcabb->SetDblParam((DCABB::DCABB_PARAM_DBL)param, value);
}

// Set int parameter
PDCABBLIB void SetIntParam(int param, int value) {
	dcabb->SetIntParam((DCABB::DCABB_PARAM_INT)param, value);
}

// Set bool parameter
PDCABBLIB void SetBoolParam(int param, bool value) {
	dcabb->SetBoolParam((DCABB::DCABB_PARAM_BOOL)param, value);
}

// Set Paralle mode
PDCABBLIB void SetParalMode(bool value) {
	dcabb->SetBoolParam(DCABB::DCABB_PARAM_BOOL::USEPARAL, value);
}

// Set Gap tolerance
PDCABBLIB double SetTolGap(double value) {
	dcabb->SetDblParam(DCABB::DCABB_PARAM_DBL::TOLGAP, value);
	return dcabb->TolGap;
}

// Set tolerance for restart DCA
PDCABBLIB double SetTolGapRestartDCA(double value) {
	dcabb->SetDblParam(DCABB::DCABB_PARAM_DBL::TOLGAPRESTARTDCA, value);
	return dcabb->TolGapRestart;
}

// Read model of lp or mps format
PDCABBLIB void ReadModel(const char* filename) {
	dcabb = new DCABB();
	// read model
	if (dcabb->ReadModel(filename) == false) { return; }
	// create a dcp
	dcabb->CreateDCP(G, H, DH, DGC);
	// set parameters
	dcabb->SetBoolParam(DCABB::DCABB_PARAM_BOOL::USEPARAL, true);
	dcabb->SetIntParam(DCABB::DCABB_PARAM_INT::NBPARALNODES, omp_get_num_procs());
	dcabb->SetDblParam(DCABB::DCABB_PARAM_DBL::TOLGAP, 1e-2);
	dcabb->SetDblParam(DCABB::DCABB_PARAM_DBL::TOLGAPRESTARTDCA, 1e-2);
	dcabb->SetBoolParam(DCABB::DCABB_PARAM_BOOL::VERBOSE, false);
}

// Get solution status, 0 for optimized, -1 for unbounded or infeasible
PDCABBLIB int GetSolutionStatus() {
	return dcabb_output.stat;
}

// Get optimal objective value
PDCABBLIB double GetObjVal() {
	return dcabb_output.fopt;
}

// Get number of variables
PDCABBLIB int GetNbVars() {
	return dcabb->nvars;
}

// Get number of integer variables
PDCABBLIB int GetNbIntVars() {
	return dcabb->nintvars;
}

// 1 for integer variable, 0 for continuous variable
PDCABBLIB int GetVarType(const int idx) {
	for (int i = 0; i < dcabb->nintvars; i++) {
		if (dcabb->idxintvars[i] == idx) { return 1; }
	}
	return 0;
}

// Get optimal solution (index from 0)
PDCABBLIB double GetX(const int idx) {
	return dcabb_output.xopt[idx][0];
}

// Get solution time
PDCABBLIB double GetCPUTime() {
	return dcabb_output.time;
}

// Get number of restarting dca
PDCABBLIB int GetNbDCA() {
	return dcabb_output.nbstartdca;
}

// Get number of nodes
PDCABBLIB int GetNbNodes() {
	return dcabb_output.nbnodes;
}

// Get number of iterations
PDCABBLIB int GetNbIter() {
	return dcabb_output.nbiter;
}

// Get variable's name (index from 0)
PDCABBLIB char * GetVarName(const int idx) {
	return dcabb->varnames[idx];
}


/******** The following functions are necessary to be defined before using dca ******/
// function G
const double G(const Matrix X) {
	if (dcabb->P.isempty())
		return (dcabb->f.trans()*X)[0][0];
	else
		return (0.5*X.trans()*dcabb->P*X + dcabb->f.trans()*X)[0][0];
};
// function H
const double H(const Matrix X) {
	double px = 0.0;
	for (auto i = 0; i < dcabb->nintvars; i++) {
		// p(x) = sum(min(xi,1-xi)), for all xi integer variables
		px += fmin(X[dcabb->idxintvars[i]][0], 1 - X[dcabb->idxintvars[i]][0]);
	}
	if (dcabb->Q.isempty())
		return -dcabb->t*px;
	else
		return (0.5*X.trans()*dcabb->Q*X)[0][0] - dcabb->t*px;
};
// function sub-differential of H
const Matrix DH(const Matrix X) {
	Matrix dpx = zeros(dcabb->nvars, 1);
	for (auto i = 0; i < dcabb->nintvars; i++) {
		dpx[dcabb->idxintvars[i]][0] = (X[dcabb->idxintvars[i]][0] > 0.5 ? -1.0 : 1.0);
	}
	if (dcabb->Q.isempty())
		return -dcabb->t*dpx;
	else
		return dcabb->Q * X - dcabb->t*dpx;
}
// function sub-differential of conjugate of G
// equivalent to solving a convex optimization
const Matrix DGC(const Matrix X) {
	CPX cpxsolver;
	OBJ subobj;
	Matrix ret;
	int status;

	// initialize obj
	subobj.Q = dcabb->P;
	subobj.c = dcabb->f - 1.0 * X;
	subobj.objsense = CPX_MIN;
	if (dcabb->P.isempty())
		status = cpxsolver.SolveLP(subobj, dcabb->C.lc);
	else
		status = cpxsolver.SolveQP(subobj, dcabb->C.lc, dcabb->C.qc);

	if (status == CPX_STAT_OPTIMAL)
		ret = cpxsolver.GetOptx();
	return ret;
}

