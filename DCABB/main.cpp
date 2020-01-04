/*********************************************
Language:	C++
Project:	test DCABB_PARAL
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Sept 2017
*********************************************/

#include "stdafx.h"

static DCABB* dcabb;

void testCPX(const char* filename);
void testDCABB(const char* filename);

// declarations for dc decompositions
const double G(const Matrix X);
const double H(const Matrix X);
const Matrix DH(const Matrix X);
const Matrix DGC(const Matrix X);

int main()
{
	//char *filename = "data/t001.lp";
	//char *filename = "data/nonconvexmiqp.lp";
	char *filename = "data/test_diff_bbdca.lp";
	//char *filename = "data/tmp2.lp";
	//char *filename = "data/p0033.mps";
	//char *filename = "data/myprob.lp";

	testCPX(filename);
	testDCABB(filename);
	return 0;
}

void testCPX(const char* filename) {
	CPX cpx;
	cpx.ReadModel((char*)filename);

	auto cpx_start = clock();
	cpx.SetOptimizationTarget(CPX_OPTIMALITYTARGET_OPTIMALGLOBAL);
	cpx.Optimize();
	cout << "The optimal solution of CPLEX : " << cpx.GetObjVal() << endl;
	cout << "CPU Time = " << TimeGap(cpx_start, clock()) << endl;
}

void testDCABB(const char* filename) {
	try {
		dcabb = new DCABB();
		// read model
		if (dcabb->ReadModel(filename) == false) { return; }
		// create a dcp
		dcabb->CreateDCP(G, H, DH, DGC);
		// set parameters
		dcabb->SetBoolParam(DCABB::DCABB_PARAM_BOOL::USEPARAL, true);
		dcabb->SetIntParam(DCABB::DCABB_PARAM_INT::NBPARALNODES, omp_get_num_procs());
		dcabb->SetDblParam(DCABB::DCABB_PARAM_DBL::TOLGAP, 1e-2);
		dcabb->SetDblParam(DCABB::DCABB_PARAM_DBL::TOLGAPRESTARTDCA, 2e-2); // INFINITY for non restart dca
		dcabb->SetBoolParam(DCABB::DCABB_PARAM_BOOL::VERBOSE, true);
		// optimize
		/*	vector<double> cputimes;
		for (auto i = 1; i < 50; i += 2) {
		dcabb->SetIntParam(DCABB::DCABB_PARAM_INT::NBPARALNODES, i);
		DCABB::OUTPUT output = dcabb->Optimize();
		dcabb->PrintResults(output);
		cputimes.push_back(output.time);
		}
		Matrix m(cputimes.size(), 1, cputimes);
		m.write("cputimes.txt");
		*/
		DCABB::OUTPUT output = dcabb->Optimize();
		// print results
		dcabb->PrintResults(output);
	}
	catch (char* errmsg) {
		cout << errmsg <<endl;
	}
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

