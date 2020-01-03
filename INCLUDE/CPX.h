/*********************************************
Language:	C++
Project:	CPX class for using cplex
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Sept 2017
*********************************************/
#pragma once
#include "Matrix.h"
#include "cplex.h"

using namespace std;

// objective
// 1/2 x'.Q.x + c'.x
struct OBJ
{
	int objsense; // CPX_MIN / CPX_MAX
	Matrix Q;
	Matrix c;
};

// quadratic constraints
// x'.Q[i].x + q[i]'.x sense[i] r[i]
struct QuadCons
{
	vector<Matrix> Q;
	vector<Matrix> q;
	vector<double> r;
	vector<char>   sense; //'L'<=,'G'>=
};

// linear constraints
// A.x<=b; Aeq.x==beq; lb<=x<=ub
struct LinCons
{
	Matrix A;
	Matrix b;
	Matrix Aeq;
	Matrix beq;
	Matrix lb;
	Matrix ub;
};

typedef struct QuadCons QC;
typedef struct LinCons LC;
typedef struct OBJ OBJ;

class CPX
{
protected:
	CPXENVptr env = NULL; // cplex env
	CPXLPptr lp = NULL; // cplex model
	int algo = 0; // algorithm, see CPX_ALG_

public:
	CPX();
	CPX(const CPX& cpx);
	~CPX();

	// Read a model from lp or mps file
	bool ReadModel(char* modelname);
	// Write Cplex model into file
	int WriteModel(char* modelname);
	// Optimize
	int Optimize();
	// Set Algorithm of Cplex (see CPX_ALG_)
	void SetAlgo(int algotype);
	// Get problem type
	int GetProbType();
	// Set double parameter
	int SetDblParam(int whichparam, double newvalue);
	// Set integer parameter
	int SetIntParam(int whichparam, int newvalue);
	// Get double parameter
	int GetDblParam(int whichparam, double& newvalue);
	// Get integer parameter
	int GetIntParam(int whichparam, int& newvalue);
	// Get solution stat
	int GetSolutionStat();
	// Get objective value
	double GetObjVal();
	// Get number of variables
	int GetNumVars();
	// Get optimal solution
	Matrix GetOptx();
	// Create linear programming problem
	bool CreateLP(OBJ obj, LC lc);
	// Create quadratic programming problem
	bool CreateQP(OBJ obj, LC lc, QC qc);
	// Create MIP
	bool CreateMIP(OBJ obj, LC lc, QC qc, char* ctype);
	// Solve linear programming problem
	int SolveLP(OBJ obj, LC lc);
	// Solve quadratic programming problem
	int SolveQP(OBJ obj, LC lc, QC qc);
	// Solve MIP
	int SolveMIP(OBJ obj, LC lc, QC qc, char* ctype);
	// Get best computed objective value for MIP
	double GetBestObj();
	// Get relative gap for mip (bestinteger - bestobjective) / (1e-10 + |bestinteger|)
	double GetMipRelgap();
	// Get number of node during MIP
	int GetNodeCnt();
	// Get number of cuts of the specified type (see CPX_CUT)
	int GetNumCuts(int cuttype);
	// Change coefficient of quadratic objective
	int SetObjQCoef(int row, int col, double coef);
	// Get coefficient of quadratic objective
	int GetObjQCoef(int row, int col, double* coef);
	// Get number of quadratic constraints
	int GetNumQcCons();
	// change bound of variable
	int SetBounds(int idxvar, char whichtype, double newval);
	// change bound of variable
	int GetBounds(int idxvar, double* lb, double* ub);
	// set coefficient of linear part
	int SetLCCoef(int row, int col, double newval);
	// get coefficient of linear part
	int GetLCCoef(int row, int col, double* val);
	// Set coefficient of linear part in objective
	int SetObjcCoef(int idx, double newval);
	// Set RHS
	int SetRhs(int idx, double newval);
	// Get RHS
	int GetRhs(int idx, double* rhs);
	// Get objsense
	int GetObjsense(int* objsense);
	// Set Objsense
	int SetObjsense(int objsense);
	// Set optimization target (see CPX_OPTIMALTARGET)
	int SetOptimizationTarget(int optimality);
	// Set sense of a linear constraint
	int SetLCSense(int idx, char sense);
	// Get sense of linear constraint
	int GetLCSense(int idx, char* sense);
	// Get variable type
	char GetVarType(int idx);
	// Set variable type
	int SetVarType(int idx, char type);
	// Get variable name
	char* GetVarName(int idx);
	// Set variable name
	int SetVarName(int idx, const char* varname);
	// extract objective
	OBJ ExtractObj();
	// Get linear objective coefficient
	double GetObjcCoef(int idx);
	// is QP
	bool isQP();
	// is LP
	bool isLP();
	// is MIP
	bool isMIP();
	// is QCP
	bool isQCP();
	// extract LC constraints
	LC ExtractLC();
	// Get number of binary variables
	int GetNumBin();
	// Get number of logical cores of the machine 
	int GetNumCors();
	// Get number of integer variables
	int GetNumInt();
	// Get number of LC
	int GetNumLC();
	// Extract variables type
	int ExtractVarType(char* ctype);
	// Get number of continuous variables
	int GetNumCont();
	// add new linear constraints
	int AddLC(const LC& lc);
	// add new QC constraints
	int AddQC(const QC& qc);
	// Delete linear constraints
	int DelLC(int begin, int end);
	// Delete QC constraints
	int DelQC(int begin, int end);
	// Delete variables
	int DelVars(int begin, int end);
	// overload =
	const CPX& operator=(const CPX& cpx);
};

