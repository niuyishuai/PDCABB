/*********************************************
Language:	C++
Project:	CPX class for using cplex
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Sept 2017
*********************************************/
#include "stdafx.h"
#include "CPX.h"

CPX::CPX()
{
	int status = 0;
	//Open CPLEX - environment
	env = CPXopenCPLEX(&status);
	if (env == NULL)
	{
		throw exception("CPX::CPX::Can't open the environement of CPLEX");
	}

	// Create the problem to solve
	lp = CPXcreateprob(env, &status, "lp");
	if (lp == NULL)
	{
		throw exception("CPX::CPX::Can't create model in CPLEX");
	}
	// activate paralle mode
	SetIntParam(CPX_PARAM_PARALLELMODE, CPX_PARALLEL_AUTO);
}

CPX::CPX(const CPX& cpx)
{
	int status = 0;
	//Open CPLEX - environment
	env = CPXopenCPLEX(&status);
	if (env == NULL)
	{
		throw exception("CPX::CPX::Can't open the environement of CPLEX");
	}

	// Copy problem
	lp = CPXcloneprob(env, cpx.lp, &status);
	if (lp == NULL)
	{
		throw exception("CPX::CPX::Can't copy model in CPLEX");
	}
	// activate paralle mode
	SetIntParam(CPX_PARAM_PARALLELMODE, CPX_PARALLEL_AUTO);
}

CPX::~CPX()
{
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
}



// Read a model from lp or mps file
bool CPX::ReadModel(char* modelname)
{
	auto status = CPXreadcopyprob(env, lp, modelname, NULL);
	if (status != 0)
	{
		return false;
	}
	return true;
}

// Write Cplex model into file
int CPX::WriteModel(char* modelname)
{
	return CPXwriteprob(env, lp, modelname, NULL);
}


// Optimize
// return CPX_STAT
int CPX::Optimize()
{
	auto probtype = GetProbType();
	switch (algo)
	{
	case CPX_ALG_PRIMAL:
		CPXprimopt(env, lp);
		break;
	case CPX_ALG_DUAL:
		CPXdualopt(env, lp);
		break;
	case CPX_ALG_BARRIER:
		CPXbaropt(env, lp);
		break;
	case CPX_ALG_MIP:
		CPXmipopt(env, lp);
		break;
	default:
		if (probtype == CPXPROB_LP) CPXlpopt(env, lp);
		else if (probtype == CPXPROB_QP || probtype == CPXPROB_QCP) CPXqpopt(env, lp);
		else if (probtype == CPXPROB_MILP || probtype == CPXPROB_MIQP || probtype == CPXPROB_MIQCP) CPXmipopt(env, lp);
	}
	return CPXgetstat(env,lp);
}


// Set Algorithm of Cplex (see CPX_ALG_)
void CPX::SetAlgo(int algotype)
{
	algo = algotype;
}


// Get problem type
int CPX::GetProbType()
{
	return CPXgetprobtype(env, lp);
}


// Set double parameter
int CPX::SetDblParam(int whichparam, double newvalue)
{
	return CPXsetdblparam(env, whichparam, newvalue);
}

// Set integer parameter
int CPX::SetIntParam(int whichparam, int newvalue)
{
	return CPXsetintparam(env, whichparam, newvalue);
}

// Get double parameter
int CPX::GetDblParam(int whichparam, double& newvalue)
{
	return CPXgetdblparam(env, whichparam, &newvalue);
}

// Set integer parameter
int CPX::GetIntParam(int whichparam, int& newvalue)
{
	return CPXgetintparam(env, whichparam, &newvalue);
}

// Get solution state
int CPX::GetSolutionStat()
{
	return CPXgetstat(env, lp);
}


// Get objective value
double CPX::GetObjVal()
{
	double obj;
	CPXgetobjval(env, lp, &obj);
	return obj;
}


// Get number of variables
int CPX::GetNumVars()
{
	return CPXgetnumcols(env, lp);
}


// Get optimal solution
Matrix CPX::GetOptx()
{
	Matrix ret;
	auto nbvars = GetNumVars();
	auto xx = new double[nbvars];
	auto stat = CPXgetx(env, lp, xx, 0, nbvars - 1);
	ret = Matrix(nbvars, 1, xx);
	delete[] xx;
	return ret;
}


// Create linear programming problem
// CPX_MIN/CPX_MAX {c'x : x in lc}
bool CPX::CreateLP(OBJ obj, LC lc)
{
	if (env == NULL || lp == NULL) { return false; }
	// check problem size consistency
	int n = obj.c.rows()*obj.c.cols();
	if (lc.A.isempty() == false && (lc.A.rows() != lc.b.rows() || lc.A.cols() != n)) throw exception("CPX::SolveLP::dimension of inequalities not match");
	if (lc.Aeq.isempty() == false && (lc.Aeq.rows() != lc.beq.rows() || lc.Aeq.cols() != n)) throw exception("CPX::SolveLP::dimension of equalities not match");
	if (lc.lb.isempty()) { lc.lb = zeros(n, 1); }
	else if (lc.lb.cols()*lc.lb.rows() != n) { throw exception("CPX::SolveLP::dimension of lower bound not match"); }
	if (lc.ub.isempty()) { lc.ub = INFINITY*ones(n, 1); }
	else if (lc.ub.cols()*lc.ub.rows() != n) { throw exception("CPX::SolveLP::dimension of upper bound not match"); }

	// Prepare CPX data
	// joint A and Aeq
	auto AA = lc.A;
	AA.add_block_by_row(lc.Aeq);
	// joint b and beq
	auto bb = lc.b;
	bb.add_block_by_row(lc.beq);
	// transform matrix into cpxmat format
	CPXMAT cpxA;
	format_cpxmat(AA, cpxA);
	// get objective coef
	auto cc = obj.c.data_col_major();
	// get constraint sense
	auto sense = new char[bb.rows()];
	for (auto i = 0; i < lc.A.rows(); i++) {
		sense[i] = 'L'; // <=
	}
	for (auto i = lc.A.rows(); i < AA.rows(); i++) {
		sense[i] = 'E'; // ==
	}

	// Copy problem into memory
	auto status = CPXcopylp(env, lp, AA.cols(), AA.rows(), obj.objsense, cc, bb.data_col_major(), sense,
		cpxA.matbeg, cpxA.matcnt, cpxA.matind, cpxA.matval, lc.lb.data_col_major(), lc.ub.data_col_major(), NULL);
	if (status != 0)
	{
		throw exception("CPX::SolveLP::Can't create CPLEX problem");
	}
	delete[] sense;
	return true;
}

// Create quadratic programming problem
// CPX_MIN/CPX_MAX {1/2 x'.Q.x + c'.x : x in lc, x in qc}
bool CPX::CreateQP(OBJ obj, LC lc, QC qc)
{
	CreateLP(obj, lc);
	
	// Create quadratic part of objective
	// Prepare cpx data
	if (obj.Q.isempty() == false) {
		cpxmat cpxQ;
		format_cpxmat(obj.Q, cpxQ);
		// copy cpxQ into cplex
		auto status = CPXcopyquad(env, lp, cpxQ.matbeg, cpxQ.matcnt, cpxQ.matind, cpxQ.matval);
		if (status != 0)
		{
			throw exception("CPX::SolveQP::Can't create quadratic objective of CPLEX problem");
		}
	}

	// Create quadratic constraints
	if (qc.Q.size() != 0) { AddQC(qc); }
	return true;
}


// Create MIP
// ctype: 'C' 'B' 'I' 'S' semi - continuous variable 'N' semi - integer variable
bool CPX::CreateMIP(OBJ obj, LC lc, QC qc, char* ctype)
{
	CreateQP(obj, lc, qc);
	// define variables type
	CPXcopyctype(env, lp, ctype);
	return true;
}


// Solve linear programming problem
// CPX_MIN/CPX_MAX {c'x : x in lc}
// return CPX_STAT
int CPX::SolveLP(OBJ obj, LC lc)
{
	CreateLP(obj, lc);
	return Optimize();
}

// Create quadratic programming problem
// CPX_MIN/CPX_MAX {1/2 x'.Q.x + c'.x : x in lc, x in qc}
// return CPX_STAT
int CPX::SolveQP(OBJ obj, LC lc, QC qc)
{
	CreateQP(obj, lc, qc);
	return Optimize();
}


// Solve MIP
// CPX_MIN/CPX_MAX {1/2 x'.Q.x + c'.x : x in lc, x in qc, x of ctype }
// ctype: 'C' 'B' 'I' 'S' semi - continuous variable 'N' semi - integer variable
// return CPXMIP_
int CPX::SolveMIP(OBJ obj, LC lc, QC qc, char* ctype)
{
	CreateMIP(obj, lc, qc, ctype);
	return Optimize();
}

// Get best computed objective value for MIP
double CPX::GetBestObj()
{
	double obj;
	CPXgetbestobjval(env, lp, &obj);
	return obj;
}


// Get relative gap for mip (bestinteger - bestobjective) / (1e-10 + |bestinteger|)
double CPX::GetMipRelgap()
{
	double relgap;
	CPXgetmiprelgap(env, lp, &relgap);
	return relgap;
}


// Get number of node during MIP
int CPX::GetNodeCnt()
{
	return CPXgetnodecnt(env, lp);
}


// Get number of cuts of the specified type (see CPX_CUT)
int CPX::GetNumCuts(int cuttype)
{
	int nbcuts;
	CPXgetnumcuts(env, lp, cuttype, &nbcuts);
	return nbcuts;
}


// Set coefficient of quadratic objective
int CPX::SetObjQCoef(int row, int col, double coef)
{
	return CPXchgqpcoef(env, lp, row, col, coef);
}


// Get coefficient of quadratic objective
int CPX::GetObjQCoef(int row, int col, double* coef)
{
	return CPXgetqpcoef(env, lp, row, col, coef);
}



// Get number of quadratic constraints
int CPX::GetNumQcCons()
{
	return CPXgetnumqconstrs(env, lp);
}


// change bound of variable
// whichtype: 'L','U','B'
int CPX::SetBounds(int idxvar, char whichtype, double newval)
{
	return CPXchgbds(env, lp, 1, &idxvar, &whichtype, &newval);
}

// get bounds of variable
int CPX::GetBounds(int idxvar, double* lb, double* ub)
{
	CPXgetlb(env, lp, lb, idxvar, idxvar);
	CPXgetub(env, lp, ub, idxvar, idxvar);
	return 0;
}

// set coefficient of linear part
int CPX::SetLCCoef(int row, int col, double newval)
{
	return CPXchgcoef(env, lp, row, col, newval);
}

// get coefficient of linear part
int CPX::GetLCCoef(int row, int col, double* val)
{
	return CPXgetcoef(env, lp, row, col, val);
}


// Set coefficient of linear part in objective
int CPX::SetObjcCoef(int idx, double newval)
{
	return CPXchgobj(env, lp, 1, &idx, &newval);
}


// Set RHS
int CPX::SetRhs(int idx, double newval)
{
	return CPXchgcoef(env, lp, idx, -1, newval);
}

// Get RHS
int CPX::GetRhs(int idx, double* rhs)
{
	return CPXgetrhs(env, lp, rhs, idx, idx);
}

// Get objsense
int CPX::GetObjsense(int* objsense)
{
	*objsense = CPXgetobjsen(env, lp);
	return 0;
}


// Set Objsense
int CPX::SetObjsense(int objsense)
{
	return CPXchgobjsen(env, lp, objsense);
}


// Set optimization target (see CPX_OPTIMALTARGET)
int CPX::SetOptimizationTarget(int optimality)
{
	return SetIntParam(CPX_PARAM_OPTIMALITYTARGET, optimality);
}


// Set sense of a linear constraint
int CPX::SetLCSense(int idx,char sense)
{
	return CPXchgsense(env, lp, 1, &idx, &sense);
}


// Get sense of linear constraint
int CPX::GetLCSense(int idx, char* sense)
{
	return CPXgetsense(env, lp, sense, idx, idx);
}


// Get variable type
char CPX::GetVarType(int idx)
{
	if (idx + 1 > GetNumVars()) { throw exception("CPX::GetVarType::index out of range"); }
	char type;
	CPXgetctype(env, lp, &type, idx, idx);
	return type;
}


// Set variable type
int CPX::SetVarType(int idx, char type)
{
	if (idx + 1 > GetNumVars()) { throw exception("CPX::GetVarType::index out of range"); }
	return CPXchgctype(env,lp,1,&idx,&type);
}



// extract objective
OBJ CPX::ExtractObj()
{
	OBJ obj;
	GetObjsense(&obj.objsense);
	int n = GetNumVars();
	obj.c.resize(n, 1);
	auto probtype = GetProbType();
	if (isQP()) {
		obj.Q.resize(n, n);
	}
	for (auto i = 0; i < n; i++) {
		obj.c[i][0] = GetObjcCoef(i);
		if (isQP())
		{
			for (auto j = i; j < n; j++) {
				GetObjQCoef(i, j, &obj.Q[i][j]);
				obj.Q[j][i] = obj.Q[i][j];
			}
		}
	}
	return obj;
}


// Get linear objective coefficient
double CPX::GetObjcCoef(int idx)
{
	double obj;
	CPXgetobj(env, lp, &obj, idx, idx);
	return obj;
}


// is QP
bool CPX::isQP()
{
	int probtype = GetProbType();
	return (probtype == CPXPROB_QP || probtype == CPXPROB_QCP || probtype == CPXPROB_MIQP || probtype == CPXPROB_MIQCP);
}


// is LP
bool CPX::isLP()
{
	int probtype = GetProbType();
	return (probtype == CPXPROB_LP);
}


// is MIP
bool CPX::isMIP()
{
	int probtype = GetProbType();
	return (probtype == CPXPROB_MILP || probtype == CPXPROB_MIQP || probtype == CPXPROB_MIQCP);
}

// is QCP
bool CPX::isQCP()
{
	int probtype = GetProbType();
	return (probtype == CPXPROB_QCP || probtype == CPXPROB_MIQCP);
}


// extract LC constraints
LC CPX::ExtractLC()
{
	LC lc;
	int m = GetNumLC();
	int n = GetNumVars();

	// count constraints type
	int cntE = 0; // number of == constraints 
	int cntL = 0; // <=
	int cntG = 0; // >=
	char* lcsense = new char[m];
	for (auto i = 0; i < m; i++) { 
		GetLCSense(i, &lcsense[i]); 
		switch (lcsense[i])
		{
		case 'E':
			cntE++;
			break;
		case 'L':
			cntL++;
			break;
		case 'G':
			cntG++;
			break;
		default:
			break;
		}
	}

	// allocate memory
	lc.A.resize(cntL + cntG, n);
	lc.b.resize(cntL + cntG, 1);
	lc.Aeq.resize(cntE, n);
	lc.beq.resize(cntE, 1);
	lc.lb.resize(n, 1);
	lc.ub.resize(n, 1);

	// get coefficients
	// get bounds
	for (auto i = 0; i < n; i++) {
		GetBounds(i, &lc.lb[i][0], &lc.ub[i][0]);
	}
	// get lc
	int idxE = 0; // index for equalities
	int idxNE = 0; // index for inequalities
	for (auto i = 0; i < m; i++) {
		// extract linear equalities ==
		if (lcsense[i] == 'E') {
			GetRhs(i, &lc.beq[idxE][0]);
			for (auto j = 0; j < n; j++) {
				GetLCCoef(i, j, &lc.Aeq[idxE][j]);
			}
			idxE++;
		}
		// extract linear inequalities <=
		if (lcsense[i] == 'L') {
			GetRhs(i, &lc.b[idxNE][0]);
			for (auto j = 0; j < n; j++) {
				GetLCCoef(i, j, &lc.A[idxNE][j]);
			}
			idxNE++;
		}
		// extract linear inequalities >=
		if (lcsense[i] == 'G') {
			GetRhs(i, &lc.b[idxNE][0]);
			lc.b[idxNE][0] = -lc.b[idxNE][0];
			for (auto j = 0; j < n; j++) {
				GetLCCoef(i, j, &lc.A[idxNE][j]);
				lc.A[idxNE][j] = -lc.A[idxNE][j];
			}
			idxNE++;
		}
	}
	delete[] lcsense;
	return lc;
}


// Get number of binary variables
int CPX::GetNumBin()
{
	return CPXgetnumbin(env, lp);
}


// Get number of logical cores of the machine 
int CPX::GetNumCors()
{
	int numcores;
	CPXgetnumcores(env, &numcores);
	return numcores;
}


// Get number of integer variables
int CPX::GetNumInt()
{
	return CPXgetnumint(env, lp) + GetNumBin();
}


// Get number of LC
int CPX::GetNumLC()
{
	return CPXgetnumrows(env, lp);
}


// Extract variables type
int CPX::ExtractVarType(char* ctype)
{
	int n = GetNumVars();
	for (int i = 0; i < n; i++) {
		ctype[i] = GetVarType(i);
	}
	return 0;
}


// Get number of continuous variables
int CPX::GetNumCont()
{
	return GetNumVars() - GetNumInt();
}


// add new linear constraints
int CPX::AddLC(const LC& lc)
{
	if (lc.b.rows() > 0) {
		auto rows = lc.b.rows();
		cpxmat cpxA;
		format_cpxmat(lc.A, cpxA);
		CPXaddrows(env, lp, 0, rows,*cpxA.matcnt,lc.b.data_col_major(),vector<char>(rows,'L').data(),cpxA.matbeg,cpxA.matind,cpxA.matval,NULL,NULL);
	}
	if (lc.beq.rows() > 0) {
		auto rows = lc.beq.rows();
		cpxmat cpxA;
		format_cpxmat(lc.Aeq, cpxA);
		CPXaddrows(env, lp, 0, rows, *cpxA.matcnt, lc.beq.data_col_major(), vector<char>(rows, 'E').data(), cpxA.matbeg, cpxA.matind, cpxA.matval, NULL, NULL);
	}
	return 0;
}


// add new QC constraints
int CPX::AddQC(const QC& qc)
{
	// Create quadratic constraints
	int nqc = (int)qc.r.size();
	cpxqcmat cpxQi;
	cpxvec cpxqi;
	for (auto i = 0; i < nqc; i++) {
		format_cpxqcmat(qc.Q[i], cpxQi);
		format_cpxvec(qc.q[i], cpxqi);
		auto status = CPXaddqconstr(env, lp, cpxqi.linnzcnt, cpxQi.quadnzcnt, qc.r[i], qc.sense[i], cpxqi.linind, cpxqi.linval,
			cpxQi.quadrow, cpxQi.quadcol, cpxQi.quadval, NULL);
		if (status != 0)
		{
			throw exception("CPX::SolveQP::Can't create quadratic constraints of CPLEX problem");
		}
	}
	return 0;
}


// Delete linear constraints
int CPX::DelLC(int begin, int end)
{
	return CPXdelrows(env, lp, begin, end);
}


// Delete QC constraints
int CPX::DelQC(int begin, int end)
{
	return CPXdelqconstrs(env, lp, begin, end);
}


// Delete variables
int CPX::DelVars(int begin, int end)
{
	return CPXdelcols(env, lp, begin, end);
}

// overload =
const CPX& CPX::operator=(const CPX& cpx) {
	int status = 0;
	//Open CPLEX - environment
	env = CPXopenCPLEX(&status);
	if (env == NULL)
	{
		throw exception("CPX::CPX::Can't open the environement of CPLEX");
	}
	// Copy problem
	lp = CPXcloneprob(env, cpx.lp, &status);
	if (lp == NULL)
	{
		throw exception("CPX::CPX::Can't copy model in CPLEX");
	}
	// activate paralle mode
	SetIntParam(CPX_PARAM_PARALLELMODE, CPX_PARALLEL_AUTO);

	return *this;
}

// Get variable name
char* CPX::GetVarName(int idx)
{
	if (idx + 1 > GetNumVars()) { throw exception("CPX::GetVarName::index out of range"); return NULL; }
	// determine how many space we need for varname
	int surplus;
	auto status = CPXgetcolname(env, lp, NULL, NULL, 0, &surplus, idx,idx);
	// allocate spaces
	int cur_colnamespace = -surplus;
	char** cur_colname = NULL;
	char* cur_colnamestore = NULL;
	if (cur_colnamespace > 0) {
		cur_colname = (char **)malloc(sizeof(char *));
		cur_colnamestore = (char *)malloc(cur_colnamespace);
		if (cur_colname == NULL ||
			cur_colnamestore == NULL) {
			throw exception("CPX::GetVarName::Failed to get memory for column names.");
			return NULL;
		}
		status = CPXgetcolname(env, lp, cur_colname, cur_colnamestore,
			cur_colnamespace, &surplus, idx, idx);
	}
	return cur_colnamestore;
}

// Set vairable name
int SetVarName(int idx, const char* varname) {
	return 0;
}
