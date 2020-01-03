/*********************************************
Language:	C++
Project:	DCA class
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Sept 2017
*********************************************/
#include "stdafx.h"
#include "DCA.h"

// Set integer parameter
int DCA::SetIntParam(enum dca_param_int paramname, int value)
{
	switch (paramname)
	{
	case DCA::DCA_PARAM_MAXITER:
		maxIter = value;
		break;
	case DCA::DCA_PARAM_VERBOSE:
		verbose = value;
		break;
	default:
		break;
	}
	return 0;
}

// Set double parameter
int DCA::SetDblParam(enum dca_param_double paramname, double value)
{
	switch (paramname)
	{
	case DCA::DCA_PARAM_EPS_OBJ:
		eps_f = value;
		break;
	case DCA::DCA_PARAM_EPS_X:
		eps_x = value;
		break;
	default:
		break;
	}
	return 0;
}

// Get integer parameter
int DCA::GetIntParam(enum dca_param_int paramname)
{
	int ret;
	switch (paramname)
	{
	case DCA::DCA_PARAM_MAXITER:
		ret = maxIter;
		break;
	case DCA::DCA_PARAM_VERBOSE:
		ret = verbose;
		break;
	default:
		ret = -1;
		break;
	}
	return ret;
}

// Get double parameter
double DCA::GetDblParam(enum dca_param_double paramname)
{
	double ret;
	switch (paramname)
	{
	case DCA::DCA_PARAM_EPS_OBJ:
		ret = eps_f;
		break;
	case DCA::DCA_PARAM_EPS_X:
		ret = eps_x;
		break;
	default:
		ret = -1;
		break;
	}
	return ret;
}

// Solve DCP with DCA
// output.status: 0 convergence, -1 infeasible or unbounded, 1 maxiter exceed
DCA::OUTPUT DCA::Solve(const DCP& dcp, const Matrix& X0)
{
	OUTPUT output;
	int numDCA = 0;
	double diff_x, diff_f;
	double timeDCA = 0.0;
	Matrix XK, XK1;
	double FK, FK1;
	char msg[256];

	clock_t     start, end;

	// start dca
	start = clock();

	XK = X0;
	// Objective value at x0
	FK = dcp.EvalDC(XK);

	if (verbose == 1) {
		printcopyright();
		sprintf_s(msg, "%8s | %15s | %15s | %15s\n", "Iter", "Obj", "Dx", "Df");
		if (callback!=NULL) callback(msg);
		printf_s(msg);
		printf_s("------------------------------------------------------------------\n");
	}
	while (numDCA < maxIter)
	{
		numDCA++;

		// Sous gradient of H
		auto YK = dcp.EvalDH(XK);

		// Sous gradient of conjugate of G
		XK1 = dcp.DGC(YK);

		if (XK1.isempty() == false)
		{
			// compute new objective value
			FK1 = dcp.EvalDC(XK1);

			// Stopping conditions
			// diff_x = norm of delta_x
			diff_x = norm(XK1 - XK, NORM_MAX);
			//diff_f = relative error of f
			diff_f = fabs(FK - FK1) / (1.0 + fmax(fabs(FK), fabs(FK1)));

			if (verbose == 1) {
				sprintf_s(msg, "%8.5d | %15.4e | %15.4e | %15.4e\n", numDCA, FK1, diff_x, diff_f);
				printf_s(msg);
				if (callback != NULL) callback(msg);
			}
			//if ((diff_x < eps_x)||(diff_f < eps_obj))
			if (diff_f < eps_f || diff_x < eps_x)
			{
				output.status = 0;
				goto TERMINATE;
			}
			else
			{
				FK = FK1;
				XK = XK1;
			}
		}
		else
		{
			if (verbose == 1) {
				sprintf_s(msg, "Error: The problem is infeasible!!!\n");
				printf_s(msg);
				callback(msg);
			}
			output.status = -1; // infeasible or unbounded
			goto TERMINATE;
		}
	}
	output.status = 1; // max iteration exceed
	if (verbose == 1) {
		sprintf_s(msg, "Max iteration of DCA exceed.\n");
		printf_s(msg);
		if (callback != NULL) callback(msg);
	}
TERMINATE:
	end = clock();
	timeDCA = (double)(end - start) / CLOCKS_PER_SEC;
	if (verbose == 1) {
		printf("Iteration of DCA: %d\n", numDCA);
		printf("Computed solution: %lf, Computing time: %lf sec.\n", FK1, timeDCA);
	}
	output.iter = numDCA;
	output.time = timeDCA;
	output.fopt = FK1;
	output.xopt = XK1;

	return output;
}

// print copyright
void DCA::printcopyright() {
	cout << "*********************************************************************" << endl;
	cout << "A C++ class of DCA for solving general dc programming problem ver 1.0" << endl;
	cout << "Author: 牛一帅 (Yi-Shuai Niu)" << endl;
	cout << "Copyright since Sept 2017 by Yi-Shuai Niu, All right reserved." << endl;
	cout << "*********************************************************************" << endl;
}