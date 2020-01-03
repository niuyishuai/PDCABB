/*********************************************
Language:	C++
Project:	DCA class
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Sept 2017
*********************************************/
#pragma once

#include "DCP.h"

/* define a real valued multivariate convex function */
typedef void(*PTR_DCACALLBACK_FUNC) (const char* msg);

class DCA
{
public:
	struct OUTPUT
	{
		Matrix xopt; // optimal solution
		double fopt; // optimal objective value
		int iter; // iteration of dca
		int status; // soluation status
		double time; // computing time
	};
	enum dca_param_int
	{
		DCA_PARAM_MAXITER,
		DCA_PARAM_VERBOSE
	};
	enum dca_param_double {
		DCA_PARAM_EPS_OBJ,
		DCA_PARAM_EPS_X
	};
private:
	// parameters
	int	maxIter = 1000; // max iter of dca
	double eps_f = 1e-6; // epsilon for objective function
	double eps_x = 1e-6; // epsilon for optimal solution
	int verbose = 1; // 1: with verbose, 0: silence mode
	PTR_DCACALLBACK_FUNC callback;

public:
	DCA() {};
	~DCA() {};
	// Set integer parameter
	int SetIntParam(enum dca_param_int paramname, int value);
	// Set double parameter
	int SetDblParam(enum dca_param_double paramname, double value);
	// Get integer parameter
	int GetIntParam(enum dca_param_int paramname);
	// Get double parameter
	double GetDblParam(enum dca_param_double paramname);

	// Solve DCP with DCA
	OUTPUT Solve(const DCP& dcp, const Matrix& X0);

	void SetCallback(PTR_DCACALLBACK_FUNC pcallback){ callback = pcallback; }
	
	// print copyright
	void printcopyright();
};
