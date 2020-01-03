/*********************************************
Language:	C++
Project:	DC Program class
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Sept 2017
*********************************************/

#pragma once
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Matrix.h"
#include "CPX.h"

/* define a real valued multivariate convex function */
typedef const double(*PTR_CVX_FUNC) (const Matrix X);
/* define a function to compute subdifferential */
typedef const Matrix(*PTR_SUBDIFF_FUNC) (const Matrix X);

// convex constraint structure
// currently support linear constraints LC and quadratic constraints QC
struct CVXC
{
	LC lc;
	QC qc;
};

class DCP
{
public:
	PTR_CVX_FUNC G;
	PTR_CVX_FUNC H;
	PTR_SUBDIFF_FUNC DH; // subdifferential of H
	PTR_SUBDIFF_FUNC DGC; // subdifferential of conjugate of G
	CVXC C;
public:
	DCP() {};
	DCP(PTR_CVX_FUNC g, PTR_CVX_FUNC h, PTR_SUBDIFF_FUNC dh, PTR_SUBDIFF_FUNC dgc, CVXC& c) : G(g), H(h), DH(dh), DGC(dgc) { C = c; };
	~DCP() {};
	// Set pointer of dc component G
	int SetFuncG(PTR_CVX_FUNC g);
	// Set pointer of dc component H
	int SetFuncH(PTR_CVX_FUNC h);
	// Set subdifferential function of H
	int SetFuncDH(PTR_SUBDIFF_FUNC dh);
	// Set subdifferential function of conjugate of G
	int SetFuncDGC(PTR_SUBDIFF_FUNC dgc);
	// Set convex constraints
	int SetConstraints(const CVXC& c);
	// Evaluate G(X)
	double EvalG(const Matrix& X) const;
	// Evaluate H(X)
	double EvalH(const Matrix& X) const;
	// Evaluate G(X)-H(X)
	double EvalDC(const Matrix& X) const;
	// Evaluate subdifferential of H at X
	Matrix EvalDH(const Matrix& X) const;
	// Evaluate subdifferential of conjugate of G at X
	Matrix EvalDGC(const Matrix& X) const;
};

