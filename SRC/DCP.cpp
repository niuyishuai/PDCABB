/*********************************************
Language:	C++
Project:	DC Program class
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Sept 2017
*********************************************/
#include "stdafx.h"
#include "DCP.h"

// Set pointer of dc component G
int DCP::SetFuncG(PTR_CVX_FUNC g)
{
	G = g;
	return 0;
}


// Set pointer of dc component H
int DCP::SetFuncH(PTR_CVX_FUNC h)
{
	H = h;
	return 0;
}


// Set subdifferential function of H
int DCP::SetFuncDH(PTR_SUBDIFF_FUNC dh)
{
	DH = dh;
	return 0;
}

// Set subdifferential function of conjugate of G
int DCP::SetFuncDGC(PTR_SUBDIFF_FUNC dgc)
{
	DGC = dgc;
	return 0;
}


// Set convex constraints
int DCP::SetConstraints(const CVXC& c)
{
	C = c;
	return 0;
}


// Evaluate G(X)
double DCP::EvalG(const Matrix& X) const
{
	return G(X);
}


// Evaluate H(X)
double DCP::EvalH(const Matrix& X) const
{
	return H(X);
}


// Evaluate G(X)-H(X)
double DCP::EvalDC(const Matrix& X) const
{
	return G(X)-H(X);
}


// Evaluate subdifferential of H at X
Matrix DCP::EvalDH(const Matrix& X) const
{
	return DH(X);
}

// Evaluate subdifferential of conjugate of G at X
Matrix DCP::EvalDGC(const Matrix& X) const
{
	return DGC(X);
}