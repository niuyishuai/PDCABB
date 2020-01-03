/*********************************************
Language:	C++
Project:	Matrix class
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Sept 2017
*********************************************/
#pragma once
#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <time.h>  
#include <stdio.h>
#include <direct.h>

#include "f2c_user.h"
#include "clapack.h"

using namespace std;

#ifndef MATRIX_ZERO
#define MATRIX_ZERO 1e-12
#endif
////////////////////////////////////////////
// Matrix class templete
// Define Basic MATRIX operations (supports all types as double, int etc.)
///////////////////////////////////////////
template <typename _T>
class MATRIX
{
protected:
	vector<vector<_T>> array;

public:
	/***** Constructors *****/

	// Constructors default
	MATRIX() : array(0) {}
	// Constructor with number of rows and cols (default zero matrix)
	MATRIX(int rows, int cols);
	// Constructor copy
	MATRIX(const MATRIX<_T>& m) { *this = m; }
	// Constructor with number of rows and cols and elements in a (col major order)
	MATRIX(int rows, int cols, _T* a);
	// Constructor with number of rows and cols and elements in a (col major order)
	MATRIX(int rows, int cols, vector<_T> a);
	// Constructor from sparse matrix
	MATRIX(int rows, int cols, int nz, int* I, int* J, _T* val);

	/***** Matrix props *****/

	// get number of rows
	int  rows() const { return (int)(array.size()); }
	// get number of cols
	int  cols() const { return rows() ? (int)(array[0].size()) : 0; }
	// is empty matrix
	bool isempty() const { return (rows() == 0 || cols() == 0); }
	// is square matrix
	bool issquare() const { return (!(isempty()) && rows() == cols()); }
	// is symmetric matrix
	bool issymmetric() const;

	/***** Data export *****/
	// export matrix as as array (col major order)
	_T* data_col_major() const;

	/***** Basic matrix operations *****/

	// resize matrix
	void resize(int rows, int cols);
	// add a row
	bool add_row(const vector<_T>& v);
	// add a col
	bool add_col(const vector<_T>& v);
	// add a matrix block by row
	bool add_block_by_row(const MATRIX<_T>& m);
	// add a matrix block by col
	bool add_block_by_col(const MATRIX<_T>& m);
	// swap two rows
	void swap_row(int row1, int row2);
	// swap two cols
	void swap_col(int col1, int col2);
	// vectorization
	const MATRIX<_T> vec() const;
	// transpose
	const MATRIX<_T> trans() const;
	// trace
	const _T trace() const;
	// diag
	const MATRIX<_T> diag() const;

	/***** IO Matrix *****/

	// Write matrix into a file
	void write(const char* filename);
	// Read a matrix from a file
	void read(const char* filename);

	/***** Overload operators *****/

	// overload +=
	const MATRIX<_T>& operator+=(const MATRIX<_T>& m);
	// overload -=
	const MATRIX<_T>& operator-=(const MATRIX<_T>& m);
	// overload *=
	const MATRIX<_T>& operator*=(const MATRIX<_T>& m);
	// overload ==
	bool operator==(const MATRIX<_T>& rhs) const;
	// overload !=
	bool operator!=(const MATRIX<_T>& rhs) const;
	// overload +
	const MATRIX<_T> operator+(const MATRIX<_T>& rhs) const;
	// overload -
	const MATRIX<_T> operator-(const MATRIX<_T>& rhs) const;
	// overload *
	const MATRIX<_T> operator*(const MATRIX<_T>& rhs) const;
	// overload | (joint block matrix)
	const MATRIX<_T> operator|(const MATRIX<_T>& rhs) const;
	// overload []
	const vector<_T>& operator[](int row) const { return array[row]; }
	vector<_T>& operator[](int row) { return array[row]; }
	// overload <<
	friend ostream& operator<< (ostream& out, const MATRIX<_T>& s) {
		for (auto i = 0; i < s.rows(); i++) {
			for (auto j = 0; j < s.cols(); j++) {
				out << "\t" << s[i][j] << "\t";
			}
			out << endl;
		}
		return out;
	}
};

/***** Constructors *****/

// Constructor with number of rows and cols (default zero matrix)
template <typename _T>
MATRIX<_T>::MATRIX(int rows, int cols) : array(rows)
{
	if (rows < 0 || cols < 0) { throw exception("MATRIX::MATRIX(int,int)::negative matrix size"); }
	for (int i = 0; i < rows; ++i)
	{
		array[i].resize(cols);
	}
}

// Constructor with number of rows and cols and elements in a (col major order)
template <typename _T>
MATRIX<_T>::MATRIX(int rows, int cols, _T* a) :array(rows) {
	if (rows < 0 || cols < 0) { throw exception("MATRIX::MATRIX(int,int)::negative matrix size"); }
	for (int i = 0; i < rows; ++i)
	{
		array[i].resize(cols);
	}
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			array[i][j] = a[i + j*rows];
		}
	}
}

// Constructor with number of rows and cols and elements in a (col major order)
template <typename _T>
MATRIX<_T>::MATRIX(int rows, int cols, vector<_T> a) :array(rows) {
	if (rows < 0 || cols < 0) { throw exception("MATRIX::MATRIX(int,int)::negative matrix size"); }
	for (int i = 0; i < rows; ++i)
	{
		array[i].resize(cols);
	}
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			array[i][j] = a[i + j*rows];
		}
	}
}

// Constructor from sparse matrix
template <typename _T>
MATRIX<_T>::MATRIX(int rows, int cols, int nz, int* I, int* J, _T* val) : MATRIX<_T>(rows, cols) {
	for (i = 0; i < nz; i++)
	{
		array[I[i]][J[i]] = val[i];
	}
}

// is a symmetric matrix (empty matrix is also symmetric)
template <typename _T>
bool MATRIX<_T>::issymmetric() const {
	if (rows()*cols() == 0 || rows() != cols()) {
		return false;
	}
	for (int i = 0; i < rows(); i++) {
		for (int j = i; j < cols(); j++) {
			if (fabs(array[i][j] - array[j][i])>MATRIX_ZERO) {
				return false;
			}
		}
	}
	return true;
}

// export matrix as as array (col major order)
template <typename _T>
_T* MATRIX<_T>::data_col_major() const {
	_T* ret;
	int row = rows();
	int col = cols();
	(ret) = (_T*)calloc(row*col, sizeof(_T));
	int s = 0;
	for (auto j = 0; j < col; j++) {
		for (auto i = 0; i < row; i++) {
			ret[s] = array[i][j];
			s++;
		}
	}
	return ret;
}

/***** Basic matrix operations *****/

// resize matrix
template <typename _T>
void MATRIX<_T>::resize(int rows, int cols)
{
	int rs = this->rows();
	int cs = this->cols();

	if (rows == rs && cols == cs)
	{
		return;
	}
	else if (rows == rs && cols != cs)
	{
		for (int i = 0; i < rows; ++i)
		{
			array[i].resize(cols);
		}
	}
	else if (rows != rs && cols == cs)
	{
		array.resize(rows);
	}
	else
	{
		array.resize(rows);
		for (int i = 0; i < rows; ++i)
		{
			array[i].resize(cols);
		}
	}
}


// add a new row
template <typename _T>
bool MATRIX<_T>::add_row(const vector<_T>& v)
{
	if (rows() == 0 || cols() == (int)v.size())
	{
		array.push_back(v);
	}
	else
	{
		return false;
	}

	return true;
}

// add a new col
template <typename _T>
bool MATRIX<_T>::add_col(const vector<_T>& v)
{
	if (cols() == 0 || rows() == (int)v.size())
	{
		int row = rows();
		int col = cols() + 1;
		resize(row, col);
		for (int i = 0; i < row; i++) {
			array[i][col - 1] = v[i];
		}
	}
	else
	{
		return false;
	}
	return true;
}

// add a matrix block by row
template <typename _T>
bool MATRIX<_T>::add_block_by_row(const MATRIX<_T>& m)
{
	if (rows() == 0 || cols() == (int)m.cols())
	{
		for (auto i = 0; i < m.rows(); i++) {
			add_row(m[i]);
		}
	}
	else
	{
		return false;
	}

	return true;
}


// add a matrix block by col
template <typename _T>
bool MATRIX<_T>::add_block_by_col(const MATRIX<_T>& m)
{
	if (cols() == 0 || rows() == (int)m.rows())
	{
		int row = m.rows();
		int oldcol = cols();
		int newcol = oldcol + m.cols();
		resize(row, newcol);
		for (auto j = oldcol; j < newcol; j++) {
			for (auto i = 0; i < row; i++) {
				array[i][j] = m[i][j - oldcol];
			}
		}
	}
	else
	{
		return false;
	}

	return true;
}

// swap two rows
template <typename _T>
void MATRIX<_T>::swap_row(int row1, int row2)
{
	if (fmax(row1, row2) > rows() || fmin(row1, row2) < 0) { throw exception("MATRIX::swap_row(int,int)::illegal input"); }
	if (row1 != row2)
	{
		vector<_T>& v1 = array[row1];
		vector<_T>& v2 = array[row2];
		vector<_T> tmp = v1;
		v1 = v2;
		v2 = tmp;
	}
}

// swap two cols
template <typename _T>
void MATRIX<_T>::swap_col(int col1, int col2)
{
	if (fmax(col1, col2) > cols() || fmin(col1, col2) < 0) { throw exception("MATRIX::swap_col(int,int)::illegal input"); }
	if (col1 != col2)
	{
		for (int i = 0; i < rows(); i++) {
			auto tmp = array[i][col1];
			array[i][col1] = array[i][col2];
			array[i][col2] = tmp;
		}
	}
}

// matrix vectorization
template <typename _T>
const MATRIX<_T> MATRIX<_T>::vec() const {
	MATRIX<_T> ret;
	if (isempty()) { return ret; }
	auto row = rows();
	auto col = cols();
	ret.resize(rows()*cols(), 1);
	int s = 0;
	for (auto j = 0; j < col; j++) {
		for (auto i = 0; i < row; i++) {
			ret[s][0] = array[i][j];
			s++;
		}
	}
	return ret;
}

// transpose matrix
template <typename _T>
const MATRIX<_T> MATRIX<_T>::trans() const
{
	MATRIX<_T> ret;
	if (isempty()) return ret;

	int row = cols();
	int col = rows();
	ret.resize(row, col);

	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			ret[i][j] = array[j][i];
		}
	}
	return ret;
}

// trace of matrix
template <typename _T>
const _T MATRIX<_T>::trace() const {
	auto minsize = (int)fmin(rows(), cols());
	_T ret = 0.0;
	for (int i = 0; i < minsize; i++) {
		ret += array[i][i];
	}
	return ret;
}

// return diagonal elements of matrix or form a diagnal matrix from vector 
template <typename _T>
const MATRIX<_T> MATRIX<_T>::diag() const {
	MATRIX<_T> ret;
	auto n = (int)fmax(cols(), rows());
	if (cols() == 1) {
		ret.resize(n, n);
		for (auto i = 0; i < n; i++) {
			ret[i][i] = array[i][0];
		}
	}
	else if (rows() == 1) {
		ret.resize(n, n);
		for (auto i = 0; i < n; i++) {
			ret[i][i] = array[0][i];
		}
	}
	else {
		auto m = (int)fmin(cols(), rows());
		ret.resize(m, 1);
		for (auto i = 0; i < m; i++) {
			ret[i][0] = array[i][i];
		}
	}
	return ret;
}

/***** IO Matrix *****/

// Write matrix into a file
template <typename _T>
void MATRIX<_T>::write(const char* filename) {
	FILE *pf;
	fopen_s(&pf, filename, "w+");
	fprintf_s(pf, "%d %d\n", rows(), cols());
	for (auto i = 0; i < rows(); i++) {
		for (auto j = 0; j < cols(); j++) {
			fprintf_s(pf, "%lf  ", array[i][j]);
		}
		fprintf_s(pf, "\n");
	}
	fclose(pf);
};

// Read a matrix from a file
template <typename _T>
void MATRIX<_T>::read(const char* filename) {
	FILE *pf;
	fopen_s(&pf, filename, "r+");
	int row, col;
	fscanf_s(pf, "%d %d", &row, &col);
	resize(row, col);
	for (auto i = 0; i < row; i++)
		for (auto j = 0; j < col; j++) fscanf_s(pf, "%lf", &array[i][j]);
	fclose(pf);
};


/***** Overload operators *****/

// Overload +=
template <typename _T>
const MATRIX<_T>& MATRIX<_T>::operator+=(const MATRIX<_T>& m)
{
	if (rows() != m.rows() || cols() != m.cols())
	{
		if (m.isempty()) { return *this; }
		else { throw exception("MATRIX::operator+=(MATRIX &)::dimension not match"); }
	}

	int r = rows();
	int c = cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			array[i][j] += m[i][j];
		}
	}

	return *this;
}

// Overload -=
template <typename _T>
const MATRIX<_T>& MATRIX<_T>::operator-=(const MATRIX<_T>& m)
{
	if (rows() != m.rows() || cols() != m.cols())
	{
		if (m.isempty()) { return *this; }
		else { throw exception("MATRIX::operator+=(MATRIX &)::dimension not match"); }
	}

	int r = rows();
	int c = cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			array[i][j] -= m[i][j];
		}
	}

	return *this;
}

// Overload *=
template <typename _T>
const MATRIX<_T>& MATRIX<_T>::operator*=(const MATRIX<_T>& m)
{
	if (cols() != m.rows())
	{
		throw exception("MATRIX::operator*=(MATRIX &)::dimension not match");
	}

	MATRIX<_T> ret(rows(), cols());

	int r = rows();
	int c = cols();

	/*	for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				double sum = 0.0;
				for (int k = 0; k < c; ++k)
				{
					sum += array[i][k] * m[k][j];
				}
				ret[i][j] = sum;
			}
		}
	*/
	for (int i = 0; i < r; ++i)
	{
		for (int k = 0; k < c; ++k)
		{
			double sum = array[i][k];
			for (int j = 0; j < c; ++j)
			{
				ret[i][j] += sum * m[k][j];
			}
		}
	}

	*this = ret;
	return *this;
}

// overload ==
template <typename _T>
bool MATRIX<_T>::operator==(const MATRIX<_T>& rhs) const
{
	if (rows() != rhs.rows() || cols() != rhs.cols())
	{
		return false;
	}

	for (int i = 0; i < rows(); ++i)
	{
		for (int j = 0; j < cols(); ++j)
		{
			if (rhs[i][j] != array[i][j])
			{
				return false;
			}
		}
	}
	return true;
}

// overload !=
template <typename _T>
bool MATRIX<_T>::operator!=(const MATRIX<_T>& rhs) const
{
	return !(*this == rhs);
}

// overload +
template <typename _T>
const MATRIX<_T> MATRIX<_T>::operator+(const MATRIX<_T>& rhs) const
{
	Matrix m;
	if (rows() != rhs.rows() || cols() != rhs.cols())
	{
		if (rhs.isempty()) {
			return *this;
		}
		throw exception("MATRIX::operator+(MATRIX &)::dimension not match");
	}
	m = *this;
	m += rhs;
	return m;
}

// overload -
template <typename _T>
const MATRIX<_T> MATRIX<_T>::operator-(const MATRIX<_T>& rhs) const
{
	Matrix m;
	if (rows() != rhs.rows() || cols() != rhs.cols())
	{
		if (rhs.isempty()) {
			return *this;
		}
		throw exception("MATRIX::operator-=(MATRIX &)::dimension not match");
	}
	m = *this;
	m -= rhs;
	return m;
}

// overload *
template <typename _T>
const MATRIX<_T> MATRIX<_T>::operator*(const MATRIX<_T>& rhs) const
{
	Matrix m;
	if (cols() != rhs.rows())
	{
		throw exception("MATRIX::operator*(MATRIX &)::dimension not match");
	}

	m.resize(rows(), rhs.cols());

	int r = m.rows();
	int c = m.cols();
	int K = cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			double sum = 0.0;
			for (int k = 0; k < K; ++k)
			{
				sum += array[i][k] * rhs[k][j];
			}
			m[i][j] = sum;
		}
	}
	return m;
}

// overload | (joint two matrices)
template <typename _T>
const MATRIX<_T> MATRIX<_T>::operator|(const MATRIX<_T>& rhs) const {
	MATRIX<_T> ret;
	if (rows() != rhs.rows()) { throw exception("MATRIX::operator|::dimension not match"); }
	ret = *this;
	ret.add_block_by_col(rhs);
	return ret;
}


//////////////////////////////////////////////////////////
// Double Matrix class
//////////////////////////////////////////////////////////
class Matrix :public MATRIX<double>
{
public:
	Matrix() :MATRIX<double>() {}
	Matrix(int c, int r) :MATRIX<double>(c, r) {}
	Matrix(const Matrix& m) { *this = m; }
	Matrix(const MATRIX<double>& m) : MATRIX<double>(m.rows(), m.cols(), m.data_col_major()) {}
	Matrix(int rows, int cols, double* a) : MATRIX<double>(rows, cols, a) {}
	Matrix(int rows, int cols, vector<double> a) : MATRIX<double>(rows, cols, a) {};

	// matrix division (inverse a matrix)
	const Matrix& operator/=(const Matrix& m);
	// is zero matrix or not
	bool iszeros();
	// count number of zeros
	int nnz() const;
	// is positive semi-definit
	bool ispsd();
	// is positive definit
	bool ispd();
};

// Double Vector type define (Matrix can not multiply vector but only matrix)
typedef vector<double> Vector;

/***** Enumerations *****/

// Matrix norm types
enum normtype {
	NORM_1,
	NORM_2,
	NORM_INF,
	NORM_MAX,
	NORM_FRO
};

// Random matrix type
enum randomtype {
	RANDOM_REAL,
	RANDOM_INTEGER,
	RANDOM_BINARY
};

/***** Double matrix operations and functions *****/

// overload /
const Matrix operator/(const Matrix& lhs, const Matrix& rhs);
// overload r*m 
const Matrix operator*(const double r, const Matrix& m);
// overload m*r
const Matrix operator*(const Matrix& m, const double r);
// determinant (algorithm 1: LU)
const double det(const Matrix& m);
// determinant (algorithm 2: definition)
const double det2(const Matrix& m);
// 计算子矩阵行列式
const double det(const Matrix& m, int start, int end);
// 计算所有元素的绝对值
const Matrix abs(const Matrix& m);
// 所有元素的最大值
const double max_m(const Matrix& m);
// 所有元素中的最大值及其下标
const double max_m(const Matrix& m, int& row, int& col);
// 所有元素的最小值
const double min_m(const Matrix& m);
// 所有元素的最小值及其下标
const double min_m(const Matrix& m, int& row, int& col);
// 返回子矩阵
const Matrix submatrix(const Matrix& m, int rb, int re, int cb, int ce);
// 删除矩阵某行和某列的元素
const Matrix remove_rowandcol(const Matrix& m, int r, int c);
// 提取矩阵块
const Matrix extractmatrix(const Matrix& m, vector<int> ridx, vector<int> cidx);
// inverse matrix
const Matrix inverse(const Matrix& m);
// LU decomposition
const Matrix LU(const Matrix& m);
// sum (as in matlab)
const Matrix sum(const Matrix& m);
// matrix norm
const double norm(const Matrix& m, enum normtype type);
const double norm_1(const Matrix& m); // norm_1
const double norm_2(const Matrix& m); // norm_2
const double norm_max(const Matrix& m); // norm_max
const double norm_inf(const Matrix& m); // norm_inf
const double norm_fro(const Matrix& m); // norm_fro
// transpose of matrix
const Matrix trans(const Matrix& m);
// trace of matrix
const double trace(const Matrix& m);
// return diagonal elements of matrix or form a diagnal matrix from vector
const Matrix diag(const Matrix& A);
// print matrix
void printMatrix(const Matrix& m, const char* msg, const char* format);
// rang of matrix
int rang(const Matrix& m);
// condition number
double cond(const Matrix& m);
// matrix power
Matrix mpower(const Matrix& m, int p);
// floor of matrix
Matrix floor(const Matrix& m);
// ciel of matrix
Matrix ceil(const Matrix& m);
// round of matrix
Matrix round(const Matrix& m);

/***** Advanced routine with lapack *****/

// SVD A = U*S*VT
void SVD(const Matrix& A, Matrix& S, Matrix& U, Matrix& VT);
// Compute eigenvalues V(Re and Im part) and left and right eigenvecteurs LE, RE
void EIG(const Matrix& A, Matrix& Vr, Matrix& Vi, Matrix& REr, Matrix& REi, Matrix& LEr, Matrix& LEi);
// Compute eigenvalues V and eigenvecteurs E for real and symmetric matrix
void eig_symmetric(const Matrix& A, Matrix& V, Matrix& E);
// Compute eigenvalues V and left and right eigenvecteurs LE, RE for nonsymmetric matrix
void eig_nonsymmetric(const Matrix& A, Matrix& Vr, Matrix& Vi, Matrix& REr, Matrix& REi, Matrix& LEr, Matrix& LEi);
// Spectral decomposition of matrix
void spdecomp(const Matrix& A, Matrix& P, Matrix& N);

/***** Sepcial matrix *****/

// identity matrix
const Matrix eye(int n);
// zero matrix
const Matrix zeros(int row, int col);
// one matrix
const Matrix ones(int row, int col);
// generate a matrix from matlab matrix expression
const Matrix genMatlabMatrix(const char* str, int rows, int cols);
//generate random matrix
const Matrix random(enum randomtype type, int m, int n, double range_min, double range_max);


/***** Matrix formats conversion *****/

// Special structures for supporting CPX
struct CPXMAT {
	int* matbeg; // 1st
	int* matcnt; // 2nd
	int* matind; // 3rd
	double* matval; //4th
};

struct CPXQCMAT {
	int quadnzcnt; // 1st
	int* quadrow; // 2nd
	int* quadcol; // 3rd
	double* quadval; // 4th
};

struct CPXVEC {
	int linnzcnt; // 1st
	int* linind; // 2nd
	double* linval; // 3rd
};

typedef CPXMAT cpxmat;
typedef CPXQCMAT cpxqcmat;
typedef CPXVEC cpxvec;

void format_cpxmat(const Matrix& m, cpxmat& cpx);
void format_cpxqcmat(const Matrix& m, cpxqcmat& cpx);
void format_cpxvec(const Matrix& m, cpxvec& cpx);

// reset random seed
void ResetRandomSeed();

#endif
