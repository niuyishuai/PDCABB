/*********************************************
Language:	C++
Project:	Matrix class
Author:		Yi-Shuai Niu
Copyright:	All rights reserved
Date:		Since Sept 2017
*********************************************/
#include "stdafx.h"
#include "Matrix.h"

// operator /= for double matrix
const Matrix& Matrix::operator/=(const Matrix& m)
{
	Matrix tmp = inverse(m);
	*this *= tmp;
	return *this;
}

// is zero matrix or not
bool Matrix::iszeros() {
	if (isempty()) { return true; }
	for (auto i = 0; i < rows(); i++) {
		for (auto j = 0; j < cols(); j++) {
			if (fabs(array[i][j]) > MATRIX_ZERO) { return false; }
		}
	}
	return true;
}

// number of zeros
int Matrix::nnz() const {
	int ret = 0;
	for (auto i = 0; i < rows(); i++) {
		for (auto j = 0; j < cols(); j++) {
			if (fabs(array[i][j]) > MATRIX_ZERO) { ret++; }
		}
	}
	return ret;
}

// transpose of matrix
const Matrix trans(const Matrix& m)
{
	return m.trans();
}


const Matrix operator*(const double r, const Matrix& m) {
	Matrix ret(m);

	int row = m.rows();
	int col = m.cols();

	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			ret[i][j] = r*m[i][j];
		}
	}
	return ret;
}

const Matrix operator*(const Matrix& m, const double r) {
	return r*m;
}

const Matrix operator/(const Matrix& lhs, const Matrix& rhs)
{
	Matrix tmp = inverse(rhs);
	Matrix m;

	if (tmp.isempty())
	{
		// should be completed
		return m;
	}
	return m = lhs * tmp;
}

inline static double LxAbs(double d)
{
	return (d >= 0) ? (d) : (-d);
}

inline
static bool isSignRev(const vector<double>& v)
{
	int p = 0;
	int sum = 0;
	int n = (int)v.size();

	for (int i = 0; i < n; ++i)
	{
		p = (int)v[i];
		if (p >= 0)
		{
			sum += p + i;
		}
	}

	if (sum % 2 == 0) // 如果是偶数，说明不变号
	{
		return false;
	}
	return true;
}

// 计算方阵行列式 (通过LU分解)
const double det(const Matrix& m)
{
	double ret = 0.0;

	if (m.isempty() || !m.issquare()) return ret;

	Matrix N = LU(m);

	if (N.isempty()) return ret;

	ret = 1.0;
	for (int i = 0; i < N.cols(); ++i)
	{
		ret *= N[i][i];
	}

	/*
	if (isSignRev(N[N.rows() - 1]))
	{
		return -ret;
	}
	*/

	return ret;
}

// 计算矩阵指定子方阵的行列式 
const double det(const Matrix& m, int start, int end)
{
	return det(submatrix(m, start, end, start, end));
}

inline bool iszero(double a) {
	return (fabs(a) < MATRIX_ZERO);
}

// 计算方阵行列式方法2 根据定义计算
const double det2(const Matrix& m)
{
	if (m.cols() == 1 && m.rows() == 1) return m[0][0];

	double ret = 0;
	for (int i = 0; i < m.cols(); i++) {
		auto sign = pow(-1.0, i);
		if (!iszero(m[0][i])) {
			ret += sign * m[0][i] * det2(remove_rowandcol(m, 0, i));
		}
	}
	return ret;
}

// 计算逆矩阵
const Matrix inverse(const Matrix& m)
{
	Matrix ret;

	if (m.isempty() || !m.issquare())
	{
		return ret;
	}

	int n = m.rows();

	ret.resize(n, n);
	Matrix A(m);

	for (int i = 0; i < n; ++i) ret[i][i] = 1.0;

	for (int j = 0; j < n; ++j)  //每一列
	{
		int p = j;
		double maxV = LxAbs(A[j][j]);
		for (int i = j + 1; i < n; ++i)  // 找到第j列中元素绝对值最大行
		{
			if (maxV < LxAbs(A[i][j]))
			{
				p = i;
				maxV = LxAbs(A[i][j]);
			}
		}

		if (maxV < 1e-20)
		{
			ret.resize(0, 0);
			return ret;
		}

		if (j != p)
		{
			A.swap_row(j, p);
			ret.swap_row(j, p);
		}

		double d = A[j][j];
		for (int i = j; i < n; ++i) A[j][i] /= d;
		for (int i = 0; i < n; ++i) ret[j][i] /= d;

		for (int i = 0; i < n; ++i)
		{
			if (i != j)
			{
				double q = A[i][j];
				for (int k = j; k < n; ++k)
				{
					A[i][k] -= q * A[j][k];
				}
				for (int k = 0; k < n; ++k)
				{
					ret[i][k] -= q * ret[j][k];
				}
			}
		}
	}

	return ret;
}

// 计算绝对值
const Matrix abs(const Matrix& m)
{
	Matrix ret;

	if (m.isempty())
	{
		return ret;
	}

	int r = m.rows();
	int c = m.cols();
	ret.resize(r, c);

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			double t = m[i][j];
			if (t < 0) ret[i][j] = -t;
			else ret[i][j] = t;
		}
	}

	return ret;
}

// 返回矩阵所有元素的最大值
const double max_m(const Matrix& m)
{
	if (m.isempty()) return 0.;

	double ret = m[0][0];
	int r = m.rows();
	int c = m.cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			if (m[i][j] > ret) ret = m[i][j];
		}
	}
	return ret;
}

// 计算矩阵最大值，并返回该元素的引用
const double max_m(const Matrix& m, int& row, int& col)
{
	if (m.isempty()) return 0.;

	double ret = m[0][0];
	row = 0;
	col = 0;

	int r = m.rows();
	int c = m.cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			if (m[i][j] > ret)
			{
				ret = m[i][j];
				row = i;
				col = j;
			}
		}
	}
	return ret;
}

// 计算矩阵所有元素最小值
const double min_m(const Matrix& m)
{
	if (m.isempty()) return 0.;

	double ret = m[0][0];
	int r = m.rows();
	int c = m.cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			if (m[i][j] > ret) ret = m[i][j];
		}
	}

	return ret;
}

// 计算矩阵最小值，并返回该元素的引用
const double min_m(const Matrix& m, int& row, int& col)
{
	if (m.isempty()) return 0.;

	double ret = m[0][0];
	row = 0;
	col = 0;
	int r = m.rows();
	int c = m.cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			if (m[i][j] > ret)
			{
				ret = m[i][j];
				row = i;
				col = j;
			}
		}
	}

	return ret;
}

// 取矩阵中指定位置的子矩阵 
const Matrix submatrix(const Matrix& m, int rb, int re, int cb, int ce)
{
	Matrix ret;
	if (m.isempty()) return ret;

	if (rb < 0 || re >= m.rows() || rb > re) return ret;
	if (cb < 0 || ce >= m.cols() || cb > ce) return ret;

	ret.resize(re - rb + 1, ce - cb + 1);

	for (int i = rb; i <= re; ++i)
	{
		for (int j = cb; j <= ce; ++j)
		{
			ret[i - rb][j - cb] = m[i][j];
		}
	}

	return ret;
}

// 删除某行和某列的元素
const Matrix remove_rowandcol(const Matrix& m, int r, int c) {
	Matrix ret;
	if (m.isempty()) return ret;

	if (r > m.rows() || c > m.cols() || r < 0 || c < 0) return ret;

	ret = submatrix(m, 0, r - 1, 0, c - 1);
	ret.add_block_by_col(submatrix(m, 0, r - 1, c + 1, m.cols() - 1));
	auto B = submatrix(m, r + 1, m.rows() - 1, 0, c - 1);
	B.add_block_by_col(submatrix(m, r + 1, m.rows() - 1, c + 1, m.cols() - 1));
	ret.add_block_by_row(B);

	return ret;
}

inline static
int max_idx(const Matrix& m, int k, int n)
{
	int p = k;
	for (int i = k + 1; i < n; ++i)
	{
		if (LxAbs(m[p][k]) < LxAbs(m[i][k]))
		{
			p = i;
		}
	}
	return p;
}

// 计算方阵 M 的 LU 分解
// 其中L为对角线元素全为1的下三角阵，U为对角元素依赖M的上三角阵
// 使得 M = LU
// 返回矩阵下三角部分存储L(对角元素除外)，上三角部分存储U(包括对角线元素)
const Matrix LU(const Matrix& m)
{
	Matrix ret;

	if (m.isempty() || !m.issquare()) return ret;

	int n = m.rows();
	ret.resize(n + 1, n);

	for (int i = 0; i < n; ++i)
	{
		ret[n][i] = -1.0;
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			ret[i][j] = m[i][j];
		}
	}

	for (int k = 0; k < n - 1; ++k)
	{
		int p = max_idx(ret, k, n);
		if (p != k)              // 进行行交换
		{
			ret.swap_row(k, p);
			ret[n][k] = (double)p; // 记录将换信息
		}

		if (ret[k][k] == 0.0)
		{
			cout << "ERROR: " << endl;
			ret.resize(0, 0);
			return ret;
		}

		for (int i = k + 1; i < n; ++i)
		{
			ret[i][k] /= ret[k][k];
			for (int j = k + 1; j < n; ++j)
			{
				ret[i][j] -= ret[i][k] * ret[k][j];
			}
		}
	}

	return ret;
}

#ifdef max
#undef max
#undef min
#endif // !max

int max(vector<int> v) {
	int maxval = 0;
	for (auto i : v) {
		if (maxval < i) {
			maxval = i;
		}
	}
	return maxval;
}

int min(vector<int> v) {
	int minval = 0;
	for (auto i : v) {
		if (minval > i) {
			minval = i;
		}
	}
	return minval;
}

// 提取矩阵块
const Matrix extractmatrix(const Matrix& m, vector<int> ridx, vector<int> cidx) {
	Matrix ret;
	if (m.isempty()) return ret;

	if (max(ridx) > m.rows() - 1 || min(ridx) < 0 || max(cidx) > m.cols() - 1 || min(cidx) < 0) return ret;

	ret.resize((int)ridx.size(), (int)cidx.size());
	for (auto i = 0; i < (int)ridx.size(); i++) {
		for (auto j = 0; j < (int)cidx.size(); j++) {
			ret[i][j] = m[ridx[i]][cidx[j]];
		}
	}
	return ret;
}

// identity matrix
const Matrix eye(int n) {
	Matrix ret(n, n);
	for (auto i = 0; i < n; i++) {
		ret[i][i] = 1.0;
	}
	return ret;
}

// zero matrix
const Matrix zeros(int row, int col) {
	Matrix ret(row, col);
	return ret;
}

// one matrix
const Matrix ones(int row, int col) {
	Matrix ret(row, col);
	for (auto i = 0; i < row; i++) {
		for (auto j = 0; j < col; j++) {
			ret[i][j] = 1.0;
		}
	}
	return ret;
}

// 计算矩阵范数
const double norm(const Matrix& m, enum normtype type) {
	double ret;
	switch (type)
	{
	case NORM_1:
		ret = norm_1(m);
		break;
	case NORM_2:
		ret = norm_2(m);
		break;
	case NORM_INF:
		ret = norm_inf(m);
		break;
	case NORM_MAX:
		ret = norm_max(m);
		break;
	case NORM_FRO:
		ret = norm_fro(m);
		break;
	default:
		break;
	}
	return ret;
}

// sum of matrix as in matlab
const Matrix sum(Matrix & m) {
	int row = m.rows();
	int col = m.cols();

	if (row == 1) {
		Matrix ret(1, 1);
		ret[0][0] = 0.0;
		for (int j = 0; j < col; j++) {
			ret[0][0] += m[0][j];
		}
		return ret;
	}
	Matrix ret(1, col);
	for (int j = 0; j < col; j++) {
		ret[0][j] = 0.0;
		for (int i = 0; i < row; i++) {
			ret[0][j] += m[i][j];
		}
	}
	return ret;
}

// compute matrix norm 1
const double norm_1(const Matrix& m) {
	auto a = abs(m);
	return max_m(sum(a));
}

// compute matrix norm 2
const double norm_2(const Matrix& m) {
	Matrix s, u, vt;
	SVD(m, s, u, vt);
	return max_m(s);
}

// compute matrix norm max
const double norm_max(const Matrix& m) {
	return max_m(abs(m));
}

// compute matrix norm inf
const double norm_inf(const Matrix& m) {
	return norm_1(trans(m));
}

// compute matrix norm fro
const double norm_fro(const Matrix& m) {
	return sqrt(trace(trans(m)*m));
}

// trace of matrix
const double trace(const Matrix& m) {
	return m.trace();
}

// rang de matrice
int rang(const Matrix& m) {
	int ret = 0;
	Matrix s, u, vt;
	SVD(m, s, u, vt);
	auto tol = fmax(m.rows(), m.cols())*max_m(s)*MATRIX_ZERO; // method get from Matlab rank function
	for (auto i = 0; i < s.cols(); i++) {
		if (fabs(s[0][i]) > tol) {
			ret++;
		}
	}
	return ret;
}

// condition number
double cond(const Matrix& m) {
	if (m.isempty()) { return 0.0; }
	Matrix s, u, vt;
	SVD(m, s, u, vt);
	return s[0][0] / s[0][s.cols() - 1];
}

// is positive semidefinit
bool Matrix::ispsd() {
	if (isempty()) { return false; }
	Matrix s, u, vt;
	SVD(*this, s, u, vt);
	return (s[0][s.cols() - 1] >= 0.0);
}

// is positive definit
bool Matrix::ispd() {
	if (isempty()) { return false; }
	Matrix s, u, vt;
	SVD(*this, s, u, vt);
	return (s[0][s.cols() - 1] > (s[0][0] - s[0][s.cols() - 1])*MATRIX_ZERO);
}

// matrix power
Matrix mpower(const Matrix& m, int p) {
	Matrix ret;
	if (!m.issquare()) { return ret; }
	if (p == 0) { return eye(m.rows()); }
	ret = m;
	for (auto i = 2; i <= p; i++) {
		ret *= m;
	}
	return ret;
}

// return diagonal elements of matrix or form a diagnal matrix from vector 
const Matrix diag(const Matrix& A) {
	return A.diag();
}


// SDV A = U*S*VT
void SVD(const Matrix& A, Matrix& S, Matrix& U, Matrix& VT) {
	auto m = (integer)A.rows();
	auto n = (integer)A.cols();
	auto a = (doublereal*)A.data_col_major();
	auto lda = m;
	auto ldu = m;
	auto ldvt = n;
	integer info;
	integer lwork;

	doublereal wkopt;
	doublereal* work;
	doublereal* s = new doublereal[n];
	doublereal* u = new doublereal[ldu*m];
	doublereal* vt = new doublereal[ldvt*n];

	/* Query and allocate the optimal workspace */
	lwork = -1;
	dgesvd_((char*)"All", (char*)"All", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info);
	lwork = (integer)wkopt;
	work = (doublereal*)malloc(lwork * sizeof(doublereal));
	/* Solve eigenproblem */
	dgesvd_((char*)"All", (char*)"All", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
	S = Matrix(1, n, s);
	U = Matrix(m, n, u);
	VT = Matrix(n, n, vt);

	free((void*)work);
	delete[] s;
	delete[] u;
	delete[] vt;
}

// compute eigenvalues V and left and right eigenvecteurs LE and RE
void EIG(const Matrix& A, Matrix& Vr, Matrix& Vi, Matrix& REr, Matrix& REi, Matrix& LEr, Matrix& LEi) {
	if (A.issymmetric()) {
		eig_symmetric(A, Vr, REr);
	}
	else if (A.rows() == A.cols()) {
		eig_nonsymmetric(A, Vr, Vi, REr, REi, LEr, LEi);
	}
}

// compute eigenvalues V and eigenvecteurs E for real and symmetric matrix
void eig_symmetric(const Matrix& A, Matrix& V, Matrix& E) {
	auto n = (integer)A.rows();
	auto lda = n;
	auto a = (doublereal*)A.data_col_major();
	integer info;
	integer lwork;

	doublereal wkopt;
	doublereal* work;
	doublereal* w = new doublereal[n];

	/* Query and allocate the optimal workspace */
	lwork = -1;
	dsyev_((char*)"Vectors", (char*)"Upper", &n, a, &lda, w, &wkopt, &lwork, &info);
	lwork = (integer)wkopt;
	work = (doublereal*)malloc(lwork * sizeof(doublereal));
	/* Solve eigenproblem */
	dsyev_((char*)"Vectors", (char*)"Upper", &n, a, &lda, w, work, &lwork, &info);
	V = Matrix(1, n, w); // Eigenvalues
	E = Matrix(n, n, a); // Eigenvectors (stored columnwise)

	free((void*)work);
	delete[] w;
}

// extract from lapack complex matrix
void extractcomplexmatrix(Matrix& Mr, Matrix& Mi, int n, double *wi, double* v, int ldv) {
	Mr = Matrix(n, n);
	Mi = Matrix(n, n);
	for (int i = 0; i < n; i++) {
		int j = 0;
		while (j < n) {
			if (wi[j] == (double)0.0) {
				Mr[i][j] = v[i + j*ldv];
				Mi[i][j] = (double)0.0;
				j++;
			}
			else {
				Mr[i][j] = v[i + j*ldv];
				Mi[i][j] = v[i + (j + 1)*ldv];
				Mr[i][j + 1] = Mr[i][j];
				Mi[i][j + 1] = -Mi[i][j];
				j += 2;
			}
		}
	}
}

// compute complex eigenvalues V and left and right eigenvecteurs LE, RE for nonsymmetric matrix
// Vr : real part of eigenvalues, Vi : image part of eigenvalues
// REr: real part of right eigenvectors, REi: image part of right eigenvectors
// LEr: real part of left eigenvectors, LEi: image part of left eigenvectors
void eig_nonsymmetric(const Matrix& A, Matrix& Vr, Matrix& Vi, Matrix& REr, Matrix& REi, Matrix& LEr, Matrix& LEi) {
	auto n = (integer)A.rows();
	auto a = (doublereal*)A.data_col_major();
	auto lda = n;
	auto ldvl = n;
	auto ldvr = n;
	integer info;
	integer lwork;

	doublereal wkopt;
	doublereal* work;
	doublereal* wr = new doublereal[n]; //  for real part eigenvalues
	doublereal* wi = new doublereal[n]; //  for image part eigenvalues
	doublereal* vl = new doublereal[ldvl*n];
	doublereal* vr = new doublereal[ldvr*n];

	/* Query and allocate the optimal workspace */
	lwork = -1;
	dgeev_((char*)"Vectors", (char*)"Vectors", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, &info);
	lwork = (integer)wkopt;
	work = (doublereal*)malloc(lwork * sizeof(doublereal));
	/* Solve eigenproblem */
	dgeev_((char*)"Vectors", (char*)"Vectors", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
	Vr = Matrix(1, n, wr);
	Vi = Matrix(1, n, wi);
	extractcomplexmatrix(REr, REi, n, wi, vr, ldvr);
	extractcomplexmatrix(LEr, LEi, n, wi, vl, ldvl);

	free((void*)work);
	delete[] wr;
	delete[] wi;
	delete[] vl;
	delete[] vr;
}

// Spectral decomposition of matrix
// decompose A into P-N where P and N are psd matrices
void spdecomp(const Matrix& A, Matrix& P, Matrix& N) {
	// only for symmetric matrix
	if (A.issymmetric() == false) return;
	Matrix V, E;
	eig_symmetric(A, V, E);
	P = V.diag();
	N = -1.0*P;
	for (auto i = 0; i < A.rows(); i++) {
		if (P[i][i] < 0) { P[i][i] = 0.0; }
		if (N[i][i] < 0) { N[i][i] = 0.0; }
	}
	auto invE = E.trans(); // E is orthogonal matrix
	P = E*P*invE;
	N = E*N*invE;
	// for some numerical issue, P and N are not tested to be symmetric, we have to symmetrize them
	P = (P + P.trans()) * 0.5;
	N = (N + N.trans())*0.5;
}


// generate a matrix from matlab matrix expression
const Matrix genMatlabMatrix(const char* str, int rows, int cols) {
	Matrix ret(rows, cols);
	string cmd(str);
	auto size = cmd.size();
	auto posbeg = cmd.find('[');
	auto posend = cmd.find(']');
	auto substrs = cmd.substr(posbeg + 1, posend - posbeg - 1);
	auto buff = substrs.data();
	char * stopbuff;

	for (auto i = 0; i < rows; i++) {
		for (auto j = 0; j < cols; j++) {
			while (buff[0] == ' ' || buff[0] == ',' || buff[0] == ';') { buff++; }
			ret[i][j] = strtod(buff, &stopbuff);
			buff = stopbuff;
		}
	}
	return ret;
}

// print a Matrix in console
void printMatrix(const Matrix& m, const char* msg, const char* format) {
	if (msg != NULL) {
		cout << msg << endl;
	}
	if (format == NULL) {
		format = "%13.4lf";
	}
	for (auto i = 0; i < m.rows(); i++) {
		for (auto j = 0; j < m.cols(); j++) {
			printf(format, m[i][j]);
		}
		printf("\n");
	}
}

// floor of matrix
Matrix floor(const Matrix& m) {
	Matrix ret(m);
	for (auto i = 0; i < m.rows(); i++) {
		for (auto j = 0; j < m.cols(); j++) {
			ret[i][j] = floor(m[i][j]);
		}
	}
	return ret;
}

// ciel of matrix
Matrix ceil(const Matrix& m) {
	Matrix ret(m);
	for (auto i = 0; i < m.rows(); i++) {
		for (auto j = 0; j < m.cols(); j++) {
			ret[i][j] = ceil(m[i][j]);
		}
	}
	return ret;
}

// round of matrix
Matrix round(const Matrix& m) {
	Matrix ret(m);
	for (auto i = 0; i < m.rows(); i++) {
		for (auto j = 0; j < m.cols(); j++) {
			ret[i][j] = round(m[i][j]);
		}
	}
	return ret;
}

// convert a matrix into cpxmat format
void format_cpxmat(const Matrix& m, cpxmat& cpx) {
	int i, j;
	int nzcols;

	cpx.matbeg = (int*)calloc(m.cols(), sizeof(int));
	cpx.matcnt = (int*)calloc(m.cols(), sizeof(int));
	auto nz = m.nnz();
	cpx.matind = (int*)calloc(nz, sizeof(int));
	cpx.matval = (double*)calloc(nz, sizeof(double));

	nz = 0;
	for (j = 0; j < m.cols(); j++)
	{
		nzcols = 0;
		cpx.matbeg[j] = nz;
		for (i = 0; i < m.rows(); i++)
		{
			if (fabs(m[i][j]) > MATRIX_ZERO)
			{
				cpx.matval[nz] = m[i][j];
				cpx.matind[nz] = i;
				if (nzcols == 0) cpx.matbeg[j] = nz;
				nz++;
				nzcols++;
			}
		}
		cpx.matcnt[j] = nzcols;
	}
}

// convert a matrix into cpxqcmat format
void format_cpxqcmat(const Matrix& m, cpxqcmat& cpx) {
	int i, j;
	int nz;

	if (!m.issquare()) { return; }

	int size = m.rows();
	cpx.quadnzcnt = m.nnz();
	cpx.quadrow = (int*)calloc(cpx.quadnzcnt, sizeof(int));
	cpx.quadcol = (int*)calloc(cpx.quadnzcnt, sizeof(int));
	cpx.quadval = (double*)calloc(cpx.quadnzcnt, sizeof(double));

	nz = 0;
	for (j = 0; j < size; j++)
	{
		for (i = 0; i < size; i++)
		{
			if (fabs(m[i][j]) > MATRIX_ZERO)
			{
				cpx.quadval[nz] = m[i][j];
				cpx.quadrow[nz] = i;
				cpx.quadcol[nz] = j;
				nz++;
			}
		}
	}
}

// convert a matrix into cpxvec format
void format_cpxvec(const Matrix& m, cpxvec& cpx) {
	int i;
	int nz;

	Matrix v = m.vec();

	cpx.linnzcnt = v.nnz();
	cpx.linind = (int*)calloc(cpx.linnzcnt, sizeof(int));
	cpx.linval = (double*)calloc(cpx.linnzcnt, sizeof(double));

	nz = 0;
	auto cols = m.cols();
	for (i = 0; i < cols; i++)
	{
		if (fabs(v[i][0]) > MATRIX_ZERO)
		{
			cpx.linval[nz] = v[i][0];
			cpx.linind[nz] = i;
			nz++;
		}
	}
}

// return a random real value in [a,b]
inline double randomreal(double a, double b)
{
	return (a + fabs(b - a) * (double)abs(rand()) / RAND_MAX);
}

// return a random integer value in [a,b]
inline int randomint(double a, double b)
{
	return (int)round(a + fabs(b - a) * (double)abs(rand()) / RAND_MAX);
}

//generate random matrix
const Matrix random(enum randomtype type, int m, int n, double range_min, double range_max) {
	Matrix ret;
	if (m < 0 || n < 0) { return ret; }
	ret.resize(m, n);
	for (auto i = 0; i < m; i++) {
		for (auto j = 0; j < n; j++) {
			switch (type)
			{
			case RANDOM_REAL:
				ret[i][j] = randomreal(range_min, range_max);
				break;
			case RANDOM_INTEGER:
				ret[i][j] = (double)randomint(range_min, range_max);
				break;
			case RANDOM_BINARY:
				ret[i][j] = (double)randomint(0, 1);
				break;
			default:
				break;
			}
		}
	}
	return ret;
}

// reset random seed
void ResetRandomSeed() { srand((unsigned)time(0)); }
