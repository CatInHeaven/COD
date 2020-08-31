#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
void Rotation_x(double Angle, double Mat[9])
{
	double C, S;

	C = cos(Angle);
	S = sin(Angle);

	*(Mat + 0 * 3 + 0) = 1.0;
	*(Mat + 0 * 3 + 1) = 0.0;
	*(Mat + 0 * 3 + 2) = 0.0;

	*(Mat + 1 * 3 + 0) = 0.0;
	*(Mat + 1 * 3 + 1) = +C;
	*(Mat + 1 * 3 + 2) = +S;

	*(Mat + 2 * 3 + 0) = 0.0;
	*(Mat + 2 * 3 + 1) = -S;
	*(Mat + 2 * 3 + 2) = +C;
}

void Rotation_y(double Angle, double Mat[9])
{
	double C, S;

	C = cos(Angle);
	S = sin(Angle);

	*(Mat + 0 * 3 + 0) = +C;
	*(Mat + 0 * 3 + 1) = 0.0;
	*(Mat + 0 * 3 + 2) = -S;

	*(Mat + 1 * 3 + 0) = 0.0;
	*(Mat + 1 * 3 + 1) = 1.0;
	*(Mat + 1 * 3 + 2) = 0.0;

	*(Mat + 2 * 3 + 0) = +S;
	*(Mat + 2 * 3 + 1) = 0.0;
	*(Mat + 2 * 3 + 2) = +C;
}

void Rotation_z(double Angle, double Mat[9])
{
	double C, S;

	C = cos(Angle);
	S = sin(Angle);

	*(Mat + 0 * 3 + 0) = +C;
	*(Mat + 0 * 3 + 1) = +S;
	*(Mat + 0 * 3 + 2) = 0.0;

	*(Mat + 1 * 3 + 0) = -S;
	*(Mat + 1 * 3 + 1) = +C;
	*(Mat + 1 * 3 + 2) = 0.0;

	*(Mat + 2 * 3 + 0) = 0.0;
	*(Mat + 2 * 3 + 1) = 0.0;
	*(Mat + 2 * 3 + 2) = 1.0;
}

void MatrixMultiply(int m1, int n1, int m2, int n2,
	const double M1[], const double M2[], double M3[])
{
	int i, j, k;
	double Sum;

	if (n1 != m2)
	{
		printf("MatrixMultiply fail: n1!=m2. \n");
		return;
	}

	for (i = 0; i < m1; i++)
	{
		for (j = 0; j < n2; j++)
		{
			Sum = 0.0;

			for (k = 0; k < n1; k++)
			{
				Sum += *(M1 + i * n1 + k) * *(M2 + k * n2 + j);
			}

			*(M3 + i * n2 + j) = Sum;
		}
	}
}

void MatrixAddition(int m, int n, const double M1[], const double M2[], double M3[])
{
	int i, j;

	if (m < 1 || n < 1)
	{
		printf("MatrixAddition fail: m<1 or n<1. \n");
		return;
	}

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			*(M3 + i * n + j) = *(M1 + i * n + j) + *(M2 + i * n + j);
		}
	}
}

void MatrixAddition2(int m, int n, const double M1[], double M2[])
{
	int i, j;

	if (m < 1 || n < 1)
	{
		printf("MatrixAddition2 fail: m<1 or n<1. \n");
		return;
	}

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			*(M2 + i * n + j) = *(M1 + i * n + j) + *(M2 + i * n + j);
		}
	}
}

void MatrixTranspose(int m, int n, const double M1[], double MT[])
{
	int i, j;

	if (m < 1 || n < 1)
	{
		printf("MatrixTranspose fail: m<1 or n<1. \n");
		return;
	}

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			*(MT + j * m + i) = *(M1 + i * n + j);
		}
	}
}

double VectDot(int m, int n, const double A[], const double B[])
{
	int i;
	double Sum;

	if (m != n || m < 1)
	{
		printf("Vector dot fail: m<1 or n<1. \n");
		return 0.0;
	}

	Sum = 0.0;

	for (i = 0; i < m; i++)
	{
		Sum = Sum + *(A + i) * *(B + i);
	}

	return (Sum);
}

void CopyArray(int n, double Dist[], const double Sour[])
{
	int i;

	if (n <= 0)
	{
		printf("CopyArray fail: n<=0. \n");
		return;
	}

	for (i = 0; i < n; i++)
	{
		Dist[i] = Sour[i];
	}
}

int CopySubMatrix(int Frow, int Fcol, int Brow, int Bcol, int row, int col, double FMat[], double SubMat[])
{
	int i, j;

	if (Brow + row > Frow || Bcol + col > Fcol)
	{
		printf("Father matrix don't contain submatrix.\n");
		return -1;
	}

	for (i = Brow; i < Brow + row; i++)
	{
		for (j = Bcol; j < Bcol + col; j++)
		{
			*(FMat + i * Fcol + j) = *(SubMat + (i - Brow) * col + j - Bcol);
		}
	}

	return 0;
}

void MatrixMultiply_MPMT(int m, int n, const double M[], const double P[], double M2[])
{
	int i, j, k;
	double Sum;
	double* Temp;

	if (m < 1 || n < 1)
	{
		printf("MatrixMultiply_MPMT fail: m<1 or n<1. \n");
		return;
	}

	//Temp = new double[m * n];
	Temp = (double*)malloc(m*n*sizeof(double));

	MatrixMultiply(m, n, n, n, M, P, Temp);   // Temp is m1*n1

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			Sum = 0.0;

			for (k = 0; k < n; k++)
			{
				Sum += *(Temp + i * n + k) * *(M + j * n + k);
			}

			*(M2 + i * n + j) = Sum;
		}
	}
	//delete[]Temp;
	free(Temp);
}

int GetSubMatrix(int Frow, int Fcol, int Brow, int Bcol, int row, int col, double FMat[], double SubMat[])
{
	int i, j;

	if (Brow + row > Frow || Bcol + col > Fcol)
	{
		printf("Father matrix don't contain GetSubMatrix.\n");
		return -1;
	}

	for (i = Brow; i < Brow + row; i++)
	{
		for (j = Bcol; j < Bcol + col; j++)
		{
			*(SubMat + (i - Brow) * col + j - Bcol) = *(FMat + i * Fcol + j);
		}
	}

	return 0;
}

int MatrixInv(int n, double a[], double b[])
{
	int i, j, k, l, u, v;   /* matrix dimension <= 5 */
	double d, p;

	int* is = (int*)malloc(sizeof(int) * n);
	int* js = (int*)malloc(sizeof(int) * n);

	for (i = 0; i < n; i++)
	{
		is[i] = 0;
		js[i] = 0;
	}

	if (n <= 0)
	{
		printf("MatrixInv fail because of n<0!\n");
		return 0;
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			b[i * n + j] = a[i * n + j];
		}
	}

	for (k = 0; k < n; k++)
	{
		d = 0.0;
		for (i = k; i < n; i++) 
		{
			for (j = k; j < n; j++)
			{
				l = n * i + j;
				p = fabs(b[l]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		}

		if (d < 0.0)
		{
			printf("MatrixInv fail!\n");
			return 0;
		}

		if (is[k] != k)
		{
			for (j = 0; j < n; j++)
			{
				u = k * n + j;
				v = is[k] * n + j;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}

		if (js[k] != k) 
		{
			for (i = 0; i < n; i++)
			{
				u = i * n + k;
				v = i * n + js[k];
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}

		l = k * n + k;
		b[l] = 1.0 / b[l]; 
		for (j = 0; j < n; j++)
		{
			if (j != k)
			{
				u = k * n + j;
				b[u] = b[u] * b[l];
			}
		}
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				for (j = 0; j < n; j++)
				{
					if (j != k)
					{
						u = i * n + j;
						b[u] = b[u] - b[i * n + k] * b[k * n + j];
					}
				}
			}
		}
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				u = i * n + k;
				b[u] = -b[u] * b[l];
			}
		}
	}

	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
		{
			for (j = 0; j < n; j++)
			{
				u = k * n + j;
				v = js[k] * n + j;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}
		if (is[k] != k)
		{
			for (i = 0; i < n; i++)
			{
				u = i * n + k;
				v = is[k] + i * n;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}
	}
	free(is);
	free(js);
	return 1;
}

void PrintSpraseMatrix(SpraseMatrix* SMatrix) {
	int m, n, sm, sn, i, j, k, l;
	SpraseMatrixNode* p;
	m = SMatrix->mu;
	n = SMatrix->nu;
	sm = SMatrix->sm;
	sn = SMatrix->sn;
	for (i = 0; i < m; i++) {
		for (k = 0; k < sm; k++) {
			p = SMatrix->rhead[i].right;
			for (j = 0; j < n; j++) {
				if (p == NULL || p->j != j) {
					for (l = 0; l < sn; l++) {
						printf("0.0 ");
					}
				}
				else {
					for (l = 0; l < sn; l++) {
						printf("%0.1f ", *(p->submatrix + k * sn + l));
					}
					p = p->right;
				}
			}
			printf("\n");
		}
	}
}

int CreateSpraseMatrix(SpraseMatrix* SMatrix, int m, int n, int sm, int sn) {
	SMatrix->mu = m;
	SMatrix->nu = n;
	SMatrix->sm = sm;
	SMatrix->sn = sn;
	if (!(SMatrix->rhead = (SpraseMatrixNode*)malloc(m * sizeof(SpraseMatrixNode))))
		return -1;
	if (!(SMatrix->chead = (SpraseMatrixNode*)malloc(n * sizeof(SpraseMatrixNode))))
		return -1;
	int i;

	for (i = 0; i < m; i++) {
		SMatrix->rhead[i].j = -1;
		SMatrix->rhead[i].right = NULL;
	}
	for (i = 0; i < n; i++) {
		SMatrix->chead[i].i = -1;
		SMatrix->chead[i].down = NULL;
	}
	return 0;
}

void EmptySpraseMatrix(SpraseMatrix* SMatrix) {
	SpraseMatrixNode* p, * pre;
	int i;
	for (i = 0; i < SMatrix->mu; i++) {
		p = SMatrix->rhead[i].right;
		SMatrix->rhead[i].right = NULL;
		while (p) {
			pre = p;
			p = p->right;
			free(pre->submatrix);
			free(pre);
		}
	}
	for (i = 0; i < SMatrix->nu; i++) {
		SMatrix->chead[i].down = NULL;
	}
}

void DestroySpraseMatrix(SpraseMatrix* SMatrix) {
	SpraseMatrixNode* p, * pre;
	int i;
	for (i = 0; i < SMatrix->mu; i++) {
		p = SMatrix->rhead[i].right;
		while (p) {
			pre = p;
			p = p->right;
			free(pre->submatrix);
			free(pre);
		}
	}
	free(SMatrix->rhead);
	free(SMatrix->chead);
	free(SMatrix);
}

int CopySubSpraseMatrix(SpraseMatrix* SMatrix, int i, int j, int m, int n, double SubMat[]) {

	if (m != SMatrix->sm || n != SMatrix->sn) {
		printf("size of submatrix wrong\n");
		return -1;
	}

	if (i >= SMatrix->mu || j >= SMatrix->nu) {
		printf("Copy submatrix to sprase matrix out of rank\n");
		return -1;
	}

	SpraseMatrixNode* SMatNode = NULL;
	double* sub = NULL;
	int k, l;

	SpraseMatrixNode* p, * pre;
	int flag = 0;

	pre = &SMatrix->rhead[i];
	for (p = pre->right; pre != NULL; p = p->right) {
		if ((p != NULL && p->j > j) || (p == NULL && pre->j < j)) {
			sub = (double*)malloc(sizeof(double) * m * n);
			SMatNode = (SpraseMatrixNode*)malloc(sizeof(SpraseMatrixNode));
			SMatNode->submatrix = sub;
			pre->right = SMatNode;
			SMatNode->right = p;
			break;
		}
		else if (p->j == j) {
			SMatNode = p;
			sub = p->submatrix;
			flag = 1;
			break;
		}
		else {
			pre = p;
		}
	}

	if (!flag) {
		pre = &SMatrix->chead[j];
		for (p = pre->down; p != NULL || pre != NULL; p = p->down) {
			if ((p != NULL && p->i > i) || (p == NULL && pre->i < i)) {
				pre->down = SMatNode;
				SMatNode->down = p;
				break;
			}
			else {
				pre = p;
			}
		}
	}

	SMatNode->i = i;
	SMatNode->j = j;
	for (k = 0; k < m; k++) {
		for (l = 0; l < n; l++) {
			*(sub + k * n + l) = *(SubMat + k * n + l);
		}
	}
	return 0;
}

int GetSubSpraseMatrix(SpraseMatrix* SMatrix, int i, int j, int m, int n, double SubMat[]) {
	if (m != SMatrix->sm || n != SMatrix->sn) {
		printf("size of submatrix wrong\n");
		return -1;
	}
	int k, l;
	if (i >= SMatrix->mu || j >= SMatrix->nu) {
		printf("Copy submatrix to sprase matrix out of rank\n");
		return -1;
	}
	SpraseMatrixNode* p;
	for (p = SMatrix->rhead[i].right; p != NULL; p = p->right) {
		if (p->j == j) {
			for (k = 0; k < m; k++) {
				for (l = 0; l < n; l++) {
					*(SubMat + k * n + l) = *(p->submatrix + k * n + l);
				}
			}
			return 0;
		}
		else if (p->j > j) {
			break;
		}
	}

	for (k = 0; k < m; k++) {
		for (l = 0; l < n; l++) {
			*(SubMat + k * n + l) = 0.0;
		}
	}

	return 0;
}

int SpraseMatrixAddition(SpraseMatrix* A, SpraseMatrix* B, SpraseMatrix* C) {
	int i;
	double* buff;
	SpraseMatrixNode* pa, * pb;
	int sm, sn;
	if (A->mu != B->mu || A->nu != B->nu || A->sm != B->sm || A->sn != B->sn) {
		printf("size of submatrix wrong\n");
		return -1;
	}
	sm = A->sm;
	sn = A->sn;

	buff = (double*)malloc(sizeof(double) * sm * sn);
	for (i = 0; i < A->mu; i++) {
		pa = A->rhead[i].right;
		pb = B->rhead[i].right;
		while (pa != NULL || pb != NULL) {
			if (pa != NULL && pb != NULL) {
				if (pa->j == pb->j) {
					MatrixAddition(sm, sn, pa->submatrix, pb->submatrix, buff);
					CopySubSpraseMatrix(C, i, pa->j, sm, sn, buff);
					pa = pa->right;
					pb = pb->right;
				}
				else if (pa->j < pb->j) {
					CopySubSpraseMatrix(C, i, pa->j, sm, sn, pa->submatrix);
					pa = pa->right;
				}
				else {
					CopySubSpraseMatrix(C, i, pb->j, sm, sn, pb->submatrix);
					pb = pb->right;
				}
			}
			else {
				if (pa != NULL) {
					while (pa) {
						CopySubSpraseMatrix(C, i, pa->j, sm, sn, pa->submatrix);
						pa = pa->right;
					}
				}
				else {
					while (pb) {
						CopySubSpraseMatrix(C, i, pb->j, sm, sn, pb->submatrix);
						pb = pb->right;
					}
				}

			}
		}
	}
	free(buff);
	return 0;
}

int SpraseMatrixAddition1(SpraseMatrix* A, SpraseMatrix* B) {
	// A + B -> B
	int i;
	double* buff;
	SpraseMatrixNode* pa, * pb;
	int sm, sn;
	if (A->mu != B->mu || A->nu != B->nu || A->sm != B->sm || A->sn != B->sn) {
		printf("size of submatrix wrong\n");
		return -1;
	}
	sm = A->sm;
	sn = A->sn;

	buff = (double*)malloc(sizeof(double) * sm * sn);
	for (i = 0; i < A->mu; i++) {
		pa = A->rhead[i].right;
		pb = B->rhead[i].right;
		while (pa != NULL || pb != NULL) {
			if (pa != NULL && pb != NULL) {
				if (pa->j == pb->j) {
					MatrixAddition(sm, sn, pa->submatrix, pb->submatrix, pb->submatrix);
					pa = pa->right;
					pb = pb->right;
				}
				else if (pa->j < pb->j) {
					CopySubSpraseMatrix(B, i, pa->j, sm, sn, pa->submatrix);
					pa = pa->right;
				}
				else {
					pb = pb->right;
				}
			}
			else {
				if (pa != NULL) {
					while (pa) {
						CopySubSpraseMatrix(B, i, pa->j, sm, sn, pa->submatrix);
						pa = pa->right;
					}
				}
				else {
					break;
				}

			}
		}
	}
	free(buff);
}

int SpraseMatrixMultiply(SpraseMatrix* A, SpraseMatrix* B, SpraseMatrix* C) {
	if (A->nu != B->mu || A->sn != B->sm) {
		printf("size of submatrix wrong\n");
		return -1;
	}
	int m, n, sm, sn;
	m = A->mu;
	n = B->nu;
	sm = A->sm;
	sn = B->sn;

	SpraseMatrixNode* pa, * pb;
	int i, j;
	double* buff, * buff_;
	buff = (double*)malloc(sizeof(double) * sm * sn);
	buff_ = (double*)malloc(sizeof(double) * sm * sn);

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			pa = A->rhead[i].right;
			pb = B->chead[j].down;
			memset(buff, 0, sm * sn * sizeof(double));
			while (pa != NULL && pb != NULL) {
				if (pa->j < pb->i) {
					pa = pa->right;
				}
				else if (pa->j > pb->i) {
					pb = pb->down;
				}
				else {
					MatrixMultiply(A->sm, A->sn, B->sm, B->sn, pa->submatrix, pb->submatrix, buff_);
					MatrixAddition(sm, sn, buff, buff_, buff);
					CopySubSpraseMatrix(C, i, j, sm, sn, buff);
					pa = pa->right;
					pb = pb->down;
				}
			}
		}
	}
	free(buff);
	free(buff_);

	return 0;
}

int SpraseMatrixMultiply_MPMT(SpraseMatrix* M, SpraseMatrix* P, SpraseMatrix* C) {
	// M is a diagonal matrix
	if (M->mu != M->nu || M->sm != M->sn) {
		printf("M is not a square matrix\n");
		return -1;
	}
	if (P->mu != P->nu || M->nu != P->mu || P->sm != P->sn || M->sn != P->sm) {
		printf("size of submatrix wrong\n");
		return -1;
	}

	int i, a, b, k;
	SpraseMatrixNode *m1, *p, *m2;
	double *buff, *buff_;
	double sum;
	int m, sm;
	m = M->mu;
	sm = M->sm;
	buff = (double*)malloc(sizeof(double) * sm * sm);
	buff_ = (double*)malloc(sizeof(double) * sm * sm);

	for (i = 0; i < m; i++) {
		p = P->rhead[i].right;
		m1 = M->rhead[i].right;
		if (m1 == NULL) continue;
		if (m1->j != i) {
			printf("M is not a diagonal matrix\n");
			return -1;
		}
		while (p) {
			m2 = M->chead[p->j].down;
			if (m2 != NULL) {
				if (m2->i != p->j) {
					printf("M is not a diagonal matrix\n");
					return -1;
				}
				MatrixMultiply(sm, sm, sm, sm, m1->submatrix, p->submatrix, buff_);
				for (a = 0; a < sm; a++) {
					for (b = 0; b < sm; b++) {
						sum = 0;
						for (k = 0; k < sm; k++) {
							sum += *(buff_ + a * sm + k) * *(m2->submatrix + b * sm + k);
						}
						*(buff + a * sm + b) = sum;
					}
				}
				CopySubSpraseMatrix(C, i, p->j, sm, sm, buff);
				
			}
			p = p->right;
			
		}
	}

	free(buff);
	free(buff_);

	return 0;
}

double SpraseVectDot(SpraseMatrix* A, SpraseMatrix* B) {
	if (A->nu != 1 || B->nu != 1 || A->sn != 1 || B->sn != 1) {
		printf("not a vector\n");
		return 0;
	}
	if (A->mu != B->mu || A->sm != B->sm) {
		printf("size of submatrix wrong\n");
		return 0;
	}
	double sum;
	sum = 0.0;
	SpraseMatrixNode *pa, *pb;
	pa = A->chead[0].down;
	pb = B->chead[0].down;
	while (pa != NULL && pb != NULL) {
		if (pa->i < pb->i) {
			pa = pa->down;
		}
		else if (pa->i > pb->i) {
			pb = pb->down;
		}
		else {
			sum += VectDot(A->sm, B->sm, pa->submatrix, pb->submatrix);
			pa = pa->down;
			pb = pb->down;
		}
	}
	return sum;
}