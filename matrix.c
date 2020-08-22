#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
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
	double* Temp;    // ����ʹ�ö�̬���䣬��������Ϊ����������ά��

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

	/* ���������ֵ���������b�������b�������棬a���󲻱� */
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
		for (i = k; i < n; i++)   /* �������½Ƿ�������Ԫ�ص�λ�� */
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

		if (d < 0.0)   /* ��Ԫ�ؽӽ���0�����󲻿��� */
		{
			printf("MatrixInv fail!\n");
			return 0;
		}

		if (is[k] != k)  /* ����Ԫ�����ڵ��������½Ƿ�������н��е��� */
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

		if (js[k] != k)  /* ����Ԫ�����ڵ��������½Ƿ�������н��е��� */
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
		b[l] = 1.0 / b[l];  /* �����б任 */
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

	for (k = n - 1; k >= 0; k--)  /* ����������е������»ָ� */
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