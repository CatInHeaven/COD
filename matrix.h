#ifndef _MATRIX_H
#define _MATRIX_H
void Rotation_x(double Angle, double Mat[]);
void Rotation_y(double Angle, double Mat[]);
void Rotation_z(double Angle, double Mat[]);

void MatrixMultiply(int m1, int n1, int m2, int n2, const double M1[], const double M2[], double M3[]);
void MatrixAddition(int m, int n, const double M1[], const double M2[], double M3[]);
void MatrixAddition2(int m, int n, const double M1[], double M2[]);
void MatrixTranspose(int m, int n, const double M1[], double MT[]);
double VectDot(int m, int n, const double A[], const double B[]);
void CopyArray(int n, double Dist[], const double Sour[]);
int CopySubMatrix(int Frow, int Fcol, int Brow, int Bcol, int row, int col, double FMat[], double SubMat[]);
void MatrixMultiply_MPMT(int m, int n, const double M[], const double P[], double M2[]);
int GetSubMatrix(int Frow, int Fcol, int Brow, int Bcol, int row, int col, double FMat[], double SubMat[]);
int MatrixInv(int n, double a[], double b[]);
#endif
