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

typedef struct SpraseMatrixNode SpraseMatrixNode;

typedef struct SpraseMatrixNode {
    int i, j;               // index of submatrix
    double *submatrix;      // submatrix
    SpraseMatrixNode *right, *down;   // next none zore submatrix;
}SpraseMatrixNode;

typedef struct SpraseMatrix {
    int mu, nu;             // size of sprase matrix
    int sm, sn;             // size of submatrix
    SpraseMatrixNode *rhead, *chead;  // row head and column head of sprase matrix
}SpraseMatrix;

void PrintSpraseMatrix(SpraseMatrix* SMatrix);
int CreateSpraseMatrix(SpraseMatrix* SMatrix, int m, int n, int sm, int sn);
void EmptySpraseMatrix(SpraseMatrix* SMatrix);
void DestroySpraseMatrix(SpraseMatrix* SMatrix);
int CopySubSpraseMatrix(SpraseMatrix* SMatrix, int i, int j, int m, int n, double SubMat[]);
int GetSubSpraseMatrix(SpraseMatrix* SMatrix, int i, int j, int m, int n, double SubMat[]);
int SpraseMatrixAddition(SpraseMatrix* A, SpraseMatrix* B, SpraseMatrix* C);
int SpraseMatrixAddition1(SpraseMatrix* A, SpraseMatrix* B);
int SpraseMatrixMultiply(SpraseMatrix* A, SpraseMatrix* B, SpraseMatrix* C);
int SpraseMatrixMultiply_MPMT(SpraseMatrix* M, SpraseMatrix* P, SpraseMatrix* C);
double SpraseVectDot(SpraseMatrix* A, SpraseMatrix* B);
#endif
