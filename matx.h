#ifndef MATX_H
#define MATX_H
#include"usef.h"

using namespace std;


template <class T> class Matrix
{
public: 

bool Hbool = true; 	// is it hermitian?

long Dim; 
double Gap;

T** matrix;		//hermitian matrix buffer

T* Eval;		//eigenvalue buffer
T** Evec;		//eigenvector buffer

T DET; 			//determinant buffer

long ERROR = 0;		//errors buffer
const string Algorithm = "DividiEtImpera";
string Name;

Matrix (){};
Matrix (long);
Matrix (T**&,long);
Matrix (long,long,long);

~Matrix() = default;


T Determinant();
Matrix Transposed ();
Matrix CConjugate ();
Matrix HConjugate ();
Matrix Inverse ();

long EigenProblem();//{ERROR = Jacobi(); return ERROR;}

void MPlot 	() 	{USEF::Plot (matrix,Dim);}
void EvecPlot 	() 	{USEF::Plot (Evec,Dim);}
void EvalPlot 	();
void MPrint	()	{MPlot();}

bool HermitianTest ();
void ALLNULL ();


Matrix operator + (Matrix&);
Matrix operator - (Matrix&);
Matrix operator * (Matrix);
const Matrix& operator = (const Matrix &); 

template <typename G> 	friend Matrix <G> 	operator *  (G,Matrix <G>);
template <typename X> 	friend Matrix 		operator >> (istream& inp,Matrix A);	// undebuggable
template <class X> 	friend ostream& 	operator << (ostream &, Matrix <X>);	// undebuggable


protected:



long Jacobi ();
void DividiEtImpera();


};
#endif
