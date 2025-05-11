#ifndef MATX_CC
#define MATX_CC
#include"matx.h"

using namespace std;



template <class T> void Matrix <T> :: ALLNULL () {
	for (int i=0;i<Dim;++i) {
		for (int j=0;j<Dim;++j) {
			matrix[i][j] =  (T) 0.;
			Evec  [i][j] =  (T) 0.;
		}
		Eval [i] = 0.;
	}
}




template <class T> long Matrix <T> :: EigenProblem () {
	if (Algorithm == "Jacobi") {
		ERROR = Jacobi();
		return ERROR;
	}
	/*if (Algorithm == "Householder") {
		ERROR = Householder();
		return ERROR;
	}*/	//maybe one day
	
    /*
	if (Algorithm == "DividiEtImpera") {
		DividiEtImpera();
		return 0;
	}
    */
	
	else {
		ERROR = Jacobi();
		return ERROR;
	}
}



template <class T> bool Matrix <T> :: HermitianTest() {
	bool hh=true;
	Matrix <T> AH (Dim);
	
	AH = *this - *this.HConjugate();
	for (int i=0;i<Dim;++i) {
		for (int j=0;j<Dim;++j) {
			if (AH.matrix [i][j] != (T) 0.){
				hh = false;
				break;
			}
			
		}	if (hh == false) break;
	}
	
	Hbool = hh;
	return hh;
}





template <class T> void Matrix <T> :: EvalPlot () {
	for (long i=0;i<Dim;++i) {
		cout << "\n" <<Eval[i];
	}	cout << endl;
}



template <class T> Matrix <T> :: Matrix  (T **&mat,long dim){

	Dim=dim;
	matrix	= new T*[Dim];
	Evec	= new T*[Dim];
	Eval 	= new T [Dim];
	for (long i=0;i<Dim;++i){
		matrix[i]=new T [Dim];
		Evec[i]=new T [Dim];
		Eval[i] = 0.;
		for (long j=0;j<Dim;++j){
			matrix[i][j]=mat[i][j];
			Evec[i][j]=0.;
		}
	}
	
	
}

template <class T> Matrix <T> :: Matrix (long dim){

	Dim=dim;
	matrix	= new T*[Dim];
	Evec	= new T*[Dim];
	Eval 	= new T [Dim];
	for (long i=0;i<Dim;++i){
		matrix[i]=new T [Dim];
		Evec[i]=new T [Dim];
		Eval[i] = 0.;
		for (long j=0;j<Dim;++j){
			matrix[i][j]=0.;
			Evec[i][j]=0.;
		}
	}
}






template <class T> Matrix <T> Matrix <T> :: operator + (Matrix& B){
	Matrix Sum(Dim);
	for (long i=0;i<Sum.Dim;++i){
		for (long j=0;j<Sum.Dim;++j){
			Sum.matrix[i][j]=matrix[i][j]+B.matrix[i][j];
		}
	}
	return Sum;
}


template <class T> Matrix <T> Matrix <T> :: operator - (Matrix& B){
	Matrix Sum(Dim);
	for (long i=0;i<Sum.Dim;++i){
		for (long j=0;j<Sum.Dim;++j){
			Sum.matrix[i][j]=matrix[i][j]-B.matrix[i][j];
		}
	}
	return Sum;
}

template <class T> Matrix <T> Matrix <T> :: operator * (Matrix B){

	Matrix Product(Dim);
	for (long i=0;i<Dim;++i){
		for (long j=0;j<Dim;++j){
			for (long n=0;n<Dim;++n){
				Product.matrix[i][j]+=matrix[i][n]*B.matrix[n][j];
			}
		}
	}
	return Product;
}

template <class T> const Matrix <T> & Matrix <T> :: operator = (const Matrix &B) {
	if (this!=&B) {
		for (long i=0;i<Dim;++i){
			for (long j=0;j<Dim;++j){
				matrix[i][j]=B.matrix[i][j];
			}
		}
	}
	return *this;
}


template <class T> Matrix <T> Matrix <T> :: Transposed () {
	Matrix Tran(Dim);
	for (long i=0;i<Dim;++i){
		for (long j=0;j<Dim;++j){
			Tran.matrix [i][j]=matrix[j][i];
		}
	}
	return Tran;
}

template <class T> Matrix <T> Matrix <T> :: CConjugate () {
	Matrix CC(Dim);
	for (long i=0;i<Dim;++i){
		for (long j=0;j<Dim;++j){
			CC.matrix [i][j]=conj(matrix[i][j]);
		}
	}
	return CC;
}

template <class T> Matrix <T> Matrix <T> :: HConjugate () {
	Matrix HC(Dim);
	HC = *this;
	return HC.Transposed().CConjugate();
}



template <class T> long Matrix <T> :: Jacobi () {
	using namespace USEF;

	T** A 	= new T*[Dim];
	T** Q	= new T*[Dim];
	T*  w	= new T [Dim];
	for (long i=0;i<Dim;++i){
		A [i]=new T [Dim];
		Q [i]=new T [Dim];
		
		for (long j=0;j<Dim;++j){
			A [i][j]= matrix[i][j];
			Q[i][j]=0;
		}
	}
		
	A = matrix;
	Q = Evec;
	w = Eval;
	
	const int MAXNOFiterations = 50;
	
	const long n = Dim;
	double sd, so;                  // Sums of diagonal resp. off-diagonal elements
	complex <double> s, t;          // sin(phi), tan(phi) and temporary storage
	double c;                       // cos(phi)
	double g, h, z;                 // Temporary storage
	double thresh;
	  
	// Initialize Q to the identitity matrix


	#ifndef EVALS_ONLY
	for (long i=0; i < n; i++){
		Q[i][i] = 1.0;
		for (long j=0; j < i; j++)
			Q[i][j] = Q[j][i] = 0.0;
	}
	#endif

	// Initialize w to diag(A)
	for (long i=0; i < n; i++)
		w[i] = real(A[i][i]); // why real tho
	
	// Calculate SQR(tr(A))  
	sd = 0.0;
	for (long i=0; i < n; i++)
		sd += abs(w[i]);
	sd = real(sd*sd);
	
	// Main iteration loop
	for (long nIter=0; nIter < MAXNOFiterations; nIter++){
		// Test for convergence 
		so = 0.0;
		for (long p=0; p < n; p++)
			for (long q=p+1; q < n; q++)
				so += abs(real(A[p][q])) + abs(imag(A[p][q]));
		if (so == 0.0)
			return 0; //When everything goes smoothly
		
		if (nIter < 4)
			thresh = real(0.2 * so / (n*n));
		else
			thresh = 0.0;
		
		// Do sweep
		for (long p=0; p < n; p++)
		for (long q=p+1; q < n; q++){
			g = 100.0 * (abs(real(A[p][q])) + abs(imag(A[p][q])));
			if (nIter > 4	&&	abs(w[p]) + g == abs(w[p])	&&	abs(w[q]) + g == abs(w[q])){
				A[p][q] = 0.0;
			}
			else if (abs(real(A[p][q])) + abs(imag(A[p][q])) > thresh){
				// Calculate Jacobi transformation
				h = real(w[q] - w[p]);
				if (abs(h) + g == abs(h))
					t = A[p][q] / h;
				else{
					if (h < 0.0)
						t = -2.0 * A[p][q] / (sqrt((h*h) + 4.0*abs(A[p][q])*abs(A[p][q]) ) - h);
					else if (h == 0.0)
						t = A[p][q] * (1.0 / abs(A[p][q]));  // A[p][q]/abs(A[p][q]) could cause overflows
					else
						t = 2.0 * A[p][q] / (sqrt((h*h) + 4.0*abs(A[p][q])*abs(A[p][q])) + h);
				}
				c = 1.0/sqrt(1.0 + abs(t)*abs(t));
				s = t * c;
				z = real(t * conj(A[p][q]));

				// Apply Jacobi transformation
				
				A[p][q] = 0.0;
				w[p] -= z;
				w[q] += z;
				for (long r=0; r < p; r++){
					t = A[r][p];
					A[r][p] = c*t - conj(s)*A[r][q];
					A[r][q] = s*t + c*A[r][q];
				}
				for (long r=p+1; r < q; r++){
					t = A[p][r];
					A[p][r] = c*t - s*conj(A[r][q]);
					A[r][q] = s*conj(t) + c*A[r][q];
				}
				for (long r=q+1; r < n; r++){
					t = A[p][r];
					A[p][r] = c*t - s*A[q][r];
					A[q][r] = conj(s)*t + c*A[q][r];
				}

					  // Update eigenvectors
				#ifndef EVALS_ONLY          
				for (long r=0; r < n; r++){
					t = Q[r][p];
					Q[r][p] = c*t - conj(s)*Q[r][q];
					Q[r][q] = s*t + c*Q[r][q];
				}
				#endif
			}
		}
	}

	return -1;//When something goes wrong, if this happens MAXNOFiterations must be increased, 50 is fine so far with L=15
}

	
// EOF Jacobi

/*
template <class T> void Matrix <T> :: DividiEtImpera () { //eigenvalues only
	using namespace arma;
	cx_mat H (Dim,Dim,fill::zeros);
	for (long i=0;i<Dim;++i) {
		for (long j=0;j<Dim;++j) {
			H(i,j) = matrix[i][j];
		}
	}
	
	vec En(Dim, fill::zeros);
	
	#ifndef EVECS_ALSO
	En = eig_sym(H);
	for (long i=0;i<Dim;++i) {
		Eval[i] = En(i);
	}
	#endif
	
	#ifndef EVALS_ONLY
	cx_mat eigvec (Dim, Dim);
	eig_sym(En,eigvec,H);
	#endif
	for (long i=0;i<Dim;++i) {
		for (long j=0;j<Dim;++j) {
			Evec[i][j] = eigvec(i,j);
		}
		Eval[i] = En(i);
	}
}
*/

#endif

