#include"matx.cc"
#include"kitaev.h"

using namespace USEF;

void BuildHpp (size_t, size_t, double, double, double, double, double, Matrix <complex<double>>);
void BuildHop (size_t, size_t, double, double, double, double, double, Matrix <complex<double>>);
void BuildHpo (size_t, size_t, double, double, double, double, double, Matrix <complex<double>>);
void BuildHoo (size_t, size_t, double, double, double, double, double, Matrix <complex<double>>); //da debuggare



/*NB:	una volta che la matrice è stata diagonalizzata con Jacobi, il suo doppio puntatore diventa inservibile*/
/*g++ -std=c++11 *.h matx.cc prova3.cc*/
int main () {
	
double  NOFsteps = 30;
double	alpha_max = 5,
	alpha_min = .2,
	alpha_step = (alpha_max - alpha_min)/ NOFsteps, //square matrix phase space
    	beta  = 100.0,
    	delta = .5,
    	t = .5,
    	mu_max = 5,
    	mu_min = -1,
    	mu_step = (mu_max - mu_min) / NOFsteps,
    	* En,
    	GAP ;	//= pow(10,4);
int	L = 30,
	Lx = L,
	Ly = L,
	Nsiti = Lx*Ly,
	NOFedge = 0;
double	NOFiter = NOFsteps*NOFsteps;
	cout << "\n L = (?) ", cin>> L;
	Lx = L;
	Ly = L;
	Nsiti = L*L;
	
	 	//edge state memory, remembers the allocation of each ES, this has been done to avoid sorting
system	("touch gap.dat");
fstream PSfile ("gap.dat");

	int error 	= 12;
	int errortot	=  0;
	En = new double [Nsiti*2];
	KitaevPoint* PhaseSpace = new KitaevPoint [(const int) NOFiter];
	
	int PScounter = 0;
	
	Matrix <complex<double>> Hpp(2*Nsiti);
	
	cout << setprecision (5) << fixed << "\n'mu' step = " << mu_step << "\t'alpha' step = " << alpha_step << "\tNOFsteps = " << NOFsteps << "\tNOFiter = " << NOFiter;
	
	// terminal plot
	//break;/*suppressed*/
		for (double A = alpha_max-alpha_step;	A > alpha_min-alpha_step;	A -= alpha_step) {
		for (double M = mu_max-mu_step;	M > mu_min-mu_step;	M -= mu_step) {
			Hpp.ALLNULL();
			BuildHpp(Lx, Ly, t, delta, M, A, beta, Hpp);
			error 	= Hpp.EigenProblem();
			
			//buffering PP
			for (long i=0;i<2*Nsiti;++i) {
				for (long j=i+1;j<2*Nsiti;++j) {
					En [i] = real(Hpp.Eval[i]);
				} 
			}
			//computing gap
			for (long i=0;i<2*Nsiti;++i) {
				if (En[i] > 0) {
					GAP = En[i]<GAP?En[i]:GAP;
				}
			}
			
			if (PScounter > NOFiter-1) break;
			errortot-=error;
			PhaseSpace[PScounter].Gap	= GAP;
			PhaseSpace[PScounter].alpha	= A;
			PhaseSpace[PScounter].mu	= M;
			PhaseSpace[PScounter].Plot();
			cout << "\terror = " <<error ;
			cout << setprecision(1) << fixed << "\tPROGRESS " << 100*(PScounter)/NOFiter << "%";
			
			PScounter ++;
			GAP = pow(10,4);	
			
		}
		
		cout << "\n'mu' step = " << mu_step << "\t'alpha' step = " << alpha_step << "\tNOFsteps = " << NOFsteps << "\tNOFiter = " << NOFiter;
	}
	
	int auxPL=0;
	for (int i=NOFiter-1;i>-1;--i) {
		if ( auxPL == NOFsteps ) {PSfile << '\n';auxPL=0;}
		PSfile << PhaseSpace[i].Gap << ' ';
		auxPL++;
	}
	cout << '\n';
	return 0;//suppressing edge state plot
	
	system	("touch evePP.dat");
fstream gappedfile ("evePP.dat");
	system	("touch evePO.dat");
fstream gaplesfile ("evePO.dat");
	double mm = 1, aa = 3;
	Matrix <complex <double>> HPP (2*Nsiti);
	Matrix <complex <double>> HPO (2*Nsiti);
	BuildHpp(Lx, Ly, t, delta, mm, aa, beta, HPP);
	BuildHpo(Lx, Ly, t, delta, mm, aa, beta, HPO);
	HPP.EigenProblem();
	HPO.EigenProblem();
	
	//buffering PP eigenvalues
	for (long i=0;i<2*Nsiti;++i) {
		for (long j=i+1;j<2*Nsiti;++j) {
			En [i] = real(HPP.Eval[i]);
		} 
	}
	

	//computing gap
	GAP = pow(10,4);
	for (long i=0;i<2*Nsiti;++i) {
		if (En[i] > 0) {
			GAP = En[i]<GAP?En[i]:GAP;
		}
	}
	
	//buffering PO eigenvalues
	for (long i=0;i<2*Nsiti;++i) {
		for (long j=i+1;j<2*Nsiti;++j) {
			En [i] = real(HPO.Eval[i]);
		} 
	}
	
	long EdgeIndex;
	double minES = pow(10,4);
	// looking for the smallest positive eigenvalue and buffering its index
	for (long i=0;i<2*Nsiti;++i) {
		if (En[i] > 0) {
			if (En[i] < GAP- (GAP/20.)) { /*to avoid ambiguities*/
				if (En[i] < minES) {
					EdgeIndex = i;
					minES = En[i];
				}
			}
		}
	}
	
	
	Matrix <double> MP (Lx*4);
	double* auxEV = new double [2*Nsiti];
	//buffering PP random eigenvector (10 th) // may cause dump!!
	for (long i = 0;i<2*Nsiti;++i) {
		auxEV[i] = abs(HPP.Evec[i][10]);
		cout << ' '<<auxEV[i];
	}
	
	
	
	// building Evector 'matrix'
	for (int i=0;i<2*Nsiti;++i) {
		MP.matrix[i%Lx][(int) i/Lx] = auxEV[i];
	}
	//writing on file
	for (int i=0;i<Lx;++i) {
		for (int j=0;j<Lx;++j) {
			gappedfile << MP.matrix[i][j] << ' ';
		}	gappedfile << '\n';
	}
	
	// repeating for PO
	
	for (long i = 0;i<2*Nsiti;++i) {
		auxEV[i] = abs(HPO.Evec[i][EdgeIndex]);
		cout << ' '<<auxEV[i];
	}
	for (int i=0;i<2*Nsiti;++i) {
		MP.matrix[i%Lx][(int) i/Lx] = auxEV[i];
	}
	for (int i=0;i<Lx;++i) {
		for (int j=0;j<Lx;++j) {
			gaplesfile << MP.matrix[i][j] << ' ';
		}	gaplesfile << '\n';
	}
	
	
	
	gappedfile.close();
	gaplesfile.close();
	
	PSfile.close();
	
	
	cout << "\nerrortot = " <<errortot;
	
	
	
cout << "\n\a";
return 0;
}//fine main
















void BuildHpp (size_t Lx, size_t Ly, double t, double delta, double mu, double alpha, double beta, Matrix <complex<double>> H) {
	long Nsiti = Lx * Ly;
	for (long j = 0; j < Nsiti; j++) {
		H.matrix[j][j] = (complex<double>) (-mu+4.*t);
		H.matrix[j+Nsiti][j+Nsiti] =(complex<double>) (mu -4.*t);
		long mj = j % Lx;
		long nj = j / Lx;
		for (long rk = 0; rk < Lx; rk++)
			for (long sk = 0; sk < Ly; sk++) {
				if ( (rk == 0 ) && (sk == 0) ) continue;
				long mk = (mj + rk) % Lx;
				long nk = (nj + sk) % Ly;
				long k  = mk + nk * Lx;
				double r = rk;
				if (rk > Lx - rk) r = -((double) Lx) + rk;
				double s = sk;
				if (sk > Ly - sk) s = -((double) Ly) + sk;
				double d = sqrt (r * r + s * s);
				double tt = - t   * exp(-beta * log (d));
				double dd = delta * exp (-alpha * log (d));
				H.matrix [j][k] = tt;
				H.matrix [j+Nsiti][k+Nsiti] = -tt;
				H.matrix [j][k+Nsiti]=dd/d*complex<double>(r,s);
				H.matrix [j+Nsiti][k]=-dd/d*complex<double>(r,-s);
		}
    	}
    	return;
}





void BuildHpo (size_t Lx, size_t Ly, double t, double delta, double mu, double alpha, double beta, Matrix <complex<double>> H){
    	long Nsiti = Lx * Ly;
	for (long j = 0; j < Nsiti; j++) {
		H.matrix[j][j] = (complex<double>) (-mu+4.*t);
		H.matrix[j+Nsiti][j+Nsiti] =(complex<double>) (mu -4.*t);
		long mj = j % Lx;
		long nj = j / Lx;
		for (long rk = 0; rk < Lx; rk++)
			for (long nk = 0; nk < Ly; nk++) {
				if ( (rk == 0) && (nk == nj) ) continue;
				long mk = (mj + rk) % Lx;
				long k  = mk + nk * Lx;
				double r = rk;
				if (rk > Lx - rk) r = - ((double) Lx) + r;
				double s = nk - nj;
				double d = sqrt (r * r + s * s);
				double tt = - t   * exp(-beta * log (d));
				double dd = delta * exp (-alpha * log (d));
				H.matrix [j][k] = tt;
				H.matrix [j+Nsiti][k+Nsiti] = -tt;
				H.matrix [j][k+Nsiti]=dd/d*complex<double>(r,s);
				H.matrix [j+Nsiti][k]=-dd/d*complex<double>(r,-s);

		}
   	}
   	return;
}







void BuildHop (size_t Lx, size_t Ly, double t, double delta, double mu, double alpha, double beta, Matrix <complex<double>> H){
	long Nsiti = Lx * Ly;
	for (long j = 0; j < Nsiti; j++) {
		H.matrix[j][j] = (complex<double>) (-mu+4.*t);
		H.matrix[j+Nsiti][j+Nsiti] =(complex<double>) (mu -4.*t);
		long mj = j % Lx;
		long nj = j / Lx;
		for (long mk = 0; mk < Lx; mk++)
		    for (long sk = 0; sk < Ly; sk++) {
			if ( (mk == mj) && (sk == 0) ) continue;
			long nk = (nj + sk) % Ly;
			long k  = mk + nk * Lx;
			double r = mk - mj;
			double s = sk;
			if (sk > Ly) s = -((double) Ly) + s;
			double d = sqrt (r * r + s * s);
			double tt = - t   * exp(-beta * log (d));
			double dd = delta * exp (-alpha * log (d));
			H.matrix [j][k] = tt;
			H.matrix [j+Nsiti][k+Nsiti] = -tt;
			H.matrix [j][k+Nsiti]=dd/d*complex<double>(r,s);
			H.matrix [j+Nsiti][k]=-dd/d*complex<double>(r,-s);
	    	
	    	}
    	}
    	return;
}






void BuildHoo (size_t Lx, size_t Ly, double t, double delta, double mu, double alpha, double beta, Matrix <complex<double>> H){ //non può funzionare neanche Hpp
	//il debug dovrebbe essere un semplice copia-incolla degli altri pezzi off-diagonal
	//appena ho tempo ci provo

    	long Nsiti = Lx * Ly;
    	for (long j = 0; j < Nsiti; j++) {	//secondo me i segni sono sbagliati
		H.matrix[j][j] = (complex<double>) (-mu+2);
		H.matrix[j+Nsiti][j+Nsiti] =(complex<double>) (mu -2);
		long mj = j % Lx; 
		long nj = j / Lx;
		for (long k = 0; k < j; k++) {
			long mk = k % Lx;
			long nk = k / Lx;
			double r = mk - mj; 
			double s = nk - nj;
			double d = sqrt (r * r + s * s);
			double tt = - t   * exp(-beta * log (d));
			double dd = delta * exp (-alpha * log (d));
			
			/*secondo me ci va 
			
				H.matrix [j][k] = tt;
				H.matrix [j+Nsiti][k+Nsiti] = -tt;
				H.matrix [j][k+Nsiti]=dd/d*complex<double>(r,s);
				H.matrix [j+Nsiti][k]=-dd/d*complex<double>(r,-s);
			
			*/
			
			/*H.matrix [j][k] = tt;
		  	H.matrix [k][j] = tt;
		  	H.matrix [j+Nsiti][k+Nsiti] = -tt;
			H.matrix [k+Nsiti][j+Nsiti] = -tt;
			H.matrix [j][k+Nsiti]=dd/d*complex<double>(r,s);
			H.matrix [k+Nsiti][j]=-dd/d*complex<double>(r,-s);
			H.matrix [j+Nsiti][k]=-dd/d*complex<double>(r,-s);
			H.matrix [k][j+Nsiti]=dd/d*complex<double>(r,s);*/
		}
		return;
    	}	
}




