#ifndef KITAEV_H
#define KITAEV_H
#include"matx.cc"

class KitaevPoint {
	public:
	double alpha;
	double beta;
	double mu;
	double t;
	double delta;
	double Gap;
	double Chern;
	int NOFedge;
	KitaevPoint () {beta = 100.; t = .5; delta = .5;} 
	
	void Print () {
		cout << setprecision(5) << fixed << "\n[a,m,GAP] = " << "[" << alpha << "," << mu << "," << Gap << "]" ;
	}
	void Plot () {Print();}
};
#endif
