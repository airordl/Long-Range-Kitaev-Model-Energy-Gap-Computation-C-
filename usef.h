#ifndef USEF_H
#define USEF_H
#include"INC.h"

namespace USEF {
	using namespace std;
	
	
	template <class x,class y> x  minabs (x a, y b)    	{return (a>(-a)?a:-a) < (b>(-b)?b:-b) ? a : b;}
	template <class x,class y> x  maxabs (x a, y b)    	{return (a>(-a)?a:-a) < (b>(-b)?b:-b) ? b : a;}

	template <typename X> void Swap (X &x,X &y)		{X t=x;x=y;y=t;}
	template <class Y> void Plot (Y** p,long d){
		cout<<endl;
		for(long i=0;i<d;++i){
			for(long j=0;j<d;++j){
				cout<<setw(8)<<setprecision(10)<<fixed<<p[i][j]<<"  ";
			}
			cout<<endl;
		}
	}

	constexpr unsigned int str2int(const char* str, int h = 0) {
		return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
	}


}
#endif
