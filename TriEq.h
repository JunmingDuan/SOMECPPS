#ifndef TRIEQ_H
#define TRIEQ_H

#include "tridiag.h"

template<class T>
class TriEq : public tridiag<T>{
	public:
		vvector<T> d;
	private:
		typedef tridiag<T> Base;
		
		using Base::n;
	public:
		TriEq(int N,vvector<T>& a1,vvector<T>& b1,vvector<T>& c1,vvector<T> d1):Base(N,a1,b1,c1),d(d1){ }

		void zhuiganfa(vvector<T> &a1,vvector<T> &b1,vvector<T> &c1,vvector<T> &u)
		{
			u[0]=d[0];
			for(int i=1;i!=n;++i)
			{
				u[i]=d[i]-b1[i-1]*u[i-1];
			}

			u[n-1]=u[n-1]/a1[n-1];
			for(int i=n-2;i!=-1;--i)
			{
				u[i]=(u[i]-u[i+1]*c1[i])/a1[i];
			}
		}

		void TriEqsolve(vvector<T>& x){//zhui gan fa
			vvector<T> a1(n,0);
			vvector<T> b1(n-1,0);
			vvector<T> c1(n-1,0);
			Base::tridiagfenjie(a1,b1,c1);
			zhuiganfa(a1,b1,c1,x);
		}

		vvector<T> TriEqSolveJacobi(T delta,int& count){
			Base::a[n-2]/=-Base::b[n-1];
			Base::c[0]/=-Base::b[0];
			for(int i=1;i!=n-1;++i){
				Base::a[i-1]/=-Base::b[i];
				Base::c[i]/=-Base::b[i];
			}
			for(int i=0;i!=n;++i){
				d[i]/=Base::b[i];
			}

			vvector<T> x_now(n,0.5);
			vvector<T> x_last(n,0);
			count=0;
			while((x_now-x_last).Norm2()>delta){
				count++;
				x_last=x_now;

				x_now[0]=Base::c[0]*x_last[1]+d[0];
				x_now[n-1]=Base::a[n-2]*x_last[n-2]+d[n-1];
				for(int i=1;i!=n-1;++i){
					x_now[i]=Base::a[i-1]*x_last[i-1]+Base::c[i]*x_last[i+1]+d[i];
				}
			}
					
			return x_now;
		}

		vvector<T> TriEqSolveGS(T delta,int& count){
			Base::a[n-2]/=-Base::b[n-1];
			Base::c[0]/=-Base::b[0];
			for(int i=1;i!=n-1;++i){
				Base::a[i-1]/=-Base::b[i];
				Base::c[i]/=-Base::b[i];
			}
			for(int i=0;i!=n;++i){
				d[i]/=Base::b[i];
			}

			vvector<T> x_now(n,0.5);
			vvector<T> x_last(n,0);
			count=0;
			while((x_now-x_last).Norm2()>delta){
				count++;
				x_last=x_now;

				x_now[0]=Base::c[0]*x_now[1]+d[0];
				for(int i=1;i!=n-1;++i){
					x_now[i]=Base::a[i-1]*x_now[i-1]+Base::c[i]*x_now[i+1]+d[i];
				}
				x_now[n-1]=Base::a[n-2]*x_now[n-2]+d[n-1];
			}
					
			return x_now;
		}

		vvector<T> TriEqSolveSOR(T delta,T w,int& count){
			Base::a[n-2]/=-Base::b[n-1];
			Base::c[0]/=-Base::b[0];
			for(int i=1;i!=n-1;++i){
				Base::a[i-1]/=-Base::b[i];
				Base::c[i]/=-Base::b[i];
			}
			for(int i=0;i!=n;++i){
				d[i]/=Base::b[i];
			}

			vvector<T> x_now(n,0.5);
			vvector<T> x_last(n,0);
			count=0;
			while((x_now-x_last).Norm2()>delta){
				count++;
				x_last=x_now;

				x_now[0]=(1-w)*x_now[0]+w*(Base::c[0]*x_now[1]+d[0]);
				for(int i=1;i!=n-1;++i){
					x_now[i]=(1-w)*x_now[i]+w*(Base::a[i-1]*x_now[i-1]+Base::c[i]*x_now[i+1]+d[i]);
				}
				x_now[n-1]=(1-w)*x_now[n-1]+w*(Base::a[n-2]*x_now[n-2]+d[n-1]);
			}
					
			return x_now;
		}

};

#endif //TRIEQU_H
