#ifndef laplace2_H
#define laplace2_H
//二维位势方程

#include <iostream>
#include "MAT.h"

typedef double (*fun2)(double x,double y);

template<class T>
class Laplace2D {
	private:
		int N;
		int n;//线性方程组维数
		T x_start;
		T x_end;
		T y_start;
		T y_end;

		fun2 f;//右端函数
		fun2 g;//u前系数
		fun2 bound;//边界条件

		T h;

	public:
		vvector<double> Exact;
		vvector<double> sol;


		Laplace2D(int N,T x0,T x1,T y0,T y1,fun2 f,fun2 g,fun2 bound):N(N),n((N-1)*(N-1)),
	        x_start(x0),x_end(x1),y_start(y0),y_end(y1),f(f),g(g),bound(bound){
			h=(x_end-x_start)/N;
			Exact.resize(n);
			sol.resize(n);
		}
		
		T xcor(int i){
			return x_start+h*i;
		}

		T ycor(int j){
			return y_start+h*j;
		}

	private:
		void CreateMat(MAT<T>& A){
			for(int j=0;j!=N-1;++j){
				for(int i=0;i!=N-1;++i){
					A[i+(N-1)*j][i+(N-1)*j]=4+h*h*g(xcor(i+1),ycor(j+1));
					if(i-1>=0) A[i+(N-1)*j][i-1+(N-1)*j]=-1;
					if(i+1<=N-2) A[i+(N-1)*j][i+1+(N-1)*j]=-1;
					if(j-1>=0) A[i+(N-1)*j][i+(N-1)*(j-1)]=-1;
					if(j+1<=N-2) A[i+(N-1)*j][i+(N-1)*(j+1)]=-1;
				}
			}
		}

		void CreateConst(vvector<T>& b)//生成二维常数项
		{
			for(int j=0;j!=N-1;++j){
				for(int i=0;i!=N-1;++i){
					b[i+(N-1)*j]=h*h*f(xcor(i+1),ycor(j+1));
					if(i-1<0) b[i+(N-1)*j]+=bound(xcor(0),ycor(j+1));
					if(i+1>N-2) b[i+(N-1)*j]+=bound(xcor(N),ycor(j+1));
					if(j-1<0) b[i+(N-1)*j]+=bound(xcor(i+1),ycor(0));
					if(j+1>N-2) b[i+(N-1)*j]+=bound(xcor(i+1),ycor(N)); 
				}
			}
		}

		void GSarrmultiply(MAT<T>& A,vvector<T>& x_now,vvector<T>& b){
			for(int j=0;j!=N-1;++j){
				for(int i=0;i!=N-1;++i){
					x_now[i+(N-1)*j]=0;
					if(i-1>=0) x_now[i+(N-1)*j]+=x_now[i-1+(N-1)*j]*A[i+(N-1)*j][i-1+(N-1)*j];
					if(i+1<=N-2) x_now[i+(N-1)*j]+=x_now[i+1+(N-1)*j]*A[i+(N-1)*j][i+1+(N-1)*j];
					if(j-1>=0) x_now[i+(N-1)*j]+=x_now[i+(N-1)*(j-1)]*A[i+(N-1)*j][i+(N-1)*(j-1)];
					if(j+1<=N-2) x_now[i+(N-1)*j]+=x_now[i+(N-1)*(j+1)]*A[i+(N-1)*j][i+(N-1)*(j+1)];

					x_now[i+(N-1)*j]+=b[i+(N-1)*j];
				}
			}
		}

		void SORarrmultiply(MAT<T>& A,vvector<T>& x_now,T w,vvector<T>& b){
			for(int j=0;j!=N-1;++j){
				for(int i=0;i!=N-1;++i){
					x_now[i+(N-1)*j]*=(1-w);
					if(i-1>=0) x_now[i+(N-1)*j]+=w*x_now[i-1+(N-1)*j]*A[i+(N-1)*j][i-1+(N-1)*j];
					if(i+1<=N-2) x_now[i+(N-1)*j]+=w*x_now[i+1+(N-1)*j]*A[i+(N-1)*j][i+1+(N-1)*j];
					if(j-1>=0) x_now[i+(N-1)*j]+=w*x_now[i+(N-1)*(j-1)]*A[i+(N-1)*j][i+(N-1)*(j-1)];
					if(j+1<=N-2) x_now[i+(N-1)*j]+=w*x_now[i+(N-1)*(j+1)]*A[i+(N-1)*j][i+(N-1)*(j+1)];

					x_now[i+(N-1)*j]+=w*b[i+(N-1)*j];
				}
			}
		}

		vvector<T> arrmultiply(MAT<T>& A,vvector<T>& x_now){
			vvector<T> temp((N-1)*(N-1));
			for(int j=0;j!=N-1;++j){
				for(int i=0;i!=N-1;++i){
					temp[i+(N-1)*j]=x_now[i+(N-1)*j]*A[i+(N-1)*j][i+(N-1)*j];
					if(i-1>=0) temp[i+(N-1)*j]+=x_now[i-1+(N-1)*j]*A[i+(N-1)*j][i-1+(N-1)*j];
					if(i+1<=N-2) temp[i+(N-1)*j]+=x_now[i+1+(N-1)*j]*A[i+(N-1)*j][i+1+(N-1)*j];
					if(j-1>=0) temp[i+(N-1)*j]+=x_now[i+(N-1)*(j-1)]*A[i+(N-1)*j][i+(N-1)*(j-1)];
					if(j+1<=N-2) temp[i+(N-1)*j]+=x_now[i+(N-1)*(j+1)]*A[i+(N-1)*j][i+(N-1)*(j+1)];
				}
			}
			return temp;
		}

	public:
		vvector<T> Laplace2DSolveGS(T delta,int& count){
			MAT<T> A(n,n,0);
			vvector<T> b(n,0);

			CreateMat(A);
			CreateConst(b);

			for(int j=0;j!=N-1;++j){//A变为Jacobi迭代矩阵B，b变为D^{-1}b.
				for(int i=0;i!=N-1;++i){
					b[i+(N-1)*j]/=A[i+(N-1)*j][i+(N-1)*j];
					if(i-1>=0) A[i+(N-1)*j][i-1+(N-1)*j]/=-A[i+(N-1)*j][i+(N-1)*j];
					if(i+1<=N-2) A[i+(N-1)*j][i+1+(N-1)*j]/=-A[i+(N-1)*j][i+(N-1)*j];
					if(j-1>=0) A[i+(N-1)*j][i+(N-1)*(j-1)]/=-A[i+(N-1)*j][i+(N-1)*j];
					if(j+1<=N-2) A[i+(N-1)*j][i+(N-1)*(j+1)]/=-A[i+(N-1)*j][i+(N-1)*j];
				}
			}

			vvector<T> x_now(n,1);
			vvector<T> x_last(n,0);
			count=0;
			while((x_now-x_last).Norm2()>delta){
				count++;
				x_last=x_now;
				GSarrmultiply(A,x_now,b);
			}

			return x_now;
		}

		vvector<T> Laplace2DSolverSOR(T delta,T w,int& count){
			MAT<T> A(n,n,0);
			vvector<T> b(n,0);

			CreateMat(A);
			CreateConst(b);

			for(int j=0;j!=N-1;++j){//A变为Jacobi迭代矩阵B，b变为D^{-1}b.
				for(int i=0;i!=N-1;++i){
					b[i+(N-1)*j]/=A[i+(N-1)*j][i+(N-1)*j];
					if(i-1>=0) A[i+(N-1)*j][i-1+(N-1)*j]/=-A[i+(N-1)*j][i+(N-1)*j];
					if(i+1<=N-2) A[i+(N-1)*j][i+1+(N-1)*j]/=-A[i+(N-1)*j][i+(N-1)*j];
					if(j-1>=0) A[i+(N-1)*j][i+(N-1)*(j-1)]/=-A[i+(N-1)*j][i+(N-1)*j];
					if(j+1<=N-2) A[i+(N-1)*j][i+(N-1)*(j+1)]/=-A[i+(N-1)*j][i+(N-1)*j];
				}
			}

			vvector<T> x_now(n,1);
			vvector<T> x_last(n,0);
			count=0;
			while((x_now-x_last).Norm2()>delta){
				count++;
				x_last=x_now;
				SORarrmultiply(A,x_now,w,b);
			}
					
			return x_now;
		}

		vvector<T> Laplace2DSolverCG(T delta,int K,int& count){
			MAT<T> A(n,n,0);
			vvector<T> b(n,0);

			CreateMat(A);
			CreateConst(b);

			vvector<T> x(b.size(),1);
			int k=0;
			T rho1,alpha,beta;

			vvector<T> r(b-arrmultiply(A,x));

			vvector<T> p(b.size(),0);
			vvector<T> w(b.size(),0);
			T rho(r.arrmultiply(r));
			while(sqrt(rho) > delta*b.Norm2() && k < K){
				k++;
				if(k==1) p=r;
				else { beta=rho/rho1; p=r+beta*p; }
				w=arrmultiply(A,p);
				alpha=rho/(p.arrmultiply(w));
				x+=alpha*p;
				r-=alpha*w;
				rho1=rho;
				rho=r.arrmultiply(r);
			}
					
			count=k;
			return x;
		}

		template<class Type> friend std::ostream& operator<<(std::ostream& out,const Laplace2D& R);
};

template<class T>std::ostream&
operator<<(std::ostream& out,const Laplace2D<T>& R)
{
	//out<<(R.sol-R.Exact).Norm2()<<" ";
	//out<<R.Sol<<std::endl;
	return out;
}
#endif
