#ifndef LinEq_H
#define LinEq_H
//线性方程组

#include "MAT.h"

template<class T>
class LinEq : public MAT<T>{
	private:
		typedef MAT<T> Base;

		typedef vvector<T> Sol;

		typedef vvector<int> SWAP;

		using Base::row;
		using Base::col;
	public:
		vvector<T> b;//常数项

		LinEq(){
		}

		LinEq(const MAT<T> &A1,const vvector<T> &b1):Base(A1),b(b1){ } 
		
		void updiagsolve(){//解上三角方程
			int N=col;
			for(int i=N-1;i!=0;--i){
				if(std::abs((*this)[i][i]-0)<pow(10,-16)) { std::cerr<<"Singular"<<std::endl; }
				b[i]=b[i]/(*this)[i][i];
				for(int j=0;j!=i;++j){
					b[j]-=b[i]*(*this)[j][i];
				}
			}
			if(std::abs((*this)[0][0]-0)<pow(10,-16)) { std::cerr<<"Singular"<<std::endl; }
			b[0]=b[0]/(*this)[0][0];
		}

//		void updiagsolve(Sol& u){//解上三角方程
//			for(int i=N-1;i!=-1;--i){
//				u[i]=b[i];
//				for(int j=N-1;j!=i;--j){
//					u[i]-=u[j]*(*this)[i][j];
//				}
//				u[i]/=(*this)[i][i];
//			}
//		}

		void lowdiagsolve(){//解下三角方程
			int N=col;
			for(int i=0;i!=N-1;++i){
				if(std::abs((*this)[i][i]-0)<pow(10,-14)) { std::cerr<<"Singular"<<std::endl; }
				b[i]=b[i]/(*this)[i][i];
				for(int j=i+1;j!=N;++j){
					b[j]-=b[i]*(*this)[j][i];
				}
			}
			if(std::abs((*this)[N-1][N-1]-0)<pow(10,-14)) { std::cerr<<"Singular"<<std::endl; }
			b[N-1]=b[N-1]/(*this)[N-1][N-1];
		}

		void unitupdiagsolve(){//解单位上三角方程
			int N=col;
			for(int i=N-1;i!=0;--i){
				for(int j=0;j!=i;++j){
					b[j]-=b[i]*(*this)[j][i];
				}
			}
		}

		void unitlowdiagsolve(){//解单位下三角方程
			int N=col;
			for(int i=0;i!=N-1;++i){
				for(int j=i+1;j!=N;++j){
					b[j]-=b[i]*(*this)[j][i];
				}
			}
		}

//		void unitlowdiagsolve(Sol& u){//解单位下三角方程
//			for(int i=0;i!=N;++i){
//				u[i]=b[i];
//				for(int j=0;j!=i;++j){
//					u[i]-=u[j]*(*this)[i][j];
//				}
//			}
//		}


		void ConstRowSwap(std::vector<int> &P)//根据形参完成常数项行交换
		{
			int N=col;
			for(int i=0;i!=N;++i){
				int j=P[i];
				if(j!=-1){
					std::swap(b[i],b[j]);
				}
			}
		}

		void ConstRowSwapTrans(std::vector<int> &P)//根据形参完成常数项行交换
		{
			int N=col;
			for(int i=N-1;i!=-1;--i){
				int j=P[i];
				if(j!=-1){
					std::swap(b[i],b[j]);
				}
			}
		}

		void LinEqsolveG(){//Gauss消去法后解方程
			Base::gauss();

			unitlowdiagsolve();
			updiagsolve();
		}

		void LinEqsolvePG(){//列主元Gauss消去法后解方程
			SWAP a;
			a=Base::colgauss();
			ConstRowSwap(a);

			unitlowdiagsolve();
			updiagsolve();
		}

		void LinEqsolveQR(){
			vvector<T> d;
			d=Base::QR();
			MAT<T> QT=(Base::QR_Q(d)).trans();
			b=QT.arrmultiply(b);

			updiagsolve();
		}

		T LSsolveQR(){//QR求解最小二乘
			vvector<T> d;
			d=Base::QR();
			MAT<T> QT=(Base::QR_Q(d)).trans();
			b=QT.arrmultiply(b);

			updiagsolve();

			T sum(0);
			for(int i=col;i!=row;++i){
				sum+=pow(b[i],2);
			}

			return sqrt(sum);
		}

		T LSsolveRe(){//正则化求解最小二乘
			vvector<T> b0(b);
			MAT<T> A0((*this));
			MAT<T> AT=Base::trans();
			b=AT.arrmultiply(b);
			Base::Regular();

			pdLinEqsolve2();

			return (A0.arrmultiply(b)-b0).Norm2();
		}

		void LinEqsolve0(){//已做列主元Gauss消去法,只解方程
			unitlowdiagsolve();
			updiagsolve();
		}

		void LinEqTranssolve0(){//已做列主元Gauss消去法,只解转置方程
			int N=col;
			for(int i=0;i!=N-1;++i){
				b[i]=b[i]/(*this)[i][i];
				for(int j=i+1;j!=N;++j){
					b[j]-=b[i]*(*this)[i][j];
				}
			}
			b[N-1]=b[N-1]/(*this)[N-1][N-1];

			for(int i=N-1;i!=0;--i){
				for(int j=0;j!=i;++j){
					b[j]-=b[i]*(*this)[i][j];
				}
			}

		}

		bool pdLinEqsolve1(){//对称正定方程用Cholesky法求解
			int N=col;
			if(Base::cholesky()==0){
				std::cerr<<"未完成cholesky分解"<<std::endl;
				return 0;
			}

			lowdiagsolve();
			for(int i=N-1;i!=0;--i){
				b[i]=b[i]/(*this)[i][i];
				for(int j=0;j!=i;++j){
					b[j]-=b[i]*(*this)[i][j];
				}
			}
			b[0]=b[0]/(*this)[0][0];
			return 1;
		}

		bool pdLinEqsolve2(){//对称正定方程用改进的Cholesky法求解
			int N=col;
			if(Base::imcholesky()==0){
				std::cerr<<"未完成改进cholesky分解"<<std::endl;
				return 0;
			}

			unitlowdiagsolve();
			for(int i=N-1;i!=-1;--i){
				b[i]=b[i]/(*this)[i][i];
			}

			for(int i=N-1;i!=0;--i){
				for(int j=0;j!=i;++j){
					b[j]-=b[i]*(*this)[i][j];
				}
			}
			return 1;
		}

		vvector<T> LinEqSolveJacobi(T delta,int& count){
			MAT<T> B(row,col,0);
			for(int i=0;i!=row;++i){
				for(int j=0;j!=col;++j){
					if(i==j){ B[i][j]=0; }
					else { B[i][j]=-(*this)[i][j]/(*this)[i][i]; }
				}
			}
			vvector<T> g(b);
			for(int i=0;i!=g.size();++i){
				g[i]/=(*this)[i][i];
			}

			vvector<T> x_now(g.size(),1);
			vvector<T> x_last(g.size(),0);
			count=0;
			while((x_now-x_last).Norm2()>delta){
				count++;
				x_last=x_now;
				x_now=B.arrmultiply(x_now);
				x_now+=g;
			}
					
			return x_now;
		}

		vvector<T> LinEqSolveGS(T delta,int& count){
			MAT<T> B(row,col,0);
			for(int i=0;i!=row;++i){
				for(int j=0;j!=col;++j){
					if(i==j){ B[i][j]=0; }
					else { B[i][j]=-(*this)[i][j]/(*this)[i][i]; }
				}
			}
			vvector<T> g(b);
			for(int i=0;i!=g.size();++i){
				g[i]/=(*this)[i][i];
			}

			vvector<T> x_now(g.size(),1);
			vvector<T> x_last(g.size(),0);
			count=0;
			while((x_now-x_last).Norm2()>delta){
				count++;
				x_last=x_now;
				for(int i=0;i!=row;++i){
					x_now[i]=B[i].arrmultiply(x_now);
					x_now[i]+=g[i];
				}
			}
					
			return x_now;
		}

		vvector<T> LinEqSolveSOR(T delta,T w,int& count){
			MAT<T> B(row,col,0);
			for(int i=0;i!=row;++i){
				for(int j=0;j!=col;++j){
					if(i==j){ B[i][j]=0; }
					else { B[i][j]=-(*this)[i][j]/(*this)[i][i]; }
				}
			}
			vvector<T> g(b);
			for(int i=0;i!=g.size();++i){
				g[i]/=(*this)[i][i];
			}

			vvector<T> x_now(g.size(),1);
			vvector<T> x_last(g.size(),0);
			count=0;
			while((x_now-x_last).Norm2()>delta){
				count++;
				x_last=x_now;
				for(int i=0;i!=row;++i){
					x_now[i]=(1-w)*x_now[i];
					x_now[i]+=w*B[i].arrmultiply(x_now);//B_{ii}=0!
					x_now[i]+=w*g[i];
				}
			}
					
			return x_now;
		}

		vvector<T> LinEqSolverCG(T delta,int K,int& count){//实用共轭梯度法
			vvector<T> x(b.size(),1);
			int k=0;
			T rho1,alpha,beta;
			vvector<T> r(b-Base::RArr(x));
			vvector<T> p(b.size(),0);
			vvector<T> w(b.size(),0);
			T rho(r.LArr(r));
			while(sqrt(rho) > delta*b.Norm2() && k < K){
				k++;
				if(k==1) p=r;
				else { beta=rho/rho1; p=r+beta*p; }
				w=Base::RArr(p);
				alpha=rho/(p.LArr(w));
				x+=alpha*p;
				r-=alpha*w;
				rho1=rho;
				rho=r.LArr(r);
			}
			count=k;

			return x;
		}

		void LinEqSolverPCG(Sol& x, MAT<T>& M, T delta, const int MAX, int& ite){//预优共轭梯度法,初值,预优矩阵,
			int feva(0), tmp;
			ite = 0;
			T rho1,alpha,beta;
			vvector<T> r(b-Base::RArr(x));
			vvector<T> p(b.size(),0);
			vvector<T> w(b.size(),0);
			vvector<T> z(b.size(),0);
			T rho;
			MAT<T> A(b.size(), b.size());
			for(int i = 0; i != b.size(); ++i){
				for(int j = 0; j != b.size(); ++j){
					A[i][j] = (*this)[i][j];
				}
			}
			for(int i = 0; i != b.size(); ++i){
				for(int j = 0; j != b.size(); ++j){
					(*this)[i][j] = M[i][j];
				}
			}
			while(sqrt(r.LArr(r)) > delta*b.Norm2() && ite < MAX){
				//解子问题
				b = r;
				z = LinEqSolverCG(delta, MAX, tmp);
				feva += tmp;
				ite++;

				//下一步
				if(ite == 1) {
					p = z;
					rho = z.LArr(r);
				}
				else { 
					rho1 = rho;
					rho = z.LArr(r);
					beta = rho/rho1; 
					p = z + beta*p;
				}

				w = A.RArr(p);
				alpha = rho/(w.LArr(p));
				x += alpha*p;
				r -= alpha*w;
				rho1=rho;
			}
		}

		T condInfiniteOptimi(){
			int N=col;
			int times(0);
			T ANorm(Base::NormInfinite());
			vvector<T> x;
			x.assign(N,1./N);
			vvector<T> w(N),z(N);

			SWAP P=Base::colgauss();
			while(1){
				times++;
				b=x;
				LinEqTranssolve0();
				ConstRowSwapTrans(P);
				w=b;

				b=w.sign();
				ConstRowSwap(P);
				LinEqsolve0();
				z=b;

				if(z.NormInfinite()<=w.Norm1()){
					return ANorm*w.Norm1();
				}
				else if(times>pow(10,1)){
					return ANorm*w.Norm1();
				}
				else {
					x.assign(N,0);
					T Norm=z.NormInfinite();
					for(int i=0;i!=z.size();++i){
						if(std::abs(z[i])==Norm) { x[i]=1; break; }
					}
				}
			}

			std::cerr<<"Cannot stop"<<std::endl;
			return 0;
		}

		T ErrorEsti(){
			int N=col;
			MAT<T> A0(*this);
			vvector<T> b0(b);
			T ANorm(A0.NormInfinite());
			T bNorm(b0.NormInfinite());

			SWAP P=Base::colgauss();
			ConstRowSwap(P);
			LinEqsolve0();

			T rNorm=(b0-A0.arrmultiply(b)).NormInfinite();

			vvector<T> x;
			x.assign(N,1./N);
			vvector<T> w(N),z(N);

			while(1){
				b=x;
				ConstRowSwap(P);
				LinEqTranssolve0();
				w=b;

				b=w.sign();
				ConstRowSwap(P);
				LinEqsolve0();
				z=b;

				if(z.NormInfinite()<=z.arrmultiply(x)){
					return rNorm*ANorm*w.Norm1()/bNorm;
				}
				else {
					x.assign(N,0);
					T Norm=z.NormInfinite();
					for(int i=0;i!=z.size();++i){
						if(std::abs(z[i])==Norm) { x[i]=1; break; }
					}
				}
			}

			std::cerr<<"Cannot stop"<<std::endl;
			return 0;
		}

		void printsol(std::ostream& out){//打印解向量
			out.precision(2);
			out<<std::showpos;
			out.setf(std::ios::scientific);

			for(int i=0;i!=b.size();++i){
				out<<i<<"	"<<b[i]<<std::endl;
			}
		}

		template<class type> friend std::ostream& operator<<(std::ostream&,const LinEq<type>&);
};

template<class T>
	std::ostream&
operator<<(std::ostream& out,const LinEq<T>& X)
{
	out<<X.b<<std::endl;
	//	for(int i=0;i!=X.row;++i){
	//		out<<X.b[i]<<"\n";
	//	}
	//	out<<std::endl;
	return out;
}
#endif
