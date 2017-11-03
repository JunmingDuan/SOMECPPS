#ifndef NLEq_H
#define NLEq_H
//非线性方程组

#include "MAT.h"
#include "LinEq.h"

template<class T>
class NLEq{
	public:
		typedef vvector<T> Sol;
		typedef vvector<T> (*func)(vvector<T>&);
		typedef MAT<T> (*H_func)(vvector<T>&);
	private:
		func f;
		Sol x0;
		int Max_Ite;
		double e1;
		double e2;
		double e3;
	public:
		NLEq(Sol &x0, func f): x0(x0), f(f), Max_Ite(50), e1(1e-8), e2(1e-8), e3(1e-8){ } 

		Sol FixedPoint(func f1){
			int ite(0);
			Sol SOL1(x0.size());
			Sol SOL2(x0);
			Sol exa(1);
			exa[0]=3.076421163795e+00;
			//Sol exa(3);
			//exa[0]=0,exa[1]=1./3,exa[2]=0;
			do{
				SOL1 = SOL2;
				SOL2 = f1(SOL2);
				ite++;
				std::cout<<(SOL2-exa).Norm2()<<std::endl;
				if (ite >= Max_Ite) break;
			} while((SOL1-SOL2).Norm2()>e1+e2*SOL2.Norm2() || f(SOL2).Norm2()>e3);

			std::cout<<"不动点,ite:"<<ite<<std::endl;
			return SOL2;
		}

		Sol Newton(H_func H){
			int ite(0);
			Sol SOL1(x0.size());
			Sol SOL2(x0);
			Sol exa(1);
			exa[0]=3.076421163795e+00;
			//Sol exa(3);
			//exa[0]=0,exa[1]=1./3,exa[2]=0;
			do{
				SOL1 = SOL2;
				LinEq<double> sub(H(SOL2), f(SOL2));
				sub.LinEqsolvePG();
				SOL2 -= sub.b;
				ite++;
				std::cout<<(SOL2-exa).Norm2()<<std::endl;
				if (ite >= Max_Ite) break;
			} while((SOL1-SOL2).Norm2()>e1+e2*SOL2.Norm2() || f(SOL2).Norm2()>e3);
		//} while((SOL2-exa).Norm2()>2.22e-16);

			std::cout<<"Newton,ite:"<<ite<<std::endl;
			return SOL2;
		}

		Sol Broyden(MAT<double> &H0){//参数是初始Jacobi的逆矩阵
			int ite(0);
			Sol SOL1(x0);
			Sol SOL2(x0 - H0.RArr(f(x0)));
			Sol temp1;
			Sol temp2;
			MAT<double> A(H0);
			Sol g;
			Sol y;
			//Sol exa(2);
			//exa[0]=0,exa[1]=1;
			do{
				y = SOL2 - SOL1;
				g = f(SOL2) - f(SOL1);
				SOL1 = SOL2;
				temp1 = A.RArr(g)-y;
				temp2 = A.LArr(y);
				A -= ArrRArrT(temp1, temp2)/(g.LArr(temp2));
				SOL2 = SOL2 - A.RArr(f(SOL2));
				//std::cout<<(SOL2-exa).Norm2()<<std::endl;
				ite++;
				if (ite >= Max_Ite) break;
			} while((SOL1-SOL2).Norm2()>e1+e2*SOL2.Norm2() || f(SOL2).Norm2()>e3);
		//} while((SOL2-exa).Norm2()>2.22e-16);

			std::cout<<"Broyden,ite:"<<ite<<std::endl;
			return SOL2;
		}

};

#endif
