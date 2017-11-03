#ifndef tridiag_H
#define tridiag_H
//追赶法解三对角方程
#include "vvector.h"

template<class T>
class tridiag{
	protected:
		int n;//维数
		vvector<T> a;//靠下对角,n-1
		vvector<T> b;//主对角,n
		vvector<T> c;//靠上对角,n-1

	public:
		tridiag(int dim,vvector<T> &a1,vvector<T> &b1,vvector<T> &c1) : n(dim), a(a1), b(b1), c(c1){ }

		int dim(){
			return n;
		}

		void tridiagfenjie(vvector<T> &a1,vvector<T> &b1,vvector<T> &c1)//b1是分解后下三角矩阵靠下对角，主对角为1,a1是上三角矩阵的主对角，c1是靠上对角
		{
			a1[0]=b[0],c1[0]=c[0];
			for(int ind=0;ind!=n-1;++ind)
			{
				c1[ind]=c[ind];
				b1[ind]=a[ind]/a1[ind];
				a1[ind+1]=b[ind+1]-b1[ind]*c1[ind];
			}
		}

};
#endif
