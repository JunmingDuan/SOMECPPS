#ifndef _MVECTOR_H
#define _MVECTOR_H

#include <valarray>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

template<class Type, int Size>
class mvector : public std::valarray<Type>
{
    public:
        enum { size = Size };
	private:
		void resize(size_t n0){ }

        typedef Type value_type;

        typedef std::valarray<Type> Base;
    public:
        mvector() : Base(Size) { }

        mvector(Type val) : Base(val,Size) { }

		template<class C> mvector& operator=(const C& CC) { Base::operator=(CC); return *this; }

		template<class C> mvector& operator+=(const C& CC) { Base::operator+=(CC); return *this; }

		template<class C> mvector& operator-=(const C& CC) { Base::operator-=(CC); return *this; }

		template<class C> mvector& operator*=(const C& CC) { Base::operator*=(CC); return *this; }

		template<class C> mvector& operator/=(const C& CC) { Base::operator/=(CC); return *this; }

		Type Norm1(){
			Type a(0);
			for(int i=0;i!=size;++i){
				a+=std::abs((*this)[i]);
			}
			return a;
		}

		Type Norm2(){
			Type a(0);
			for(int i=0;i!=size;++i){
				a+=pow((*this)[i],2);
			}
			return sqrt(a);
		}

		Type NormInfinite(){
			Type a(0);
			for(int i=0;i!=size;++i){
				if(std::abs((*this)[i])>a){
					a=std::abs((*this)[i]);
				}
			}
			return a;
		}

		template<class T,int _Size> friend std::istream& operator>>(std::istream&, mvector&);

		template<class T,int _Size> friend std::ostream& operator<<(std::ostream&, mvector&);
};

template<class T,int Size1>
mvector<T,Size1> operator+(const mvector<T,Size1>& v1, const mvector<T,Size1>& v2)
{
    mvector<T,Size1> v3(v1);
    v3 += v2;
    return v3;
}

template<class T,int Size1>
mvector<T,Size1> operator-(const mvector<T,Size1>& v1, const mvector<T,Size1>& v2)
{
    mvector<T,Size1> v3(v1);
    v3 -= v2;
    return v3;
}

template<class T1,int Size,class T2>
mvector<T1,Size> operator*(const mvector<T1,Size>& v1, const T2& a)
{
    mvector<T1,Size> v3(v1);
    v3 *= a;
    return v3;
}

template<class T1,int Size,class T2>
mvector<T1,Size> operator*(const T2& a,const mvector<T1,Size>& v1)
{
    mvector<T1,Size> v3(v1);
    v3 *= a;
    return v3;
}

template<class T1,int Size,class T2>
mvector<T1,Size> operator/(const mvector<T1,Size>& v1, const T2& a)
{
	if(a==0){
		std::cerr<<"Cannot be divided by 0."<<std::endl;
		return v1;
	}
    mvector<T1,Size> v3(v1);
    v3 /= a;
    return v3;
}

template<class T,int Size> 
std::istream& operator>>(std::istream& in, mvector<T,Size>& mv){
	std::cout<<"Please input numbers"<<std::endl;
	T value;
	for(int i=0;i!=Size && in>>value;++i){
		mv[i]=value;
	}
	return in;
}

template<class T,int Size>
std::ostream& operator<<(std::ostream& os, const mvector<T,Size>& mv){
	os.precision(6);
	os<<std::showpos;
	os<<std::fixed;
	//os.setf(std::ios::scientific);
	for(int i=0;i!=mv.size;++i){
		os<<mv[i]<<"\t";
	}
	return os;
}

#endif //_mvector_H
