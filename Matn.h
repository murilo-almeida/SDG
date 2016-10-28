# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;


// **********************  Classe Matn ***********************
// OBSERVACAO: Nunca devem ser passadas por valor;
// Devem ser passadas sempre por referencia!!!
// ***********************************************************
template<class T, int N> class Matn
{
public:
    Matn();
    Matn(int I);
    ~Matn();
    void set (int I);
    Matn<T,N-1>& operator[](int idx){return contents[idx];};
    const Matn<T,N-1>& operator[](int idx) const {return contents[idx];};
  
private:
    int i0;
    Matn<T,N-1> * contents;
};

template<class T, int N>  Matn<T,N>::Matn()
{
    i0=0;
};

template<class T, int N>  Matn<T,N>::Matn(int I)
{
  assert(I>0);
  i0=I;
  contents = new Matn<T,N-1> [i0];

};

template<class T, int N> void Matn<T,N>::set(int I)
{
  assert(I>0);
  if(i0==0) {
    i0=I;
		contents = new Matn<T,N-1> [i0];
	}
	else {
		cout << "Tentativa de reiniciar Matn<T," << N << "> já iniciado\n";
		exit(1);
	}
};

template<class T, int N>  Matn<T,N>::~Matn()
{
 // cout<<"Destructor Matn<T,"<<N<<">\n";
    if(i0>0){
			for(int i=0;i<i0;++i){
				contents[i].~Matn();
			}
			delete [] contents; contents = nullptr;
			i0=0;
		}
 // cout<<"Terminou Matn<T,"<<N<<">\n";
};

// **********************  Classe Matn Specialization ********************
template<class T>
class Matn<T,1>
{
public:
    Matn();
    Matn(int I);
    ~Matn();
    void set (int I);
    T& operator[](int idx){return contents[idx];};
    const T& operator[](int idx) const {return contents[idx];};
    
private:
    int i0;
    T * contents;
};

template<class T>  Matn<T,1>::Matn()
{
    i0=0;
};

template<class T>  Matn<T,1>::Matn(int I)
{
  assert(I>0);
  i0=I;
  contents = new T [i0];
};

template<class T> void Matn<T,1>::set(int I)
{
  assert(I>0);
	if(i0==0) {
    i0=I;
    contents = new T [i0];
	}
	else {
		cout << "Tentativa de reiniciar Matn<T,1> já iniciada\n";
		exit(1);
	}
};

template<class T>  Matn<T,1>::~Matn()
{
   // cout<<"Destructor Matn<T,1>\n";
    if(i0>0) {
			delete [] contents; contents = nullptr;
			i0=0;
		}
   // cout<<"Terminou Matn<T,1>\n";
};
