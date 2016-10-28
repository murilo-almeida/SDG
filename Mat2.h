# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;
// **********************  Classe Mat1 ***********************
// OBSERVACAO: Nunca devem ser passadas por valor;
// Devem ser passadas sempre por referencia!!!
// ***********************************************************
template<typename T> class Mat1 
{
public:
	Mat1();
  Mat1(int I);
  ~Mat1();
  void set (int I);
  T& operator[](int idx){return contents[idx];};
  const T& operator[](int idx) const {return contents[idx];};
	T * _contents() {return contents;};
  const int size() const {return i0;};
  
private:
  int i0;
  T * contents;
};

template<typename T>  Mat1<T>::Mat1()
{
  i0=0;
};

template<typename T>  Mat1<T>::Mat1(int I)
{
	i0=0;
	set(I);
};

template<typename T> void Mat1<T>::set(int I)
{
	if(i0==0 && I > 0) {
			i0=I;
			contents = new T [i0];
		}
	else {
		cout << "Mat1: contents não necessitou ser iniciado\n";
		exit(0);
	}
};

template<typename T>  Mat1<T>::~Mat1()
{
  //cout<<"Inicio Terminou Mat1\n";
  if(i0>0) {
		delete [] contents; contents = nullptr;
		i0=0;
	}
  //cout<<"Destruiu Mat1\n";
};

// **********************  Classe Mat2 ***********************
// OBSERVACAO: Nunca devem ser passadas por valor;
// Devem ser passadas sempre por referencia!!!
// ***********************************************************
template<typename T> class Mat2 
{
public:
	Mat2();
  Mat2(int I,int J);
  ~Mat2();
  void set (int I,int J);
  
  Mat1<T>& operator[](int idx){return contents[idx];};
  const Mat1<T>& operator[](int idx) const {return contents[idx];};
  Mat1<T> * _contents() {return contents;}
  
private:
  int i0;
  Mat1<T> * contents;
};

template<typename T>  Mat2<T>::Mat2()
{
  i0=0;
};

template<typename T>  Mat2<T>::Mat2(int I, int J)
{
	i0=0;
	set(I,J);
};

template<typename T> void Mat2<T>::set(int I,int J)
{
	if(i0 == 0  && I > 0) {
		i0=I;
		contents = new Mat1<T> [i0];
		for(int i=0;i<i0;++i) contents[i].set(J);
	}
	else {
		cout << "Mat2 já está iniciada\n";
		exit(0);
	}
};

template<typename T>  Mat2<T>::~Mat2()
{
	//cout<<"Inicio Terminou Mat2\n";
	if(i0>0) {
		//cout<<"Dentro de Mat2\n";
		for(int i=0;i<i0;++i){
			contents[i].~Mat1();
		}
		delete [] contents; contents = nullptr;
		i0=0;
	}
  //cout<<"Destruiu Mat2\n";
};

// **********************  Classe Mat3 ***********************
// OBSERVACAO: Nunca devem ser passadas por valor;
// Devem ser passadas sempre por referencia!!!
// ***********************************************************
template<typename T> class Mat3 
{
 public:
	Mat3(){i0=0;};
  Mat3(int I,int J,int K);
  ~Mat3();
  void set (int I,int J,int K);
  
  Mat2<T>& operator[](int idx){return contents[idx];};
  const Mat2<T>& operator[](int idx) const {return contents[idx];};
  Mat2<T> * _contents() {return contents;}
  
 private:
  int i0;
  Mat2<T> * contents;
};

template<typename T>  Mat3<T>::Mat3(int I, int J, int K)
{
	i0=0;
	set(I,J,K);
		
};

template<typename T> void Mat3<T>::set(int I,int J,int K)
{
	if(i0 == 0 && (I*J*K) > 0) {
		i0=I;
		contents = new Mat2<T>  [i0];
		for(int i=0; i<i0;++i) contents[i].set(J,K);
	}
	else {
		cout << "Mat3 já está iniciada\n";
		exit(0);
	}
		
};

template<typename T>  Mat3<T>::~Mat3()
{
	//cout<<"Inicio Terminou Mat3\n";
	if(i0>0) {
		//cout<<"Dentro de Mat3\n";
		for(int i=0;i<i0;++i){
			contents[i].~Mat2();
		}
		delete [] contents; contents = nullptr;
		i0=0;
	}
  //cout<<"Terminou Mat3\n";
};
