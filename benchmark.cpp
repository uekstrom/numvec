#include <iostream>
#include <ctime>
#include <cassert>
using namespace std;
#include "numvec.hpp"


// The is the formula of the popular Lee-Yang-Parr correlation functional.
template<class T>
void lypc(T &out,
	  const T &a,
	  const T &b,
	  const T &gaa,
	  const T &gnn,
	  const T &gbb)
{
  const double A = 0.04918;
  const double B = 0.132;
  const double C = 0.2533;
  const double Dd = 0.349;
  const double CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);
  T n = a+b;
  T icbrtn = pow(n,-1.0/3.0);
  T P = 1/(1+Dd*icbrtn);
  T omega = exp(-C*icbrtn)*P*pow(n,-11.0/3.0);
  T delta = icbrtn*(C+Dd*P);
  T n2 = n*n;
  out =
    -A*(4*a*b*P/n +
	B*omega*(a*b*(pow(2,11.0/3.0)*CF*(pow(a,8.0/3.0)+pow(b,8.0/3.0))
			  +(47.0 - 7.0*delta)*gnn/18.0
			  -(2.5 - delta/18.0)*(gaa + gbb)
			  -(delta-11.0)/9.0*(a*gaa + b*gbb)/n)
		 - 2.0/3.0*n2*gnn 
		 + (2.0/3.0*n2 - a*a)*gbb
		 + (2.0/3.0*n2 - b*b)*gaa));
}

template<class T>
void lypc_disciplined(T &out,
		      const T &a,
		      const T &b,
		      const T &gaa,
		      const T &gnn,
		      const T &gbb)
{
  const double A = 0.04918;
  const double B = 0.132;
  const double C = 0.2533;
  const double Dd = 0.349;
  const double CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);
  T n = a+b;
  T icbrtn = pow(n,-1.0/3.0);
  T tmp = Dd*icbrtn;
  tmp += 1;
  T P = 1/tmp;
  tmp = Dd*P;
  tmp += C;
  T delta = icbrtn*tmp;
  icbrtn *= -C;
  T omega = exp(icbrtn);
  omega *= P;
  omega *= pow(n,-11.0/3.0);
  T n2 = n*n;
  out = 4*a;
  out *= b;
  out *= P/n;
  tmp = pow(a,8.0/3.0);
  tmp += pow(b,8.0/3.0);
  tmp *= pow(2,11.0/3.0)*CF;
  T tmp2 = 47.0;
  tmp2 += -7.0*delta;
  tmp2 *= gnn/18.0;
  tmp += tmp2;
  
  tmp2 = -2.5;
  tmp2 += delta/18.0;
  tmp2 *= gaa + gbb;
  tmp += tmp2;

  tmp2 = 11.0 - delta;
  tmp2 /= 9.0*n;
  T tmp3 = a*gaa;
  tmp3 += b*gbb;
  tmp += tmp2*tmp3;

  tmp *= a*b;

  tmp2 = (-2.0/3.0)*n2;
  tmp2 *= gnn;
  tmp += tmp2;

  tmp2 = (-2.0/3.0)*n2;
  tmp2 += a*a;
  tmp2 *= gbb;
  tmp -= tmp2;

  tmp2 = (-2.0/3.0)*n2;
  tmp2 += b*b;
  tmp2 *= gaa;
  tmp -= tmp2;

  tmp *= B*omega;
  out += tmp;
  out *= -A;
}

void bench_double(double *dst, const double *src, size_t len)
{
  for (size_t i = 0; i < len/5; i++)
    lypc(dst[i],
	 src[i*5+0],
	 src[i*5+1],
	 src[i*5+2],
	 src[i*5+3],
	 src[i*5+4]);
}
void bench_double_disciplined(double *dst, const double *src, size_t len)
{
  for (size_t i = 0; i < len/5; i++)
    lypc_disciplined(dst[i],
		     src[i*5+0],
		     src[i*5+1],
		     src[i*5+2],
		     src[i*5+3],
		     src[i*5+4]);
}

// At blocksize = 4 this is already faster than the double version.
// with just a cos() the speed is nearly the same.
template<int blocksize>
void bench_numvec(double *dst, const double *src, size_t len)
{
  assert(len % blocksize == 0);
  const numvec<double,blocksize> *ps = reinterpret_cast< const numvec<double,blocksize>* > (src);
  numvec<double,blocksize> *pd = reinterpret_cast< numvec<double,blocksize>* > (dst);
  for (size_t i = 0; i < (len/blocksize)/5; i++)
    {
      lypc_disciplined(pd[i],
		       ps[i*5+0],
		       ps[i*5+1],
		       ps[i*5+2],
		       ps[i*5+3],
		       ps[i*5+4]);
    }
}


int main(int argc, const char *argv[])
{
  const int len = 1<<25;
  const int blocksize = 128;
  double *x = new double[len];
  double *y = new double[len];  
  for (int i=0;i<len;i++)
    x[i] = i % 10 + 0.1;
  for (int i=0;i<len;i++)
    y[i] = 0;
  cout << "Main loop starting" << endl;
  for (int i=0;i<10;i++)
    {
      clock_t tick = clock();
      bench_double(y,x,len);
      clock_t tock = clock();
      cout << y[len-1] << " Standard loop Time    : " << (tock - tick)/(double)CLOCKS_PER_SEC << endl;
      tick = clock();
      //bench_double_disciplined(y,x,len);
      bench_numvec<blocksize>(y,x,len);
      tock = clock();
      cout << y[len-1] << " Block<" << blocksize << "> vector Time: " << (tock - tick)/(double)CLOCKS_PER_SEC << endl;
    }
  delete[] x;
  delete[] y;
  return 0;
}
