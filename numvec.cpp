#include <iostream>
using namespace std;
#include "numvec.hpp"


template<typename T>
void f1(T &res, const T &x)
{
  T u,y = 2*x;
  y = 7*x;
  y -= 1.2;
  u = x + y;
  y = y/x;
  u = 3/x;
  u += -y;
  u += 1.1*x;
  res = sin(u);
}

void bench_double(double *dst, const double *src, size_t len)
{
  for (size_t i = 0; i < len; i++)
    dst[i] = sin(cos(src[i]));
}

// At blocksize = 4 this is already faster than the double version.
// with just a cos() the speed is nearly the same.
template<int blocksize>
void bench_numvec(double *dst, const double *src, size_t len)
{
  numvec<double,blocksize> tmp;
  assert(len % blocksize == 0);
  for (size_t i = 0; i < len/blocksize; i++)
    {
      const numvec<double,blocksize> *ps = reinterpret_cast< const numvec<double,blocksize>* > (src + i*blocksize);
      numvec<double,blocksize> *pd = reinterpret_cast< numvec<double,blocksize>* > (dst + i*blocksize);
      tmp = cos(*ps);
      *pd = sin(tmp);
    }
}


int main(int argc, const char *argv[])
{
  const int len = 1<<27;
  double *x = new double[len];
  double *y = new double[len];
  for (int i=0;i<len;i++)
    x[i] = i % 10 + 0.1;
  bench_double(y,x,len);
  //bench_numvec<4>(y,x,len);
  //  bench_numvec<>(y,x,len);
  cout <<  y[len-1] << endl;
  delete[] x;
  delete[] y;
  return 0;
}
