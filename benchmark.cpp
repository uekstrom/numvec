#include <iostream>
#include <ctime>
#include <cassert>
using namespace std;
#include "numvec.hpp"

void bench_double(double *dst, const double *src, size_t len)
{
  for (size_t i = 0; i < len; i++)
    dst[i] = exp(cos(src[i])) + 1.8*src[i];
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
      *pd = exp(tmp);
      *pd += 1.8*(*ps);
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
      bench_numvec<blocksize>(y,x,len);
      tock = clock();
      cout << y[len-1] << " Block<" << blocksize << "> vector Time: " << (tock - tick)/(double)CLOCKS_PER_SEC << endl;
    }
  delete[] x;
  delete[] y;
  return 0;
}
