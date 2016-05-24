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
  y = sin(u);
  res = pow(x,y);
}

void test_correctness()
{
  const int len = 4;
  numvec<double, len> inp, out;
  for (int i=0;i<len;i++)
    inp[i] = i+0.1;
  f1(out,inp);
  double sumres = 0;
  for (int i=0;i<len;i++)
    {
      double outd;
      f1(outd,inp[i]);
      sumres += fabs(out[i] - outd);
      cout << outd << " " << out[i] << endl;
    }
  cout.precision(15);
  cout << "Sum of absolute differences " << sumres << endl;
}


int main(int argc, const char *argv[])
{
  test_correctness();
  return 0;
}
