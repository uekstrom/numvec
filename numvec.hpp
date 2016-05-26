#ifndef NUMVEC_HPP
#define NUMVEC_HPP
#include <cmath>
#ifdef NUMVEC_USE_VML
#include <mkl.h>
#endif

// Light weight fixed length vectors with some delayed operations
// By Ulf Ekstrom (uekstrom@gmail.com) 2016.
// The idea is to be able to use a diciplined subset of the C++
// math features without creating any temporaries. For example
// u += 12*v;
// can be evaluated without creating a temporary vector. 
// TODO: Complete the list of math functions
// TODO: Full coverage testing
// TODO: Refactor the code to use better name space separation for
//       the auxiliary templates and types.

template<typename S, typename OP> 
struct numvec_delayed;

struct numvec_op_neg;

struct numvec_op_add_sv;
struct numvec_op_add_vv;

struct numvec_op_sub_sv;
struct numvec_op_sub_vs;
struct numvec_op_sub_vv;

struct numvec_op_mul_sv;
struct numvec_op_mul_vv;

struct numvec_op_div_sv;
struct numvec_op_div_vs;
struct numvec_op_div_vv;

template<typename T, int len>
class numvec
{
public:
  T c[len];
  numvec() {}
  numvec(const T &val) { *this = val; }
  int size() const { return len; }
  T &operator[](int i)
  {
    return c[i];
  }
  const T &operator[](int i) const
  {
    return c[i];
  }
  numvec &operator=(const T &scalar)  { for (int i=0;i<len;i++) c[i] = scalar;  return *this; }
  numvec &operator=(const numvec &v)  { for (int i=0;i<len;i++) c[i] = v.c[i];  return *this; }
  numvec &operator+=(const T &scalar) { for (int i=0;i<len;i++) c[i] += scalar; return *this; }
  numvec &operator+=(const numvec &v) { for (int i=0;i<len;i++) c[i] += v.c[i]; return *this; }
  numvec &operator-=(const T &scalar) { for (int i=0;i<len;i++) c[i] -= scalar; return *this; }
  numvec &operator-=(const numvec &v) { for (int i=0;i<len;i++) c[i] -= v.c[i]; return *this; }
  numvec &operator*=(const T &scalar) { for (int i=0;i<len;i++) c[i] *= scalar; return *this; }
  numvec &operator*=(const numvec &v) { for (int i=0;i<len;i++) c[i] *= v.c[i]; return *this; }
  numvec &operator/=(const T &scalar) { for (int i=0;i<len;i++) c[i] /= scalar; return *this; }
  numvec &operator/=(const numvec &v) { for (int i=0;i<len;i++) c[i] /= v.c[i]; return *this; }
  // Operations that need to be delayed for efficiency.
  numvec_delayed<numvec<T,len>,numvec_op_add_sv> operator+(const T &scalar) const;
  numvec_delayed<numvec<T,len>,numvec_op_add_vv> operator+(const numvec<T,len> &v) const;
  numvec_delayed<numvec<T,len>,numvec_op_sub_vs> operator-(const T &scalar) const;
  numvec_delayed<numvec<T,len>,numvec_op_sub_vv> operator-(const numvec<T,len> &v) const;
  numvec_delayed<numvec<T,len>,numvec_op_mul_sv> operator*(const T &scalar) const;
  numvec_delayed<numvec<T,len>,numvec_op_mul_vv> operator*(const numvec<T,len> &v) const;
  numvec_delayed<numvec<T,len>,numvec_op_div_vs> operator/(const T &scalar) const;
  numvec_delayed<numvec<T,len>,numvec_op_div_vv> operator/(const numvec<T,len> &v) const;
  numvec_delayed<numvec<T,len>,numvec_op_neg> operator-() const;

  // Construct or assign from delayed operations
  template<typename OP>
  numvec(const numvec_delayed<numvec<T,len>,OP> &from)
  {
    from.apply(*this);
  }
  template<typename OP>
  numvec &operator=(const numvec_delayed<numvec<T,len>,OP> &from)
  {
    from.apply(*this);
    return *this;
  }
  template<typename OP>
  numvec &operator+=(const numvec_delayed<numvec<T,len>,OP> &from)
  {
    from.apply_addto(*this);
    return *this;
  }
  template<typename OP>
  numvec &operator*=(const numvec_delayed<numvec<T,len>,OP> &from)
  {
    from.apply_multo(*this);
    return *this;
  }
};

// NEGATION
template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_neg>	
{
  numvec_delayed(const numvec<T,len> &right_) : right(right_) {}
  const numvec<T,len> &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = -right[i];
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] -= right[i];
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= -right[i];
  }
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_neg> numvec<T,len>::operator-() const
{
  return numvec_delayed<numvec<T,len>,numvec_op_neg>(*this);
}


// ADDITION

template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_add_sv>	
{
  numvec_delayed(const T &left_, const numvec<T,len> &right_) : left(left_), right(right_) {}
  const T left;
  const numvec<T,len> &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = left + right[i];
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += left + right[i];
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= left + right[i];
  }
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_add_sv> numvec<T,len>::operator+(const T &scalar) const
{
  return numvec_delayed<numvec<T,len>,numvec_op_add_sv>(scalar, *this);
}
template<typename T, int len, typename S>
numvec_delayed<numvec<T,len>,numvec_op_add_sv> operator+(const S &left, const numvec<T,len> &right)
{
  return numvec_delayed<numvec<T,len>,numvec_op_add_sv>(left, right);
}


template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_add_vv>	
{
  numvec_delayed(const numvec<T,len> &left_, const numvec<T,len> &right_) : left(left_), right(right_) {}
  const numvec<T,len> &left, &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = left[i] + right[i];
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += left[i] + right[i];
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= left[i] + right[i];
  }
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_add_vv> numvec<T,len>::operator+(const numvec<T,len> &right) const
{
  return numvec_delayed<numvec<T,len>,numvec_op_add_vv>(*this,right);
}


// SUBTRACTION

template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_sub_sv>	
{
  numvec_delayed(const T &left_, const numvec<T,len> &right_) : left(left_), right(right_) {}
  const T left;
  const numvec<T,len> &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = left - right[i];
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += left - right[i];
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= left - right[i];
  }
};
template<typename T, int len, typename S>
numvec_delayed<numvec<T,len>,numvec_op_sub_sv> operator-(const S &left, const numvec<T,len> &right)
{
  return numvec_delayed<numvec<T,len>,numvec_op_sub_sv>(left, right);
}


template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_sub_vs>	
{
  numvec_delayed(const numvec<T,len> &left_, const T &right_) : left(left_), right(right_) {}
  const numvec<T,len> &left;
  const T &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = left[i] - right;
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += left[i] - right;
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= left[i] - right;
  }
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_sub_vs> numvec<T,len>::operator-(const T &scalar) const
{
  return numvec_delayed<numvec<T,len>,numvec_op_sub_vs>(*this, scalar);
}


template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_sub_vv>	
{
  numvec_delayed(const numvec<T,len> &left_, const numvec<T,len> &right_) : left(left_), right(right_) {}
  const numvec<T,len> &left, &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = left[i] - right[i];
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += left[i] - right[i];
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= left[i] - right[i];
  }
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_sub_vv> numvec<T,len>::operator-(const numvec<T,len> &right) const
{
  return numvec_delayed<numvec<T,len>,numvec_op_sub_vv>(*this, right);
}

// MULTIPLICATION

template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_mul_sv>	
{
  numvec_delayed(const T &left_, const numvec<T,len> &right_) : left(left_), right(right_) {}
  const T left;
  const numvec<T,len> &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = left*right[i];
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += left*right[i];
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= left*right[i];
  }
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_mul_sv> numvec<T,len>::operator*(const T &scalar) const
{
  return numvec_delayed<numvec<T,len>,numvec_op_mul_sv>(scalar, *this);
}
template<typename T, int len, typename S>
numvec_delayed<numvec<T,len>,numvec_op_mul_sv> operator*(const S &left, const numvec<T,len> &right)
{
  return numvec_delayed<numvec<T,len>,numvec_op_mul_sv>(left, right);
}


template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_mul_vv>	
{
  numvec_delayed(const numvec<T,len> &left_, const numvec<T,len> &right_) : left(left_), right(right_) {}
  const numvec<T,len> &left, &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = left[i]*right[i];
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += left[i]*right[i];
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= left[i]*right[i];
  }
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_mul_vv> numvec<T,len>::operator*(const numvec<T,len> &right) const
{
  return numvec_delayed<numvec<T,len>,numvec_op_mul_vv>(*this,right);
}

// DIVISION

template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_div_sv>	
{
  numvec_delayed(const T &left_, const numvec<T,len> &right_) : left(left_), right(right_) {}
  const T left;
  const numvec<T,len> &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = left/right[i];
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += left/right[i];
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= left/right[i];
  }
};
template<typename T, int len, typename S>
numvec_delayed<numvec<T,len>,numvec_op_div_sv> operator/(const S &left, const numvec<T,len> &right)
{
  return numvec_delayed<numvec<T,len>,numvec_op_div_sv>(left, right);
}


template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_div_vs>	
{
  numvec_delayed(const numvec<T,len> &left_, const T &right_) : left(left_), right(right_) {}
  const numvec<T,len> &left;
  const T &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = left[i]/right;
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += left[i]/right;
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= left[i]/right;
  }
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_div_vs> numvec<T,len>::operator/(const T &scalar) const
{
  return numvec_delayed<numvec<T,len>,numvec_op_div_vs>(*this, scalar);
}


template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_div_vv>	
{
  numvec_delayed(const numvec<T,len> &left_, const numvec<T,len> &right_) : left(left_), right(right_) {}
  const numvec<T,len> &left, &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = left[i]/right[i];
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += left[i]/right[i];
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= left[i]/right[i];
  }
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_div_vv> numvec<T,len>::operator/(const numvec<T,len> &right) const
{
  return numvec_delayed<numvec<T,len>,numvec_op_div_vv>(*this, right);
}

/// MATH FUNCTIONS BELOW THIS LINE

// pow (only binary math function supported so far)

struct numvec_op_pow_sv;

template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_pow_sv>	
{
  numvec_delayed(const T &left_, const numvec<T,len> &right_) : left(left_), right(right_) {}
  const T left;
  const numvec<T,len> &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = pow(left,right[i]);
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += pow(left,right[i]);
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= pow(left,right[i]);
  }
};
template<typename T, int len, typename S>
numvec_delayed<numvec<T,len>,numvec_op_pow_sv> pow(const S &left, const numvec<T,len> &right)
{
  return numvec_delayed<numvec<T,len>,numvec_op_pow_sv>(left, right);
}

struct numvec_op_pow_vs;

template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_pow_vs>	
{
  numvec_delayed(const numvec<T,len> &left_, const T &right_) : left(left_), right(right_) {}
  const numvec<T,len> &left;
  const T &right;
#ifndef NUMVEC_USE_VML
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = pow(left[i],right);
  }
#else
  void apply(numvec<T,len> &arg) const
  {
    vdPowx(arg.size(),left.c, right, arg.c);
  }
#endif
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += pow(left[i],right);
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= pow(left[i],right);
  }
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_pow_vs> pow(const numvec<T,len> &left, const T &scalar)
{
  return numvec_delayed<numvec<T,len>,numvec_op_pow_vs>(left, scalar);
}

struct numvec_op_pow_vi;

template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_pow_vi>	
{
  numvec_delayed(const numvec<T,len> &left_, const T &right_) : left(left_), right(right_) {}
  const numvec<T,len> &left;
  int right;
#ifndef NUMVEC_USE_VML
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = pow(left[i],right);
  }
#else
  void apply(numvec<T,len> &arg) const
  {
    vdPowx(arg.size(),left.c, right, arg.c);
  }
#endif
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += pow(left[i],right);
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= pow(left[i],right);
  }
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_pow_vi> pow(const numvec<T,len> &left, int scalar)
{
  return numvec_delayed<numvec<T,len>,numvec_op_pow_vi>(left, scalar);
}

struct numvec_op_pow_vv;

template<typename T, int len>
struct numvec_delayed<numvec<T,len>,numvec_op_pow_vv>	
{
  numvec_delayed(const numvec<T,len> &left_, const numvec<T,len> &right_) : left(left_), right(right_) {}
  const numvec<T,len> &left, &right;
  void apply(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] = pow(left[i],right[i]);
  }
  void apply_addto(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] += pow(left[i],right[i]);
  }
  void apply_multo(numvec<T,len> &arg) const
  {
    for (int i=0;i<arg.size();i++)
      arg[i] *= pow(left[i],right[i]);
  }
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_pow_vv> pow(const numvec<T,len> &left, const numvec<T,len> &right)
{
  return numvec_delayed<numvec<T,len>,numvec_op_pow_vv>(left, right);
}


// UNARY MATH FUNCTIONS

#define NUMVEC_UNARY(FUN)\
struct numvec_op_##FUN;\
template<typename S>\
struct numvec_delayed<S,numvec_op_##FUN>\
{\
  numvec_delayed(const S &val_) : val(val_) {}\
  const S &val;\
  void apply(S &arg) const\
  {\
    for (int i=0;i<arg.size();i++)\
      arg[i] = FUN(val[i]);\
  }\
  void apply_addto(S &arg) const\
  {\
    for (int i=0;i<arg.size();i++)\
      arg[i] += FUN(val[i]);\
  }\
  void apply_multo(S &arg) const\
  {\
    for (int i=0;i<arg.size();i++)\
      arg[i] *= FUN(val[i]);\
  }\
};\
template<typename T, int len>\
numvec_delayed<numvec<T,len>,numvec_op_##FUN> FUN(const numvec<T,len> &arg)\
{\
  return numvec_delayed<numvec<T,len>,numvec_op_##FUN>(arg);\
}

NUMVEC_UNARY(exp)
NUMVEC_UNARY(log)
NUMVEC_UNARY(sin)
NUMVEC_UNARY(cos)
NUMVEC_UNARY(tan)
NUMVEC_UNARY(asin)
NUMVEC_UNARY(acos)
NUMVEC_UNARY(atan)
NUMVEC_UNARY(sinh)
NUMVEC_UNARY(cosh)
NUMVEC_UNARY(tanh)
NUMVEC_UNARY(asinh)
NUMVEC_UNARY(acosh)
NUMVEC_UNARY(atanh)

#endif
