#ifndef NUMVEC_HPP
#define NUMVEC_HPP
#include <cmath>

// Light weight fixed length vectors with some delayed operations
// By Ulf Ekstrom (uekstrom@gmail.com) 2016.
// The idea is to be able to use a diciplined subset of the C++
// math features without creating any temporaries. For example
// u += 12*v; 
// can be evaluated without creating a temporary vector. 
// TODO: Complete the list of math functions
// TODO: Implement pow(,)
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

private:
  T c[len];
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
};
template<typename T, int len>
numvec_delayed<numvec<T,len>,numvec_op_div_vv> numvec<T,len>::operator/(const numvec<T,len> &right) const
{
  return numvec_delayed<numvec<T,len>,numvec_op_div_vv>(*this, right);
}

/// MATH FUNCTIONS BELOW THIS LINE

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
