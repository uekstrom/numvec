# numvec
Vectorization helper class supporting a disciplined subset of C++ math functions without temporary vectors.
The idea is to make it possible to write moderately complicated expressions that will still be vectorized
efficiently. All decisions on memory allocation must be handled by the programmer. An example

    numvec<double,32> x, tmp, result;
    for (int i=0;i<x.size();i++)
       x[i] = i; // Initialize the vector elements

Assume we want to compute

    // result = pow(x,2)/7 + 3.14*x; 

but we have to do it without vector temporaries

    result = pow(x,2);
    result /= 7;
    result += 3.14*x;

It could be done without any use of temporary vectors, but if we want

    // result = pow(x,2)/7 + 3.14*exp(x);

then one temporary is needed:

    result = pow(x,2);
    result /= 7;
    tmp = exp(x);
    result += 3.14*tmp;

The Numvec library only intends to provide the most basic delayed
operations: Arithmetic expressions with two arguments, and standard
math functions. Three-operand versions could also be implemented if
needed for performance reasons.
