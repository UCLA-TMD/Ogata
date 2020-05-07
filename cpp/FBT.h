#ifndef __FBT_H
#define __FBT_H
#include <algorithm>    // std::copy
#include <vector>       // std::vector
#include <functional>   // std::bind

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/tools/minima.hpp>

class FBT {
private:
  double nu; // nu is Bessel function order
  int N;     // N is number of function calls
  double Q;  // a rough estimate where the function x*f(x) has maximum x = Q
  constexpr static double nu_def = 0.0;
  const static int N_def = 10;
  constexpr static double Q_def = 1.; // a rough estimate where the maximum of the function f(x) is
  std::vector<double> jn_zeros0;
  void acknowledgement();
  double ogatat(double (*f)(double), double q, double h);
  double ogatau(double (*f)(double), double q, double h);
  double get_hu(double (*f)(double), double q);
  double get_ht(double hu);

public:
  FBT(double _nu = nu_def, int _N = N_def, double _Q = Q_def); // Constructor
  ~FBT(); // Deconstructor

  double fbtu(double (*g)(double), double q); // unmodified original Ogata
  double fbt(double (*g)(double), double q);  // modified Ogata

};

#endif // #ifndef __FBT_H
