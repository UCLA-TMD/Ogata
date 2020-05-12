// bessel.cpp
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <chrono>
#include "FBT.h"


double test( double x, void* data ){ return x*exp(-x);} // test function to transform data allows to send anything else to the function
double exact( double qT){ return pow((1.+qT*qT),-1.5);} // test function exact


double test1( double x, void* data ){ return x*exp(-x*x);} // test function to transform
double exact1( double qT ){ return exp(-qT*qT/4.)/2.;} // test function exact

int main( void )
{
  FBT ogata = FBT(); // Fourier Transform with Jnu, nu=0.0 and N=10
  double qT = 1.;
  double data = 1;

  auto begin = std::chrono::high_resolution_clock::now();
  double res = ogata.fbt(test,&data,qT);
  auto end = std::chrono::high_resolution_clock::now();

  std::cout << std::setprecision(30) << " FT( J0(x*qT) x*exp(-x) ) at qT= " << qT << std::endl;
  std::cout << std::setprecision(30) << "Numerical transformed = " << res << std::endl;
  auto overhead=std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout<<"Calc time: "<<overhead<<" nanoseconds\n";


  begin = std::chrono::high_resolution_clock::now();
  res = ogata.fbtu(test,&data,qT);
  end = std::chrono::high_resolution_clock::now();

  std::cout << std::setprecision(30) << "Numerical untransformed = " << res << std::endl;
  overhead=std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout<<"Calc time: "<<overhead<<" nanoseconds\n";

  std::cout << std::setprecision(30) << "Exact = " << exact(qT) << std::endl;


  begin = std::chrono::high_resolution_clock::now();
  res = ogata.fbt(test1,&data,qT);
  end = std::chrono::high_resolution_clock::now();

  std::cout << std::setprecision(30) << " FT( J0(x*qT) x*exp(-x**2) ) at qT= " << qT << std::endl;
  std::cout << std::setprecision(30) << "Numerical transformed = " << res << std::endl;
  overhead=std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout<<"Calc time: "<<overhead<<" nanoseconds\n";


  begin = std::chrono::high_resolution_clock::now();
  res = ogata.fbtu(test1,&data,qT);
  end = std::chrono::high_resolution_clock::now();

  std::cout << std::setprecision(30) << "Numerical untransformed = " << res << std::endl;
  overhead=std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout<<"Calc time: "<<overhead<<" nanoseconds\n";

  std::cout << std::setprecision(30) << "Exact = " << exact1(qT) << std::endl;

}
