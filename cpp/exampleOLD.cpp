// bessel.cpp
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/airy.hpp>
#include "FBT.h"


double test( double x){ return x*exp(-x);} // test function to transform

int main( void )
{
   double x = 2.387;
   int n = 3, c;

   printf( "Bessel functions for x = %f:\n", x );
   printf( "   Kind   Order  Function     Result\n\n" );
   printf( "   First  0      j0( x )     %f\n", j0( x ) );
   printf( "   First  1      j1( x )     %f\n", j1( x ) );
   for( c = 2; c < 5; c++ )
      printf( "   First  %d      jn( %d, x )  %f\n", c, c, jn( c, x ) );
   printf( "   Second 0      y0( x )     %f\n", y0( x ) );
   printf( "   Second 1      y1( x )     %f\n", y1( x ) );
   for( c = 2; c < 5; c++ )
      printf( "   Second %d      yn( %d, x )  %f\n", c, c, yn( c, x ) );


  try
  { // Try a zero order v.
    float dodgy_root = boost::math::cyl_bessel_j_zero(0.F, 1);
    std::cout << "boost::math::cyl_bessel_j_zero(0.F, 1) " << dodgy_root << std::endl;
    // Thrown exception Error in function boost::math::cyl_bessel_j_zero<double>(double, int):
    // Requested the 0'th zero of J0, but the rank must be > 0 !
  }
  catch (std::exception& ex)
  {
    std::cout << "Thrown exception " << ex.what() << std::endl;
  }


  FBT ogata = FBT(0); // Fourier Transform with J0
  std::cout << ogata.fbt(test,1.) << std::endl;
}
