      program example
      implicit none
      real*8 qt,Q,z,fuu,exact,test_fun
      external test_fun,fbt
      integer nu
      
      qt = 1d0
      Q=2d0
      nu=1
      z=1d0
      exact = 0.0844046546397287
      call fbt(test_fun,qt,Q,nu,z,fuu)
      print *, 'ogata results is', fuu
      print *, 'analytic results is', exact
      print *, 'error is', (fuu-exact)/exact
      
      end program

      real*8 function test_fun(b)
      implicit none
      real*8 b
      
      test_fun = b**2d0*dexp(-b)
      
      end function
