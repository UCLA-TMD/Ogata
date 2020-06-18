      program example
      implicit none
      real*8 qt,Q,z,fuu,ex,exact,test_fun
      external test_fun,fbt
      integer nu,n
      integer i
      
      open(unit = 1, file = "output.dat")
      write(1, *) 'qT fbt exact'

      do i = 1, 100
      qt = 0.01*i       
      Q=1d0             ! The initial guess for the peak of the test function
      nu=0              ! nu the order of the Bessel function
      z=1d0             ! z momentum fraction for fragmentation, set equal to one for other applications
      n = 10            ! number of nodes
      ex = exact(qT)    ! analytic result
      call fbt(test_fun,qt,Q,nu,z,n,fuu)
      write(1,*) qT,fuu,ex
      enddo
      close(1)


      end program

      real*8 function test_fun(b)
      implicit none
      real*8 b
      
      test_fun = b*dexp(-b)
      
      end function

      real*8 function exact(qT)
      real*8 qT
      real*8 pi

      pi = datan(1d0)*4d0

      exact = (1d0+qT*qT)**(-1.5d0)/2d0/pi

      end
