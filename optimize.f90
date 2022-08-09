PROGRAM optimize

!-----------------------------------------------------------------------
!----             This program minimize a two-dimension             ----
!----               function, using Steepest Descent                ----
!----                      and Newton-Raphson                       ---- 
!----                                                               ----
!----                         ~VCastor 2020                         ----
!-----------------------------------------------------------------------
! There is a pdf with documentation about the math theory necessary to
! understand properly the code. opt_explanation.pdf

IMPLICIT NONE
INTEGER                       :: iter, i, j
INTEGER, PARAMETER            :: max_iter = 40 
REAL*8                        :: x_0, y_0, z_0, alpha, e, er, tmp_z_0
REAL*8, DIMENSION(2)          :: g
REAL*8, DIMENSION(max_iter,2) :: coord
REAL*8, DIMENSION(2,2)        :: hess, invhess

alpha = 0.3        !step size
e     = 1.d-8      !stop condition

!-----------------------------------------------------------------------
!---- Starting point

WRITE(*,*) 'Where start?'
WRITE(*,FMT='(A)',ADVANCE='NO') 'x_0 = '
READ(*,*) x_0
WRITE(*,FMT='(A)',ADVANCE='NO') 'y_0 = '
READ(*,*) y_0

!-----------------------------------------------------------------------
!---- Steepest Descent

z_0     = f(x_0,y_0)
tmp_z_0 = z_0         !value starting point

!---- table format
WRITE(*,*) '**************************************************'
WRITE(*,*) '*****************Steepest Descent*****************'
WRITE(*,*) 'iteration  x       y      f(x,y)       gradient'

iter = 1; er = 1.d0
DO WHILE ( er .GT. e .AND. iter .LE. max_iter )
  coord(iter,1) = x_0
  coord(iter,2) = y_0
  g(:) = gradient(coord(iter,1),coord(iter,2))
  g(:) = g(:)/NORM2(g)
  WRITE(*,FMT='(I9,3F8.3,A4,F6.3,A1,F6.3,A1)') iter-1, x_0, y_0, z_0, &
                                                '   (',g(1),',',g(2),')'
  x_0     = x_0 - alpha*g(1)
  y_0     = y_0 - alpha*g(2)
  z_0     = f(x_0,y_0)
  er      = ABS(z_0 - tmp_z_0)
  tmp_z_0 = z_0
  iter    = iter + 1
ENDDO

!---- Final result
WRITE(*,*)
WRITE(*,FMT='(A,I2,A)') 'The compute is finished with: ', iter-1, &
                                                           ' iterations'
WRITE(*,*) 'in the point: '
WRITE(*,*) '  x        y        f(x,y)   '
WRITE(*,FMT='(3F9.4)') x_0, y_0, z_0

!---- File with coordiantes
OPEN(11,FILE='coordinates_Steepest_Descent.out')
  DO i = 1, iter-1
    WRITE(11,*) coord(i,1), coord(i,2)
  ENDDO
CLOSE(11)

!-----------------------------------------------------------------------
!---- Newton-Raphson

!coord(1,:) = coord(iter-1,:)   !the 1st one will be the resoult of SD
DO i = 2, max_iter              !restart coordinates, only stay with 1st
  coord(i,1) = 0.d0
  coord(i,2) = 0.d0
ENDDO

er = 1.d0; iter=1

WRITE(*,*)
WRITE(*,*) '**************************************************'
WRITE(*,*) '******************Newton-Raphson******************'
WRITE(*,*) 'iteration  x       y      f(x,y)       &
                                   &gradient                    hessian'
DO WHILE ( er .GT. e .AND. iter .LT. max_iter )
  g               = gradient(coord(iter,1),coord(iter,2))
  hess            = hessian(coord(iter,1),coord(iter,2))
  invhess         = inverse(hess,2)
  coord(iter+1,:) = coord(iter,:) - MATMUL(invhess,g)
  x_0 = coord(iter,1)
  y_0 = coord(iter,2)
  z_0 = f(coord(iter,1),coord(iter,2))
  WRITE(*,FMT='(I9,3F8.3,A4,F6.3,A1,F7.3,A5,4F7.3,A1)') iter-1, x_0,&
                         y_0, z_0, '  (',g(1),',',g(2),')   (', hess,')'
  er   = NORM2(g)
  iter = iter + 1
ENDDO

!---- Final result
WRITE(*,*)
WRITE(*,FMT='(A,I2,A)') 'The compute is finished with: ', iter-1, &
                                                           ' iterations'
WRITE(*,*) 'with the next values: '
WRITE(*,*) '  x        y        f(x,y)   gradient_x  gradient_y'
WRITE(*,FMT='(3F9.4,2F11.4)') x_0, y_0, z_0,g
WRITE(*,*) 'hessian'
DO i = 1, 2
  WRITE(*,FMT='(2F11.4)') (hess(i,j), j=1,2)
ENDDO

!---- File with coordiantes
OPEN(11,FILE='coordinates_Newton-Raphson.out')
  DO i = 1, iter-1
    WRITE(11,*) coord(i,1), coord(i,2)
  ENDDO
CLOSE(11)

!-----------------------------------------------------------------------
                                 CONTAINS
!-----------------------------------------------------------------------
    FUNCTION f(x,y) RESULT (val)
    IMPLICIT NONE
    REAL*8 :: x, y, val
   
    val = DSIN(x+y) + (x-y)**2 - 1.5d0*x + 3.5d0*y +3.d0
    ENDFUNCTION f
!-----------------------------------------------------------------------
    FUNCTION f_prime(d,x,y) RESULT (prime)
    IMPLICIT NONE
    INTEGER :: d
    REAL*8  :: h, x_p, y_p, x, y, prime
    ! d=1   :: \frac{\partial}{\partial x}
    ! d=2   :: \frac{\partial}{\partial y}
    h = 0.00001d0
    x_p = x;      y_p = y
    IF (d .EQ. 1) x_p = x_p + h
    IF (d .EQ. 2) y_p = y_p + h
    prime = ( f(x_p,y_p) - f(x,y) ) / h
    ENDFUNCTION f_prime
!-----------------------------------------------------------------------
    FUNCTION f_dprime(d,x,y) RESULT (dprime)
    IMPLICIT NONE
    INTEGER :: d
    REAL*8  :: h, x, y, dprime, fl, fr
    ! d=1   :: \frac{\partial^2}{\partial x^2}
    ! d=2   :: \frac{\partial^2}{\partial x \partial y}
    ! d=3   :: \frac{\partial^2}{\partial y \partial x}
    ! d=4   :: \frac{\partial^2}{\partial y^2}
    h = 0.00001d0
    SELECTCASE(d)
      CASE(1)
        fl = f_prime(1,x-h,y)
        fr = f_prime(1,x+h,y)
      CASE(2)
        fl = f_prime(1,x,y-h)
        fr = f_prime(1,x,y+h)
      CASE(3)
        fl = f_prime(2,x-h,y)
        fr = f_prime(2,x+h,y)
      CASE(4)
        fl = f_prime(2,x,y-h)
        fr = f_prime(2,x,y+h)
    ENDSELECT
    dprime = (fr-fl)/(2.d0*h)
    ENDFUNCTION
!-----------------------------------------------------------------------
    FUNCTION gradient(x,y) RESULT(g)
    IMPLICIT NONE
    REAL*8 :: x, y
    REAL*8, DIMENSION(2) :: g

    g(1) = f_prime(1,x,y)
    g(2) = f_prime(2,x,y)
    ENDFUNCTION
!-----------------------------------------------------------------------
    FUNCTION hessian(x,y) RESULT(hess)
    IMPLICIT NONE
    REAL*8 :: x, y
    REAL*8, DIMENSION(2,2) :: hess

    hess(1,1) = f_dprime(1,x,y)
    hess(1,2) = f_dprime(2,x,y)
    hess(2,1) = f_dprime(3,x,y)
    hess(2,2) = f_dprime(4,x,y)
    ENDFUNCTION
!-----------------------------------------------------------------------
    FUNCTION inverse(M,n) RESULT (inv)
    IMPLICIT NONE
    INTEGER                  :: i, j, k, n, l, aux3
    REAL*8                   :: aux, aux2
    REAL*8, DIMENSION(n,n)   :: M, inv
    REAL*8, DIMENSION(2*n)   :: V
    REAL*8, DIMENSION(n,2*n) :: S
    !---- First, add the identity
    S(:,:) = 0.d0
    DO i = 1, n
      aux3 = n+i            !"diagonal" on the right side
      S(i,aux3) = 1.d0
    ENDDO
    DO i = 1, n
      DO j = 1, n
        S(i,j) = M(i,j)     !data in the left side
      ENDDO
    ENDDO
    !---- 2nd up diagonal
    DO i = 1, n-1
      aux = S(i,i)          !diagonal elemnt
      DO j = 1, 2*n
        V(j) = S(i,j)/aux   !divided by the diagonal elemnt
      ENDDO
      DO k = i+1, n
        aux2 = S(k,i)
        DO j = 1, 2*n       !Actualize the matrix
          S(k,j) = S(k,j) - V(j)*aux2
        ENDDO
      ENDDO
    ENDDO
    !---- 3rd diagonal, zeros in the up part of matrix
    DO i = 1, n-1
      DO l = i+1, n
        DO k = l, n
          aux = S(k,k)
          DO j = 1, 2*n
            V(j) = S(k,j)/aux
          ENDDO
          aux2 = S(i,k)
          DO j = 1, 2*n
            S(i,j) = S(i,j) - V(j)*aux2
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !---- 4th diagonal equals to one
    DO i = 1, n
      aux = S(i,i)
      DO j = 1, 2*n
        S(i,j) = S(i,j)/aux
      ENDDO
    ENDDO
    !---- only the inverse
    DO i = 1, n
      DO j = 1, n
        aux3 = n+j
        inv(i,j) = S(i,aux3)
      ENDDO
    ENDDO
    ENDFUNCTION
ENDPROGRAM optimize
