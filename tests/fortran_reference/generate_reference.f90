!******************************************************************************
! generate_reference.f90
! 
! Generates reference values using MPFUN2020-Fort for validating the C++ port.
! 
! Compile:
!   gfortran -O2 -o generate_reference generate_reference.f90 \
!       mpfuna.f90 mpfunb.f90 mpfunc.f90 mpfund.f90 mpfune.f90 mpfunf.f90 \
!       mpfung1.f90 mpfunh1.f90 mpmodule.f90
!
! Run:
!   ./generate_reference > reference_output.txt
!
! The output can be parsed to update reference_values.hpp
!******************************************************************************

program generate_reference
    use mpmodule
    implicit none
    
    type(mp_real) :: a, b, c, d
    type(mp_real) :: pi_val, e_val, ln2_val
    integer :: i
    integer, parameter :: ndp = 100  ! Number of decimal places
    
    ! Initialize with 100+ digit precision
    call mpinit(ndp + 10)
    
    write(*,'(A)') '# MPFUN2020-Fort Reference Values'
    write(*,'(A)') '# Generated for Kokkos-MPFloat validation'
    write(*,'(A)') '#'
    write(*,'(A,I0,A)') '# Precision: ', ndp, ' decimal digits'
    write(*,'(A)') '#'
    write(*,'(A)') ''
    
    ! =========================================================================
    ! Constants
    ! =========================================================================
    
    write(*,'(A)') '=== CONSTANTS ==='
    
    ! Pi
    pi_val = mppi()
    write(*,'(A)') 'PI:'
    call mpwrite(6, ndp + 5, ndp, pi_val)
    
    ! e (Euler's number)
    a = mppic('1.0')
    e_val = exp(a)
    write(*,'(A)') 'E:'
    call mpwrite(6, ndp + 5, ndp, e_val)
    
    ! log(2)
    a = mppic('2.0')
    ln2_val = log(a)
    write(*,'(A)') 'LN2:'
    call mpwrite(6, ndp + 5, ndp, ln2_val)
    
    write(*,'(A)') ''
    
    ! =========================================================================
    ! Arithmetic Operations
    ! =========================================================================
    
    write(*,'(A)') '=== ARITHMETIC ==='
    
    ! Test case 1: Simple addition
    a = mppic('3.14159265358979323846')
    b = mppic('2.71828182845904523536')
    
    write(*,'(A)') 'A:'
    call mpwrite(6, ndp + 5, ndp, a)
    write(*,'(A)') 'B:'
    call mpwrite(6, ndp + 5, ndp, b)
    
    c = a + b
    write(*,'(A)') 'A + B:'
    call mpwrite(6, ndp + 5, ndp, c)
    
    c = a - b
    write(*,'(A)') 'A - B:'
    call mpwrite(6, ndp + 5, ndp, c)
    
    c = a * b
    write(*,'(A)') 'A * B:'
    call mpwrite(6, ndp + 5, ndp, c)
    
    c = a / b
    write(*,'(A)') 'A / B:'
    call mpwrite(6, ndp + 5, ndp, c)
    
    write(*,'(A)') ''
    
    ! =========================================================================
    ! Square Root
    ! =========================================================================
    
    write(*,'(A)') '=== SQRT ==='
    
    a = mppic('2.0')
    c = sqrt(a)
    write(*,'(A)') 'SQRT(2):'
    call mpwrite(6, ndp + 5, ndp, c)
    
    a = mppic('3.0')
    c = sqrt(a)
    write(*,'(A)') 'SQRT(3):'
    call mpwrite(6, ndp + 5, ndp, c)
    
    write(*,'(A)') ''
    
    ! =========================================================================
    ! Exponential and Logarithm
    ! =========================================================================
    
    write(*,'(A)') '=== EXP/LOG ==='
    
    a = mppic('1.0')
    c = exp(a)
    write(*,'(A)') 'EXP(1):'
    call mpwrite(6, ndp + 5, ndp, c)
    
    a = mppic('2.0')
    c = exp(a)
    write(*,'(A)') 'EXP(2):'
    call mpwrite(6, ndp + 5, ndp, c)
    
    a = mppic('10.0')
    c = log(a)
    write(*,'(A)') 'LOG(10):'
    call mpwrite(6, ndp + 5, ndp, c)
    
    write(*,'(A)') ''
    
    ! =========================================================================
    ! Trigonometric Functions
    ! =========================================================================
    
    write(*,'(A)') '=== TRIG ==='
    
    ! sin(pi/6) = 0.5
    a = mppi() / mppic('6.0')
    write(*,'(A)') 'PI/6:'
    call mpwrite(6, ndp + 5, ndp, a)
    c = sin(a)
    write(*,'(A)') 'SIN(PI/6):'
    call mpwrite(6, ndp + 5, ndp, c)
    
    ! cos(pi/3) = 0.5
    a = mppi() / mppic('3.0')
    c = cos(a)
    write(*,'(A)') 'COS(PI/3):'
    call mpwrite(6, ndp + 5, ndp, c)
    
    ! sin(pi/4) = cos(pi/4) = 1/sqrt(2)
    a = mppi() / mppic('4.0')
    c = sin(a)
    write(*,'(A)') 'SIN(PI/4):'
    call mpwrite(6, ndp + 5, ndp, c)
    
    ! tan(pi/4) = 1
    c = tan(a)
    write(*,'(A)') 'TAN(PI/4):'
    call mpwrite(6, ndp + 5, ndp, c)
    
    write(*,'(A)') ''
    
    ! =========================================================================
    ! Identity Checks
    ! =========================================================================
    
    write(*,'(A)') '=== IDENTITIES ==='
    
    ! sin^2 + cos^2 = 1
    a = mppic('0.7')
    b = sin(a)
    c = cos(a)
    d = b*b + c*c
    write(*,'(A)') 'SIN(0.7)^2 + COS(0.7)^2:'
    call mpwrite(6, ndp + 5, ndp, d)
    
    ! exp(log(x)) = x
    a = mppic('12.345')
    c = exp(log(a))
    write(*,'(A)') 'EXP(LOG(12.345)):'
    call mpwrite(6, ndp + 5, ndp, c)
    
    write(*,'(A)') ''
    write(*,'(A)') '# End of reference values'
    
end program generate_reference
