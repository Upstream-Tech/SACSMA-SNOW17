compile the fortran code
f2py -c sacsma.pyf sacsma.f
f2py sacsma.f -m primes -h sacsma.pyf
