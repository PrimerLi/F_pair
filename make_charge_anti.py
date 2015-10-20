#!/usr/bin/env python

def main():
    import os
    compiler = "ifort"
    os.system(compiler + " -c F_charge_anti.f90")
    os.system(compiler + " -c lev.f90")
    os.system(compiler + " -c inverse.f90")
    os.system(compiler + " -c eigenvalues.f90")
    os.system(compiler + " -o F_charge_anti.out F_charge_anti.o lev.o inverse.o eigenvalues.o -L/usr/local/packages/lapack/3.4.2/INTEL-140-MVAPICH2-2.0/lib -llapack -lrefblas")
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
