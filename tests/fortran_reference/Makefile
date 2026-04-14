# Makefile for generating MPFUN reference values
#
# Requires: gfortran and MPFUN2020-Fort source files
#
# Usage:
#   make MPFUN_SRC=/path/to/mpfun20-fort-v33/fortran-var1
#   ./generate_reference > reference_output.txt

FC = gfortran
FFLAGS = -O2 -Wall

# Path to MPFUN2020-Fort source
MPFUN_SRC ?= ../../../../mpfun20-fort-v33/fortran-var1

# MPFUN source files in dependency order
MPFUN_SRCS = \
    $(MPFUN_SRC)/mpfunf.f90 \
    $(MPFUN_SRC)/mpfuna.f90 \
    $(MPFUN_SRC)/mpfunb.f90 \
    $(MPFUN_SRC)/mpfunc.f90 \
    $(MPFUN_SRC)/mpfund.f90 \
    $(MPFUN_SRC)/mpfune.f90 \
    $(MPFUN_SRC)/mpfung1.f90 \
    $(MPFUN_SRC)/mpfunh1.f90 \
    $(MPFUN_SRC)/mpmodule.f90

TARGET = generate_reference

.PHONY: all clean run

all: $(TARGET)

$(TARGET): generate_reference.f90 $(MPFUN_SRCS)
	$(FC) $(FFLAGS) -o $@ $^

run: $(TARGET)
	./$(TARGET) > reference_output.txt
	@echo "Reference values written to reference_output.txt"

clean:
	rm -f $(TARGET) *.o *.mod reference_output.txt
