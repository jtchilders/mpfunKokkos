/**
 * Minimal test to debug the trig functions
 */

#include <cstdio>
#include <cstdint>
#include <cmath>
#include <cstring>

#include "mpfun/core/representation.hpp"
#include "mpfun/core/add.hpp"
#include "mpfun/core/mul.hpp"
#include "mpfun/core/div.hpp"
#include "mpfun/core/sqrt.hpp"
#include "mpfun/mp_float.hpp"
#include "mpfun/transcendental/trig.hpp"

using namespace mpfun;

constexpr int WORDS = 8;
using MP = MPFloat<WORDS>;

int main() {
    printf("Starting minimal trig test...\n");
    fflush(stdout);
    
    printf("Creating MP value x = 0.0...\n");
    fflush(stdout);
    MP x(0.0);
    
    printf("x created, calling sin(x)...\n");
    fflush(stdout);
    
    MP s = sin(x);
    
    printf("sin(0) = %.15g\n", s.to_double());
    
    printf("Test complete!\n");
    return 0;
}
