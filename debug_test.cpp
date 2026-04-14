/**
 * Debug test to trace division issues
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "include/mpfun/mpfun.hpp"

using namespace mpfun;

void print_mp(const char* name, const MPFloat100& x) {
    const int64_t* d = x.data();
    std::cout << name << ":" << std::endl;
    std::cout << "  double approx: " << std::setprecision(17) << x.to_double() << std::endl;
    std::cout << "  sign_length: " << d[2] << std::endl;
    std::cout << "  exponent: " << d[3] << std::endl;
    std::cout << "  mantissa[0]: " << d[4] << std::endl;
    std::cout << "  mantissa[1]: " << d[5] << std::endl;
    std::cout << "  is_zero: " << x.is_zero() << std::endl;
    std::cout << "  is_nan: " << x.is_nan() << std::endl;
    std::cout << std::endl;
}

int main() {
    std::cout << std::setprecision(17);
    
    std::cout << "=== Test 1: 7 * (1/7) ===" << std::endl;
    {
        MPFloat100 seven(7.0);
        print_mp("seven", seven);
        
        MPFloat100 one(1.0);
        print_mp("one", one);
        
        MPFloat100 inv7 = one / seven;
        print_mp("1/7", inv7);
        
        MPFloat100 product = seven * inv7;
        print_mp("7 * (1/7)", product);
        
        MPFloat100 error = product - one;
        print_mp("error", error);
        
        std::cout << "Error double: " << error.to_double() << std::endl;
    }
    
    std::cout << "\n=== Test 2: 3^2 + 4^2 and sqrt ===" << std::endl;
    {
        MPFloat100 three(3.0);
        MPFloat100 four(4.0);
        
        MPFloat100 nine = three * three;
        print_mp("3*3", nine);
        
        MPFloat100 sixteen = four * four;
        print_mp("4*4", sixteen);
        
        MPFloat100 sum = nine + sixteen;
        print_mp("9 + 16", sum);
        
        MPFloat100 result = sqrt(sum);
        print_mp("sqrt(25)", result);
    }
    
    std::cout << "\n=== Test 3: 1e20 / 1e10 ===" << std::endl;
    {
        MPFloat100 large(1e20);
        print_mp("1e20", large);
        
        MPFloat100 small(1e10);
        print_mp("1e10", small);
        
        MPFloat100 ratio = large / small;
        print_mp("1e20 / 1e10", ratio);
    }
    
    return 0;
}
