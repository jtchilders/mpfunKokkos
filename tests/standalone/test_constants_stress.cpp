/**
 * Quick stress test for constants at various precisions
 */
#include <cstdio>
#include <cmath>
#include "mpfun/mp_float.hpp"
#include "mpfun/transcendental/constants.hpp"

using namespace mpfun;

int main() {
    printf("Stress testing constants at various precisions...\n\n");
    
    // Test pi at W=2,4,8,16 
    {
        printf("Pi consistency test:\n");
        double pi2 = pi<2>().to_double();
        double pi4 = pi<4>().to_double();
        double pi8 = pi<8>().to_double();
        double pi16 = pi<16>().to_double();
        printf("  W=2:  %.15g\n", pi2);
        printf("  W=4:  %.15g\n", pi4);
        printf("  W=8:  %.15g\n", pi8);
        printf("  W=16: %.15g\n", pi16);
        
        double max_err = 0;
        max_err = std::max(max_err, std::fabs(pi2 - M_PI));
        max_err = std::max(max_err, std::fabs(pi4 - M_PI));
        max_err = std::max(max_err, std::fabs(pi8 - M_PI));
        max_err = std::max(max_err, std::fabs(pi16 - M_PI));
        printf("  Max error: %.2e (should be < 1e-14)\n\n", max_err);
    }
    
    // Test e
    {
        printf("e consistency test:\n");
        double e2 = e<2>().to_double();
        double e4 = e<4>().to_double();
        double e8 = e<8>().to_double();
        printf("  W=2: %.15g\n", e2);
        printf("  W=4: %.15g\n", e4);
        printf("  W=8: %.15g\n", e8);
        
        double max_err = 0;
        max_err = std::max(max_err, std::fabs(e2 - M_E));
        max_err = std::max(max_err, std::fabs(e4 - M_E));
        max_err = std::max(max_err, std::fabs(e8 - M_E));
        printf("  Max error: %.2e (should be < 1e-14)\n\n", max_err);
    }
    
    // Test ln2
    {
        printf("ln(2) consistency test:\n");
        double ln2_2 = ln2<2>().to_double();
        double ln2_4 = ln2<4>().to_double();
        double ln2_8 = ln2<8>().to_double();
        double ref = std::log(2.0);
        printf("  W=2: %.15g\n", ln2_2);
        printf("  W=4: %.15g\n", ln2_4);
        printf("  W=8: %.15g\n", ln2_8);
        
        double max_err = 0;
        max_err = std::max(max_err, std::fabs(ln2_2 - ref));
        max_err = std::max(max_err, std::fabs(ln2_4 - ref));
        max_err = std::max(max_err, std::fabs(ln2_8 - ref));
        printf("  Max error: %.2e (should be < 1e-14)\n\n", max_err);
    }
    
    // Test that pi * 2 == two_pi
    {
        printf("Identity: 2*pi == two_pi:\n");
        auto p = pi<6>();
        auto two = MPFloat<6>(2.0);
        auto prod = p * two;
        auto tp = two_pi<6>();
        double err = std::fabs(prod.to_double() - tp.to_double());
        printf("  2*pi:    %.15g\n", prod.to_double());
        printf("  two_pi:  %.15g\n", tp.to_double());
        printf("  Error: %.2e\n\n", err);
    }
    
    printf("All stress tests completed.\n");
    return 0;
}
