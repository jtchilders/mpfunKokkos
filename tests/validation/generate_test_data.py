#!/usr/bin/env python3
"""
Generate test data for MPFloat validation.

Uses mpmath for high-precision reference values, converts to
the MPFUN20-Fort format (60-bit words, radix 2^60).

Output format: Binary file with header + test cases
"""

import argparse
import struct
import sys
from pathlib import Path

try:
    from mpmath import mp, mpf, sqrt as mp_sqrt, exp as mp_exp, log as mp_log
    from mpmath import sin as mp_sin, cos as mp_cos, tan as mp_tan, atan as mp_atan
except ImportError:
    print("Error: mpmath is required. Install with: pip install mpmath")
    sys.exit(1)


# Constants matching representation.hpp
BITS_PER_WORD = 60
RADIX = 1 << BITS_PER_WORD
MASK = RADIX - 1
DIGITS_PER_WORD = 18.061799739838871713

# Operation codes
OP_ADD = 0
OP_SUB = 1
OP_MUL = 2
OP_DIV = 3
OP_SQRT = 4
OP_EXP = 5
OP_LOG = 6
OP_SIN = 7
OP_COS = 8
OP_TAN = 9
OP_ATAN = 10

OP_NAMES = {
    OP_ADD: "add",
    OP_SUB: "sub",
    OP_MUL: "mul",
    OP_DIV: "div",
    OP_SQRT: "sqrt",
    OP_EXP: "exp",
    OP_LOG: "log",
    OP_SIN: "sin",
    OP_COS: "cos",
    OP_TAN: "tan",
    OP_ATAN: "atan",
}


def mpf_to_mpfloat_words(x: mpf, words_per_number: int) -> list[int]:
    """
    Convert mpmath mpf to MPFloat word representation.
    
    Returns list of int64 values:
      [sign_and_length, exponent, mantissa[0], mantissa[1], ...]
    
    The mantissa is stored as 60-bit words, most significant first.
    Total words = words_per_number (includes sign_length and exponent).
    """
    mantissa_words = words_per_number - 2  # Reserve 2 for metadata
    
    if x == 0:
        return [0, 0] + [0] * mantissa_words
    
    # mpmath internal representation: x = sign * man * 2^exp
    # where man is an arbitrary-precision integer
    sign = 1 if x > 0 else -1
    man = abs(int(x.man))  # mantissa as Python int
    exp = int(x.exp)        # binary exponent
    
    # Compute the total number of bits in mantissa
    if man == 0:
        return [0, 0] + [0] * mantissa_words
    
    bit_length = man.bit_length()
    
    # We want to normalize so the mantissa, when interpreted as words,
    # has the most significant word in the range [1, RADIX).
    # 
    # Total value = sign * man * 2^exp
    # We want: sign * (w0 + w1/RADIX + w2/RADIX^2 + ...) * RADIX^E
    # where w0 is in [1, RADIX)
    
    # First, shift mantissa so it aligns to 60-bit word boundaries
    # Number of 60-bit words needed for the full mantissa
    full_words_needed = (bit_length + BITS_PER_WORD - 1) // BITS_PER_WORD
    
    # Pad mantissa to align to word boundary (shift left)
    total_bits = full_words_needed * BITS_PER_WORD
    shift_amount = total_bits - bit_length
    shifted_man = man << shift_amount
    
    # The effective exponent adjustment
    # Original: value = man * 2^exp
    # After shift: value = shifted_man * 2^(exp - shift_amount)
    # In terms of RADIX: value = (shifted_man / RADIX^full_words_needed) * 2^(exp - shift_amount + total_bits)
    #                          = (shifted_man / RADIX^full_words_needed) * RADIX^((exp - shift_amount + total_bits) / 60)
    
    # The binary exponent for the result
    adjusted_exp = exp + bit_length  # exponent if normalized to [1, 2)
    
    # Convert to RADIX exponent
    # We want: value = mantissa_frac * RADIX^radix_exp
    # where mantissa_frac = w0 + w1/RADIX + ... is in [1, RADIX)
    
    # radix_exp = ceil(adjusted_exp / 60) would give us the number of RADIX powers
    # Actually: we have value = man * 2^exp = (man * 2^bit_length) * 2^(exp - bit_length)
    #                        = (normalized_man) * 2^(exp + bit_length - bit_length)
    # Let's think differently...
    
    # radix_exp such that: value / RADIX^radix_exp is in [1, RADIX)
    # log_RADIX(value) = (exp + log2(man)) / 60 = (exp + bit_length) / 60
    radix_exp = (exp + bit_length + BITS_PER_WORD - 1) // BITS_PER_WORD
    
    # Now scale mantissa to fit with this exponent
    # value = man * 2^exp = normalized_value * RADIX^radix_exp
    # normalized_value = man * 2^exp / RADIX^radix_exp
    #                  = man * 2^(exp - radix_exp * 60)
    shift_for_radix = radix_exp * BITS_PER_WORD - exp - bit_length + BITS_PER_WORD
    
    if shift_for_radix >= 0:
        normalized_man = man << shift_for_radix
    else:
        normalized_man = man >> (-shift_for_radix)
    
    # Extract words (most significant first)
    words = []
    temp = normalized_man
    for _ in range(mantissa_words):
        # Extract from most significant end
        pass
    
    # Actually, let's extract words in a cleaner way
    # We have normalized_man, and we want to extract mantissa_words 60-bit chunks
    # starting from the most significant
    
    # Find total bits in normalized_man
    if normalized_man == 0:
        return [0, 0] + [0] * mantissa_words
        
    norm_bits = normalized_man.bit_length()
    total_word_bits = mantissa_words * BITS_PER_WORD
    
    # Pad to total_word_bits
    if norm_bits < total_word_bits:
        normalized_man <<= (total_word_bits - norm_bits)
    elif norm_bits > total_word_bits:
        normalized_man >>= (norm_bits - total_word_bits)
    
    # Extract words most-significant first
    words = []
    for i in range(mantissa_words - 1, -1, -1):
        word = (normalized_man >> (i * BITS_PER_WORD)) & MASK
        words.append(word)
    
    # Count actual non-zero words (from end)
    length = mantissa_words
    while length > 0 and words[length - 1] == 0:
        length -= 1
    if length == 0:
        length = 1  # At least one word
    
    # Verify first word is in valid range
    if words[0] == 0 and length > 0:
        # Need to shift left and adjust exponent
        shift = 0
        while shift < len(words) and words[shift] == 0:
            shift += 1
        if shift < len(words):
            words = words[shift:] + [0] * shift
            radix_exp -= shift
            length = mantissa_words - shift
            while length > 0 and words[length - 1] == 0:
                length -= 1
            if length == 0:
                length = 1
    
    sign_and_length = sign * length
    
    # Ensure we have exactly mantissa_words entries
    while len(words) < mantissa_words:
        words.append(0)
    words = words[:mantissa_words]
    
    return [sign_and_length, radix_exp] + words


def mpf_to_mpfloat_words_v2(x: mpf, words_per_number: int) -> list[int]:
    """
    Simpler, more direct conversion from mpmath to MPFloat format.
    
    Strategy: Convert the mpf value to our representation by repeated
    division/multiplication by RADIX.
    """
    mantissa_words = words_per_number - 2
    
    if x == 0:
        return [0, 0] + [0] * mantissa_words
    
    sign = 1 if x > 0 else -1
    x_abs = abs(x)
    
    # Find the exponent: largest E such that x_abs >= RADIX^E
    # Use high-precision logarithm
    if x_abs >= 1:
        radix_exp = 0
        temp = x_abs
        while temp >= RADIX:
            temp /= RADIX
            radix_exp += 1
    else:
        radix_exp = 0
        temp = x_abs
        while temp < 1:
            temp *= RADIX
            radix_exp -= 1
    
    # Now temp should be in [1, RADIX), and x_abs = temp * RADIX^radix_exp
    # Actually let's be more careful - re-normalize
    radix_exp += 1  # We want x_abs = temp * RADIX^radix_exp with temp in [1/RADIX, 1)
    # No, let's match MPFUN: value = mantissa * RADIX^exponent where mantissa[0] is in [1, RADIX)
    
    # Reset and compute properly
    radix_exp = 0
    temp = x_abs
    
    # Normalize so temp is in [1, RADIX)
    radix_f = mpf(RADIX)
    while temp >= radix_f:
        temp /= radix_f
        radix_exp += 1
    while temp < 1 and temp != 0:
        temp *= radix_f
        radix_exp -= 1
    
    # Now extract mantissa words
    words = []
    for _ in range(mantissa_words):
        word_val = int(temp)
        words.append(word_val)
        temp = (temp - word_val) * radix_f
    
    # Count actual used words
    length = mantissa_words
    while length > 1 and words[length - 1] == 0:
        length -= 1
    
    sign_and_length = sign * length
    
    return [sign_and_length, radix_exp] + words


def generate_binary_ops_tests(precision_digits: int, words_per_number: int, 
                               num_cases: int = 100) -> list[tuple]:
    """Generate test cases for binary operations."""
    mp.dps = precision_digits + 50  # Extra precision for intermediate calculations
    
    cases = []
    
    # Test values: various magnitudes and special cases
    test_values = [
        mpf("1.0"),
        mpf("2.0"),
        mpf("0.5"),
        mpf("3.14159265358979323846264338327950288419716939937510"),
        mpf("2.71828182845904523536028747135266249775724709369995"),
        mpf("1e10"),
        mpf("1e-10"),
        mpf("1e50"),
        mpf("1e-50"),
        mpf("0.1"),
        mpf("0.333333333333333333333333333333333333333333333333333"),
        mpf("1.4142135623730950488016887242096980785696718753769"),  # sqrt(2)
        mpf("1.7320508075688772935274463415058723669428052538104"),  # sqrt(3)
        mpf("123456789.987654321"),
        mpf("9999999999999999999.0000000000000000001"),
    ]
    
    # Generate random values with good coverage
    import random
    random.seed(42)  # Reproducible
    for _ in range(20):
        # Random mantissa
        digits = ''.join(str(random.randint(0, 9)) for _ in range(precision_digits))
        exp = random.randint(-30, 30)
        val = mpf(f"0.{digits}e{exp}")
        test_values.append(val)
    
    # Binary operations
    for op in [OP_ADD, OP_SUB, OP_MUL, OP_DIV]:
        for i, a in enumerate(test_values[:min(len(test_values), num_cases // 4)]):
            for b in test_values[:min(len(test_values), num_cases // 4 // len(test_values) + 1)]:
                if op == OP_DIV and b == 0:
                    continue
                
                if op == OP_ADD:
                    result = a + b
                elif op == OP_SUB:
                    result = a - b
                elif op == OP_MUL:
                    result = a * b
                elif op == OP_DIV:
                    result = a / b
                
                a_words = mpf_to_mpfloat_words_v2(a, words_per_number)
                b_words = mpf_to_mpfloat_words_v2(b, words_per_number)
                result_words = mpf_to_mpfloat_words_v2(result, words_per_number)
                
                cases.append((op, [a_words, b_words], result_words))
    
    return cases


def generate_unary_ops_tests(precision_digits: int, words_per_number: int,
                              num_cases: int = 50) -> list[tuple]:
    """Generate test cases for unary operations."""
    mp.dps = precision_digits + 50
    
    cases = []
    
    # SQRT test values (positive only)
    sqrt_values = [
        mpf("2.0"),
        mpf("3.0"),
        mpf("0.5"),
        mpf("10.0"),
        mpf("100.0"),
        mpf("1e10"),
        mpf("1e-10"),
        mpf("2.5"),
        mpf("123456789.0"),
    ]
    
    for x in sqrt_values[:num_cases]:
        result = mp_sqrt(x)
        x_words = mpf_to_mpfloat_words_v2(x, words_per_number)
        result_words = mpf_to_mpfloat_words_v2(result, words_per_number)
        cases.append((OP_SQRT, [x_words], result_words))
    
    # EXP test values (small enough to not overflow)
    exp_values = [
        mpf("1.0"),
        mpf("0.5"),
        mpf("0.1"),
        mpf("-1.0"),
        mpf("-0.5"),
        mpf("2.0"),
        mpf("10.0"),
        mpf("-10.0"),
    ]
    
    for x in exp_values[:num_cases // 2]:
        result = mp_exp(x)
        x_words = mpf_to_mpfloat_words_v2(x, words_per_number)
        result_words = mpf_to_mpfloat_words_v2(result, words_per_number)
        cases.append((OP_EXP, [x_words], result_words))
    
    # LOG test values (positive only)
    log_values = [
        mpf("2.0"),
        mpf("10.0"),
        mpf("0.5"),
        mpf("2.718281828459045"),
        mpf("100.0"),
        mpf("1e10"),
        mpf("1e-5"),
    ]
    
    for x in log_values[:num_cases // 2]:
        result = mp_log(x)
        x_words = mpf_to_mpfloat_words_v2(x, words_per_number)
        result_words = mpf_to_mpfloat_words_v2(result, words_per_number)
        cases.append((OP_LOG, [x_words], result_words))
    
    # Trig test values
    trig_values = [
        mpf("0.0"),
        mpf("0.5"),
        mpf("1.0"),
        mpf("-1.0"),
        mpf("3.14159265358979323846264338327950288419716939937510") / 4,  # pi/4
        mpf("3.14159265358979323846264338327950288419716939937510") / 6,  # pi/6
        mpf("3.14159265358979323846264338327950288419716939937510") / 3,  # pi/3
    ]
    
    for x in trig_values[:num_cases // 4]:
        for op, func in [(OP_SIN, mp_sin), (OP_COS, mp_cos), (OP_TAN, mp_tan)]:
            if op == OP_TAN and abs(x - mpf("1.5707963267948966")) < 0.01:
                continue  # Skip near pi/2
            result = func(x)
            x_words = mpf_to_mpfloat_words_v2(x, words_per_number)
            result_words = mpf_to_mpfloat_words_v2(result, words_per_number)
            cases.append((op, [x_words], result_words))
    
    # ATAN test values
    atan_values = [
        mpf("0.0"),
        mpf("1.0"),
        mpf("-1.0"),
        mpf("0.5"),
        mpf("2.0"),
        mpf("10.0"),
    ]
    
    for x in atan_values[:num_cases // 4]:
        result = mp_atan(x)
        x_words = mpf_to_mpfloat_words_v2(x, words_per_number)
        result_words = mpf_to_mpfloat_words_v2(result, words_per_number)
        cases.append((OP_ATAN, [x_words], result_words))
    
    return cases


def write_test_file(output_path: Path, precision_digits: int, words_per_number: int,
                    cases: list[tuple]) -> None:
    """Write test cases to binary file."""
    
    # Compute precision in bits (approximate)
    precision_bits = int(precision_digits * 3.32193)  # log2(10)
    
    with open(output_path, 'wb') as f:
        # Header (32 bytes)
        # Format: magic(4) + version(4) + precision_bits(4) + num_cases(4) + words_per_number(4) + reserved(12)
        header = struct.pack(
            '<4sIIII',  # little-endian: magic, version, precision_bits, num_cases, words_per_number
            b'MPFT',
            1,  # version
            precision_bits,
            len(cases),
            words_per_number
        )
        # Add 12 bytes of reserved/padding to reach 32 bytes
        header += b'\x00' * 12
        assert len(header) == 32, f"Header size mismatch: {len(header)}"
        f.write(header)
        
        # Write each test case
        for op, inputs, expected in cases:
            # Case header: operation (1 byte), num_inputs (1 byte), padding (2 bytes)
            case_header = struct.pack('<BBxx', op, len(inputs))
            f.write(case_header)
            
            # Write inputs
            for inp in inputs:
                for word in inp:
                    # Convert to signed int64 for struct.pack
                    f.write(struct.pack('<q', word))
            
            # Write expected output
            for word in expected:
                f.write(struct.pack('<q', word))
    
    print(f"Wrote {len(cases)} test cases to {output_path}")
    print(f"  Precision: {precision_digits} decimal digits ({precision_bits} bits)")
    print(f"  Words per number: {words_per_number}")
    print(f"  File size: {output_path.stat().st_size} bytes")


def main():
    parser = argparse.ArgumentParser(description='Generate MPFloat test data')
    parser.add_argument('--precision', '-p', type=int, default=100,
                        help='Precision in decimal digits (default: 100)')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output file path')
    parser.add_argument('--num-cases', '-n', type=int, default=100,
                        help='Number of test cases per operation type')
    parser.add_argument('--ops', type=str, default='all',
                        help='Operations to test: all, basic (add/sub/mul/div/sqrt), transcendental')
    
    args = parser.parse_args()
    
    # Calculate words needed for the precision
    # Each word provides ~18.06 decimal digits
    # Add a few extra words for safety
    mantissa_words = int(args.precision / DIGITS_PER_WORD) + 2
    words_per_number = mantissa_words + 2  # +2 for sign_length and exponent
    
    print(f"Generating test data:")
    print(f"  Precision: {args.precision} decimal digits")
    print(f"  Mantissa words: {mantissa_words}")
    print(f"  Words per number: {words_per_number}")
    
    # Generate test cases
    cases = []
    
    if args.ops in ('all', 'basic'):
        cases.extend(generate_binary_ops_tests(args.precision, words_per_number, args.num_cases))
    
    if args.ops in ('all', 'basic'):
        cases.extend(generate_unary_ops_tests(args.precision, words_per_number, args.num_cases))
    
    # Default output path
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = Path(__file__).parent / f"test_data_{args.precision}d.bin"
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_test_file(output_path, args.precision, words_per_number, cases)
    
    # Print summary by operation
    op_counts = {}
    for op, _, _ in cases:
        op_counts[op] = op_counts.get(op, 0) + 1
    
    print("\nTest cases by operation:")
    for op, count in sorted(op_counts.items()):
        print(f"  {OP_NAMES.get(op, f'op_{op}')}: {count}")


if __name__ == '__main__':
    main()
