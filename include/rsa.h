#pragma once

#include <stdint.h>
#include <math.h>
#include <limits>
#include <iostream>
#include <random>
#include <chrono> 

// InfInt library
#include "lib/InfInt.h"

struct Timer 
{
    std::chrono::time_point<std::chrono::high_resolution_clock> s;

    void set_time() { s = std::chrono::high_resolution_clock::now(); }
    void print_time() { std::cout << (std::chrono::high_resolution_clock::now() - s).count() << std::endl; }
    auto get_time() { return (std::chrono::high_resolution_clock::now() - s).count(); }
};

std::vector<int> first_primes(
{ 
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31,  37,  41,  43, 47, 53, 59, 61, 67,
    71,  73,  79,  83,  89, 97, 101, 103,
    107, 109, 113, 127, 131, 137, 139,
    149, 151, 157, 163, 167, 173, 179,
    181, 191, 193, 197, 199, 211, 223,
    227, 229, 233, 239, 241, 251, 257,
    263, 269, 271, 277, 281, 283, 293,
    307, 311, 313, 317, 331, 337, 347, 349 });

InfInt li_pow(InfInt base, InfInt exp)
{
    InfInt output = 1;
    for (InfInt i = 0; i < exp; i++)
    {
        output *= base;
    }

    return output;
}

// Generate large random numbers
InfInt large_rng(int digits)
{
    std::random_device rand;
    InfInt ret_val;

    size_t cv = 4294967295;

    for (size_t i = 0; i * 10 < digits; i++)
    {
        ret_val += ((InfInt) rand()) * li_pow(cv, i);
    } 

    return ret_val;
}

// Function to return gcd of a and b (not extended)
InfInt gcd(InfInt a, InfInt b)
{
    if (a ==  0) return b;
    return gcd(b % a, a);
}

// Function to return an lcm of a and b
InfInt lcm(InfInt a, InfInt b)
{
    return a * b / gcd(a, b);
}

// Function for extended Euclidean Algorithm
InfInt ext_gcd(InfInt a, InfInt b, InfInt* x, InfInt* y)
{
    // Base Case
    if (a ==  0) 
    {
        *x = 0, *y = 1;
        return b;
    }

    // To store results of recursive call
    InfInt x1, y1;
    InfInt gcd = ext_gcd(b % a, a, &x1, &y1);

    // Update x and y using results of recursive call
    *x = y1 - (b / a) * x1;
    *y = x1;

    return gcd;
}

// Square and multiply for modular exponentation
InfInt sq_and_mul(InfInt base, InfInt exp, InfInt mod)
{
    InfInt t = 1;
    while (exp >  0)
    {
        // for cases where exponent
        // is not an even value
        if (exp %  2 !=  0) t = (t * base) % mod;
        base = (base * base) % mod;
        exp /= 2;
    }
    return t % mod;
}

// Function to find modulo inverse of a
InfInt mod_inverse(InfInt a, InfInt m)
{
    InfInt x, y;
    InfInt g = ext_gcd(a, m, &x, &y);
    if (g !=  1) return 0;
    // m is added to handle negative x
    else return (x % m + m) % m;
}

// This function is called for all k trials. It returns
// false if n is composite and returns true if n is
// probably prime.
// d is an odd number such that  d*2^r = n-1
// for some r >= 1
bool mr_test(InfInt d, InfInt n, InfInt a = 0)
{
    // Pick a random number in [2..n-2]
    // Corner cases make sure that n > 4
    std::random_device rand;
    if (a == 0) a = (InfInt) 2 + (InfInt) large_rng(n.numberOfDigits() - rand() % (n.numberOfDigits() / 5)) % (n - 4);
 
    // Compute a^d % n
    InfInt x = sq_and_mul(a, d, n);
 
    if (x == 1 || x == n - 1) return true;
 
    // Keep squaring x while one of the following doesn't
    // happen
    // (i)   d does not reach n-1
    // (ii)  (x^2) % n is not 1
    // (iii) (x^2) % n is not n-1
    while (d != n-1)
    {
        x = (x * x) % n;
        d *= 2;
 
        if (x == 1) return false;
        if (x == n-1) return true;
    }
 
    // Return composite
    return false;
}
 
// It returns false if n is composite and returns true if n is probably prime. 
// k is an input parameter that determines accuracy level.
// Higher value of k indicates more accuracy.
bool is_prime_mr(InfInt n, int k, InfInt a = 0)
{
    // Corner cases
    if (n <= 1 || n == 4) return false;
    if (n <= 3) return true;
 
    // Find r such that n = 2^d * r + 1 for some r >= 1
    InfInt d = n - 1;
    while (d % 2 == 0) d /= 2;
 
    // Iterate given number of 'k' times
    for (int i = 0; i < k; i++)
         if (!mr_test(d, n, a))
              return false;
 
    return true;
}
 
InfInt jacobi(InfInt a, InfInt m) 
{
    // assumes a an integer and
    // m an odd positive integer
    a = a % m;
    InfInt t = 1;
    while (a != 0)
    {
        InfInt z = (m % 8 == 3 || m % 8 == 4 || m % 8 == 5) ? -1 : 1;
        while (a % 2 == 0)
        {
            a /= 2;
            t *= z;            
        }
        if (a % 4 == 3 && m % 4 == 3) t = -t;
        InfInt tmp = a;
        a = m % a;
        m = tmp;
        return m == 1 ? t : 0;
    }
}
 
std::pair<InfInt, std::pair<InfInt, InfInt>> selfridge(InfInt n)
{
    InfInt d = 5, s = 1;
    while (true)
    {
        InfInt ds = d * s;
        if (gcd(ds, n) > 1) return {ds, {0, 0}};
        if (jacobi(ds, n) == -1) return {ds, {1, ((InfInt) 1 - ds) / 4}};
        d += 2;
        s *= -1;
    }
}
 
InfInt lucasPQ(InfInt p, InfInt q, InfInt m, InfInt n) 
{
    // nth element of lucas sequence with
    // parameters p and q (mod m);
    auto half = [m](InfInt x) -> InfInt 
    {
        if (x % 2 == 1) x = x + m;
        return x / 2 % m;
    };
        
    InfInt un = 1, vn = p, qn = q;
    InfInt u = n % 2 == 0 ? 0 : 1;
    InfInt v = n % 2 == 0 ? 2 : p;
    InfInt k = n % 2 == 0 ? 1 : q;
    n /= 2;
    InfInt d = p * p - (InfInt) 4 * q;

    while (n > 0)
    {
        InfInt u2 = un * vn % m;
        InfInt v2 = (vn * vn - (InfInt) 2 * qn) % m;
        InfInt q2 = qn * qn % m;
        InfInt n2 = n / 2;
        if (n % 2 == 1)
        {
            InfInt uu = half(u * v2 + u2 * v);
            v = half(v * v2 + d * u * u2);
            u = uu;
            k *= q2;
        }
            
        un = u2;
        vn = v2;
        qn = q2;
        n = n2;
    }

    return u;
}

bool lucas_test(InfInt n)
{
    auto sfr = selfridge(n);
    if (sfr.second.first == 0) return n == sfr.first;
    return  (InfInt) 0 == lucasPQ(sfr.second.first, sfr.second.second, n, n + 1);
}

bool ballie_psw(InfInt n)
{
    return mr_test(n, 1, 2) && mr_test(n, 1, 3) && lucas_test(n);
}

// Check if a prime is low level
bool is_low_level(InfInt n)
{
    for (int divisor : first_primes)
    {
        if (n < divisor) break;
        if (n % divisor == 0) return false;
    }

    return true;
}

bool is_prime(InfInt n, bool use_mr)
{
    if (is_low_level(n))
    {
        if (use_mr) return is_prime_mr(n, 20);
        else return ballie_psw(n); //Better and hopefully faster??;
    }
    return false;
}

InfInt if_pow(InfInt base, InfInt exp)
{
    if (exp == 0) return 1;
    base *= if_pow(base, exp - 1);
    return base;
}

// Randomly generate a prime
InfInt find_prime(size_t digits, bool mr)
{
    // Generate random number between 2^(n - 1) and 2^n, with a random offset
    InfInt output = large_rng(digits);
    if ((output % 2) == 0) output += 1;
    
    while (!is_prime(output, mr)) 
    {
        output += 2;
    }
    // std::cout << "find_p" << std::endl;
    // t.print_time();

    // The loop will have increased it by 1
    return output;
}

class rsa_handler
{
public:
    // Size of the primes
    size_t digits = 0;
    // The primes
    InfInt p;
    InfInt q;
    // The product of the primes
    InfInt n;
    // lcm(p - 1, q - 1)
    InfInt tot_n;
    // EEEEE
    int e = 65537;
    // Math stuff here
    InfInt d;
public:
    rsa_handler()
    {
        digits = 150;
    }

    rsa_handler(size_t digits)
    {
        this->digits = digits;
    }

    void generate_keys(bool mr = 1)
    {
        std::random_device rand;
        // Make sure E is valid
        do 
        {
            // Generate 2 different primes for p and q
            Timer t;
            size_t digit_sub = rand() % 10;
            p = find_prime(digits - digit_sub, mr);
            q = find_prime(digits + digit_sub, mr);

            // Calculate n and tot_n
            n = p * q;
            tot_n = lcm(p - 1, q - 1);
        } while (InfInt(e) > tot_n || gcd(e, tot_n) != 1);
        d = mod_inverse(e, tot_n);
    }

    InfInt encrypt(InfInt c)
    {
        return sq_and_mul(c,  e, n);
    }

    InfInt decrypt(InfInt i)
    {
        return sq_and_mul(i, d, n);
    }
};