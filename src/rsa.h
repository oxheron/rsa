#pragma once

#include <stdint.h>
#include <math.h>
#include <limits>
#include <iostream>
#include <random>

// InfInt library
#include "../lib/InfInt.h"

std::array<int> small_primes({
    
});

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
bool miiller_test(InfInt d, InfInt n)
{
    // Pick a random number in [2..n-2]
    // Corner cases make sure that n > 4
    std::random_device rand;
    InfInt a = 2 + rand() % (n - 4);
 
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
bool is_prime_miiller(InfInt n, int k)
{
    // Corner cases
    if (n <= 1 || n == 4)  return false;
    if (n <= 3) return true;
 
    // Find r such that n = 2^d * r + 1 for some r >= 1
    int d = n - 1;
    while (d % 2 == 0) d /= 2;
 
    // Iterate given number of 'k' times
    for (int i = 0; i < k; i++)
         if (!miillerTest(d, n))
              return false;
 
    return true;
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

template <typename T>
class rsa_handler
{
public:
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
    void generate_keys()
    {
        // Make sure E is valid
        do 
        {
            // Generate 2 different primes for p and q
            p = find_prime();
            q = find_prime();

            std::cout << "p: " << p << "q: " << q << std::endl;

            // Test 
            while (q == p) q = find_prime();

            // Calculate n and tot_n
            n = p * q;
            tot_n = lcm(p - 1, q - 1);
        } while (InfInt(e) > tot_n || gcd(e, tot_n) !=  1);

        std::cout << "here" << std::endl;

        d = mod_inverse(e, tot_n);
    }

    InfInt encrypt(InfInt c)
    {
        return sq_and_mul(c,  e, n);
    }

    // Not cryptographically safe
    InfInt decrypt(InfInt i)
    {
        return sq_and_mul(i, d, n);
    }
private: 
    // Randomly generate a prime
    InfInt find_prime()
    {
        // Start output from 2 to the power of half of the bits in T
        std::numeric_limits<T> limit;
        InfInt output = limit.max() / 8;

        // Generate random number between 1 and 10000
        std::random_device rand;
        short num_primes = rand() % (10 * (int) log(limit.max())) + 1;
        
        // Loop through, generating as many primes as the random number
        for (InfInt i = 0; i < 5; i++)
        {
            while (!is_prime(output)) 
            {
                std::cout << "entered loop" << std::endl;
                output += 1;
            }
            output += 1;
        }
            std::cout << "output" << std::endl;


        std::cout << is_prime(output) << std::endl;

        // The loop will have increased it by 1
        return output - 1;
    }
};