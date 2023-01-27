#pragma once

#include <stdint.h>
#include <math.h>
#include <limits>
#include <iostream>
#include <random>

// Function for extended Euclidean Algorithm
long long ext_gcd(long long a, long long b, long long* x, long long* y)
{
    // Base Case
    if (a == 0) 
    {
        *x = 0, *y = 1;
        return b;
    }

    // To store results of recursive call
    long long x1, y1;
    long long gcd = ext_gcd(b % a, a, &x1, &y1);

    // Update x and y using results of recursive
    // call
    *x = y1 - (b / a) * x1;
    *y = x1;

    return gcd;
}

// Function to find modulo inverse of a
size_t mod_inverse(size_t a, long long m)
{
    long long x, y;
    size_t g = ext_gcd(a, m, &x, &y);
    if (g != 1) return 0;
    // m is added to handle negative x
    else return (x % m + m) % m;
}

// Check if a number is prime
bool is_prime(size_t n)
{
    if (n == 2 || n == 3) return true;
    if (n <= 1 || n % 2 == 0 || n % 3 == 0) return false;
    for (int i = 5; i * i <= n; i += 6) if (n % i == 0 || n % (i + 2) == 0) return false;
    return true;
}

// Function to return gcd of a and b (not extended)
size_t gcd(size_t a, size_t b)
{
    if (a == 0) return b;
    return gcd(b % a, a);
}

// Function to return an lcm of a and b
size_t lcm(size_t a, size_t b)
{
    return a * b / gcd(a, b);
}

template <typename T>
class rsa_handler
{
public:
    // The primes
    T p;
    T q;
    // The product of the primes
    size_t n;
    // lcm(p - 1, q - 1)
    size_t tot_n;
    // EEEEE
    // uint32_t e = 65537;
    uint32_t e = 17;
    // Math stuff here
    size_t d;
public:
    void generate_keys()
    {
        // Make sure E is valid
        do 
        {
            // Generate 2 different primes for p and q
            p = find_prime();
            q = find_prime();
            // Test 
            while (q == p) q = find_prime();

            // Calculate n and tot_n
            n = p * q;
            tot_n = lcm(p - 1, q - 1);
        } while (e > tot_n || gcd(e, tot_n) != 1);

        d = mod_inverse(e, tot_n);
    }

private: 
    // Randomly generate a prime
    T find_prime()
    {
        // Start output from 2 to the power of half of the bits in T
        std::numeric_limits<T> limit;
        T output = limit.max() / 8;

        // Generate random number between 1 and 10000
        std::random_device rand;
        short num_primes = rand() % (10 * (int) log(limit.max())) + 1;
        
        // Loop through, generating as many primes as the random number
        for (size_t i = 0; i < num_primes; i++)
        {
            while (!is_prime(output)) output += 1;
            output += 1;
        }

        // The loop will have increased it by 1
        return output - 1;
    }
};