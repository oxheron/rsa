#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <math.h>

template <typename T>
class rsa_handler
{
private:
    // The primes
    T p;
    T q;
    // The product of the primes
    T n;
    // lcm(p - 1, q - 1)
    T tot_n;
    // EEEEE
    uint32_t e = 65537;
    // Math stuff here
    T d;
public:
    void generate_keys()
    {
        // Generate 2 different primes for p and q
        p = find_prime();
        q = find_prime();
        while (q == p) q = find_prime();

        // Calculate n and tot_n
        n = p * q;
        T tot_n = lcm(p - 1, q - 1);
    }

private: 
    // Randomly generate a prime
    T find_prime()
    {
        // Start output from 2 to the power of half of the bits in T
        T output = pow(2, sizeof(T) / 2);
        // Generate random number between 1 and 10000
        int random = rand() % 9999 + 1;
        
        // Loop through, generating as many primes as the random number
        for (size_t i = 0; i < random; i++)
        {
            while (!is_prime(output)) output += 1;
            output += 1;
        }

        // The loop will have increased it by 1
        return output - 1;
    }

    // Check if a number is prime
    bool is_prime(T n)
    {
        if (n == 2 || n == 3) return true;
        if (n <= 1 || n % 2 == 0 || n % 3 == 0) return false;
        for (int i = 5; i * i <= n; i += 6) if (n % i == 0 || n % (i + 2) == 0) return false;
        return true;
    }

    // Function to return gcd of a and b
    T gcd(T a, T b)
    {
        if (a == 0) return b;
        return gcd(b % a, a);
    }

    // Function to return an lcm of a and b
    T lcm(T a, T b)
    {
        return abs(a * b) / gcd(a, b);
    }
};