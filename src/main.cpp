#include "rsa.h"

#include <iostream>

int main(void)
{
    rsa_handler<short> handler;
    handler.generate_keys();
    std::cout << "p: " << handler.p << "\nq: " << handler.q << "\nn: " << handler.n << "\ntot_n: " << handler.tot_n << std::endl;
    std::cout << "d: " << handler.d << std::endl;
}