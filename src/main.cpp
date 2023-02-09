#include "rsa.h"

#include <iostream>
#include <map> 
#include <future>
#include <thread>

int main(void)
{
    rsa_handler handler(300);
    handler.generate_keys();

    std::cout << "handler.d: " << handler.d << std::endl;
}
