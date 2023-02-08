#include "rsa.h"

#include <iostream>
#include <map> 
#include <future>
#include <thread>

int main(void)
{
    rsa_handler handler(200);

    auto fn = [handler](int offset) mutable -> size_t 
    {
        Timer t;
        t.set_time();
        for (size_t i = 0; i < 200; i++)
        {
            handler.generate_keys(offset);
        }
        return t.get_time();
    };

    auto two = std::async(std::launch::async, fn, 2);
    auto four = std::async(std::launch::async, fn, 4);
    auto eight = std::async(std::launch::async, fn, 8);
    auto sixteen = std::async(std::launch::async, fn, 16);

    two.wait();
    four.wait();
    eight.wait();
    sixteen.wait();

    std::cout << "Two " << two.get() << "\nFour " << four.get() << "\nEight " << eight.get() << "\nSixteen " << sixteen.get() << std::endl;
}
