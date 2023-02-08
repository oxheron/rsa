#include "rsa.h"

#include <iostream>
#include <unordered_map> 

std::unordered_map<size_t, size_t> num_totime;

int main(void)
{
    rsa_handler handler(30);
    for (size_t i = 1; i < 17; i++)
    {
        Timer t;
        t.set_time();
        for (size_t j = 0; j < 20; j++)
        {
            handler.generate_keys(pow(2, i), 1);
        }
        num_totime[i] = t.get_time();
    }

    for (auto [n, t] : num_totime)
    {
        std::cout << "Num: " << n << " Time " << t << std::endl;
    }
}