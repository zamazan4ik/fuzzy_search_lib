#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>

#include "Randl_fuzzy_search.hpp"
#include "Bitap.hpp"
#include "libflasm.h"


int main()
{
    const std::string data = "i live in fer Mensk and it's really beautiful place for living.";
    const std::string search = "Minsk";

auto result2 = fuzzy_search(data.begin(), data.end(), search.begin(), search.end(), 1, true);

    auto result = bitap_fuzzy_bitwise_search_new(data, search, 1);

    for(; result.first != result.second; result.first++)
    {
        std::cout << *result.first;
    }


    auto start = std::chrono::system_clock::now();

    const int times = 1000;
    for(size_t i = 0; i < times; ++i)
    auto result3 = libflasm::flasm_ed(reinterpret_cast<unsigned char*>(const_cast<char*>(data.c_str())), data.length(),
                                      reinterpret_cast<unsigned char*>(const_cast<char*>(search.c_str())), search.length(),
                                      search.length(), 1, false);
    //auto result = bitap_fuzzy_bitwise_search_new(data, search, 1);

    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - start).count() << std::endl;


libflasm::ResultTupleSetIterator it;
    for ( it = result3.begin(); it != result3.end(); ++it )
    {
        libflasm::ResultTuple res = *it;
        //fprintf( out_fd, "(%u,%u,%u)\n", res.pos_t, res.pos_x, res.error );
        //std::cout << res.pos_t << " " << res.pos_x << " " << res.error << std::endl;
        std::cout << data.substr(res.pos_t - search.length() + 1, search.length()) << std::endl;
    }

    std::cout << std::endl;
    return 0;
}

