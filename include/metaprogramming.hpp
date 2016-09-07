/*!
  \file metaprogramming.hpp

  \brief Metapgrogramming bits
*/

#ifndef FWDPY_METAPROGRAMMING_HPP
#define FWDPY_METAPROGRAMMING_HPP

#include <future>
#include <thread>
#include <type_traits>
#include <vector>
namespace fwdpy
{
    namespace meta
    {
        template <typename T>
        inline typename std::enable_if<std::is_void<T>::value, void>::type
        return_sampler_futures(std::vector<std::future<T>> &futures)
        {
            for (std::size_t i = 0; i < futures.size(); ++i)
                futures[i].get();
            return;
        }

        template <typename T>
        inline typename std::enable_if<!std::is_void<T>::value,
                                       std::vector<T>>::type
        return_sampler_futures(std::vector<std::future<T>> &futures)
        {
            std::vector<T> rv(futures.size());
            for (std::size_t i = 0; i < futures.size(); ++i)
                {
                    rv[i] = futures[i].get();
                }
            return rv;
        }
    }
}

#endif
