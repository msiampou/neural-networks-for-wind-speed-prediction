#ifndef HASH_FUNCTION
#define HASH_FUNCTION

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iterator>
#include <random>
#include <stdlib.h>
#include <string>
#include <vector>

#include "../../core/utils/utils.h"

namespace hash {

  template <typename T>
  class HashFunction {
    private:
      std::default_random_engine generator;
      std::uniform_real_distribution<double> distribution;
      const uint16_t D;
      const uint32_t m;
      const uint32_t M;
      const double w;
      std::vector<double> s;
      std::vector<int> a;
    public:
      /** \brief HashFunction class constructor
        This class illustrates the following hash function:
        h(x)= (a_d−1 + m*a_d−2 +···+ m^(d−1)*a_0) modM , m > max a_i
        where a_i = floor((x_i - s_i) / w)
        @par D - Space dimension
        @par m - parameter m in the hash function
        @par M - parameter M in the hash function
        @par w - window size
      */
      HashFunction(const uint16_t D, const uint32_t m, const uint32_t M,
        const double w): D(D), m(m), M(M), w(w), distribution(0,w), s(D), a(D),
        generator(std::chrono::system_clock::now().time_since_epoch().count()) {

        /* Initialize s vector of dimension D using uniform_real_distribution */
        for (size_t i = 0; i < D; ++i) {
          s[i] = distribution(generator);
        }
      };
      /**
        \brief HashFunction class default destructor
      */
      ~HashFunction() = default;
      /*
      /** \brief Hash point as follows:
        1) Compute a_i = floor((x_i - s_i) / w) for i = 0...D-1
        2) Compute h(x) = (a_d−1 + m*a_d−2 +···+ m^(d−1)*a_0) modM
      */
      uint32_t Hash(const std::vector<T> &points, int offset) {
        uint32_t hash_value{};
        /* Computing a_i */
        for (size_t i = 0; i < D; ++i) {
          a[i] = floor((points[offset * D + i] - s[i]) / w);
        }
        /* Reverse vector a */
        std::reverse(a.begin(),a.end());
        /* Computing h(x) */
        for (size_t i = 0; i < D; ++i) {
          hash_value += (utils::mod(a[i],M) * utils::mod_exp(m,i,M)) % M;
        }
        return hash_value % M;
      };
  };

  template <typename T>
  class AmplifiedHashFunction {
    private:
      std::vector<HashFunction<T>> h;
      const uint8_t K;
      const uint16_t D;
      const uint32_t m;
      const uint32_t M;
      const double w;
    public:
      /** \brief AmplifiedHashFunction class constructor
        This class illustrates g(x) = [h1(x)|h2(x)| · · · |hk (x)].
        @par K - Number of Hash Functions selected uniformly
        @par D - Space dimension
        @par m - parameter m in the hash function
        @par M - parameter M in the hash function
        @par w - window size
      */
      AmplifiedHashFunction(const uint8_t K, const uint16_t D, const uint32_t m,
        const uint32_t M, const double w): K(K), D(D), m(m), M(M), w(w) {
          /* Select uniformly K hash functions */
          for (size_t i = 0; i < K; ++i) {
            h.push_back(HashFunction<T>(D, m, M, w));
          }
        }
      /**
        \brief AmplifiedHashFunction class default destructor
      */
      ~AmplifiedHashFunction() = default;
      /** Hash point as follows:
        1) Hashing using h_i for i = 1..K
        2) Concat h_i and modulo with table_size
      */
      uint64_t Hash(const std::vector<T> &points, int offset) {
        std::string str_value{};
        for (size_t i = 0; i < K; ++i) {
          str_value += std::to_string(h[i].Hash(points,offset));
        }
        // convert str_value to uint64_t
        char *p_end;
        uint64_t hash_value = strtoull(str_value.c_str(), &p_end, 10);
        return hash_value;
      }
  };

}


#endif
