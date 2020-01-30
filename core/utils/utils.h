#ifndef UTILS
#define UTILS

#include <cmath>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <queue>

namespace utils {
  /** \brief Computes modular exponentiation
    @par base - base
    @par exp - expotent
    @par mod - modulo divisor
  */
  inline uint64_t mod_exp(uint32_t base, uint16_t exp, uint32_t mod) {
    uint64_t res = 1;
    while (exp > 0) {
      if (exp % 2 == 1)
        res = (res * base) % mod;
      exp = exp >> 1;
      base = (base * base) % mod;
    }
    return res;
  };
  /**
    \brief Split a string given a delimiter
  */
  inline std::vector<std::string> Split(const std::string& s, char delimiter) {
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter)) {
      tokens.push_back(token);
   }
   return tokens;
}
  /** \brief Computes both negative and positive modulos
    @par int a - Dividend
    @par int b - modulo divisor
  */
  constexpr int mod(int a, int b) { return (a % b + b) % b; }

  /**
    \brief Compute mean of the points given in R^d
  */
  template <typename T>
  double ComputeMean(const std::vector<T>& points, const uint16_t& points_dim,
    const uint32_t& no_points) {

    /* Compute mean point coordinate by coordinate */
    std::vector<double> mean_point(points_dim);
    for (size_t i = 0; i < points_dim; ++i) {
      for (size_t j = 0; j < no_points; ++j) {
        mean_point[i] += points[j * points_dim + i];
      }
      mean_point[i] /= no_points;
    }
    /* Compute the mean coordinate */
    double mean = 0.0;
    for (size_t i = 0; i < points_dim; ++i) {
      mean += mean_point[i];
    }
    mean /= points_dim;
    // Return result
    return mean;
  }
  /** \brief Compute delta as the average of the euclidian distance of
    concecutive points for all curves
  */
  template <typename T>
  double ComputeDelta(const std::vector<std::pair<T,T>>& curves,
    const std::vector<int>& lengths, const std::vector<int>& offsets) {

    // Get number of curves in the input file
    size_t N = lengths.size();
    // vector to store average euclidian distance of points for each curve
    std::vector<double> avg_eucl_dists_of_curves_points(N);
    // Repeat for all curves in the dataset
    for (size_t i = 0; i < N; ++i) {
      T sum_eucl_dists{};
      for (size_t j = 0; j < lengths[i] - 1; ++j) {
        std::pair<T,T> p_1 = curves[offsets[i] + j];
        std::pair<T,T> p_2 = curves[offsets[i] + j + 1];
        T x_diff = std::abs((std::get<0>(p_1)-std::get<0>(p_2)));
        T y_diff = std::abs((std::get<1>(p_1)-std::get<1>(p_2)));
        sum_eucl_dists += x_diff + y_diff;
      }
      avg_eucl_dists_of_curves_points[i] = (double) sum_eucl_dists / lengths[i];
    }
    // calculate total sum to average with number of curves
    double total_sum{};
    for (size_t i = 0; i < N; ++i) {
      total_sum += avg_eucl_dists_of_curves_points[i];
    }
    // return averaged result
    return total_sum / N;
  }
  /** \brief Compute parameter R as the average of the distances of each curve
    to its nearest neighbor
  */
  template <typename T, typename U>
  double ComputeParameterR(const std::vector<std::tuple<T,U,double>>& exact) {

    /* Get number of points */
    int N = exact.size();
    /* Sum up distances from the nearest neighbor */
    double distance_to_nn{};
    for (size_t i = 0; i < N; ++i) {
      distance_to_nn += std::get<0>(exact[i]);
    }
    /* Return its average */
    return distance_to_nn / N;
  }
  /**
    \brief Gets a std::string to convert in the specified type T
  */
  template <typename T>
  T convert_to (const std::string& str) {
    std::istringstream ss(str);
    T num;
    ss >> num;
    return num;
  }
  /**
    \brief Variadic templated min function
  */
  template <typename T>
  T min(T&& t) {
    return std::forward<T>(t);
  }

  template <typename T0, typename T1, typename... Ts>
  typename std::common_type<T0, T1, Ts... >
    ::type min(T0&& val1, T1&& val2, Ts&&... vs) {
    if (val2 < val1) {
      return min(val2, std::forward<Ts>(vs)...);
    } else {
      return min(val1, std::forward<Ts>(vs)...);
    }
  }
  /**
    \brief Variadic templated max function
  */
  template <typename T>
  T max(T&& t) {
    return std::forward<T>(t);
  }

  template <typename T0, typename T1, typename... Ts>
  typename std::common_type<T0, T1, Ts... >
    ::type max(T0&& val1, T1&& val2, Ts&&... vs) {
    if (val2 > val1) {
      return max(val2, std::forward<Ts>(vs)...);
    } else {
      return max(val1, std::forward<Ts>(vs)...);
    }
  }
  /** \brief Returns a vector of k indices with smallest values
    @par[in] const std::vector<T>& input - Given array of values
    @par[in] k - number of indices we need
  */
  template <typename T>
  std::vector<int> ArgMin(const std::vector<T>& input, const int& k = 1) {
    // Define a priority_queue data structure
    std::priority_queue<std::pair<T,int>, std::vector<std::pair<T,int>>,
                        std::greater <std::pair<T,int>>> q;
    for (size_t i = 0; i < input.size(); ++i) {
      if(q.size() < k) {
        q.push(std::pair<T,int>(input[i], i));
      }
      else if(q.top().first < input[i]) {
        q.pop();
        q.push(std::pair<T,int>(input[i], i));
      }
    }
    // Define result vector to be returned
    std::vector<int> result(k);
    for (size_t i = 0; i < k; ++i) {
      result[k - i - 1] = q.top().second;
      q.pop();
    }
    return result;
  }
}

#endif
