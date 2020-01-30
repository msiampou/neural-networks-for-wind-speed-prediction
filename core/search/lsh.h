#ifndef LSH_H_
#define LSH_H_

#include <iterator>
#include <random>
#include <tuple>
#include <unordered_map>

#include "../../core/hash/hash_function.h"
#include "../../core/metric/metric.h"
#include "../../core/utils/utils.h"

using namespace std::chrono;

namespace search {
  /** \brief General LSH parameters
    @par table_size - a virtual size of unordered_map to have more colissions
      in each bucket of LSH
    @par m - parameter used by Hash Functions
    @par M - parameter used by Hash Functions
    @par w - window size used by Hash Functions
    @par K - number of h hash functions
    @par L - number of LSH hashtables and g Amplified Hash Functions
    @par D - points' dimension
    @par N - number of dataset points
    @par R - average of exact NN distances calculated using brute force search
  */
  /**
    \brief LSH class for vectors
  */
  namespace vectors {

    template <typename T, typename U>
    class LSH {
      private:
        uint32_t table_size;
        uint32_t m;
        uint32_t M;
        double w;
        const uint8_t K;
        const uint8_t L;
        const uint16_t D;
        const uint32_t N;
        const double R;

        const std::vector<T>& feature_vector;
        const std::vector<U>& feature_vector_ids;

        std::vector<hash::AmplifiedHashFunction<T>> hash_functions;
        std::vector<std::unordered_map<int,std::vector<int>>> hash_tables;
      public:
        /**
          \brief class LSH constructor
        */
        LSH(const uint8_t K, const uint8_t L, const uint16_t D,
            const uint32_t N, const double r, const std::vector<T>& points,
            const std::vector<U>& ids) : K(K), L(L), D(D), N(N), R(r),
            feature_vector(points), feature_vector_ids(ids) {

            w = 2 * r;
            m = (1ULL << 32) - 5;
            M = 1ULL << (32 / K);
            table_size = N / 32;
            // Preprocess step
            // 1) Randomly select L amplified hash functions g1 , . . . , gL .
            for (size_t i = 0; i < L; ++i) {
              hash_functions.push_back(hash::AmplifiedHashFunction<T>(K,D,m,M,w));
            }
            // 2) Initialize L hash-tables, hash all points to all tables using g
            for (size_t i = 0; i < L; ++i) {
              std::unordered_map<int,std::vector<int>> ht;
              for (size_t j = 0; j < N; ++j) {
                ht[hash_functions[i].Hash(feature_vector,j) % table_size].push_back(j);
              }
              hash_tables.push_back(ht);
            }
        };
        /**
          \brief class LSH default destructor
        */
        ~LSH() = default;

        /** \brief Executes approximate Nearest tNeighbor
          @par const std::vector<T>& query_points - Pass by reference query points
          @par const int offset - Offset to get correspodent point
        */
        std::tuple<T,U,double> NearestNeighbor(const std::vector<T>& query_points,
          const int offset) {

          auto start = high_resolution_clock::now();
          /* Initialize min_dist to max value of type T */
          T min_dist = std::numeric_limits<T>::max();
          /* Initialize correspodent min_id using the C++11 way */
          U min_id{};
          for (size_t i = 0; i < L; ++i) {
            // get i_th hashtable
            std::unordered_map<int,std::vector<int>>& ht_i = hash_tables[i];
            // get all points in the same bucket
            std::vector<int>& bucket = ht_i[hash_functions[i].Hash(query_points,offset) % table_size];
            // iterate over all points in the buck
            for (auto const& fv_offset: bucket) {
              T dist = metric::ManhattanDistance<T>(
                std::next(feature_vector.begin(), fv_offset * D),
                std::next(query_points.begin(), offset * D),
                std::next(query_points.begin(), offset * D + D));
              if (dist < min_dist) {
                min_dist = dist;
                min_id = feature_vector_ids[fv_offset];
              }
            }
          }
          auto stop = high_resolution_clock::now();
          duration <double> total_time = duration_cast<duration<double>>(stop - start);
          /* return result as a tuple of min_dist, min_id and total_time */
          return std::make_tuple(min_dist,min_id,total_time.count());
        };

        /** \brief Executes approximate Radius Nearest tNeighbor
          @par const std::vector<T>& query_points - Pass by reference query points
          @par const int offset - Offset to get correspodent point
        */
        std::vector<std::pair<T,U>> RadiusNearestNeighbor(
          const std::vector<T>& query_points, const int offset,
          const int radius) {

          /* Define result as a vector of pairs of min_dist and min_id */
          std::vector<std::pair<T,U>> result;
          /* Initialize min_dist to max value of type T */
          T min_dist = std::numeric_limits<T>::max();
          /* Initialize correspodent min_id using the C++11 way */
          U min_id{};
          for (size_t i = 0; i < L; ++i) {
            // get i_th hashtable
            std::unordered_map<int,std::vector<int>>& ht_i = hash_tables[i];
            // get all points in the same bucket
            std::vector<int> &bucket = ht_i[hash_functions[i].Hash(query_points,offset) % table_size];
            // iterate over all points in the buck
            for (auto const& fv_offset: bucket) {
              T dist = metric::ManhattanDistance<T>(
                std::next(feature_vector.begin(), fv_offset * D),
                std::next(query_points.begin(), offset * D),
                std::next(query_points.begin(), offset * D + D));
              if (dist <= radius) {
                result.push_back(std::make_pair(dist,feature_vector_ids[fv_offset]));
              }
            }
          }
          return result;
        };
    };
  }
  /**
    \brief LSH class for curves
  */
  namespace curves {

    template <typename T, typename U>
    class LSH {
      private:
        uint32_t table_size;
        uint32_t m;
        uint32_t M;
        double w;
        const uint8_t K;
        const uint8_t L;
        const uint16_t D;
        const uint32_t N;
        const double R;

        const std::vector<double> &feature_vector;
        const std::vector<std::pair<T,T>>& input_curves;
        const std::vector<U> &input_curves_ids;
        const std::vector<int>& input_curves_lengths;
        const std::vector<int>& input_curves_offsets;

        std::vector<hash::AmplifiedHashFunction<double>> hash_functions;
        std::vector<std::unordered_map<int,std::vector<int>>> hash_tables;
      public:
        /** \brief class LSH constructor
        */
        LSH(const uint8_t K, const uint8_t L, const uint16_t D, const uint32_t N,
            const double r, const std::vector<std::pair<T,T>>& curves,
            const std::vector<U> &ids, const std::vector<int>& lengths,
            const std::vector<int>& offsets, const std::vector<T> &points) :
            K(K), L(L), D(D), N(N), R(r), input_curves(curves),
            input_curves_ids(ids), input_curves_lengths(lengths),
            input_curves_offsets(offsets), feature_vector(points) {

            w = 40 * R;
            m = (1ULL << 32) - 5;
            M = pow(2, 32 / K);
            table_size = N / 8;
            // Preprocess step
            // 1) Randomly select L amplified hash functions g1 , . . . , gL .
            for (size_t i = 0; i < L; ++i) {
              hash_functions.push_back(hash::AmplifiedHashFunction<T>(K,D,m,M,w));
            }
            // 2) Initialize L hash-tables, hash all points to all tables using g
            for (size_t i = 0; i < L; ++i) {
              std::unordered_map<int,std::vector<int>> ht;
              for (size_t j = 0; j < N; ++j) {
                ht[hash_functions[i].Hash(feature_vector,j) % table_size].push_back(j);
              }
              hash_tables.push_back(ht);
            }
        };
        /**
          \brief class LSH default construct
        */
        ~LSH() = default;
        /** \brief Executes approximate Nearest tNeighbor
          @par const std::vector<T>& query_points - Pass by reference query points
          @par const int offset - Offset to get correspodent point
        */
        std::pair<T,U> NearestNeighbor(const std::vector<T>& query_points,
          const int offset, const std::vector<std::pair<T,T>>& query_curves,
          const std::vector<int>& query_curves_lengths,
          const std::vector<int>& query_curves_offsets) {

          /* Initialize min_dist to max value of type T */
          T min_dist = std::numeric_limits<T>::max();
          /* Initialize correspodent min_id using the C++11 way */
          U min_id{};
          for (size_t i = 0; i < L; ++i) {
            // get i_th hashtable
            std::unordered_map<int,std::vector<int>> &ht_i = hash_tables[i];
            // get all curves in the same bucket
            std::vector<int> &bucket = ht_i[hash_functions[i].Hash(query_points,offset) % table_size];
            //iterate over all curves in the bucket
            for (auto const& fv_offset: bucket) {
              T dist =  metric::DTWDistance<T> (
                std::next(input_curves.begin(),input_curves_offsets[fv_offset]),
                std::next(input_curves.begin(),
                          input_curves_offsets[fv_offset] + input_curves_lengths[fv_offset]),
                std::next(query_curves.begin(),query_curves_offsets[offset]),
                std::next(query_curves.begin(),
                          query_curves_offsets[offset] + query_curves_lengths[offset])
              );
              if (dist < min_dist) {
                min_dist = dist;
                min_id = input_curves_ids[fv_offset];
              }
            }
          }
          /* return result as a tuple of min_dist and min_id */
          return std::make_pair(min_dist,min_id);
        };

        std::pair<T,U> NearestNeighbor(const std::vector<T>& query_points,
          const int offset, const std::vector<std::pair<T,T>>& query_curves,
          const std::vector<int>& query_curves_lengths,
          const std::vector<int>& query_curves_offsets,
          const U idx) {

          /* Initialize min_dist to max value of type T */
          T min_dist = std::numeric_limits<T>::max();
          /* Initialize correspodent min_id using the C++11 way */
          U min_id{};
          for (size_t i = 0; i < L; ++i) {
            // get i_th hashtable
            std::unordered_map<int,std::vector<int>> &ht_i = hash_tables[i];
            // get all curves in the same bucket
            std::vector<int> &bucket = ht_i[hash_functions[i].Hash(query_points,offset) % table_size];
            //iterate over all curves in the bucket
            //auto const& fv_offset:bucket[0];
            for (auto const& fv_offset: bucket) {
              int of = fv_offset % input_curves_offsets.size(); //fv_offset % input_curves_offsets.size();
              //std::cout << fv_offset << std::endl;
              T dist =  metric::DTWDistance<T> (
                std::next(input_curves.begin(),input_curves_offsets[of]),
                std::next(input_curves.begin(),
                          input_curves_offsets[of] + input_curves_lengths[of]),
                std::next(query_curves.begin(),query_curves_offsets[offset]),
                std::next(query_curves.begin(),
                          query_curves_offsets[offset] + query_curves_lengths[offset])
              );
              if (dist < min_dist) {
                min_dist = dist;
                min_id = input_curves_ids[of];
              }
            }
          }
          /* return result as a tuple of min_dist and min_id */
          return std::make_pair(min_dist,min_id);
        };
    };
  }
}

#endif
