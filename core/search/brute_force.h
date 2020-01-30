#ifndef BRUTE_FORCE
#define BRUTE_FORCE

#include <limits>
#include <tuple>

#include "../../core/metric/metric.h"

using namespace std::chrono;

namespace search {
  /**
    \brief BruteForce class for vectors
  */
  namespace vectors {

    template <typename T, typename U>
    class BruteForce {
      private:
        const uint32_t N;
        const uint16_t D;
        const std::vector<T>& feature_vector;
        const std::vector<U>& feature_vector_ids;
      public:
        /** \brief class BruteForce constructor
          @par const std::vector<T>& points - Pass by reference given points
          @par const std::vector<T>& ids - Pass by reference points' ids
          @par const int N  - Number of points
          @par const int D - Points' dimension
        */
        BruteForce(const uint32_t N, const uint16_t D,
            const std::vector<T>& points, const std::vector<U>& ids)
          : N(N), D(D), feature_vector(points), feature_vector_ids(ids) {};
        /**
          \brief class BruteForce default construct
        */
        ~BruteForce() = default;
        /** \brief Executes exact Nearest tNeighbor
          @par const std::vector<T>& query_points - Pass by reference query points
          @par const int idx - idx to get correspodent point
        */
        std::tuple<T,U,double> NearestNeighbor(const std::vector<T>& query_points,
          const int idx) {

          auto start = high_resolution_clock::now();
          /* Initialize min_dist to max value of type T */
          T min_dist = std::numeric_limits<T>::max();
          /* Initialize correspodent min_id using the C++11 way */
          U min_id{};
          /* Run NearestNeighbor for all points in the dataset */
          for (size_t i = 0; i < N; ++i) {
            T dist = metric::ManhattanDistance<T>(
              std::next(feature_vector.begin(), i * D),
              std::next(query_points.begin(), idx * D),
              std::next(query_points.begin(), idx * D + D));
            if (dist < min_dist) {
              min_dist = dist;
              min_id = feature_vector_ids[i];
            }
          }
          auto stop = high_resolution_clock::now();
          duration <double> total_time = duration_cast<duration<double>>(stop - start);
          /* return result as a tuple of min_dist, min_id and total_time */
          return std::make_tuple(min_dist,min_id,total_time.count());
        };
        /** \brief Executes (r,c)-Nearest tNeighbor
          @par const std::vector<T> &query_points - Pass by reference query points
          @par const int idx - idx to get correspodent point
        */
        std::vector<std::pair<T,U>> RadiusNearestNeighbor(
          const std::vector<T>& query_points,
          const int idx, const double R) {

          /* Define result vector */
          std::vector<std::pair<T,U>> result;
          /* Run (r,c)-NearestNeighbor for all points in the dataset */
          for (size_t i = 0; i < N; ++i) {
            T dist = metric::ManhattanDistance<T>(
              std::next(feature_vector.begin(), i * D),
              std::next(query_points.begin(), idx * D),
              std::next(query_points.begin(), idx * D + D));
            if (dist <= R) {
              result.push_back(std::make_pair(dist,feature_vector_ids[i]));
            }
          }
          return result;
        };
    };
  }
  /**
    \brief BruteForce class for curves
  */
  namespace curves {

    template <typename T, typename U>
    class BruteForce {
      private:
        const std::vector<std::pair<T,T>>& input_curves;
        const std::vector<U>& input_curves_ids;
        const std::vector<int>& input_curves_lengths;
        const std::vector<int>& input_curves_offsets;
      public:
        /** \brief class BruteForce constructor
          @par const std::vector<T>& curves - Pass by reference given curves
          @par const std::vector<U>& ids - Pass by reference curves' ids
          @par const std::vector<int>& lengths - Pass by reference curves' lengths
          @par const std::vector<int>& curves_offsets - Pass by reference offsets
            to given curves to have access on them
        */
        BruteForce(const std::vector<std::pair<T,T>>& curves,
          const std::vector<U>& ids, const std::vector<int>& lengths,
          const std::vector<int>& offsets)
          : input_curves(curves), input_curves_ids(ids),
            input_curves_lengths(lengths), input_curves_offsets(offsets) {};
        /**
          \brief class BruteForce default construct
        */
        ~BruteForce() = default;
        /** \brief Executes exact Nearest Neighbor for curves
          @par const std::vector<std::pair<T,T>>& query_curves - Pass by reference
            query curves
          @par const std::vector<int>& query_lengths - Pass by reference query
            curves' lengths
          @par const std::vector<int> query_offsets - Pass by reference offsets
            to given query curves to have access on them
          @par const int idx - Index to current query curve
        */
        std::tuple<T,U,double> NearestNeighbor(
          const std::vector<std::pair<T,T>>& query_curves,
          const std::vector<int>& query_lengths,
          const std::vector<int>& query_offsets, const int idx) {

          auto start = high_resolution_clock::now();
          /* Initialize min_dist to max value of type T */
          T min_dist = std::numeric_limits<T>::max();
          /* Initialize correspodent min_id using the C++11 way */
          U min_id{};
          /* Run NearestNeighbor for all points in the input */
          for (size_t i = 0; i < input_curves_ids.size(); ++i) {
            T dist = metric::DTWDistance<T>(
              std::next(input_curves.begin(),input_curves_offsets[i]),
              std::next(input_curves.begin(),
                        input_curves_offsets[i] + input_curves_lengths[i]),
              std::next(query_curves.begin(),query_offsets[idx]),
              std::next(query_curves.begin(),
                        query_offsets[idx] + query_lengths[idx]));
            if (dist < min_dist) {
              min_dist = dist;
              min_id = input_curves_ids[i];
            }
          }
          auto stop = high_resolution_clock::now();
          duration <double> total_time = duration_cast<duration<double>>(stop - start);
          /* return result as a tuple of min_dist, min_id and total_time */
          return std::make_tuple(min_dist,min_id,total_time.count());
        };
    };
  }
}

#endif
