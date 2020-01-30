#ifndef INITIALIZATION
#define INITIALIZATION

#include <algorithm>
#include <tuple>
#include <random>

#include "../../../core/metric/metric.h"
#include "../../../core/utils/utils.h"

namespace cluster {
  namespace initialization {
    /* Initialization Algorithms for vectors */
    namespace vectors {
      /** \brief Initilization of cluster centroids
        Details:
          (1) Calculate symmetric n × n distance matrix of all objects,
          i.e. all distances d_ij from every object i = 1,...,n to every other
          object j = 1,...,n,i != j.
          (2) For object i compute:
            v_i = sum from j = 1 to n (d_ij / sum from t = 1 to n d_jt)
            for i = 1,...,n
          (3) Return the k objects with k smallest vi values.
          Algorithm proposed in [Park-Jun’09]
        @par[in] const std::vector<T>& dataset_vectors : vectors given from dataset
        @par[in] const int& no_vectors : number of vectors
        @par[in] const int& vectors_dim : vectors' dimensions (all vectors are
                                          dimensionally equal)
        @par[in] const int& no_clusters : number of clusters
      */
      template <typename T>
      std::vector<T> ParkJunInit(const std::vector<T>& dataset_vectors,
        const int& no_vectors, const int& vectors_dim, const int& no_clusters) {
        // Declare 2D d array which holds all distances between the vectors
        T** d_array = new T*[no_vectors];
        for(int i = 0; i < no_vectors; ++i)
          d_array[i] = new T[no_vectors];
        // Calculate all distances using ManhattanDistance
        for (size_t i = 0; i < no_vectors; ++i) {
          for (size_t j = 0; j < no_vectors; ++j) {
            if (i == j) {
              d_array[i][j] = {};
            } else {
              d_array[i][j] = metric::ManhattanDistance<T>(
                std::next(dataset_vectors.begin(), i * vectors_dim),
                std::next(dataset_vectors.begin(), j * vectors_dim),
                std::next(dataset_vectors.begin(), j * vectors_dim + vectors_dim));
            }
          }
        }
        // Declare 1D v array and computes its values as the algorithm proposes
        std::vector<T> v(no_vectors);
        for (size_t i = 0; i < no_vectors; ++i) {
          T o_sum{};
          for (size_t j = 0; j < no_vectors; ++j) {
            T i_sum{};
            for (size_t t = 0; t < no_vectors; ++t) {
              i_sum += d_array[j][t];
            }
            o_sum += d_array[i][j] / i_sum;
          }
          v[i] = o_sum;
        }
        // Find indices of no_clusters smallest v_array values
        std::vector<int> idx = utils::ArgMin<T>(v, no_clusters);
        // Initialize vector centroids
        std::vector<T> centroids(no_clusters * vectors_dim);
        for (size_t i = 0; i < no_clusters; ++i) {
          for (size_t j = 0; j < vectors_dim; ++j) {
            centroids[i * vectors_dim + j] = dataset_vectors[idx[i] * vectors_dim + j];
          }
        }
        for(int i = 0; i < no_vectors; ++i)
          delete[] d_array[i];
        delete[] d_array;
        // Return Initialized centroids
        return centroids;
      }

      /** \brief Initilization of cluster centroids.
        Details:
          (1) Create a vector containing numbers in range 1-no_vectors.
          (2) Suffle vector and keep the k - first cells.
          (3) Match each number with vectors' offsets.
          (4) Return selected vectors.
          @par[in] const std::vector<T>& dataset_vectors : vectors given from dataset
          @par[in] const int& no_vectors : number of vectors
          @par[in] const int& vectors_dim : vectors' dimensions (all vectors are
                                            dimensionally equal)
          @par[in] const int& no_clusters : number of clusters
      */
      template <typename T>
      std::vector<T> RandomInit(const std::vector<T>& dataset_vectors,
        const int& no_vectors, const int& vectors_dim, const int& no_clusters) {

          // Create a vector with numbers in range(0,no_clusters)
          std::vector<size_t> rand_vec(no_vectors);
          std::vector<T> centers(no_clusters * vectors_dim);
          for (size_t i = 0; i < no_vectors; i++) {
            rand_vec[i] = i;
          }
          // Suffle vector
          std::srand(std::time(0));
          std::random_shuffle(rand_vec.begin(), rand_vec.end());
          // Match offsets with vectors
          for (size_t i = 0; i < no_clusters; ++i) {
            for (size_t j = 0; j < vectors_dim; ++j) {
              centers[i * vectors_dim + j] = dataset_vectors[rand_vec[i] * vectors_dim + j];
            }
          }
          // Return Initialized centers
          return centers;
      }

    }
    /* Initialization Algorithms for curves */
    namespace curves {
      /** \brief Initilization of cluster centroids using Park-Jun Algorithm
      @par[in] const std::vector<std::pair<T,T>>& dataset_curves :
                                      curves givern from dataset
      @par[in] const int& dataset_curves_lengths : curves' lengths
      @par[in] const int& dataset_curves_offsets : curves' offsets
      @par[in] const int& no_curves : number of curves
      @par[in] const int& no_clusters : number of clusters
        return A tuple of:
          1) std::vector<std::pair<T,T> which represents the centroids curves
          2) std::vector<int> which stores centroids curves lengths
          3) std::vector<int> which stores centroids curves offsets
      */
      template <typename T>
      std::tuple<std::vector<std::pair<T,T>>,std::vector<int>,std::vector<int>>
        ParkJunInit(const std::vector<std::pair<T,T>>& dataset_curves,
        const std::vector<int>& dataset_curves_lengths,
        const std::vector<int>& dataset_curves_offsets,
        const int& no_curves, const int& no_clusters) {

        // Declare 2D d array which holds all distances between the curves
        //T d_array[no_curves][no_curves];
        T** d_array = new T*[no_curves];
        for(int i = 0; i < no_curves; ++i)
          d_array[i] = new T[no_curves];
        // Calculate all distances using ManhattanDistance
        for (size_t i = 0; i < no_curves; ++i) {
          for (size_t j = 0; j < no_curves; ++j) {
            if (i == j) {
              d_array[i][j] = {};
            } else {
              d_array[i][j] = metric::DTWDistance<T> (
                std::next(dataset_curves.begin(),dataset_curves_offsets[i]),
                std::next(dataset_curves.begin(),
                          dataset_curves_offsets[i] + dataset_curves_lengths[i]),
                std::next(dataset_curves.begin(),dataset_curves_offsets[j]),
                std::next(dataset_curves.begin(),
                          dataset_curves_offsets[j] + dataset_curves_lengths[j]));
            }
          }
        }
        // Declare 1D v array and computes its values as the algorithm proposes
        std::vector<T> v(no_curves);
        for (size_t i = 0; i < no_curves; ++i) {
          T o_sum{};
          for (size_t j = 0; j < no_curves; ++j) {
            T i_sum{};
            for (size_t t = 0; t < no_curves; ++t) {
              i_sum += d_array[j][t];
            }
            o_sum += d_array[i][j] / i_sum;
          }
          v[i] = o_sum;
        }

        // Find indices of no_clusters smallest v_array values
        std::vector<int> idx = utils::ArgMin<T>(v, no_clusters);
        // Initialize curves centroids
        std::vector<std::pair<T,T>> centroids;
        std::vector<int> centroids_lengths(no_clusters);
        std::vector<int> centroids_offsets(no_clusters);
        int offset{};
        for (size_t i = 0; i < idx.size(); ++i) {
          centroids_lengths[i] = dataset_curves_lengths[idx[i]];
          for (size_t j = 0; j < centroids_lengths[i]; ++j) {
            centroids.push_back(dataset_curves[dataset_curves_offsets[idx[i]] + j]);
          }
          centroids_offsets[i] = offset;
          offset += centroids_lengths[i];
        }
        // Return Initialized centroids
        for(int i = 0; i < no_curves; ++i)
          delete[] d_array[i];
        delete[] d_array;
        return std::make_tuple(centroids,centroids_lengths,centroids_offsets);
      }

      /** \brief Initilization of cluster centroids.
        Details:
          (1) Create a vector containing numbers in range 1-no_curves.
          (2) Suffle vector and keep the k - first cells.
          (3) Match each number with curves' offsets.
          (4) Calculate curves' new offsets.
          (4) Return centers as a tuple containg also their legths and offsets.
          @par[in] const std::vector<std::pair<T,T>>& dataset_curves :
                                          curves givern from dataset
          @par[in] const int& dataset_curves_lengths : curves' lengths
          @par[in] const int& dataset_curves_offsets : curves' offsets
          @par[in] const int& no_curves : number of curves
          @par[in] const int& no_clusters : number of clusters
      */
      template <typename T>
      std::tuple<std::vector<std::pair<T,T>>,std::vector<int>,std::vector<int>>
      RandomInit(const std::vector<std::pair<T,T>>& dataset_curves,
        const std::vector<int>& dataset_curves_lengths,
        const std::vector<int>& dataset_curves_offsets,
        const int& no_curves, const int& no_clusters) {

        // Create a vector with numbers in range (0,no_clusters)
        std::vector<size_t> rand_vec(no_curves);
        for (size_t i = 0; i < no_curves; i++) {
          rand_vec[i] = i;
        }
        // Suffle vector
        std::srand(std::time(0));
        std::random_shuffle(rand_vec.begin(), rand_vec.end());
        // Match ids with vectors
        int offset{};
        std::vector<std::pair<T,T>> centers;
        std::vector<int> centers_offsets(no_clusters);
        std::vector<int> centers_lengths(no_clusters);
        for (size_t i = 0; i < no_clusters; ++i) {
          centers_lengths[i] = dataset_curves_lengths[rand_vec[i]];
          for (size_t j = 0; j < centers_lengths[i]; ++j) {
            centers.push_back(dataset_curves[dataset_curves_offsets[rand_vec[i]] + j]);
          }
          centers_offsets[i] = offset;
          offset += centers_lengths[i];
        }
        // Return Initialized centroids
        return std::make_tuple(centers,centers_lengths,centers_offsets);
      }
    }
  }
}


#endif
