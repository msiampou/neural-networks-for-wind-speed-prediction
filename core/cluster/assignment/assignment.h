#ifndef ASSIGNMENT
#define ASSIGNMENT

#include <algorithm>
#include <map>
#include <tuple>

#include "../../../core/search/lsh.h"
#include "../../../core/metric/metric.h"
#include "../../../core/utils/utils.h"

namespace cluster {
  namespace assignment {
    /* Assignment Algorithms for vectors */
    namespace vectors {
      /** \brief Lloyd's Approach for Assignment.
       * Details:
        (1) For each vector calculate distances from each centroid.
        (2) Assign each vector to its closest cluster.
        (3) Store assignment costs.
        (4) Return clusters' offsets and total assignment costs.
        @par[in] const std::vector<T>& dataset_vectors : vectors given from dataset
        @par[in] const std::vector<T>& centroids : selected centroids of each cluster
        @par[in] const int& no_vectors : number of vectors
        @par[in] const int& vectors_dim : vectors' dimensions (all vectors are
          dimensionally equal)
        @par[in] const int& no_clusters : number of clusters
       * We use Euclidean Distance to calculate minimum distances.
      **/
      template <typename T>
      std::tuple<std::vector<std::vector<size_t>>, std::vector<T>> LloydsAssignment (
        const std::vector<T>& dataset_vectors, const std::vector<T>& centroids,
        const int& no_vectors, const int& vectors_dim, const int& no_clusters) {

        // Declare 2D vector which holds all vectors
        // to their assigned cluster.
        std::vector<std::vector<size_t>> d_array(no_clusters);
        std::vector<T> assign_costs(no_clusters, 0);
        // Calculate all distances using Euclidean Distance
        T total_cost{};
        for (size_t i = 0; i < no_vectors; i++) {
          T min_dist = std::numeric_limits<T>::max();
          size_t assigned_cluster{};
          for (size_t j = 0; j < no_clusters; j++) {
            T dist = metric::SquaredEuclidianDistance<T>(
              std::next(dataset_vectors.begin(), i * vectors_dim),
              std::next(centroids.begin(), j * vectors_dim),
              std::next(centroids.begin(), j * vectors_dim + vectors_dim));
            // Store minimum dist from clusters
            if (dist < min_dist) {
              min_dist = dist;
              assigned_cluster = j;
            }
          }
          // Store to closest cluster and
          // compute assignment cost :
          // Assignment Cost (aka Obj. Function)
          // equals to the summary of distances
          // of each vector to its cluster's centroid
          d_array[assigned_cluster].push_back(i);
          assign_costs[assigned_cluster] += min_dist;
        }
        return std::make_tuple(d_array,assign_costs);
      }
    }
    /* Assignment Algorithms for curves */
    namespace curves {
      /** \brief Lloyd's Approach for Assignment.
       * Details:
        (1) For each vector calculate distances from each centroid.
        (2) Assign each vector to its closest cluster.
        (3) Store assignment costs.
        (4) Return clusters' offsets and total assignment costs.
        @par[in] const std::vector<std::pair<T,T>>& dataset_curves :
                                          curves givern from dataset
        @par[in] std::tuple<std::vector<std::pair<T,T>>,
                 std::vector<int>, std::vector<int>>& centroids :
                                        centroids assigned so far.
        @par[in] const int& dataset_curves_lengths : curves' lengths
        @par[in] const int& dataset_curves_offsets : curves' offsets
        @par[in] const int& no_curves : number of curves
        @par[in] const int& no_clusters : number of clusters
       * We use Euclidean Distance to calculate minimum distances.
      **/
      template <typename T>
      std::tuple<std::vector<std::vector<size_t>>, std::vector<T>>
        LloydsAssignment (const std::vector<std::pair<T,T>>& dataset_curves,
        std::tuple<std::vector<std::pair<T,T>>, std::vector<int>,
        std::vector<int>>& centroids,
        const std::vector<int>& dataset_curves_lengths,
        const std::vector<int>& dataset_curves_offsets,
        const int& no_curves, const int& no_clusters) {

        // Declare 2D vector which holds all curves
        // to their assigned cluster
        std::vector<std::vector<size_t>> d_array(no_clusters);
        std::vector<T> assign_costs(no_clusters, 0);
        // Get lengths, offsets and points of curves
        // Warning: DTW does not work without const&
        const std::vector<std::pair<T,T>>& centroid_curves = std::get<0>(centroids);
        const std::vector<int>& centroid_lengths = std::get<1>(centroids);
        const std::vector<int>& centroid_offsets = std::get<2>(centroids);
        // Calculate all distances using DTW Distance
        T total_cost{};
        for (size_t i = 0; i < no_curves; i++) {
          T min_dist = std::numeric_limits<T>::max();
          size_t assigned_cluster{};
          for (size_t j = 0; j < no_clusters; j++) {
            T dist = metric::DTWDistance<T> (
              std::next(dataset_curves.begin(),dataset_curves_offsets[i]),
              std::next(dataset_curves.begin(),
                        dataset_curves_offsets[i] + dataset_curves_lengths[i]),
              std::next(centroid_curves.begin(),centroid_offsets[j]),
              std::next(centroid_curves.begin(),
                        centroid_offsets[j] + centroid_lengths[j]));
            // Store minimum dist from clusters
            if (dist < min_dist) {
              min_dist = dist;
              assigned_cluster = j;
            }
          }
          // Store id to closest cluster
          d_array[assigned_cluster].push_back(i);
          assign_costs[assigned_cluster] += min_dist;
        }

        return std::make_tuple(d_array,assign_costs);
      }
    }
    /** \brief ReverseAssignment assignes vectors to clusters using range LSH
      @par[in] dataset_vectors : vectors given from dataset
      @par[in] dataset_vectors_ids : vectors' ids given from the dataset
      @par[in] centroids : selected centroids of each cluster
      @par[in] no_vectors : total number of vectors
      @par[in] vectors_dim : vectors' dimensions (all vectors are
        dimensionally equal)
      @par[in] lsh_structure : lsh_structure to execute range search queries
      @par[in] map_id_to_index : a map which associates each vector id
        to each index in the dataset
    */
    template <typename T, typename U>
    std::tuple<std::vector<std::vector<size_t>>, std::vector<T>>
    ReverseAssignment(const std::vector<T>& dataset_vectors,
      const std::vector<U>& dataset_vectors_ids, const std::vector<T>& centroids,
      const int& no_vectors, const int& vectors_dim, const int& no_clusters,
      search::vectors::LSH<T,U> *lsh, std::map<U,int>& map_id_to_index) {

      /* At first create a map vector_id --> centroid_id. Unassigned vectors
        are mapped to -1. All vectors are initialized to -1 */
      std::map<int,int> assigned_vector_to_centroid;
      for (size_t i = 0; i < no_vectors; ++i) {
        assigned_vector_to_centroid[i] = -1;
      }
      /* Assign 95% of the dataset using range search. For the rest compate
        each vector distance to each centroid and pick this one which
        minimized the corresponding distance */
      int no_vectors_assigned = 0;
      /** Alogirthm key steps:
        1) At each iteration, for each centroid c,
          range/ball queries centered at c.
        2) Mark assigned points: either move them at end of LSH buckets
          (and insert "barrier", or mark them using "flag" field).
        3) Multiply radius by 2, start with min(dist between centers)/2,
        until most balls get no new point (centers mapped to buckets once)
      */
      /**
        Initialize radius compoting the distance of each centroid ot each
        other centroid
      */
      T min_dist = std::numeric_limits<T>::max();
      T max_dist = std::numeric_limits<T>::min();
      for (size_t i = 0; i < no_clusters; ++i) {
        for (size_t j = i + 1; j < no_clusters; ++j) {
          T dist = metric::SquaredEuclidianDistance<T>(
            std::next(centroids.begin(), i * vectors_dim),
            std::next(centroids.begin(), j * vectors_dim),
            std::next(centroids.begin(), j * vectors_dim + vectors_dim));
          if (dist < min_dist && !std::isinf(dist)) {
            min_dist = dist;
          }
          if (dist > max_dist && !std::isinf(dist)) {
            max_dist = dist;
          }
        }
      }
      double radius = (double) (min_dist / 2) + (max_dist / 10);
      // Declare 2D vector which holds all vectors
      // to their assigned cluster.
      std::vector<std::vector<size_t>> d_array(no_clusters);
      std::vector<T> assign_costs(no_clusters, 0);
      /* Execute step 1 */
      while (radius <= max_dist) {
        /* Execute range search for each centroid */
        for (size_t i = 0; i < no_clusters; ++i) {
          auto range_results = lsh->RadiusNearestNeighbor(centroids, i, radius);
          for (const auto& res : range_results) {
            T v_dist = std::get<0>(res);
            U v_id = std::get<1>(res);
            int v_index = map_id_to_index[v_id];
            // Unassigned vector
            if (assigned_vector_to_centroid[v_index] == -1) {
              // Store to closest cluster and
              // compute assignment cost :
              // Assignment Cost (aka Obj. Function)
              // equals to the summary of distances
              // of each vector to its cluster's centroid
              d_array[i].push_back(v_index);
              assign_costs[i] += v_dist;
              // Update structures
              assigned_vector_to_centroid[v_index] = i;
              no_vectors_assigned++;
            } else {
              // Vector has assigned to some other centroid before
              // In this case compare its old distance with the new one
              int centroid_assigned = assigned_vector_to_centroid[v_index];
              T prev_dist = metric::SquaredEuclidianDistance<T>(
                std::next(dataset_vectors.begin(), v_index * vectors_dim),
                std::next(centroids.begin(), centroid_assigned * vectors_dim),
                std::next(centroids.begin(), centroid_assigned * vectors_dim + vectors_dim));
              if (v_dist < prev_dist) {
                d_array[i].push_back(v_index);
                // Substuct previous dist and add the one
                assign_costs[i] += v_dist;
                // Update structures
                assigned_vector_to_centroid[v_index] = i;
                no_vectors_assigned++;
              }
            }
          }
        }
        /* Increase radius (Multiply by 2) */
        radius *= 2;
      }
      /* For each unassigned vector compute its distance to each centroid */
      for (size_t i = 0; i < no_vectors; ++i) {
        if (assigned_vector_to_centroid[i] == -1) {
          T min_dist = std::numeric_limits<T>::max();
          int min_centroid_index{};
          for (size_t j = 0; j < no_clusters; ++j) {
            T dist = metric::SquaredEuclidianDistance<T>(
              std::next(dataset_vectors.begin(), i * vectors_dim),
              std::next(centroids.begin(), j * vectors_dim),
              std::next(centroids.begin(), j * vectors_dim + vectors_dim));
            if (dist < min_dist) {
              min_dist = dist;
              min_centroid_index = j;
            }
          }
          if (!std::isinf(min_dist)) {
            d_array[min_centroid_index].push_back(i);
            assign_costs[min_centroid_index] += min_dist;
            // Update structures
            assigned_vector_to_centroid[i] = min_centroid_index;
            no_vectors_assigned++;
          } else {
            d_array[min_centroid_index].push_back(i);
            // Update structures
            assigned_vector_to_centroid[i] = min_centroid_index;
            no_vectors_assigned++;
          }
        }
      }
      // return a tuple of the result
      return std::make_tuple(d_array,assign_costs);
    }
  }
}

#endif
