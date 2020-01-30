#ifndef METRIC
#define METRIC

#include <cmath>
#include <iterator>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "../utils/utils.h"

namespace metric {
  /** \brief Computes Manhattan Distance of 2 points in R^2
    @par const std::pair<T,T>& p - first point
    @par const std::pair<T,T>& q - second point
  */
  template <typename T>
  T _2DManhattanDistance(const std::pair<T,T>& p, const std::pair<T,T>& q) {
    T x_diff = std::abs((std::get<0>(p)-std::get<0>(q)));
    T y_diff = std::abs((std::get<1>(p)-std::get<1>(q)));
    return x_diff + y_diff;
  }
  /** \brief Computes Manhattan Distance of 2 points in R^d
    @par iterator p - iterator of the dataset point
    @par iterator q - iterator of the query point
    @par iterator q_end - end iterator of query point
      (each point is of dimension D so end of query point signals the
       final coordinate of dataset point as well)
    returns: manhattan_distance type of T
  */
  template <typename T, typename iterator>
  T ManhattanDistance(iterator p, iterator q, iterator q_end) {
    T manhattan_distance{};
    for (; q < q_end; ++q, ++p) {
      manhattan_distance += std::abs(*q - *p);
    }
    return manhattan_distance;
  }
  /** \brief Computes Euclidian Distance of 2 points in R^2
    @par const std::pair<T,T>& p - first point
    @par const std::pair<T,T>& q - second point
  */
  template <typename T>
  T _2DEuclidianDistance(const std::pair<T,T>& p, const std::pair<T,T>& q) {
    T x_diff = (std::get<0>(p)-std::get<0>(q)) * (std::get<0>(p)-std::get<0>(q));
    T y_diff = (std::get<1>(p)-std::get<1>(q)) * (std::get<1>(p)-std::get<1>(q));
    return sqrt(x_diff + y_diff);
  }
  /** \brief Computes Euclidian Distance of 2 points in R^d
    @par iterator p - iterator of the dataset point
    @par iterator q - iterator of the query point
    @par iterator q_end - end iterator of query point
      (each point is of dimension D so end of query point signals the
       final coordinate of dataset point as well)
    returns: euclidian_distance type of T
  */
  template <typename T, typename iterator>
  T SquaredEuclidianDistance(iterator p, iterator q, iterator q_end) {
    T euclidean_distance{};
    for (; q < q_end; ++q, ++p) {
      euclidean_distance += (*q - *p) * (*q - *p);
    }
    return euclidean_distance;
  }

  /** \brief Computes Dynamic Time Warping between two curves
    @par iterator p - iterator of the dataset curve
    @par iterator p_end - end iterator of dataset curve
    @par iterator q - iterator of the query curve
    @par iterator q_end - end iterator of query curve
      (each curve represents as a vector of pairs in which a pair
      represents a 2D-point)
  */
  template <typename T, typename iterator>
  T DTWDistance(iterator p, iterator p_end, iterator q, iterator q_end) {
    T dtw_distance{};
    // save up start iterator of query curve
    iterator q_start = q;
    // Get correspodent lengths of the two curves
    size_t N = std::distance(p,p_end);
    size_t M = std::distance(q,q_end);
    size_t dtw_size = (N + 1) * (M + 1);
    // Initialize the dtw 2D array using 1D notation
    auto dtw_array = new T [dtw_size];
    dtw_array[0] = 0;
    for (size_t i = 1; i < N + 1; ++i) {
      dtw_array[i * (M + 1)] = std::numeric_limits<T>::max();
    }
    for (size_t i = 1; i < M + 1; ++i) {
      dtw_array[i] = std::numeric_limits<T>::max();
    }
    // Compute the dtw distance using dynamic programming
    for (size_t i = 1; i < N + 1, p < p_end; ++i, ++p) {
      for (size_t j = 1; j < M + 1, q < q_end; ++j, ++q) {
        T dist = _2DEuclidianDistance(*p, *q);
        dtw_array[i * (M + 1) + j] = dist +
          utils::min(dtw_array[(i - 1) * (M + 1) + j],      // increment
            dtw_array[i * (M + 1) + j - 1],                 // deletion
            dtw_array[(i - 1) * (M + 1) + j - 1]);          // match
      }
      // reset query curve iterator
      q = q_start;
    }

    // Set result to a new variable and free memory allocated
    dtw_distance = dtw_array[dtw_size - 1];
    delete[] dtw_array;
    // Return dynamic time warping distance
    return dtw_distance;
  }

  template <typename T, typename iterator>
  std::vector<std::pair<T,T>> BestTraversal(iterator p, iterator p_end,
                                            iterator q, iterator q_end) {
    T dtw_distance{};
    // save up start iterator of query curve
    iterator q_start = q;
    // Get correspodent lengths of the two curves
    size_t N = std::distance(p,p_end);
    size_t M = std::distance(q,q_end);
    size_t dtw_size = (N + 1) * (M + 1);
    // Initialize the dtw 2D array using 1D notation
    auto dtw_array = new T [dtw_size];
    dtw_array[0] = 0;
    for (size_t i = 1; i < N + 1; ++i) {
      dtw_array[i * (M + 1)] = std::numeric_limits<T>::max();
    }
    for (size_t i = 1; i < M + 1; ++i) {
      dtw_array[i] = std::numeric_limits<T>::max();
    }
    // Compute the dtw distance using dynamic programming
    for (size_t i = 1; i < N + 1, p < p_end; ++i, ++p) {
      for (size_t j = 1; j < M + 1, q < q_end; ++j, ++q) {
        T dist = _2DEuclidianDistance(*p, *q);
        dtw_array[i * (M + 1) + j] = dist +
          utils::min(dtw_array[(i - 1) * (M + 1) + j],      // increment
            dtw_array[i * (M + 1) + j - 1],                 // deletion
            dtw_array[(i - 1) * (M + 1) + j - 1]);          // match
      }
      // reset query curve iterator
      q = q_start;
    }
    // iterate dtw array to get optimal path
    size_t i = N;
    size_t j = M;
    std::vector<std::pair<T,T>> ret_vec;
    while(i>1 && j>1) {
      if (i == 1) {
        j--;
      } else if (j == 1) {
        i--;
      } else {
        T dist = utils::min(dtw_array[(i - 1) * (M + 1) + j],
                dtw_array[i * (M + 1) + j - 1],
                dtw_array[(i - 1) * (M + 1) + j - 1]);
        if (dtw_array[(i - 1) * (M + 1) + j] == dist) {
          i--;
        } else if (dtw_array[i * (M + 1) + j - 1] == dist) {
          j--;
        } else {
          i--;
          j--;
        }
      }
      ret_vec.push_back(std::make_pair(i-1,j-1));
    }
    ret_vec.push_back(std::make_pair(0,0));
    delete[] dtw_array;
    std::reverse(ret_vec.begin(), ret_vec.end());
    // for(auto& i:ret_vec) {
    //   std::cout << i.first << " ";
    // }
    // std::cout << std::endl;
    return ret_vec;
  }

  /** \brief Computes average and max distance ratio appox_dists / exact_dists
    Each tuple consists of the nearest distance found, the id of the point/curve
    with minimum distance and the time taken to be computed
    @par const std::vector<std::tuple<T,U,double>> &exact - Pass by reference
      the vector of results produced by exact Nearest Neighbor for each query.
    @par const std::vector<std::tuple<T,U,double>> &approx - Pass by reference
      the vector of results produced by approximate Nearest Neighbor for
      each query.
  */
  template <typename T, typename U>
  std::tuple<double,double,int>
    EvaluationMetric(const std::vector<std::tuple<T,U,double>> &exact,
      const std::vector<std::tuple<T,U,double>> &approx) {

    // number of queries failed to be ansewered by A-NN
    int cnt_not_found{};
    int found{};
    // distance difference between exact and approximate NN
    double distance_error;
    // the sum of distance_error for each query
    double af{};
    // the averaged af by number of queries ansewered successfully
    double avg_af{};
    // number of queries
    size_t Q = exact.size();

    double max_af = std::numeric_limits<double>::min();
    for (size_t i = 0; i < Q; ++i) {
      // check the case for zero distance found (math warning by divsion with 0)
      if (std::get<0>(approx[i]) != 0 && std::get<0>(exact[i]) != 0) {
        // max distance found, so query failed to be ansewered
        if (std::get<0>(approx[i]) == std::numeric_limits<T>::max()) {
          cnt_not_found++;
          continue;
        }
        distance_error = (double) std::get<0>(approx[i]) / std::get<0>(exact[i]);
        if (distance_error > max_af) {
          max_af = distance_error;
        }
        af += distance_error;
      }
      if (std::get<0>(approx[i]) == std::get<0>(exact[i])) {
        found++;
      }
    }
    // For total number of queries discard the ones who failed to be ansewered
    avg_af = af / (Q - cnt_not_found) ;
    // return result as tuple of max Af, average Af and total #queries failed
    return std::make_tuple(max_af,avg_af,cnt_not_found);
  }
  /* Implementation of Silhouette metric for vectors */
  namespace vectors {
    /** \brief Computes Silhouette metric
      @par[in] dataset_vectors - dataset vectors coordinates
      @par[in] no_vectors - total number of vectors in the dataset
      @par[in] vectors_dim - vectors dimension
      @par[in] clusters - vector of vectors where each dataset_vector is assigned
      @par[in] centroids - centroids coordinates stored in 1D vector
      @par[in] mapped_vectors - each vector is mapped to a cluster index
      return : a pair of average s(p) of points in cluster i and
        stotal = average s(p) of points in dataset
    */
    template <typename T>
    std::pair<std::vector<double>,double> Silhouette(
      const std::vector<T>& dataset_vectors,
      const int& no_vectors, const int& vectors_dim,
      const std::vector<std::vector<size_t>>& clusters,
      const std::vector<T>& centroids,
      const std::map<int,int>& mapped_vectors) {
      /*  s(i) = (b(i) − a(i)) / max{a(i), b(i)}  */
      std::vector<double> s(no_vectors);
      std::vector<double> a(no_vectors);
      std::vector<double> b(no_vectors);
      /*
        Precompute closest centroid to each other centroid using
        manhattan distance
      */
      std::map<int,int> closest_centroid;
      for (size_t i = 0; i < clusters.size(); ++i) {
        T min_dist = std::numeric_limits<T>::max();
        int min_id;
        // for every centroid compute manhattan distance to each other centroid
        for (size_t j = 0; j < clusters.size(); ++j) {
          if (i == j) continue;
          T dist = ManhattanDistance<T>(
            std::next(centroids.begin(), i * vectors_dim),
            std::next(centroids.begin(), j * vectors_dim),
            std::next(centroids.begin(), j * vectors_dim + vectors_dim));
          // Pick this one with the minimum manhattan distance
          if (dist < min_dist) {
            min_dist = dist;
            min_id = j;
          }
        }
        closest_centroid[i] = min_id;
      }
      /* Iterate over every vector to compute its a_i value */
      for (auto it = mapped_vectors.cbegin(); it != mapped_vectors.cend(); ++it) {
        /* Get all vectors' indexes in the same cluster */
        std::vector<size_t> cluster_vectors = clusters[it->second];
        T a_total_dist{};
        for (const auto& i: cluster_vectors) {
          // skip itself
          if (i == it->first) continue;
          a_total_dist += ManhattanDistance<T>(
            std::next(dataset_vectors.begin(), it->first * vectors_dim),
            std::next(dataset_vectors.begin(), i * vectors_dim),
            std::next(dataset_vectors.begin(), i * vectors_dim + vectors_dim));
        }
        if (cluster_vectors.size() != 1) {
          a[it->first] = (double) a_total_dist / (cluster_vectors.size() - 1);
        } else {
          a[it->first] = (double) a_total_dist;
        }
      }
      /* Iterate over every vector to compute its b_i value */
      for (auto it = mapped_vectors.cbegin(); it != mapped_vectors.cend(); ++it) {
        /* Get all vectors' indexes in the closest centroid cluster */
        std::vector<size_t> cluster_vectors = clusters[closest_centroid[it->second]];
        T b_total_dist{};
        for (const auto& i: cluster_vectors) {
          b_total_dist += ManhattanDistance<T>(
            std::next(dataset_vectors.begin(), it->first * vectors_dim),
            std::next(dataset_vectors.begin(), i * vectors_dim),
            std::next(dataset_vectors.begin(), i * vectors_dim + vectors_dim));
        }
        if (cluster_vectors.size() != 0) {
          b[it->first] = (double) b_total_dist / cluster_vectors.size();
        } else {
          b[it->first] = (double) b_total_dist;
        }
      }
      /* Iterate over every vector to compute its Silhouette value */
      for (size_t i = 0; i < s.size(); ++i) {
        s[i] = (b[i] - a[i]) / utils::max(a[i],b[i]);
      }
      /* Compute average s(p) of points in cluster i */
      std::vector<double> s_avg(clusters.size(), 0);
      for (size_t i = 0; i < clusters.size(); ++i) {
        int counter = 0;
        for (auto it = mapped_vectors.cbegin(); it != mapped_vectors.cend(); ++it) {
          if (it->second == i) {
            counter++;
            s_avg[i] += s[it->first];
          }
        }
        if (counter != 0) {
          s_avg[i] /= counter;
        }
      }
      /* Compute stotal = average s(p) of points in dataset */
      double s_total{};
      for (size_t i = 0; i < clusters.size(); ++i) {
        s_total += s_avg[i];
      }
      s_total /= clusters.size();
      // Return result as pair
      return std::make_pair(s_avg,s_total);
    }
  }
  /* Implementation of Silhouette metric for curves */
  namespace curves {
    /** \brief Computes Silhouette metric
      @par[in] dataset_curves - dataset curves coordinates
      @par[in] dataset_curves_lengths
      @par[in] dataset_curves_offsets
      @par[in] no_curves - total number of curves in the dataset
      @par[in] clusters - vector of vectors where each dataset_curve is assigned
      @par[in] centroids_curves - centroids coordinates
      @par[in] centroid_curves_lengths
      @par[in] centroid_curves_offsets
      @par[in] mapped_curves - each vector is mapped to a cluster index
      return : a pair of average s(p) of curves in cluster i and
        stotal = average s(p) of curves in dataset
    */
    template <typename T>
    std::pair<std::vector<double>,double> Silhouette(
      const std::vector<std::pair<T,T>>& dataset_curves,
      const std::vector<int>& dataset_curves_lengths,
      const std::vector<int>& dataset_curves_offsets,
      const int& no_curves, const std::vector<std::vector<size_t>>& clusters,
      const std::vector<std::pair<T,T>>& centroid_curves,
      const std::vector<int>& centroid_curves_lengths,
      const std::vector<int>& centroid_curves_offsets,
      const std::map<int,int>& mapped_curves) {
      /*  s(i) = (b(i) − a(i)) / max{a(i), b(i)}  */
      std::vector<double> s(no_curves);
      std::vector<double> a(no_curves);
      std::vector<double> b(no_curves);
      /*
        Precompute closest centroid to each other centroid using
        manhattan distance
      */
      std::map<int,int> closest_centroid;
      for (size_t i = 0; i < clusters.size(); ++i) {
        T min_dist = std::numeric_limits<T>::max();
        int min_id;
        // for every centroid compute manhattan distance to each other centroid
        for (size_t j = 0; j < clusters.size(); ++j) {
          if (i == j) continue;
          T dist = DTWDistance<T>(
          std::next(centroid_curves.begin(),centroid_curves_offsets[i]),
          std::next(centroid_curves.begin(),
                    centroid_curves_offsets[i] + centroid_curves_lengths[i]),
          std::next(centroid_curves.begin(),centroid_curves_offsets[j]),
          std::next(centroid_curves.begin(),
                    centroid_curves_offsets[j] + centroid_curves_lengths[j]));
          // Pick this one with the minimum manhattan distance
          if (dist < min_dist) {
            min_dist = dist;
            min_id = j;
          }
        }
        closest_centroid[i] = min_id;
      }

      /* Iterate over every vector to compute its a_i value */
      for (auto it = mapped_curves.cbegin(); it != mapped_curves.cend(); ++it) {
        /* Get all curves' indexes in the same cluster */
        std::vector<size_t> cluster_curves = clusters[it->second];
        T a_total_dist{};
        for (const auto& i: cluster_curves) {
          // skip itself
          if (i == it->first) continue;
          a_total_dist += DTWDistance<T>(
          std::next(dataset_curves.begin(),dataset_curves_offsets[it->first]),
          std::next(dataset_curves.begin(),
                    dataset_curves_offsets[it->first] + dataset_curves_lengths[it->first]),
          std::next(dataset_curves.begin(),dataset_curves_offsets[i]),
          std::next(dataset_curves.begin(),
                    dataset_curves_offsets[i] + dataset_curves_lengths[i]));
        }
        if (cluster_curves.size() != 1 && cluster_curves.size() != 0) {
          a[it->first] = (double) a_total_dist / (cluster_curves.size() - 1);
        } {
          a[it->first] = (double) a_total_dist;
        }
      }
      /* Iterate over every vector to compute its b_i value */
      for (auto it = mapped_curves.cbegin(); it != mapped_curves.cend(); ++it) {
        /* Get all vectors' indexes in the closest centroid cluster */
        std::vector<size_t> cluster_curves = clusters[closest_centroid[it->second]];
        T b_total_dist{};
        for (const auto& i: cluster_curves) {
          b_total_dist += DTWDistance<T>(
          std::next(dataset_curves.begin(),dataset_curves_offsets[it->first]),
          std::next(dataset_curves.begin(),
                    dataset_curves_offsets[it->first] + dataset_curves_lengths[it->first]),
          std::next(dataset_curves.begin(),dataset_curves_offsets[i]),
          std::next(dataset_curves.begin(),
                    dataset_curves_offsets[i] + dataset_curves_lengths[i]));
        }
        if (cluster_curves.size() != 0) {
          b[it->first] = (double) b_total_dist / cluster_curves.size();
        } else {
          b[it->first] = (double) b_total_dist;
        }
      }
      /* Iterate over every vector to compute its Silhouette value */
      for (size_t i = 0; i < s.size(); ++i) {
        s[i] = (a[i] - b[i]) / utils::max(a[i],b[i]);
      }
      /* Compute average s(p) of points in cluster i */
      std::vector<double> s_avg(clusters.size(), 0);
      for (size_t i = 0; i < clusters.size(); ++i) {
        int counter = 0;
        for (auto it = mapped_curves.cbegin(); it != mapped_curves.cend(); ++it) {
          if (it->second == i) {
            counter++;
            s_avg[i] += s[it->first];
          }
        }
        if (counter != 0) {
          s_avg[i] /= counter;
        }
      }
      /* Compute stotal = average s(p) of points in dataset */
      double s_total{};
      for (size_t i = 0; i < clusters.size(); ++i) {
        s_total += s_avg[i];
      }
      s_total /= clusters.size();
      // Return result as pair
      return std::make_pair(s_avg,s_total);
    }
  }
}


#endif
