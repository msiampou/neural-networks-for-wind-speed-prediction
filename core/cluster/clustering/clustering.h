#ifndef CLUSTERING
#define CLUSTERING

#include <chrono>
#include <map>
#include <tuple>

#include "../../search/lsh.h"
#include "../../utils/utils.h"
#include "../../vectorization/vectorization.h"
#include "../initialization/initialization.h"
#include "../assignment/assignment.h"
#include "../update/update.h"

using namespace std::chrono;

namespace cluster {
  namespace vectors {
    using namespace cluster::initialization::vectors;
    using namespace cluster::assignment::vectors;
    using namespace cluster::update::vectors;
    /**
      \brief Class Cluster representing k-means/k-medoids algorithm
      for vectors
    */
    template <typename T, typename U>
    class Cluster {
      private:
        /* class Cluster parameters */
        int no_clusters;
        int max_iter;
        std::string init;
        std::string assign;
        std::string update;
        /* Range-lsh parameters */
        search::vectors::LSH<T,U> *lsh_structure;
        double window;
        uint8_t k;  // number of hash functions
        uint8_t L;  // number of hash tables
        /* Dataset info */
        std::vector<T> dataset_vectors;
        std::vector<U> dataset_vectors_ids;
        uint16_t vectors_dim;
        uint32_t no_vectors;
        std::map<U,int> map_id_to_index;
      public:
        /** \brief class Cluster constructor
          @par no_clusters: int, optional, default: 8
            The number of clusters to form as well as the number of
            centroids to generate.
          @par max_iter : int, default: 300
            Maximum number of iterations of the k-means/k-medoids algorithm for
            a single run.
          @par init : {‘k-means++’, ‘random’}
            Method for initialization, defaults to ‘k-means++’:
          @par assign : {‘lloyds’, ‘lsh’}
            Method for assignment, defaults to ‘lloyds’:
          @par update : {‘mean’, ‘pam’}
            Method for update, defaults to 'mean'
        */
        Cluster(int no_clusters = 8, int max_iter = 1,
          std::string init = "random", std::string assign = "lloyd",
          std::string update = "mean", uint8_t no_hf = 3, uint8_t no_ht = 5)
          : no_clusters(no_clusters), max_iter(max_iter), init(init),
          assign(assign), update(update), k(no_hf), L(no_ht) {}
        /**
          \brief default Class Cluster destructor
        */
        ~Cluster() {
          if (assign == "range-lsh") delete lsh_structure;
        };
        /** \brief Fit method stores dataset info for clustering
          @par[in] dv : vectors given from dataset
          @par[in] no_v : number of vectors
          @par[in] v_dim : vectors' dimensions (all vectors are dimensionally equal)
        */
        void Fit(const std::vector<T>& dv, const std::vector<U>& dv_ids,
          const uint32_t& no_v, const uint16_t& v_dim) {
          dataset_vectors = dv;
          dataset_vectors_ids = dv_ids;
          no_vectors = no_v;
          vectors_dim = v_dim;
          /**
            If assignment step is executed using lsh range search, map dataset
            to corresponding lsh structures
          */
          if (assign == "range-lsh") {
            window = utils::ComputeMean<T>(dataset_vectors,
                                           vectors_dim, no_vectors);
            /* Index no_vectors points into L hashtables */
            lsh_structure = new search::vectors::LSH<T,U>(k, L, vectors_dim,
                                                          no_vectors, window,
                                                          dataset_vectors,
                                                          dataset_vectors_ids);
            /* Map each id from dataset_vectors_ids to its index */
            for (size_t i = 0; i < no_vectors; ++i) {
              map_id_to_index[dataset_vectors_ids[i]] = i;
            }
          }
        }
        /**
          \brief Predict the closest cluster each sample in dataset belongs to
          returns a vector of centroids coordinates and a vector of vectors
          which stores in each position the indexes of dataset_vectors assigned
          in this cluster
        */
        std::tuple<std::vector<T>,std::vector<std::vector<size_t>>,double>
          Predict(void) {
          /* Start time measuring */
          auto start = high_resolution_clock::now();
          /* Declare types */
          std::vector<T> centroids;
          std::tuple<std::vector<std::vector<size_t>>,std::vector<T>> clusters;
          /* At first initialize centroids */

					if (init == "random") {
						centroids = RandomInit(dataset_vectors, no_vectors,
																	 vectors_dim, no_clusters);
					} else if (init == "k-means++") {
						centroids = ParkJunInit(dataset_vectors, no_vectors,
																	  vectors_dim, no_clusters);
					}

          /* Calculate clusters and update centroids max_iter times */
          for (size_t i = 0; i < max_iter; ++i) {
            /* Assigment step */
            if (assign == "lloyd") {
							clusters = LloydsAssignment(dataset_vectors, centroids,
																					no_vectors, vectors_dim, no_clusters);
						} else if (assign == "range-lsh") {
              clusters = cluster::assignment::
                          ReverseAssignment<T,U>(dataset_vectors, dataset_vectors_ids,
                                                centroids, no_vectors, vectors_dim,
                                                no_clusters, lsh_structure,
                                                map_id_to_index);
						}
            /* Update step */
            auto prev_centroids = centroids;
						if (update == "mean") {
							centroids = LloydsUpdate(dataset_vectors, centroids,
																		 	 no_vectors, vectors_dim,
																			 no_clusters, std::get<0>(clusters));
						} else if (update == "pam") {
							centroids = PAMUpdate(dataset_vectors, centroids,
																	  no_vectors, vectors_dim, no_clusters,
																	  std::get<0>(clusters), std::get<1>(clusters));
						}
            if(prev_centroids == centroids) break;
          }
          /* End time measuring */
          auto stop = high_resolution_clock::now();
          duration <double> total_time = duration_cast<duration<double>>(stop - start);
          // Return result in the form of a tuple
          return std::make_tuple(centroids,std::get<0>(clusters),total_time.count());
        }
        /** \brief Map each vector index to the corresponding cluster
          @par[in] clusters - std::vector<std::vector<size_t>> returned by Predict
          return: a map of each index to the the corresponding cluster
        */
        std::map<int,int> MapToClusters(std::vector<std::vector<size_t>> clusters) {
          std::map<int,int> map_result;
          int cluster_idx = 0;
          for (auto const& cluster: clusters) {
            for (auto const& object: cluster) {
              map_result[object] = cluster_idx;
            }
            cluster_idx++;
          }
          return map_result;
        }
    };
  }
  namespace curves {
    using namespace cluster::initialization::curves;
    using namespace cluster::assignment::curves;
    using namespace cluster::update::curves;
    /**
      \brief Class Cluster representing k-means/k-medoids algorithm
      for curves
    */
    template <typename T, typename U>
    class Cluster {
      private:
        /* class Cluster parameters */
        int no_clusters;
        int max_iter;
        std::string init;
        std::string assign;
        std::string update;
        /* Range-LSH and Grid parameters for vectorization */
        std::vector<vectorization::Grid<T>> grids;
        std::vector<std::vector<double>> L_grid_dataset_vectors;
        std::vector<search::vectors::LSH<T,U>> lsh_structures;
        uint8_t k;
        uint8_t L;
        uint8_t L_grid;
        uint16_t vectors_dim;
        double *window;
        double delta;
        const uint8_t factor = 10;
        /* Dataset info */
        std::vector<std::pair<T,T>> dataset_curves;
        std::vector<U> dataset_curves_ids;
        std::vector<int> dataset_curves_lengths;
        std::vector<int> dataset_curves_offsets;
        int no_curves;
        std::map<U,int> map_id_to_index;
      public:
        /** \brief class Cluster constructor
          @par no_clusters: int, optional, default: 8
            The number of clusters to form as well as the number of
            centroids to generate.
          @par max_iter : int, default: 300
            Maximum number of iterations of the k-means/k-medoids algorithm for
            a single run.
          @par init : {‘k-means++’, ‘random’}
            Method for initialization, defaults to ‘k-means++’:
          @par assign : {‘lloyds’, ‘lsh’}
            Method for assignment, defaults to ‘lloyds’:
          @par update : {‘mean’, ‘pam’}
            Method for update, defaults to 'mean'
        */
        Cluster(int no_clusters = 8, int max_iter = 300,
          std::string init = "k-means++", std::string assign = "lloyds",
          std::string update = "pam", uint8_t no_hf = 3, uint8_t no_ht = 5,
          uint8_t no_grids = 5) : no_clusters(no_clusters), max_iter(max_iter),
          init(init), assign(assign), update(update), k(no_hf),  L(no_ht),
          L_grid(no_grids) {}
        /**
          \brief  Class Cluster default destructor
        */
        ~Cluster() {
          if (assign == "range-lsh") delete[] window;
        }
        /** \brief Fit method stores dataset info for clustering
          @par[in] dc : curves given from dataset
					@par[in] dcl: a vector which store the length of each curve in the dataset
					@par[in] dco: a vector which store the offset of each curve in the dc
          @par[in] no_c : number of curves
        */
        void Fit(const std::vector<std::pair<T,T>>& dc, const std::vector<U> dc_ids,
          const std::vector<int>& dcl, const std::vector<int>& dco, const int& no_c) {

          dataset_curves = dc;
          dataset_curves_ids = dc_ids;
          dataset_curves_lengths = dcl;
          dataset_curves_offsets = dco;
          no_curves = no_c;
          /**
            If assignment step is executed using lsh range search, map dataset
            to corresponding lsh structures
          */
          if (assign == "range-lsh") {
            /* Compute delta parameter */
            delta = utils::ComputeDelta(dataset_curves, dataset_curves_lengths,
                                        dataset_curves_offsets);
            vectors_dim = 2 * *max_element(std::begin(dataset_curves_lengths),
                                           std::end(dataset_curves_lengths));
            /* Building no_grids Grids */
            for (size_t i = 0; i < L_grid; ++i) {
              grids.push_back(vectorization::Grid<T>(dataset_curves,
                                            dataset_curves_lengths,
                                            dataset_curves_offsets,
                                            no_curves, vectors_dim, factor * delta));
            }
            /* Vectorizing dataset */
            for (size_t i = 0; i < L_grid; ++i) {
              L_grid_dataset_vectors.push_back(grids[i].Vectorize());
            }
            /* Compute window */
            window = new double [L_grid];
            for (size_t i = 0; i < L_grid; ++i) {
              window[i] = utils::ComputeMean<T>(L_grid_dataset_vectors[i],
                                                vectors_dim, no_curves);
            }
            /* Building LSH structures */
            for (size_t i = 0; i < L_grid; ++i) {
              lsh_structures.
                push_back(search::vectors::LSH<T,U>(k, L, vectors_dim,
                                                    no_curves, window[i],
                                                    L_grid_dataset_vectors[i],
                                                    dataset_curves_ids));
            }
            /* Map each id from dataset_curves_ids to its index */
            for (size_t i = 0; i < no_curves; ++i) {
              map_id_to_index[dataset_curves_ids[i]] = i;
            }
          }
        }
        /**
          \brief Predict the closest cluster each sample in dataset belongs to
          returns a vector of centroids coordinates and a vector of vectors
          which stores in each position the indexes of dataset_curves assigned
          in this cluster
        */
        std::tuple<std::tuple<std::vector<std::pair<T,T>>,std::vector<int>,std::vector<int>>,
        std::vector<std::vector<size_t>>,double>
          Predict(void) {
          /* Start time measuring */
          auto start = high_resolution_clock::now();
          /* Declare types */
          std::tuple<std::vector<std::pair<T,T>>,
          std::vector<int>,std::vector<int>> centroids;
          std::tuple<std::vector<std::vector<size_t>>, std::vector<T>> clusters;
          /* At first initialize centroids */
					if (init == "random") {
						centroids = RandomInit(dataset_curves, dataset_curves_lengths,
																	 dataset_curves_offsets, no_curves, no_clusters);
					} else if (init == "k-means++") {
						centroids = ParkJunInit(dataset_curves, dataset_curves_lengths,
																		dataset_curves_offsets, no_curves, no_clusters);
					}
          /* Calculate clusters and update centroids max_iter times */
          for (size_t i = 0; i < max_iter; ++i) {
            /* Assigment step */
						if (assign == "lloyd") {
							clusters = LloydsAssignment(dataset_curves, centroids,
																				  dataset_curves_lengths,
																					dataset_curves_offsets,
																					no_curves, no_clusters);
						} else if (assign == "range-lsh") {
              /* Vectorize centroids */
              std::vector<std::pair<T,T>> centroids_curves = std::get<0>(centroids);
              std::vector<int> centroids_curves_lengths = std::get<1>(centroids);
              std::vector<int> centroids_curves_offsets = std::get<2>(centroids);
              std::vector<std::vector<double>> L_grid_centroids_vectors(L_grid);
              for (size_t i = 0; i < L_grid; ++i) {
                L_grid_centroids_vectors[i] = grids[i].Vectorize(no_clusters,
                                                centroids_curves,
                                                centroids_curves_lengths,
                                                centroids_curves_offsets);
              }
              std::vector<std::tuple<std::vector<std::vector<size_t>>,std::vector<T>>>
              clusters_array(L_grid);
              /* Execute reverse assignment */
              for (size_t i = 0; i < L_grid; ++i) {
                clusters_array[i] = cluster::assignment::
                                      ReverseAssignment<T,U>(
                                        L_grid_dataset_vectors[i],
                                        dataset_curves_ids,
                                        L_grid_centroids_vectors[i],
                                        no_curves, vectors_dim,
                                        no_clusters, &lsh_structures[i],
                                        map_id_to_index);
              }
              clusters = clusters_array[L_grid - 1];
						}
            /* Update step */
						if (update == "mean") {
							centroids = LloydsUpdate(dataset_curves, centroids, dataset_curves_lengths,
  																	 	 dataset_curves_offsets, no_curves, no_clusters,
  																	 	 std::get<0>(clusters), std::get<1>(clusters));
						} else if (update == "pam") {
							centroids = PAMUpdate(dataset_curves, centroids, dataset_curves_lengths,
																	 	dataset_curves_offsets, no_curves, no_clusters,
																	 	std::get<0>(clusters), std::get<1>(clusters));
						}
          }
          /* End time measuring */
          auto stop = high_resolution_clock::now();
          duration <double> total_time = duration_cast<duration<double>>(stop - start);
          // Return result in the form of a tuple
          return std::make_tuple(centroids,std::get<0>(clusters),total_time.count());
        }
        /** \brief Map each curve index to the corresponding cluster
          @par[in] clusters - std::vector<std::vector<size_t>> returned by Predict
          return: a map of each index to the the corresponding cluster
        */
        std::map<int,int> MapToClusters(std::vector<std::vector<size_t>> clusters) {
          std::map<int,int> map_result;
          int cluster_idx = 0;
          for (auto const& cluster: clusters) {
            for (auto const& object: cluster) {
              map_result[object] = cluster_idx;
            }
            cluster_idx++;
          }
          return map_result;
        }
    };
  }
}

#endif
