#ifndef VECTORIZATION
#define VECTORIZATION

#include <set>
#include <cmath>
#include <tuple>
#include <vector>
#include <random>
#include <iostream>

#include "../lib/map_hash.h"

namespace vectorization {

  template <typename T>
  class Grid {
    private:
      uint32_t N;
      uint32_t D;

      const double delta;
      std::pair<double,double> t;

      const std::vector<std::pair<T,T>>& input_curves;
      const std::vector<int>& input_curves_lengths;
      const std::vector<int>& input_curves_offsets;

      std::default_random_engine generator;
      std::uniform_real_distribution<double> distribution;
    public:
      /** \brief Grid class constructor
        Initializing private members
      */
      Grid(const std::vector<std::pair<T,T>>& curves,
        const std::vector<int>& lengths, const std::vector<int>& offsets,
        const uint32_t N, const uint32_t D, double delta) :
          N(N), D(D), delta(delta), distribution(0,delta),
          generator(std::chrono::system_clock::now().time_since_epoch().count()),
          input_curves(curves), input_curves_lengths(lengths),
          input_curves_offsets(offsets) {
          // t is selected uniformly between 0 - delta
          t.first = distribution(generator);
          t.second = distribution(generator);
      }
      /**
        \brief class Grid default constructor
      */
      ~Grid() = default;
      /** \brief Given the input curves transform each one to an equivalent
        vector. Each vector is computed as follows:
        1) Gδ = {(a1δ,...,adδ)|∀ai ∈ Z}.
        2) Define family of shifted grids parameterized by t ∈ Rd :
            Gδt ={p+t|p∈Gδ}.
        3) Choose real vector t uniformly ∈ [0,delta)^d.
            For curve P = p1, . . . , pm, hash-function hδt (P) is:
              - foreach pi set pi′ = argminq∈Gδt||pi−q||∈Gδt(snap),
              - removeconsecutiveduplicatesinp1′,...,pm′
        4) Then hδt (P) = the resulting polygonal curve
      */
      std::vector<double> Vectorize(void) {
        /**
          Vectorize each input curve and store the corresponding vector to
          an 1D array. Each vector is of dimension D (max curve length) and
          we have in total_time N curves, and so N vectors.
        */
        std::vector<double> result(D * N);
        /* Iterave over each curve to compute its vector */
        for (size_t i = 0; i < N; ++i) {
          size_t idx = 0;
          for (size_t j = 0; j < input_curves_lengths[i]; ++j) {
            double x_1 = std::get<0>(input_curves[input_curves_offsets[i] + j]);
            double t_1 = std::get<0>(t);
            double x_2 = std::get<1>(input_curves[input_curves_offsets[i] + j]);
            double t_2 = std::get<1>(t);
            double a_1 = round((x_1 - t_1) / delta);
            double a_2 = round((x_2 - t_2) / delta);
            double s_1 = a_1 * delta + t_1;
            double s_2 = a_1 * delta + t_1;
            result[i * D + (idx++)] = s_1;
            result[i * D + (idx++)] = s_2;
          }
          // remove duplicates
          std::vector<double>::iterator ip;
          ip = std::unique(result.begin() + i * D, result.begin() + i * D + idx);
          result.resize(std::distance(result.begin() + i * D, ip));
          // Fill with pading coordinates to have equal length vectors
          size_t start = std::distance(result.begin() + i * D, ip);
          for (size_t j = start; j < D; ++j) {
            result[i * D + j] = std::numeric_limits<T>::max();
          }
        }
        return result;
      }
      /**
        \brief Given the query curves perform vectorization as above
      */
      std::vector<double> Vectorize(const int Q,
        const std::vector<std::pair<T,T>>& query_curves,
        const std::vector<int>& query_curves_lengths,
        const std::vector<int>& query_curves_offsets) {
        /**
          Vectorize each query curve and store the corresponding vector to
          an 1D array. Each vector is of dimension D (max curve length) and
          we have in total_time Q query curves, and so Q vectors.
        */
        std::vector<double> result(D * Q);
        /* Iterave over each curve to compute its vector */
        for (size_t i = 0; i < Q; ++i) {
          size_t idx = 0;
          for (size_t j = 0; j < query_curves_lengths[i]; ++j) {
            double x_1 = std::get<0>(query_curves[query_curves_offsets[i] + j]);
            double t_1 = std::get<0>(t);
            double x_2 = std::get<1>(query_curves[query_curves_offsets[i] + j]);
            double t_2 = std::get<1>(t);
            double a_1 = round((x_1 - t_1) / delta);
            double a_2 = round((x_2 - t_2) / delta);
            double s_1 = a_1 * delta + t_1;
            double s_2 = a_1 * delta + t_1;
            result[i * D + (idx++)] = s_1;
            result[i * D + (idx++)] = s_2;
          }
          // remove duplicates
          std::vector<double>::iterator ip;
          ip = std::unique(result.begin() + i * D, result.begin() + i * D + idx);
          result.resize(std::distance(result.begin() + i * D, ip));
          // Fill with pading coordinates to have equal length vectors
          size_t start = std::distance(result.begin() + i * D, ip);
          // Fill with pading coordinates to have eqaal length vectors
          for (size_t j = start; j < D; ++j) {
            result[i * D + j] = std::numeric_limits<T>::max();
          }
        }
        return result;
      }
  };

  template <typename T, typename U>
  class Projection  {
    private:
      int d;
      uint32_t M;
      uint32_t N;
      float eps;
      int K;
      std::vector<double> G;
      const std::vector<std::pair<T,T>>& input_curves;
      const std::vector<int>& input_curves_lengths;
      const std::vector<int>& input_curves_offsets;
      const std::vector<U>& input_curves_ids;
      /* M * M array containing all possible paths from (0,0) to each cell */
      std::vector<std::vector<std::vector<std::pair<T,T>>>> relevant_traversals;
      /* Storing datasets vectors' info per traversal */
      std::unordered_map<std::tuple<int,int,int>,std::vector<double>> vectors;
      std::unordered_map<std::tuple<int,int,int>,std::vector<int>> vectors_lengths;
      std::unordered_map<std::tuple<int,int,int>,std::vector<int>> vectors_offsets;
      std::unordered_map<std::tuple<int,int,int>,std::vector<U>> vectors_ids;
      /* Storing queries vectors' info per traversal */
      std::unordered_map<U,std::vector<double>> qvectors;
      std::unordered_map<U,std::vector<int>> qvectors_lengths;
      std::unordered_map<U,std::vector<int>> qvectors_offsets;
      std::unordered_map<U,std::vector<U>> qvectors_ids;

      std::default_random_engine generator;
      std::uniform_real_distribution<double> distribution;
    public:
      /**
        \brief Just a constructor
      */
      Projection(std::vector<std::pair<T,T>>& dataset_curves,
        std::vector<int>& dataset_offsets, std::vector<int>& dataset_lengths,
        std::vector<U>& dataset_ids, uint32_t N, const int K) :
          d(2), N(N), K(K), distribution(0,1),
          generator(std::chrono::system_clock::now().time_since_epoch().count()),
          input_curves(dataset_curves), input_curves_lengths(dataset_lengths),
          input_curves_offsets(dataset_offsets), input_curves_ids(dataset_ids) {

        // Get max length from all curves and
        // store all relevant traversals at an M * M array
        M = *max_element(std::begin(dataset_lengths), std::end(dataset_lengths));

        // Computing all relevant traversals
        relevant_traversals = std::vector<std::vector<std::vector<std::pair<T,T>>>> (M * M);
        int count = 0;
        for (size_t i = 0; i < M; ++i) {
          for (size_t j = 0; j < M; ++j) {
            if (abs(i-j < 4)) {
              relevant_traversals[i * M + j] = RelevantTraversals(i,j);
              count += relevant_traversals[i * M + j].size();
            }
          }
        }

        /* Generating G matrix with
         * random values ~N(0,1) */
        size_t size_G = K * d;
        G = std::vector<double> (size_G);
        for (size_t i = 0; i < size_G; ++i) {
          G[i] = distribution(generator);
        }
      };

      /**
       * For every dataset curve replace 1st pair pointer (ui)
       * with actual coordinate. For every value multiply with G
       * matrix, then concat: x = [G*u1|..|G*un]. Store each
       * vector and its info to map with key the corresponding traversal.
       * Each traversal corresponds to a NN structure.
       */
      void Vectorize (void) {
        for (size_t i = 0; i < N; ++i) {
          size_t length = input_curves_lengths[i] - 1;
          for (size_t j = 0; j < M; ++j) {
            size_t i_traversal = 0;
            for (const auto& tr:relevant_traversals[length*M+j]) {
              std::vector<std::pair<T,T>> rep_curve;
              i_traversal++;
              for (const auto& pair:tr) {
                T pos = pair.first;
                rep_curve.push_back(std::make_pair(input_curves[input_curves_offsets[i] + pos].first,pair.second));
              }
              std::tuple<size_t,size_t,size_t> key = std::make_tuple(length,j,i_traversal);
              std::vector<T> value = CreateVector(rep_curve);
              for (const auto& ivalue:value) {
                vectors[key].push_back(ivalue);
              }
              vectors_lengths[key].push_back(input_curves_lengths[i]);
              vectors_offsets[key].push_back(input_curves_offsets[i]);
              vectors_ids[key].push_back(input_curves_ids[i]);
            }
          }
        }
      };

      /**
        \brief Perfom exactly the opossite procedure as above
      */
      void Vectorize(const int Q,
        const std::vector<std::pair<T,T>>& query_curves,
        const std::vector<int>& query_curves_lengths,
        const std::vector<int>& query_curves_offsets,
        const std::vector<U>& query_curves_ids) {
          for (size_t i = 0; i < Q; ++i) {
            size_t length = query_curves_lengths[i] - 1;
            for (size_t j = 0; j < M; ++j) {
              size_t i_traversal = 0;
              for (const auto& tr:relevant_traversals[j*M+length]) {
                std::vector<std::pair<T,T>> rep_curve;
                i_traversal++;
                for (const auto& pair:tr) {
                  T pos = pair.second;
                  rep_curve.push_back(std::make_pair(pair.first,input_curves[input_curves_offsets[i] + pos].second));
                }
                U key = query_curves_ids[i];
                std::vector<T> value = qCreateVector(rep_curve);
                for (const auto& ivalue:value) {
                  qvectors[key].push_back(ivalue);
                }
                qvectors_lengths[key].push_back(query_curves_lengths[i]);
                qvectors_offsets[key].push_back(query_curves_offsets[i]);
                qvectors_ids[key].push_back(query_curves_ids[i]);
              }
            }
          }
      };

      const std::unordered_map<std::tuple<int,int,int>,
                              std::vector<double>>& GetVectors() {
                              return vectors;
                        };

      const std::unordered_map<std::tuple<int,int,int>,
                              std::vector<int>>& GetVectorsLengths() {
                              return vectors_lengths;
                        };

      const std::unordered_map<std::tuple<int,int,int>,
                              std::vector<int>>& GetVectorsOffsets() {
                              return vectors_offsets;
                        };

      const std::unordered_map<std::tuple<int,int,int>,
                              std::vector<U>>& GetVectorsIds() {
                              return vectors_ids;
                        };

      const std::unordered_map<U,
                              std::vector<double>>& qGetVectors() {
                              return qvectors;
                        };

      const std::unordered_map<U,
                              std::vector<int>>& qGetVectorsLengths() {
                              return qvectors_lengths;
                        };

      const std::unordered_map<U,
                              std::vector<int>>& qGetVectorsOffsets() {
                              return qvectors_offsets;
                        };

      const std::unordered_map<U,
                              std::vector<U>>& qGetVectorsIds() {
                              return qvectors_ids;
                        };

      std::vector<T> CreateVector(std::vector<std::pair<T,T>>& traversal) {
        std::vector<std::vector<T>> v;
        std::vector<T> x (K*d);
        for (const auto& ui:traversal) {
          std::vector<T> xi;
          for (size_t i = 0; i < K * d; ++i) {
            xi.push_back(ui.first * G[i]);
          }
          v.push_back(xi);
        }
        for (size_t i = 0; i < v.size(); ++i) {
          for (size_t j = 0; j < v[i].size(); ++j) {
            x[j] += v[i][j];
          }
        }
        return x;
      };

      std::vector<T> qCreateVector(std::vector<std::pair<T,T>>& traversal) {
        std::vector<std::vector<T>> v;
        std::vector<T> x (K*d);
        for (const auto& ui:traversal) {
          std::vector<T> xi;
          for (size_t i = 0; i < K * d; ++i) {
            xi.push_back(ui.second * G[i]);
          }
          v.push_back(xi);
        }
        for (size_t i = 0; i < v.size(); ++i) {
          for (size_t j = 0; j < v[i].size(); ++j) {
            x[j] += v[i][j];
          }
        }
        return x;
      };

      std::vector<std::vector<std::pair<T,T>>> RelevantTraversals(int i, int j) {
        std::set<std::pair<T,T>> diag_cells;
        std::set<std::pair<T,T>> rel_cells;
        diag_cells = DrawLineSegment(i, j, i+1, j+1);
        rel_cells = FindNeighbors(diag_cells,i,j);
        std::vector<std::pair<T,T>> path;
        std::vector<std::vector<std::pair<T,T>>> paths;
        FindTraversals(path, paths, rel_cells, 0, 0, i+1, j+1);
        return paths;
      }

      /** \brief
          Find cells that are crossed by main diagonal line segment
      */
      std::set<std::pair<T,T>> DrawLineSegment(int x, int y, int m, int n) {
        std::set<std::pair<T,T>> diag_cells;
        for (int xi = 0; xi < m; xi++) {
          int yi = round((double) n/m * (double)xi);
          diag_cells.insert(std::make_pair(xi,yi));
        }
        for (int yi = 0; yi < n; yi++) {
          int xi = round((double) m/n * (double)yi);
          diag_cells.insert(std::make_pair(xi,yi));
        }
        return diag_cells;
      }

      /** \brief
          Find cells that are crossed by main diagonal line segment
      */
      std::set<std::pair<T,T>> FindNeighbors(std::set<std::pair<T,T>>& diag_cells,
          int m, int n) {

        std::set<std::pair<T,T>> rel_cells (diag_cells);
        for (const auto& p:diag_cells) {
          if (p.first != 0) {
            rel_cells.insert(std::make_pair(p.first-1,p.second));
          }
        }
        return rel_cells;
      }

      /** \brief
          A path is relevant when it consists of cells that are either on
          main diagonal line segment or have a distance of 1.
      */
      const bool isRelevant(std::set<std::pair<T,T>>& s, int i, int j, int m, int n) {
        return s.find(std::make_pair(i,j)) != s.end();
      }

      void FindTraversals(std::vector<std::pair<T,T>>& path, std::vector<std::vector<std::pair<T,T>>>& paths,
        std::set<std::pair<T,T>>& s, int i, int j, int m, int n) {
        //destination point reached
        if ((i == m-1) && (j == n-1)) {
          path.push_back(std::make_pair(i,j));
          paths.push_back(path);
          path.pop_back();
          return;
        }
        //add curr cell to path
        path.push_back(std::make_pair(i,j));
        //move right
        if (isRelevant(s, i+1, j, m, n)) {
          FindTraversals(path, paths, s, i+1, j, m, n);
        }
        //move left
        if (isRelevant(s, i, j+1, m, n)) {
          FindTraversals(path, paths, s, i, j+1, m, n);
        }
        if (isRelevant(s, i+1, j+1, m, n)) {
          FindTraversals(path, paths, s, i+1, j+1, m, n);
        }
        //delete last pair
        path.pop_back();
      }
  };
}

#endif
