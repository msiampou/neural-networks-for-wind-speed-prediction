#ifndef CLUSTER_IO_UTILS
#define CLUSTER_IO_UTILS

#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <algorithm>

#include "../../core/utils/utils.h"


namespace utils {
  namespace io {
    /** \brief GetFirstLine - Reads first line of the input file
      @par[in] file_name - Relative path to the dataset
      @par[out] first_line - First line string to be returned
      @par[out] ExitCode &statues - enumerated ExitCode provided from namespace utils
      return: SUCCESS or FAIL
    */
    int GetFirstLine(std::string &file_name, std::string &first_line,
      utils::ExitCode &status);
    /** \brief ReadConfig - Reads configuration file
      @par[in] file_name - Relative path to the dataset
      @par[out] no_clusters - Total number of clusters
      @par[out] no_grids - Total number of grids used by grid LSh
      @par[out] no_hf - Total number of LSH hash function
      @par[out] no_ht - Total number of LSH hash tables
      @par[out] ExitCode &statues - enumerated ExitCode provided from namespace utils
      return: SUCCESS or FAIL
    */
    int ReadConfig(std::string &file_name, uint8_t& no_clusters,
      uint8_t& no_grids, uint8_t& no_hf, uint8_t& no_ht, utils::ExitCode &status);
    /* Vectors I/O utils */
    namespace vectors {
      /** \brief ReadFile - Reads file provided by user
        @par[in] std::string file_name - Pass by reference the path to the input file
        @par[out] std::vector<T> &vectors - Pass by reference a vector type T which
             represent the N vectors of dimension D
        @par[out] std::vector<K> &ids - Pass by reference a vector type K which stores
             vectors' ids
        @par[in] const int no_vectors - Number of point in file
        @par[in] const int dim - Vectors' dimension
        @par[out] ExitCode &statues - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      template <typename T, typename U>
      int ReadFile(const std::string& file_name, const int& no_vectors,
        const int& dim, std::vector<T>& vectors, std::vector<U>& ids,
        utils::ExitCode &status) {
        // Open file
        std::ifstream infile;
        infile.open(file_name);
        // Check if file exists
        if (!infile) {
          status = INVALID_DATASET;
          return FAIL;
        }
        // Read file
        for (size_t i = 0; i < no_vectors && infile; ++i) {
          infile >> ids[i];
          for (size_t j = 0; j < dim; ++j) {
            infile >> vectors[i * dim + j];
          }
        }
        // close the file
        infile.close();
        return SUCCESS;
      }
      /** \brief WriteFile - Output prorgam results to the given output file
        @par[in] file_name -  Relative path to the dataset
        @par[in] init - init method
        @par[in] assign - assignment method
        @par[in] update - update method
        @par[in] cluster_res - tuple of centroids, assigned vectors to clusters
          and total clustering time
        @par[in] silhouette_res - pair of Silhouette of each vector and the total one
        @par[in] complete - complete parameter determines either to print
          clusters' objects or not
        @par[in] ids - dataset vectors ids
        @par[in] vectors_dim - all vectors and centroids have the same dimension
        @par[out] status - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      template <typename T, typename U>
      int WriteFile(const std::string& file_name, const std::string& init,
        const std::string& assign, const std::string& update,
        std::tuple<std::vector<T>,std::vector<std::vector<size_t>>,double> cluster_res,
        std::pair<std::vector<double>,double> silhouette_res, bool complete,
        const std::vector<U>& ids, const int& vectors_dim, utils::ExitCode &status) {

        // Open file
        std::ofstream outfile;
        outfile.open(file_name);
        // Check if file is opened
        if (outfile.is_open()) {
          /* Display clustering methods */
          if (init == "k-means++") {
            outfile << "Initialization: K-Means++ | ";
          } else if (init == "random") {
            outfile << "Initialization: Random | ";
          }
          if (assign == "lloyd") {
            outfile << "Assignment: Lloyd's Assignment | ";
          } else if (assign == "range-lsh") {
            outfile << "Assignment: By Range LSH | ";
          }
          if (update == "pam") {
            outfile << "Update: PAM" << std::endl;
          } else if (update == "mean") {
            outfile << "Update: Mean vector" << std::endl;
          }
          /* Extract result info */
          std::vector<T> centroids = std::get<0>(cluster_res);
          auto clusters = std::get<1>(cluster_res);
          double clustering_time = std::get<2>(cluster_res);
          /* Print for each cluster its data info */
          int cl_idx = 0;
          for (const auto& cluster: clusters) {
            outfile << "CLUSTER-" << cl_idx + 1 << " {size: " << cluster.size();
            outfile << " centroid: ";
            for (size_t i = 0; i < vectors_dim; ++i) {
              outfile << centroids[cl_idx * vectors_dim + i] << " ";
            }
            outfile << "}" << std::endl;
            cl_idx++;
          }
          /* Print total clustering time as well as the Silhouette metrics */
          outfile << "clustering_time: " << clustering_time << " seconds"
                  << std::endl;
          std::vector<double> s = std::get<0>(silhouette_res);
          double s_total = std::get<1>(silhouette_res);
          outfile << "Silhouette: {";
          for (size_t i = 0; i < clusters.size(); ++i) {
            outfile << "s" << i + 1 << ": " << s[i] << ", ";
          }
          outfile << "stotal: " << s_total << "}"<< std::endl;
          /* If complete parameter was given print cluster explicitly */
          if (complete == true) {
            cl_idx = 0;
            for (const auto& cluster: clusters) {
              outfile << "CLUSTER-" << cl_idx + 1 << " {";
              bool first = true;
              for (const auto& object_idx: cluster) {
                if (first == true) {
                  outfile << ids[object_idx];
                  first = false;
                } else {
                  outfile <<  ", " << ids[object_idx];
                }
              }
              outfile << "}" << std::endl;
              cl_idx++;
            }
          }
        } else {
          status = INVALID_OUTPUT;
          return FAIL;
        }
        // close file
        outfile.close();
        return SUCCESS;
      }
      /** \brief GetNoDatasetVectors - Get the number of vectors in the dataset
        @par[in] file_name - Relative path to the dataset
        @par[out] no_vectors - Total number of vectors in the dataset to be returned
        @par[out] ExitCode &statues - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      int GetNoDataVectors(std::string &file_name, uint32_t &no_vectors,
         ExitCode &status);
      /** \brief GetVectorsDim - Get the dimension of the vectors in the dataset
        @par[in] std::string &file_name - Path to the dataset
        @par[out] int &dim - Pass by reference the vectors' dimension
        @par[out] ExitCode &statues - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      int GetVectorsDim(std::string &file_name, uint16_t &dim, ExitCode &status);
    }
    /* Curves I/O utils */
    namespace curves {
      /** \brief ReadFile - Reads file provided by user
        @par[in] std::string& file_name - Pass by reference the path to the input file
        @par[out] std::vector<T>& curves - Pass by reference a vector type T which
             represent the N curves of length m_i, i = 1..N
        @par[out] std::vector<K>& ids - Pass by reference a vector type K which stores
             curves' ids
        @par[out] std::vector<T>& lengths - Pass by reference a vector type int which
             stores curves' lenths
        @par[out] std::vector<T>& offsets - Pass by reference a vector type int which
             stores offesets to the std::vector<std::pair<T,T>>& curves
        @par[in] const int no_curves - Number of point in file
        @par[out] ExitCode& status - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      template <typename T, typename U>
      int ReadFile(const std::string& file_name, const int& no_curves,
        std::vector<std::pair<T,T>>& curves, std::vector<U>& ids,
        std::vector<int>& lengths, std::vector<int>& offsets,
        utils::ExitCode& status) {

        // Open file
        std::ifstream infile;
        infile.open(file_name);
        // Check if file exists
        if (!infile) {
          status = INVALID_DATASET;
          return FAIL;
        }
        // skip first line
        std::string temp_fl;
        infile >> temp_fl;
        // Read file
        std::string buffer_x;
        std::string buffer_y;
        int offset{};
        for (size_t i = 0; i < no_curves && infile; ++i) {
          infile >> ids[i];
          infile >> lengths[i];
          for (size_t j = 0; j < lengths[i]; ++j) {
            infile >> buffer_x;
            infile >> buffer_y;
            // remove noise characters such as '(', to convert successuflly to T
            buffer_x.erase(std::remove(buffer_x.begin(),buffer_x.end(),'('),
                           buffer_x.end());
            T point_x = convert_to<T>(buffer_x);
            T point_y = convert_to<T>(buffer_y);
            curves.push_back(std::make_pair(point_x,point_y));
          }
          offsets[i] = offset;
          offset += lengths[i];
        }
        // close the file
        infile.close();
        return SUCCESS;
      }
      /** \brief WriteFile - Output prorgam results to the given output file
        @par[in] file_name -  Relative path to the dataset
        @par[in] init - init method
        @par[in] assign - assignment method
        @par[in] update - update method
        @par[in] cluster_res - tuple of centroids, assigned vectors to clusters
          and total clustering time
        @par[in] silhouette_res - pair of Silhouette of each vector and the total one
        @par[in] complete - complete parameter determines either to print
          clusters' objects or not
        @par[in] ids - dataset curves ids
        @par[out] status - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      template <typename T, typename U>
      int WriteFile(const std::string& file_name, const std::string& init,
        const std::string& assign, const std::string& update,
        std::tuple<std::tuple<
        std::vector<std::pair<T,T>>,std::vector<int>,std::vector<int>>,
        std::vector<std::vector<size_t>>,double> cluster_res,
        std::pair<std::vector<double>,double> silhouette_res, bool complete,
        const std::vector<U>& ids, utils::ExitCode &status) {

        // Open file
        std::ofstream outfile;
        outfile.open(file_name);
        // Check if file is opened
        if (outfile.is_open()) {
          /* Display clustering methods */
          if (init == "k-means++") {
            outfile << "Initialization: K-Means++ | ";
          } else if (init == "random") {
            outfile << "Initialization: Random | ";
          }
          if (assign == "lloyd") {
            outfile << "Assignment: Lloyd's Assignment | ";
          } else if (assign == "range-lsh") {
            outfile << "Assignment: By Range LSH | ";
          }
          if (update == "pam") {
            outfile << "Update: PAM" << std::endl;
          } else if (update == "mean") {
            outfile << "Update: Mean vector" << std::endl;
          }
          /* Extract result info */
          auto centroids = std::get<0>(cluster_res);
          auto clusters = std::get<1>(cluster_res);
          double clustering_time = std::get<2>(cluster_res);
          // Break centroids to its componenets
          std::vector<std::pair<T,T>> centroids_curves = std::get<0>(centroids);
          std::vector<int> centroids_lengths = std::get<1>(centroids);
          std::vector<int> centroids_offsets = std::get<2>(centroids);
          /* Print for each cluster its data info */
          int cl_idx = 0;
          for (const auto& cluster: clusters) {
            outfile << "CLUSTER-" << cl_idx + 1 << " {size: " << cluster.size();
            outfile << " centroid: ";
            for (size_t i = 0; i < centroids_lengths[cl_idx]; ++i) {
              std::pair<T,T> point = centroids_curves[centroids_offsets[cl_idx] + i];
              outfile << "(" << std::get<0>(point) << ", " << std::get<1>(point) << ")";
            }
            outfile << "}" << std::endl;
            cl_idx++;
          }
          /* Print total clustering time as well as the Silhouette metrics */
          outfile << "clustering_time: " << clustering_time << " seconds"
                  << std::endl;
          std::vector<double> s = std::get<0>(silhouette_res);
          double s_total = std::get<1>(silhouette_res);
          outfile << "Silhouette: {";
          for (size_t i = 0; i < clusters.size(); ++i) {
            outfile << "s" << i + 1 << ": " << s[i] << ", ";
          }
          outfile << "stotal: " << s_total << "}"<< std::endl;
          /* If complete parameter was given print cluster explicitly */
          if (complete == true) {
            cl_idx = 0;
            for (const auto& cluster: clusters) {
              outfile << "CLUSTER-" << cl_idx + 1 << " {";
              bool first = true;
              for (const auto& object_idx: cluster) {
                if (first == true) {
                  outfile << ids[object_idx];
                  first = false;
                } else {
                  outfile <<  ", " << ids[object_idx];
                }
              }
              outfile << "}" << std::endl;
              cl_idx++;
            }
          }
        } else {
          status = INVALID_OUTPUT;
          return FAIL;
        }
        // close file
        outfile.close();
        return SUCCESS;
      }
      /** \brief GetNoDataCurves - Get the number of curves in the dataset
        @par[in] file_name - Relative path to the dataset
        @par[out] no_vectors - Total number of curvess in the dataset to be returned
        @par[out] ExitCode &statues - enumerated ExitCode provided from namespace utils
        return: SUCCESS or FAIL
      */
      int GetNoDataCurves(std::string& file_name, uint32_t& no_curves,
        ExitCode& status);
    }
  }
}

#endif
