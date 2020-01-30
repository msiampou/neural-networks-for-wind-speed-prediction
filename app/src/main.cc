#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "../../core/cluster/assignment/assignment.h"
#include "../../core/cluster/initialization/initialization.h"
#include "../../core/cluster/clustering/clustering.h"
#include "../../core/cluster/update/update.h"
#include "../../core/metric/metric.h"
#include "../../core/utils/utils.h"

#include "../includes/args_utils.h"
#include "../includes/io_utils.h"
#include "../includes/report_utils.h"
#include "../includes/utils.h"

using namespace std::chrono;

// NOTE These values are harcoded independently for each expirement
#define MAX_ITER 30
#define U std::string
#define T double

int main(int argc, char **argv) {
  utils::InputInfo input_info;
  utils::ExitCode status;
  std::string input_buffer, clustering_object = "vectors";
  int exit_code;

  /* Get arguments */
  exit_code = utils::args::ReadArguments(argc, argv, input_info, status);
  switch (exit_code) {
    case utils::SUCCESS:
      std::cout << "\nArguments provided correctly" << std::endl;
      break;
    case utils::FAIL:
      if (status == utils::NO_ARGS) {
        std::cout << "\nNo arguments provided" << std::endl;
        std::cout << "Proceding to input them.." << std::endl;
        exit_code = utils::args::ScanArguments(input_info, status);
        switch (exit_code) {
          case utils::SUCCESS:
            std::cout << "Arguments provided correctly" << std::endl;
            break;
          case utils::FAIL:
            utils::report::ReportError(status);
            break;
          default:
            break;
        }
      } else {
        utils::report::ReportError(status);
      }
      break;
    default:
      break;
  }

  /* Read configuration file and pass parameters to the input_info struct */
  auto start = high_resolution_clock::now();
  std::cout << "\nReading configuration file.." << std::endl;
  exit_code = utils::io::ReadConfig(input_info.config_file, input_info.K,
                                    input_info.grids, input_info.L,
                                    input_info.k, status);
  if (exit_code != utils::SUCCESS) {
    utils::report::ReportError(status);
  }
  auto stop = high_resolution_clock::now();
  auto total_time = duration_cast<duration<double>>(stop - start);
  std::cout << "Reading configuration file completed successfully." << std::endl;
  std::cout << "Time elapsed: " << total_time.count() << " seconds"
            << std::endl;

  /**
   Preprocessing input file to get number of dataset vectors
   and their dimension
  */
  start = high_resolution_clock::now();
  std::cout << "\nGetting number of dataset vectors.." << std::endl;
  exit_code = utils::io::vectors::GetNoDataVectors(input_info.input_file,
                                                   input_info.N, status);
  if (exit_code != utils::SUCCESS) {
    utils::report::ReportError(status);
  }
  stop = high_resolution_clock::now();
  total_time = duration_cast<duration<double>>(stop - start);
  std::cout << "Getting number of dataset vectors completed successfully."
            << std::endl;
  std::cout << "Time elapsed: " << total_time.count() << " seconds"
            << std::endl;

  start = high_resolution_clock::now();
  std::cout << "\nGetting dataset vectors' dimension.." << std::endl;
  exit_code = utils::io::vectors::GetVectorsDim(input_info.input_file,
                                                input_info.D, status);
  if (exit_code != utils::SUCCESS) {
    utils::report::ReportError(status);
  }
  stop = high_resolution_clock::now();
  total_time = duration_cast<duration<double>>(stop - start);
  std::cout << "Getting dataset vectors' dimension completed successfully."
            << std::endl;
  std::cout << "Time elapsed: " << total_time.count() << " seconds"
            << std::endl;

  /**
    Read dataset and create 1D vector which represents the d-dimensional
    vectors of N vectors. Also create 1D vector that stores vectors' ids.
    1D vector of points representation support cache efficiency and as a result
    faster computations
  */
  start = high_resolution_clock::now();
  std::cout << "\nReading dataset.." << std::endl;
  std::vector<T> dataset_vectors(input_info.N * input_info.D);
  std::vector<U> dataset_vectors_ids(input_info.N);
  exit_code = utils::io::vectors::ReadFile<T,U>(input_info.input_file,
                                                input_info.N, input_info.D,
                                                dataset_vectors,
                                                dataset_vectors_ids, status);
  if (exit_code != utils::SUCCESS) {
    utils::report::ReportError(status);
  }
  stop = high_resolution_clock::now();
  total_time = duration_cast<duration<double>>(stop - start);
  std::cout << "Reading dataset completed successfully." << std::endl;
  std::cout << "Time elapsed: " << total_time.count() << " seconds"
            << std::endl;

  /* Print info */
  input_info.Print(clustering_object);
  /**
    Create a clustering object to perform dataset clustering using a
    specific algorithm, specified by init, assign and update string paramters
   */
  start = high_resolution_clock::now();
  std::cout << "\nBuilding Cluster class.." << std::endl;
  cluster::vectors::Cluster<T,U> *cl = new cluster::vectors::Cluster<T,U>
                                    {input_info.K, MAX_ITER, input_info.init,
                                    input_info.assign, input_info.update,
                                    input_info.k, input_info.L};
  stop = high_resolution_clock::now();
  total_time = duration_cast<duration<double>>(stop - start);
  std::cout << "Building Cluster class completed successfully." << std::endl;
  std::cout << "Time elapsed: " << total_time.count() << " seconds"
            << std::endl;

  /**
    Fit dataset into the Cluster class
    Note: In case of range LSH assignment, fitting may take more time as
    LSH structures will be constructed
  */
  start = high_resolution_clock::now();
  std::cout << "\nFitting dataset.." << std::endl;
  cl->Fit(dataset_vectors, dataset_vectors_ids, input_info.N, input_info.D);
  stop = high_resolution_clock::now();
  total_time = duration_cast<duration<double>>(stop - start);
  std::cout << "Fitting dataset completed successfully." << std::endl;
  std::cout << "Time elapsed: " << total_time.count() << " seconds"
            << std::endl;

  /* Compute clusters running specified algorithm */
  start = high_resolution_clock::now();
  std::cout << "\nComputing clusters.." << std::endl;
  auto clusters_res = cl->Predict();
  stop = high_resolution_clock::now();
  total_time = duration_cast<duration<double>>(stop - start);
  std::cout << "Computing clusters completed successfully." << std::endl;
  std::cout << "Time elapsed: " << total_time.count() << " seconds"
            << std::endl;

  /* Extract info */
  start = high_resolution_clock::now();
  std::cout << "\nExtracting cluster info.." << std::endl;
  std::vector<T> centroids = std::get<0>(clusters_res);
  std::vector<std::vector<size_t>> clusters = std::get<1>(clusters_res);
  std::map<int,int> mapped_vectors = cl->MapToClusters(clusters);
  stop = high_resolution_clock::now();
  total_time = duration_cast<duration<double>>(stop - start);
  std::cout << "Extracting cluster info completed successfully." << std::endl;
  std::cout << "Time elapsed: " << total_time.count() << " seconds"
            << std::endl;

  /* Compute Silhouette */
  std::pair<std::vector<double>,double> silhouette_res;
  start = high_resolution_clock::now();
  std::cout << "\nComputing Silhouette.." << std::endl;
  silhouette_res = metric::vectors::Silhouette<T>(dataset_vectors, input_info.N,
                                                  input_info.D, clusters,
                                                  centroids, mapped_vectors);
  stop = high_resolution_clock::now();
  total_time = duration_cast<duration<double>>(stop - start);
  std::cout << "Computing Silhouette completed successfully." << std::endl;
  std::cout << "Time elapsed: " << total_time.count() << " seconds"
            << std::endl;

  /* Writing results to the output file */
  start = high_resolution_clock::now();
  std::cout << "\nWriting results to the output file.." << std::endl;
  exit_code = utils::io::vectors::WriteFile<T,U>(input_info.output_file,
    input_info.init, input_info.assign, input_info.update, clusters_res,
    silhouette_res, input_info.complete, dataset_vectors_ids, input_info.D,
    status);
  if (exit_code != utils::SUCCESS) {
    utils::report::ReportError(status);
  }
  stop = high_resolution_clock::now();
  total_time = duration_cast<duration<double>>(stop - start);
  std::cout << "Writing results to the output file completed successfully."
            << std::endl;
  std::cout << "Time elapsed: " << total_time.count() << " seconds"
            << std::endl;
  delete cl;

  return EXIT_SUCCESS;
}
