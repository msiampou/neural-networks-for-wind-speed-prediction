#ifndef CLUSTER_UTILS
#define CLUSTER_UTILS

#include <tuple>
#include <vector>

namespace utils {
  /* enumerated general exit codes */
  typedef enum ExitCode{
    SUCCESS,
    FAIL,
    NO_ARGS,
    INVALID_L,
    INVALID_k,
    INVALID_grids,
    INVALID_K,
    INVALID_PARARAMETERS,
    INVALID_DATASET,
    INVALID_CONFIG,
    INVALID_OUTPUT,
    DATASET_ERROR,
    CONFIG_ERROR,
    MEMORY_ERROR
  } ExitCode;
  /**
    InputInfo - Group all input parameters of LSH in a struct
  */
  struct InputInfo {
    std::string input_file;      // name of the relative path to the input file
    std::string config_file;     // name of the relative path to the configuration file
    std::string output_file;     // name of the relative path to the output file
    uint8_t L = 3;               // number of vector hash tables
    uint8_t k = 4;               // number of vector hash functions
    uint8_t grids = 2;           // number of grids
    uint8_t K;                   // number of clusters
    uint32_t N;                  // number of dataset points
    uint16_t D;                  // dimension of dataset points
    bool complete = false;       // flag to ouput each cluster analytically or not
    /* Set default string names of clustering step algortithms */
    std::string init = "k-means++";
    std::string assign = "lloyd";
    std::string update = "mean"; 
    /* Print method of struct InputInfo */
    void Print(const std::string);
  };
  /** \brief ShowUsage - Prints the usage of the program
    @par const std::string &name - Pass by reference the name of the program
    @par const struct InputInfo &input_info - Pass by reference the input parameters
  */
  void ShowUsage(const std::string &name, const struct InputInfo &input_info);
}

#endif
