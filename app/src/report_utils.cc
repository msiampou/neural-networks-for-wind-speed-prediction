#include <iostream>

#include "../includes/utils.h"
#include "../includes/report_utils.h"

void utils::report::ReportError(const utils::ExitCode &code) {
  switch (code) {
    case INVALID_L:
      std::cout << "\n[ERROR]: " << "Invalid value of vector hash tables"
                << std::endl;
      break;
    case INVALID_k:
      std::cout << "\n[ERROR]: " << "Invalid value of vector hash functions"
                << std::endl;
      break;
    case INVALID_grids:
      std::cout << "\n[ERROR]: " << "Invalid value of grids" << std::endl;
      break;
    case INVALID_K:
      std::cout << "\n[ERROR]: " << "Invalid value of clusters" << std::endl;
      break;
    case INVALID_PARARAMETERS:
      std::cout << "\n[ERROR]: " << "Invalid parameters given" << std::endl;
      break;
    case INVALID_DATASET:
      std::cout << "\n[ERROR]: " << "Invalid dataset file name" << std::endl;
      break;
    case INVALID_CONFIG:
      std::cout << "\n[ERROR]: " << "Invalid configuration file name"
                << std::endl;
      break;
    case INVALID_OUTPUT:
      std::cout << "\n[ERROR]: " << "Invalid output file name" << std::endl;
      break;
    case DATASET_ERROR:
      std::cout << "\n[ERROR]: " << "Invalid dataset file format" << std::endl;
      break;
    case CONFIG_ERROR:
      std::cout << "\n[ERROR]: " << "Invalid configuration file format"
                << std::endl;
      break;
    case MEMORY_ERROR:
      std::cout << "\n[ERROR]: " << "System ran out of memory" << std::endl;
      break;
    default:
      abort();
  }
  exit(EXIT_FAILURE);
}
