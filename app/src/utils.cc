#include <iostream>

#include "../includes/utils.h"

void utils::ShowUsage(const std::string &program_name,
  const struct InputInfo &input_info) {
  std::cerr << "usage: " << program_name << " [--help] [-i INPUT FILE]"
            << " [-c CONFIGURATION FILE] [-o OUTPUT FILE] [--complete COMPLETE]"
            << "\n\n"
            << "optional arguments:\n"
            << "  --help\tshow this help message and exit\n"
            << "  -d\tdefine the input file\n"
            << "  -c\tdefine the configuration file\n"
            << "  --complete\toutput each cluster analytically (default = false)\n"
            << "  -o\tdefine the output file\n"
            << std::endl;
  exit(EXIT_FAILURE);
}

void utils::InputInfo::Print(const std::string idf) {
  std::cout << std::endl;
  std::cout << "Input file: " << input_file << std::endl;
  std::cout << "Configuration file: " << config_file << std::endl;
  std::cout << "Output file: " << output_file << std::endl;
  if (idf == "vectors") {
    std::cout << "Number of dataset vectors: "
              << static_cast<unsigned int>(N) << std::endl;
    std::cout << "Dataset vectors dimension: "
              << static_cast<unsigned int>(D) << std::endl;
  } else if (idf == "curves") {
    std::cout << "Number of dataset curves: "
              << static_cast<unsigned int>(N) << std::endl;
  }
  std::cout << "Number of clusters: "
            << static_cast<unsigned int>(K) << std::endl;
  std::cout << "Number of grids: "
            << static_cast<unsigned int>(grids) << std::endl;
  std::cout << "Number of vector hash tables: "
            << static_cast<unsigned int>(L) << std::endl;
  std::cout << "Number of vector hash functions: "
            << static_cast<unsigned int>(k) << std::endl;
  std::cout << "Initialization algorithm: " << init << std::endl;
  std::cout << "Assignment algorithm: " << assign << std::endl;
  std::cout << "Update algorithm: " << update << std::endl;


}
