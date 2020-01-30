#include <exception>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>

#include <stdio.h>
#include <stdlib.h>

#include "../includes/utils.h"
#include "../includes/io_utils.h"
#include "../includes/report_utils.h"

#include "../../core/utils/utils.h"

int utils::io::GetFirstLine(std::string &file_name,
  std::string &first_line, utils::ExitCode &status) {

  // Open file
  std::ifstream infile;
  infile.open(file_name);
  // Check if file exists
  if (!infile) {
    status = INVALID_DATASET;
    return FAIL;
  }
  infile >> first_line;
  // if (first_line != "vectors" && first_line != "curves") {
  //   status = DATASET_ERROR;
  //   return FAIL;
  // }
  // close the file
  infile.close();
  return SUCCESS;
}

int utils::io::ReadConfig(std::string &file_name, uint8_t& no_clusters,
  uint8_t& no_grids, uint8_t& no_hf, uint8_t& no_ht, utils::ExitCode &status) {

  // Open file
  std::ifstream infile;
  infile.open(file_name);
  // Check if file exists
  if (!infile) {
    status = INVALID_CONFIG;
    return FAIL;
  }

  std::string line;
  while (getline(infile, line)) {
    // using printf() in all tests for consistency
    std::vector<std::string> tokens = Split(line,' ');
    // Each line should has exactly two tokens
    if (tokens.size() != 2) {
      status = CONFIG_ERROR;
      return FAIL;
    }
    // Extract info
    if (tokens[0] == "number_of_clusters:") {
      try {
        no_clusters = std::stoi(tokens[1]);
      } catch (std::exception& e) {
        std::cout << tokens[0] << " " << "Standard exception: "
                  << e.what() << std::endl;
      }
    } else if (tokens[0] == "number_of_grids:") {
      try {
        no_grids = std::stoi(tokens[1]);
      } catch (std::exception& e) {
        std::cout << tokens[0] << " " << "Standard exception: "
                  << e.what() << std::endl;      }
    } else if (tokens[0] == "number_of_hash_functions:") {
      try {
        no_hf = std::stoi(tokens[1]);
      } catch (std::exception& e) {
        std::cout << tokens[0] << " " << "Standard exception: "
                  << e.what() << std::endl;      }
    } else if (tokens[0] == "number_of_hash_tables:") {
      try {
        no_ht = std::stoi(tokens[1]);
      } catch (std::exception& e) {
        std::cout << tokens[0] << " " << "Standard exception: "
                  << e.what() << std::endl;      }
    } else {
      std::cout << "Unknown token: " << tokens[0] << std::endl;
    }
  }
  // close the file
  infile.close();
  return SUCCESS;
}

int utils::io::vectors::GetNoDataVectors(std::string &file_name,
  uint32_t &no_vectors, utils::ExitCode &status) {

  FILE *fp;                 // To opem the file for reading
  /* Initialize count to -1 to skip first line which contains the identifier */
  unsigned int count = -1;  // Line counter (result)
  char c;                   // To store a character read from file
  // Open the file
  fp = fopen(file_name.c_str(), "r");
  // Check if file exists
  if (!fp) {
    status = INVALID_DATASET;
    return FAIL;
  }
  // Extract characters from file and store in character c
  for (c = getc(fp); c != EOF; c = getc(fp)) {
    if (c == '\n') { // Increment count if this character is newline
      count = count + 1;
    }
  }
  // Pass info
  no_vectors = count;
  // Close the file
  fclose(fp);
  // everything ok, return SUCCESS
  return SUCCESS;
}

int utils::io::vectors::GetVectorsDim(std::string &file_name, uint16_t &dim,
  utils::ExitCode &status) {

  std::ifstream file;
  file.open(file_name);
  std::string line;
  // Check if file exists
  if (!file) {
    status = INVALID_DATASET;
    return FAIL;
  }
  else {
    // Skip first line which contains the object identifier
    std::string temp;
    getline(file, temp);
    // Get the next one line
    if (getline(file, line)) {
      // Get total number of elements
      int count =  std::distance(std::istream_iterator<std::string>(
                    std::istringstream(line) >> std::ws),
                    std::istream_iterator<std::string>());
      // The first one is the vector's id. The remaining determines the dimension
      dim = count - 1;
    }
  }
  // close the file
  file.close();
  // everything ok, return SUCCESS
  return SUCCESS;
}

int utils::io::curves::GetNoDataCurves(std::string &file_name,
  uint32_t &no_curves, utils::ExitCode &status) {

  FILE *fp;                 // To opem the file for reading
  /* Initialize count to -1 to skip first line which contains the identifier */
  unsigned int count = -1;  // Line counter (result)
  char c;                   // To store a character read from file
  // Open the file
  fp = fopen(file_name.c_str(), "r");
  // Check if file exists
  if (!fp) {
    status = INVALID_DATASET;
    return FAIL;
  }
  // Extract characters from file and store in character c
  for (c = getc(fp); c != EOF; c = getc(fp)) {
    if (c == '\n') { // Increment count if this character is newline
      count = count + 1;
    }
  }
  // Pass info
  no_curves = count;
  // Close the file
  fclose(fp);
  // everything ok, return SUCCESS
  return SUCCESS;
}
