#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <getopt.h>
#include <unistd.h>
#include <string>
#include <string.h>

#include "../includes/utils.h"
#include "../includes/args_utils.h"

int utils::args::ScanArguments(struct InputInfo &input_info,
  utils::ExitCode &status) {

  std::string input_buffer;

  std::cout << "Provide the relative path for the input file: ";
  std::cin >> input_buffer;
  input_info.input_file = input_buffer;

  std::cout << "Provide the relative path for the configuration file: ";
  std::cin >> input_buffer;
  input_info.config_file = input_buffer;

  std::cout << "Provide the relative path for the output file: ";
  std::cin >> input_buffer;
  input_info.output_file = input_buffer;

  do {
    std::cout << "Do you want to output each cluster analytically (y/n)? : ";
    std::cin >> input_buffer;
    if (input_buffer != "y" && input_buffer != "n") {
      std::cout << "Wrong input! Try again." << std::endl;
    }
  } while (input_buffer != "y" && input_buffer != "n");

  if (input_buffer != "n") {
    input_info.complete = false;
  } else {
    input_info.complete = true;
  }

  do {
    std::cout << "Do you want to declare algorithms to be used for clustering (y/n)? : ";
    std::cin >> input_buffer;
    if (input_buffer != "y" && input_buffer != "n") {
      std::cout << "Wrong input! Try again." << std::endl;
    }
  } while (input_buffer != "y" && input_buffer != "n");

  if (input_buffer != "n") {
    do {
      std::cout << "Provide initialization algorithm {k-means++,random}: ";
      std::cin >> input_buffer;
      if (input_buffer != "k-means++" && input_buffer != "random") {
        std::cout << "Wrong input! Try again." << std::endl;
      }
    } while (input_buffer != "k-means++" && input_buffer != "random");
    input_info.init = input_buffer;

    do {
      std::cout << "Provide assignment algorithm {lloyd,range-lsh}: ";
      std::cin >> input_buffer;
      if (input_buffer != "lloyd" && input_buffer != "range-lsh") {
        std::cout << "Wrong input! Try again." << std::endl;
      }
    } while (input_buffer != "lloyd" && input_buffer != "range-lsh");
    input_info.assign = input_buffer;

    do {
      std::cout << "Provide update algorithm {pam,mean}: ";
      std::cin >> input_buffer;
      if (input_buffer != "pam" && input_buffer != "mean") {
        std::cout << "Wrong input! Try again." << std::endl;
      }
    } while (input_buffer != "pam" && input_buffer != "mean");
    input_info.update = input_buffer;
  }

  return SUCCESS;
}

int utils::args::ReadArguments(int argc, char **argv,
  struct InputInfo &input_info, utils::ExitCode &status) {

  if (argc == 1) {
    status = NO_ARGS;
    return FAIL;
  }

  if (argc == 2) {
    if (!strcmp(argv[1],"--help")) {
      ShowUsage(argv[0], input_info);
    }
  }

  const char * const short_opts = "i:c:o:0:1:2:3:";
  const option long_opts[] = {
           {"input", required_argument, nullptr, 'i'},
           {"config", required_argument, nullptr, 'c'},
           {"output", required_argument, nullptr, 'o'},
           {"complete", no_argument, nullptr, 0},
           {"init", required_argument, nullptr, 1},
           {"assign", required_argument, nullptr, 2},
           {"update", required_argument, nullptr, 3},
           {nullptr, no_argument, nullptr, 0}
   };

   while (true) {

     const auto opt = getopt_long(argc, argv, short_opts, long_opts, NULL);

     if (-1 == opt) {
       break;
     }

     switch (opt) {
       case 'i': {
        input_info.input_file = optarg;
        break;
      }
      case 'c' : {
        input_info.config_file = optarg;
        break;
      }
      case 'o': {
        input_info.output_file = optarg;
        break;
      }
      case 0: {
        input_info.complete = true;
        break;
      }
      case 1: {
        input_info.init = optarg;
        break;
      }
      case 2: {
        input_info.assign = optarg;
        break;
      }
      case 3: {
        input_info.update = optarg;
        break;
      }
      case '?':
        break;
      default:
        abort();
    }
  }
  return SUCCESS;
}
