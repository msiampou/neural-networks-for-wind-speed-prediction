#ifndef CLUSTER_ARGS_UTILS
#define CLUSTER_ARGS_UTILS

#include "./utils.h"

namespace utils {
  namespace args {
    /** \brief ReadArguments - Reads arguments given by user in the command line
      @par int argc - The number of arguments
      @par char **argv - The array of the arguments' values
      @par ExitCode &statues - enumerated ExitCode provided from namespace utils
      return: SUCCESS or FAIL
    */
    int ReadArguments(int argc, char **argv, struct InputInfo &input_info,
                      ExitCode &status);
    /** \brief ScanArguments - Reads arguments given by user in the stdin
      @par const ExitCode &statues - enumerated ExitCode provided from
        namespace utils
      return: SUCCESS or FAIL
    */
    int ScanArguments(struct InputInfo &input_info, ExitCode &status);
  }
}

#endif
