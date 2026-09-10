#include <iostream>

#include <catch2/catch_session.hpp>
#include <mpi.h>

int main(int argc, char* argv[])
{
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    std::cerr << "Failed to initialize MPI." << std::endl;
    return -1;
  }

  Catch::Session session;

  // Ensure tests are run in declared order to avoid randomness across ranks
  session.configData().runOrder = Catch::TestRunOrder::Declared;

  int returnCode = session.applyCommandLine(argc, argv);
  if (returnCode != 0) { // command line error
    MPI_Finalize();
    return returnCode;
  }

  int result = session.run();

  if (MPI_Finalize() != MPI_SUCCESS) {
    std::cerr << "Failed to finalize MPI." << std::endl;
    return -1;
  }
  return result;
}
