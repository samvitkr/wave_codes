#include "mpi_reporter.h"
#include <common/runtime/manager.h>

#include <catch2/catch_session.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <mpi.h>

class runtimeSetupTeardownListener : public Catch::EventListenerBase
{
 public:
  using Catch::EventListenerBase::EventListenerBase;

  void testRunStarting(Catch::TestRunInfo const& /*unused*/) override
  {
    alps::RuntimeManager::instance().init_runtimes();
  }
};
CATCH_REGISTER_LISTENER(runtimeSetupTeardownListener);

CATCH_REGISTER_REPORTER("mpi", MpiReporter)

int main(int argc, char* argv[])
{
  int provided_thread_support = 0;
  int mpi_init_result         = MPI_Init_thread(
    &argc, &argv, MPI_THREAD_MULTIPLE, &provided_thread_support);
  if (mpi_init_result != MPI_SUCCESS) return mpi_init_result;
  if (provided_thread_support < MPI_THREAD_MULTIPLE) {
    MPI_Finalize();
    return 1;
  }

  Catch::Session session;

  // Broadcast default RNG seed to all processes
  // This ensures that even if the tests are run with randomized order, all
  // processes run the same tests
  auto seed = Catch::generateRandomSeed(Catch::GenerateFrom::Default);
  static_assert(std::is_same_v<decltype(seed), uint32_t>);
  MPI_Bcast(&seed, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
  session.configData().rngSeed = seed;
  session.configData().reporterSpecifications.clear();
  session.configData().reporterSpecifications.push_back(
    std::move(*Catch::parseReporterSpec("mpi")));

  int returnCode = session.applyCommandLine(argc, argv);
  if (returnCode != 0) { // command line error
    MPI_Finalize();
    return returnCode;
  }

  int test_result = session.run();

  alps::RuntimeManager::instance().finalize(); // also finalize MPI
  return test_result;
}
