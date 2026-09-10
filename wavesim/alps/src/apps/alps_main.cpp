//
// Created by xuanx004 on 3/2/24.
//

#include <common/runtime/manager.h>

#include <iostream>

int alps_main(int argc, char** argv); // NOLINT(misc-use-internal-linkage)

int main(int argc, char** argv)
{
  try {
    return alps::RuntimeManager::execute_main(alps_main, argc, argv);
  } catch (std::exception const& e) {
    std::cerr << "Unhandled exception in runtime: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception caught in runtime" << std::endl;
    return 1;
  }
}
