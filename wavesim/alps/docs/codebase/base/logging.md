# Logging

Logging allows one to print messages with different logging levels, which can be used to control the verbosity of the output. This could be useful for debugging. One may specify the logging level by setting the environment variable `SPDLOG_LEVEL` to `debug` or `trace` (default is `info`) to obtain more information about the execution.

The logging is implemented using [`spdlog`](https://github.com/gabime/spdlog). The default logger is modified to include the MPI rank in the message.

```{note}
It is strongly recommended to use `spdlog::info`, `spdlog::warn`, `spdlog::error`, `spdlog::debug` and `spdlog::trace` instead of `std::cout` to print messages.
```

```{note}
`spdlog` comes with [`fmt`](https://fmt.dev/) library, which is used to format the messages. The format string is similar to Python's `str.format` method. See [here](https://fmt.dev/latest/syntax.html) for the syntax.

It is strongly recommended to use `fmt::format` instead of `std::sprintf` or `std::stringstream` to format strings.
```

```cpp
#include <common/base/logging.h>

default_logger()->info("Hello, world!");
default_logger()->debug("Hello, world! {}", 42); // prints when SPDLOG_LEVEL is set to debug
spdlog::debug("Hello, world! {}", 42); // equivalent to the above
```