#pragma once

#include <cstdint>
#include <memory>
#include <sstream>

#include <mpi.h>

#include <catch2/reporters/catch_reporter_streaming_base.hpp>

namespace Catch {
class IStream;
class ColourImpl;
class ConsoleReporter;
} // namespace Catch

/**
 * MPI-aware console reporter.
 *
 * This reporter behaves the same as a console reporter, except:
 * - only one process at a time will print assertion failures
 * - add the rank to the assertion output
 * - only the root process will print the final result, which is gathered from
 *   all processes
 */
class MpiReporter final : public Catch::StreamingReporterBase
{
  using StringRef      = Catch::StringRef;
  using AssertionStats = Catch::AssertionStats;
  using SectionInfo    = Catch::SectionInfo;
  using TestCaseInfo   = Catch::TestCaseInfo;
  using TestRunInfo    = Catch::TestRunInfo;
  using SectionStats   = Catch::SectionStats;
  using TestCaseStats  = Catch::TestCaseStats;
  using TestRunStats   = Catch::TestRunStats;

 public:
  MpiReporter(Catch::ReporterConfig&& config);

  ~MpiReporter() override;

  static std::string getDescription();

  void noMatchingTestCases(StringRef unmatchedSpec) override;
  void reportInvalidTestSpec(StringRef arg) override;

  void assertionEnded(AssertionStats const& stats) override;

  void sectionStarting(SectionInfo const& sectionInfo) override;
  void sectionEnded(SectionStats const& stats) override;

  void testCaseEnded(TestCaseStats const& stats) override;
  void testCaseStarting(TestCaseInfo const& testInfo) override;
  void testRunEnded(TestRunStats const& stats) override;
  void testRunStarting(TestRunInfo const& testRunInfo) override;

  void listTests(std::vector<Catch::TestCaseHandle> const& tests) override;
  void listTags(std::vector<Catch::TagInfo> const& tags) override;

 private:
  MpiReporter(Catch::IConfig const*                     full_config,
              Catch::ColourMode                         colourMode,
              std::map<std::string, std::string>        customOptions,
              Catch::Detail::unique_ptr<Catch::IStream> _stream);

  static constexpr int ROOT = 0;

  void lockStream();
  void unlockStream();
  void outputWithRankPrefix(std::string const& str);
  void handleUniformTestCountsSummary(TestRunStats const& stats);
  void handleNonUniformTestCountsSummary(TestRunStats const& stats);
  bool is_root() const;

  // RMA window for ticket-based lock
  MPI_Win lock_win_{MPI_WIN_NULL};

  // Ticket lock structure (only allocated on ROOT)
  // Using uint64_t to prevent overflow
  struct
  {
    uint64_t next_ticket;
    uint64_t now_serving;
  } ticket_lock_{0, 0};

  MPI_Comm comm_;
  int      rank_;
  int      comm_size_;

  std::unique_ptr<Catch::ConsoleReporter> reporter_;
  std::ostringstream*                     r_stream_{nullptr};
};
