#include "mpi_reporter.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <ostream>
#include <sstream>
#include <string_view>
#include <thread>

#include <catch2/catch_test_case_info.hpp>
#include <catch2/interfaces/catch_interfaces_config.hpp>
#include <catch2/internal/catch_console_colour.hpp>
#include <catch2/internal/catch_context.hpp>
#include <catch2/internal/catch_istream.hpp>
#include <catch2/internal/catch_random_number_generator.hpp>
#include <catch2/internal/catch_test_registry.hpp>
#include <catch2/reporters/catch_reporter_console.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/reporters/catch_reporter_streaming_base.hpp>
#include <fmt/format.h>
#include <mpi.h>

namespace {
class StringStream : public Catch::IStream
{
  std::ostringstream m_oss;
  bool               is_console{false};

 public:
  explicit StringStream(bool _is_console)
    : is_console(_is_console)
  {}

  std::ostream& stream() override { return m_oss; }

  bool isConsole() const override { return is_console; }
};

void clear_stream(std::ostringstream& oss)
{
  oss.str("");
  oss.clear();
}

bool is_same_on_all_ranks(uint64_t value, MPI_Comm comm)
{
  uint64_t global_min{};
  uint64_t global_max{};

  MPI_Allreduce(&value, &global_min, 1, MPI_UINT64_T, MPI_MIN, comm);
  MPI_Allreduce(&value, &global_max, 1, MPI_UINT64_T, MPI_MAX, comm);

  return (global_min == global_max);
}

class SummaryColumn
{
 public:
  SummaryColumn(std::string suffix, Catch::Colour::Code colour)
    : m_suffix(CATCH_MOVE(suffix))
    , m_colour(colour)
  {}

  SummaryColumn&& addRow(std::uint64_t count) &&
  {
    std::string row = std::to_string(count);
    m_width         = std::max(m_width, row.size());
    m_rows.push_back(row);
    return std::move(*this);
  }

  std::string const&  getSuffix() const { return m_suffix; }
  Catch::Colour::Code getColour() const { return m_colour; }
  std::string         getRow(std::size_t index) const
  {
    return fmt::format("{:>{}}", m_rows[index], m_width);
  }

 private:
  std::string              m_suffix;
  Catch::Colour::Code      m_colour;
  std::size_t              m_width = 0;
  std::vector<std::string> m_rows;
};

void printSummaryRow(std::ostream&                     stream,
                     Catch::ColourImpl&                colour,
                     std::string_view                  label,
                     std::vector<SummaryColumn> const& cols,
                     std::size_t                       row)
{
  for (auto const& col : cols) {
    auto const& value  = col.getRow(row);
    auto const& suffix = col.getSuffix();
    if (suffix.empty()) {
      stream << label << ": ";
      if (value != "0") {
        stream << value;
      } else {
        stream << colour.guardColour(Catch::Colour::Warning) << "- none -";
      }
    } else if (value != "0") {
      stream << colour.guardColour(Catch::Colour::LightGrey) << " | "
             << colour.guardColour(col.getColour()) << value << ' ' << suffix;
    }
  }
  stream << '\n';
}

void printTestRunTotalsWithRank(std::ostream&        stream,
                                Catch::ColourImpl&   streamColour,
                                Catch::Totals const& totals,
                                int                  rank,
                                int                  comm_size)
{
  std::size_t rank_width = std::to_string(comm_size - 1).size();
  std::string rank_str =
    rank >= 0 ? fmt::format("[rk {:>{}}]", rank, rank_width) : "[rk all]";

  if (totals.testCases.total() == 0) {
    stream << streamColour.guardColour(Catch::Colour::Warning) << rank_str
           << " No tests ran\n";
    return;
  }

  if (totals.assertions.total() > 0 && totals.testCases.allPassed()) {
    auto plural = [](std::uint64_t count, std::string_view singular) {
      return fmt::format("{} {}{}", count, singular, (count == 1 ? "" : "s"));
    };
    stream << streamColour.guardColour(Catch::Colour::ResultSuccess) << rank_str
           << " All tests passed";
    if (rank >= 0) {
      // Per-rank output
      stream << fmt::format(" ({}; {})\n",
                            plural(totals.testCases.passed, "test case"),
                            plural(totals.assertions.passed, "assertion"));
    } else {
      // Aggregate output
      stream << fmt::format(" ({} per rank; {} across {})\n",
                            plural(totals.testCases.passed, "test case"),
                            plural(totals.assertions.passed, "assertion"),
                            plural(comm_size, "rank"));
    }
    return;
  }

  std::vector<SummaryColumn> columns;
  // Don't include "skipped assertions" in total count
  const auto totalAssertionCount =
    totals.assertions.total() - totals.assertions.skipped;
  columns.push_back(SummaryColumn("", Catch::Colour::None)
                      .addRow(totals.testCases.total())
                      .addRow(totalAssertionCount));
  columns.push_back(SummaryColumn("passed", Catch::Colour::Success)
                      .addRow(totals.testCases.passed)
                      .addRow(totals.assertions.passed));
  columns.push_back(SummaryColumn("failed", Catch::Colour::ResultError)
                      .addRow(totals.testCases.failed)
                      .addRow(totals.assertions.failed));
  columns.push_back(SummaryColumn("skipped", Catch::Colour::Skip)
                      .addRow(totals.testCases.skipped)
                      // Don't print "skipped assertions"
                      .addRow(0));
  columns.push_back(
    SummaryColumn("failed as expected", Catch::Colour::ResultExpectedFailure)
      .addRow(totals.testCases.failedButOk)
      .addRow(totals.assertions.failedButOk));

  printSummaryRow(stream, streamColour, rank_str + " test cases", columns, 0);
  printSummaryRow(stream, streamColour, rank_str + " assertions", columns, 1);
}
} // anonymous namespace

MpiReporter::MpiReporter(Catch::ReporterConfig&& config)
  : MpiReporter(config.fullConfig(),
                config.colourMode(),
                config.customOptions(),
                CATCH_MOVE(config).takeStream() /* .takeStream() */)
{}

MpiReporter::MpiReporter(Catch::IConfig const*              full_config,
                         Catch::ColourMode                  color_mode,
                         std::map<std::string, std::string> customOptions,
                         Catch::Detail::unique_ptr<Catch::IStream> _stream)
  : Catch::StreamingReporterBase(Catch::ReporterConfig{full_config,
                                                       CATCH_MOVE(_stream),
                                                       color_mode,
                                                       customOptions})
  , comm_([] {
    MPI_Comm new_comm{MPI_COMM_NULL};
    MPI_Comm_dup(MPI_COMM_WORLD, &new_comm);
    return new_comm;
  }())
  , rank_([] {
    int rank{};
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
  }())
  , comm_size_([] {
    int size{};
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
  }())
{
  auto string_stream = Catch::Detail::make_unique<StringStream>(
    this->m_wrapped_stream->isConsole());
  // r_stream_ owned by reporter_ via StringStream passed to ConsoleReporter.
  // Lifetime: valid as long as reporter_ exists and its stream is unchanged.
  r_stream_ = dynamic_cast<std::ostringstream*>(&string_stream->stream());
  assert(r_stream_ != nullptr
         && "Failed to get ostringstream from StringStream");
  reporter_ = std::make_unique<Catch::ConsoleReporter>(Catch::ReporterConfig{
    full_config, CATCH_MOVE(string_stream), color_mode, customOptions});

  // Create RMA window for ticket-based lock
  if (comm_size_ > 1) {
    if (is_root()) {
      // ROOT rank: Create window with ticket lock structure
      // Using uint64_t (8 bytes each) for overflow protection
      MPI_Win_create(&ticket_lock_,
                     sizeof(ticket_lock_), // 16 bytes
                     sizeof(uint64_t),     // displacement unit = 8 bytes
                     MPI_INFO_NULL,
                     comm_,
                     &lock_win_);
    } else {
      // Non-root ranks: Create window with zero size
      MPI_Win_create(nullptr,
                     0,                // size = 0
                     sizeof(uint64_t), // displacement unit (must match ROOT)
                     MPI_INFO_NULL,
                     comm_,
                     &lock_win_);
    }

    // Ensure window is ready across all ranks
    MPI_Barrier(comm_);
  }
}

MpiReporter::~MpiReporter()
{
  // Free RMA window if it was created
  if (lock_win_ != MPI_WIN_NULL) {
    MPI_Win_free(&lock_win_);
  }

  // Free the duplicated communicator
  MPI_Comm_free(&comm_);
}

void MpiReporter::lockStream()
{
  if (comm_size_ == 1) return;

  // Step 1: Atomically get a ticket number
  uint64_t my_ticket{};
  uint64_t one = 1;

  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ROOT, 0, lock_win_);
  MPI_Fetch_and_op(&one,
                   &my_ticket,
                   MPI_UINT64_T,
                   ROOT,
                   0, // offset 0: next_ticket
                   MPI_SUM,
                   lock_win_);
  MPI_Win_unlock(ROOT, lock_win_);

  // Step 2: Wait for our turn with adaptive backoff
  // Initialize to ensure loop entry
  constexpr int64_t max_backoff_us = 10000; // 10 millisecond
  int64_t           backoff_us     = 1;

  for (uint64_t current_serving = my_ticket + 1; current_serving != my_ticket;
       backoff_us               = std::min(backoff_us * 2, max_backoff_us)) {
    MPI_Win_lock(MPI_LOCK_SHARED, ROOT, 0, lock_win_);
    MPI_Get(&current_serving,
            1,
            MPI_UINT64_T,
            ROOT,
            1, // offset 1: now_serving
            1,
            MPI_UINT64_T,
            lock_win_);
    MPI_Win_unlock(ROOT, lock_win_);

    if (current_serving != my_ticket) {
      // Intentionally invoke an MPI call to progress MPI communication
      int flag{};
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &flag, MPI_STATUS_IGNORE);
      if (backoff_us < 100) {
        std::this_thread::yield();
      } else {
        std::this_thread::sleep_for(std::chrono::microseconds(backoff_us));
      }
    }
  }
}

void MpiReporter::unlockStream()
{
  if (comm_size_ == 1) return;

  // Step 1: Flush output
  m_stream << std::flush;

  // Step 2: Increment now_serving to release lock
  uint64_t one = 1;
  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ROOT, 0, lock_win_);
  MPI_Accumulate(&one,
                 1,
                 MPI_UINT64_T,
                 ROOT,
                 1, // offset 1: now_serving
                 1,
                 MPI_UINT64_T,
                 MPI_SUM,
                 lock_win_);
  MPI_Win_unlock(ROOT, lock_win_);
}

void MpiReporter::outputWithRankPrefix(const std::string& str)
{
  if (str.empty()) {
    return;
  }

  std::istringstream iss(str);
  std::string        line;

  while (std::getline(iss, line)) {
    m_stream << m_colour->guardColour(Catch::Colour::None)
             << fmt::format("[rk{}] {}", rank_, line);
    if (!iss.eof()) {
      m_stream << '\n';
    }
  }
}

std::string MpiReporter::getDescription()
{
  return "MPI aware console reporter";
}

void MpiReporter::noMatchingTestCases(StringRef unmatchedSpec)
{
  clear_stream(*r_stream_);
  reporter_->noMatchingTestCases(unmatchedSpec);
  if (is_root()) {
    lockStream();
    m_stream << r_stream_->str();
    unlockStream();
  }
  clear_stream(*r_stream_);
}

void MpiReporter::reportInvalidTestSpec(StringRef arg)
{
  clear_stream(*r_stream_);
  reporter_->reportInvalidTestSpec(arg);
  if (is_root()) {
    lockStream();
    m_stream << r_stream_->str();
    unlockStream();
  }
  clear_stream(*r_stream_);
}

void MpiReporter::assertionEnded(AssertionStats const& stats)
{
  clear_stream(*r_stream_);
  reporter_->assertionEnded(stats);
  if (r_stream_->str().empty()) {
    return;
  }
  // Acquire lock for serialized output
  lockStream();

  // Write assertion output
  outputWithRankPrefix(r_stream_->str());

  // Release lock
  unlockStream();
  clear_stream(*r_stream_);
}

void MpiReporter::sectionStarting(SectionInfo const& sectionInfo)
{
  clear_stream(*r_stream_);
  reporter_->sectionStarting(sectionInfo);
  if (r_stream_->str().empty()) {
    return;
  }

  // Acquire lock for serialized output
  lockStream();

  outputWithRankPrefix(r_stream_->str());

  // Release lock
  unlockStream();
  clear_stream(*r_stream_);
}

void MpiReporter::sectionEnded(SectionStats const& stats)
{
  clear_stream(*r_stream_);
  reporter_->sectionEnded(stats);

  if (r_stream_->str().empty()) {
    return;
  }

  // Acquire lock for serialized output
  lockStream();

  outputWithRankPrefix(r_stream_->str());

  // Release lock
  unlockStream();
  clear_stream(*r_stream_);
}

void MpiReporter::testCaseStarting(TestCaseInfo const& testInfo)
{
  reporter_->testCaseStarting(testInfo);
}

void MpiReporter::testCaseEnded(TestCaseStats const& stats)
{
  clear_stream(*r_stream_); // Clear previous contents
  reporter_->testCaseEnded(stats);

  if (r_stream_->str().empty()) {
    return;
  }

  // Acquire lock for serialized output
  lockStream();

  outputWithRankPrefix(r_stream_->str());

  // Release lock
  unlockStream();
  clear_stream(*r_stream_);
}

void MpiReporter::testRunEnded(TestRunStats const& stats)
{
  clear_stream(*r_stream_);

  static_assert(
    std::is_same_v<decltype(stats.totals.assertions.passed), uint64_t>);
  bool uniform_test_counts =
    is_same_on_all_ranks(stats.totals.testCases.passed, comm_)
    && is_same_on_all_ranks(stats.totals.testCases.failed, comm_)
    && is_same_on_all_ranks(stats.totals.testCases.failedButOk, comm_)
    && is_same_on_all_ranks(stats.totals.testCases.skipped, comm_);

  if (uniform_test_counts) {
    handleUniformTestCountsSummary(stats);
  } else {
    handleNonUniformTestCountsSummary(stats);
  }
  // still call the underlying report to ensure proper finalization
  reporter_->testRunEnded(stats);
  Catch::StreamingReporterBase::testRunEnded(stats);

  clear_stream(*r_stream_);
}

void MpiReporter::handleUniformTestCountsSummary(TestRunStats const& stats)
{
  auto globalStats = stats;
  MPI_Reduce(&stats.aborting,
             &globalStats.aborting,
             1,
             MPI_CXX_BOOL,
             MPI_LOR,
             ROOT,
             comm_);
  auto reduce_count = [comm = this->comm_, root = ROOT](uint64_t local) {
    uint64_t global{};
    MPI_Reduce(&local, &global, 1, MPI_UINT64_T, MPI_SUM, root, comm);
    return global;
  };

  auto const& local_assertions  = stats.totals.assertions;
  auto&       global_assertions = globalStats.totals.assertions;
  global_assertions.passed      = reduce_count(local_assertions.passed);
  global_assertions.failed      = reduce_count(local_assertions.failed);
  global_assertions.failedButOk = reduce_count(local_assertions.failedButOk);
  global_assertions.skipped     = reduce_count(local_assertions.skipped);

  if (!is_root()) return;

  lockStream();
  m_stream << std::string(CATCH_CONFIG_CONSOLE_WIDTH - 1, '=') << '\n';
  printTestRunTotalsWithRank(
    m_stream, *m_colour, globalStats.totals, -1, comm_size_);
  unlockStream();
}

void MpiReporter::handleNonUniformTestCountsSummary(TestRunStats const& stats)
{
  constexpr int                        kCountsPerRank = 8;
  std::array<uint64_t, kCountsPerRank> local_counts   = {
    stats.totals.testCases.passed,
    stats.totals.testCases.failed,
    stats.totals.testCases.failedButOk,
    stats.totals.testCases.skipped,
    stats.totals.assertions.passed,
    stats.totals.assertions.failed,
    stats.totals.assertions.failedButOk,
    stats.totals.assertions.skipped};
  std::vector<uint64_t> all_counts((size_t)comm_size_ * kCountsPerRank);
  MPI_Gather(local_counts.data(),
             kCountsPerRank,
             MPI_UINT64_T,
             all_counts.data(),
             kCountsPerRank,
             MPI_UINT64_T,
             ROOT,
             comm_);
  if (!is_root()) return;

  // Determine the width needed for rank number
  lockStream();
  m_stream << std::string(CATCH_CONFIG_CONSOLE_WIDTH - 1, '=') << '\n';
  for (int rank = 0; rank < comm_size_; ++rank) {
    const std::size_t base = (size_t)rank * kCountsPerRank;
    Catch::Totals     totals{};
    totals.testCases.passed       = all_counts[base + 0];
    totals.testCases.failed       = all_counts[base + 1];
    totals.testCases.failedButOk  = all_counts[base + 2];
    totals.testCases.skipped      = all_counts[base + 3];
    totals.assertions.passed      = all_counts[base + 4];
    totals.assertions.failed      = all_counts[base + 5];
    totals.assertions.failedButOk = all_counts[base + 6];
    totals.assertions.skipped     = all_counts[base + 7];

    printTestRunTotalsWithRank(m_stream, *m_colour, totals, rank, comm_size_);
  }
  unlockStream();
}

void MpiReporter::testRunStarting(TestRunInfo const& testRunInfo)
{
  clear_stream(*r_stream_);
  Catch::StreamingReporterBase::testRunStarting(testRunInfo);
  reporter_->testRunStarting(testRunInfo);

  if (is_root()) {
    lockStream();
    m_stream << r_stream_->str();
    unlockStream();
  }
  clear_stream(*r_stream_);
}

void MpiReporter::listTests(std::vector<Catch::TestCaseHandle> const& tests)
{
  if (!is_root()) return;

  lockStream();
  Catch::StreamingReporterBase::listTests(tests);
  unlockStream();
}

void MpiReporter::listTags(std::vector<Catch::TagInfo> const& tags)
{
  if (!is_root()) return;

  lockStream();
  Catch::StreamingReporterBase::listTags(tags);
  unlockStream();
}

bool MpiReporter::is_root() const
{
  return rank_ == ROOT;
}
