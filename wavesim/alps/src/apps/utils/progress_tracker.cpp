//
// Created by xuanx004 on 3/31/23.
//

#include "progress_tracker.h"

#include <fmt/chrono.h>
#include <fmt/format.h>

#include <chrono>

namespace alps::apps {

namespace {
/// @brief Format a duration to a string, with days included if the duration is
/// longer than 1 day
template<typename Rep, typename Period>
std::string
duration_to_string(std::chrono::duration<Rep, Period> const& duration_)
{
  using namespace std::chrono;
  auto const days = duration_cast<duration<int, std::ratio<86400>>>(duration_);
  if (days.count() > 0) {
    return fmt::format("{:d}d {:.1%T}", days.count(), duration_ - days);
  }
  return fmt::format("{:.1%T}", duration_);
}
} // namespace

ProgressTracker::ProgressTracker(int    total,
                                 double start_solver_time,
                                 double max_solver_time)
  : total_{total}
  , max_solver_time_{max_solver_time}
  , current_solver_time_{start_solver_time}
  , start_time_{std::chrono::high_resolution_clock::now()}
  , current_time_{std::chrono::high_resolution_clock::now()}
  , last_output_time_{std::chrono::high_resolution_clock::now()}
{
  times_.reserve(total_);
}

ProgressTracker& ProgressTracker::set_smoothing_factor(double smoothing_factor)
{
  ema_alpha_ = smoothing_factor;
  return *this;
}

ProgressTracker& ProgressTracker::set_output_step_interval(int n_steps)
{
  output_step_interval_ = n_steps;
  return *this;
}

ProgressTracker& ProgressTracker::set_output_time_interval(int seconds)
{
  output_time_interval_ = std::chrono::seconds(seconds);
  return *this;
}

bool ProgressTracker::is_complete() const
{
  return (current_ >= total_) || (current_solver_time_ >= max_solver_time_);
}

void ProgressTracker::update(double current_solver_time)
{
  current_solver_time_ = current_solver_time;
  ++current_;
  ++steps_since_last_output_;

  auto const last_time = current_time_;
  current_time_        = std::chrono::high_resolution_clock::now();
  decltype(time_per_step_ema_) const dt = current_time_ - last_time;

  // perform exponential moving average on the time per step
  if (times_.empty()) {
    time_per_step_ema_ = dt;
  } else {
    time_per_step_ema_ =
      ema_alpha_ * dt + (1 - ema_alpha_) * time_per_step_ema_;
  }
  times_.emplace_back(dt.count());
}

std::optional<std::string> ProgressTracker::output()
{
  if (!should_output()) {
    return std::nullopt;
  }

  steps_since_last_output_ = 0;
  last_output_time_        = std::chrono::high_resolution_clock::now();
  return fmt::format(
    "{}/{} [{:.0f}%] in {:} ({:.2f}ms/step, eta: {:})",
    current_,
    total_,
    static_cast<double>(current_) / total_ * 100,
    duration_to_string(elapsed()),
    std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
      time_per_step_ema_)
      .count(),
    duration_to_string(ETA()));
}

std::string ProgressTracker::complete_output()
{
  return fmt::format(
    "Done. {} steps in {:} ({:.2f}ms/step)",
    current_,
    duration_to_string(elapsed()),
    std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
      elapsed() / static_cast<double>(current_))
      .count());
}

std::chrono::duration<double> ProgressTracker::elapsed() const
{
  return current_time_ - start_time_;
}

std::chrono::duration<double> ProgressTracker::ETA() const
{
  return static_cast<double>(total_ - current_) * time_per_step_ema_;
}

bool ProgressTracker::should_output() const
{
  if (steps_since_last_output_ >= output_step_interval_) {
    return true;
  }
  if (auto const current_time = std::chrono::high_resolution_clock::now();
      current_time - last_output_time_ >= output_time_interval_) {
    return true;
  }
  return false;
}

} // namespace alps::apps
