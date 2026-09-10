//
// Created by xuanx004 on 3/31/23.
//

#pragma once

#include <chrono>
#include <optional>
#include <string>
#include <vector>

namespace alps::apps {

///@brief A class for tracking and reporting progress of a long running process.
/** This class provides functionality for tracking the progress of a long
 * running process and reporting the progress to the user. The progress can be
 * reported in terms of completed steps, time elapsed, and ETA. The
 * time-per-step is smoothed to avoid fluctuations in ETA.
 * To use this class, first create an instance of ProgressTracker by passing the
 * total number of steps in the process to the constructor. Then, optionally set
 * the smoothing factor, output step interval, and/or output time interval. To
 * update the progress tracker, call the update() method each time a step in the
 * process is completed.
 */
class ProgressTracker
{
 public:
  explicit ProgressTracker(
    int    total,
    double start_solver_time = 0,
    double max_solver_time   = std::numeric_limits<double>::max());

  /// @brief Set the smoothing factor for exponential moving average
  ProgressTracker& set_smoothing_factor(double smoothing_factor);

  /// @brief Set the output step interval
  ProgressTracker& set_output_step_interval(int n_steps);

  /// @brief Set the output time interval
  ProgressTracker& set_output_time_interval(int seconds);

  /// @brief Get the current step (1...total)
  int current() const { return current_ + 1; }

  /// @brief Get the total number of steps
  int total() const { return total_; }

  /// @brief Query whether the process is complete
  bool is_complete() const;

  /// @brief Indicate one step completed in the progress
  void update(double current_solver_time = 0);

  /// @brief Report the progress
  std::optional<std::string> output();

  /// @brief Report the completion of the process
  std::string complete_output();

  std::chrono::duration<double> elapsed() const;

  std::chrono::duration<double> ETA() const;

 private:
  /// @brief Query whether the output condition is met (i.e., reached a certain
  /// number of steps or a certain amount of time since last output)
  bool should_output() const;

  int    total_{0};
  double max_solver_time_{};
  double ema_alpha_{0.3}; /// smoothing factor for exponential moving average
  int    output_step_interval_{100};
  std::chrono::seconds output_time_interval_{std::chrono::seconds(60)};

  int    current_{0};
  double current_solver_time_{};

  std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
  std::chrono::time_point<std::chrono::high_resolution_clock> current_time_;
  std::chrono::duration<double, std::milli> time_per_step_ema_{
    0}; /// averaged time per step

  int steps_since_last_output_{0};
  std::chrono::time_point<std::chrono::high_resolution_clock> last_output_time_;

  std::vector<double> times_;
};

} // namespace alps::apps
