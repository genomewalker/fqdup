#ifndef FQDUP_LOGGER_HPP
#define FQDUP_LOGGER_HPP

#include <atomic>
#include <string>
#include <vector>

// Simple thread-safe logger + periodic progress reporter used by the CLI.
// Usage:
//   init_logger(optional_logfile);
//   start_progress_report(counters, interval_seconds);
//   // run work that updates counters (atomic<size_t>)
//   stop_progress_report();

void init_logger(const std::string& logfile = "");
void shutdown_logger();

// counters: vector of (name, pointer-to-atomic)
// counters: vector of (name, pointer-to-atomic)
// optional stage_names: if provided and a counter named "stage" is present, the numeric
// stage id will be mapped to these names for human-readable output.
void start_progress_report(const std::vector<std::pair<std::string, const std::atomic<size_t>*>>& counters, int interval_seconds = 5, const std::vector<std::string>& stage_names = {});
void stop_progress_report();

// immediate logging
void log_info(const std::string& msg);
void log_warn(const std::string& msg);
// Progress-style update: emit a single-line, carriage-return-updating message to
// stderr (no newline) for high-frequency progress updates. This intentionally
// does not write to the logfile to avoid huge files when reporting millions of
// updates. If you need the progress persisted, call log_info() occasionally.
void log_progress(const std::string& msg);

// Notify the reporter that progress counters changed and it should redraw
// the live status line. This is a lightweight wake-up; callers should use it
// at low frequency (or in a batched/debounced manner) to avoid excessive
// wakeups. It's safe to call from multiple threads.
void notify_progress_report();

#endif // FQDUP_LOGGER_HPP
