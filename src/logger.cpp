#include "fqdup/logger.hpp"

#include <chrono>
#include <cmath>
#include <condition_variable>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <iomanip>

// Small helpers for pretty-printing numbers and rates
static std::string human_int(size_t v) {
    // add thousands separators
    std::string s = std::to_string(v);
    int insertPosition = (int)s.length() - 3;
    while (insertPosition > 0) {
        s.insert(insertPosition, ",");
        insertPosition -= 3;
    }
    return s;
}

static std::string human_rate(double r) {
    // format r as SI with 1 decimal when >=1k
    const char* units[] = {"","K","M","G","T"};
    int ui = 0;
    double val = r;
    while (val >= 1000.0 && ui < 4) { val /= 1000.0; ++ui; }
    std::ostringstream ss; ss << std::fixed;
    if (ui == 0) ss << std::setprecision(0) << val << units[ui];
    else ss << std::setprecision(2) << val << units[ui];
    return ss.str();
}

static std::mutex g_log_mutex;
static std::ofstream g_log_file;
// Optional pointer to a stage counter (set by start_progress_report if provided)
static const std::atomic<size_t>* g_stage_counter = nullptr;
// Reporter counters and stage names (declared early so log_impl can use stage names)
static std::vector<std::pair<std::string, const std::atomic<size_t>*>> g_counters;
static std::vector<std::string> g_stage_names;
// Flag used by the reporter to indicate it is running (declared early so log_impl can check)
static std::atomic<bool> g_report_running{false};

// Reporter thread control primitives (declared early so log_impl can notify)
static std::thread g_report_thread;
static std::condition_variable g_report_cv;
static std::mutex g_report_mutex;

void init_logger(const std::string& logfile) {
    std::lock_guard<std::mutex> lk(g_log_mutex);
    if (!logfile.empty()) {
        g_log_file.open(logfile, std::ios::app);
    }
}

void shutdown_logger() {
    std::lock_guard<std::mutex> lk(g_log_mutex);
    if (g_log_file.is_open()) g_log_file.close();
}

void log_impl(const std::string& level, const std::string& msg) {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    // Include a stage label when available to make logs stage-specific
    std::string stage_label;
    if (g_stage_counter) {
        size_t s = g_stage_counter->load();
        if (!g_stage_names.empty() && s < g_stage_names.size()) stage_label = g_stage_names[s];
        else stage_label = std::to_string(s);
    }
    // Sanitize message to a single line: replace any internal newlines/carriage
    // returns with spaces so each log entry is exactly one line on disk.
    std::string safe_msg; safe_msg.reserve(msg.size());
    for (char c : msg) {
        if (c == '\n' || c == '\r') safe_msg.push_back(' ');
        else safe_msg.push_back(c);
    }

    std::ostringstream ss;
    ss << "[" << buf << "] " << level;
    if (!stage_label.empty()) ss << " [" << stage_label << "]";
    ss << ": " << safe_msg << '\n';
    std::string out = ss.str();
    {
        std::lock_guard<std::mutex> lk(g_log_mutex);
        // If a live progress reporter is running on the TTY, it may have emitted a
        // carriage-return-updated line without a trailing newline. To avoid the
        // timestamped log line clobbering that live line we first emit a newline
        // so the normal log line always starts on a fresh line in interactive mode.
        bool interactive = isatty(fileno(stderr));
        if (interactive && g_report_running.load()) std::cerr << std::endl;
        std::cerr << out;
        // If interactive and reporter running, re-draw the progress line after
        // emitting the timestamped log line so the user still sees the live status
        // at the bottom of the terminal.
        if (interactive && g_report_running.load() && !g_counters.empty()) {
            // invoke the reporter build to print a fresh status (non-flushing)
            // We'll reuse the reporter_loop machinery by notifying the report thread
            g_report_cv.notify_all();
        }
        if (g_log_file.is_open()) g_log_file << out;
    }
}

void log_info(const std::string& msg) { log_impl("INFO", msg); }
void log_warn(const std::string& msg) { log_impl("WARN", msg); }

void log_progress(const std::string& msg) {
    bool interactive = isatty(fileno(stderr));
    // sanitize like other messages
    std::string safe_msg; safe_msg.reserve(msg.size());
    for (char c : msg) {
        if (c == '\n' || c == '\r') safe_msg.push_back(' ');
        else safe_msg.push_back(c);
    }
    // Build a short prefix like the standard logger (timestamp + level)
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    std::string line = std::string("[") + buf + "] INFO: " + safe_msg;
    if (interactive) {
        std::lock_guard<std::mutex> lk(g_log_mutex);
        std::cerr << '\r' << line << std::flush;
    } else {
        // Non-interactive: fall back to a normal log line so output is captured
        log_impl("INFO", safe_msg);
    }
}

// Progress reporter
void reporter_loop(int interval_seconds) {
    bool interactive = isatty(fileno(stderr));
    size_t prev_len = 0;
    // store previous counter value and timestamp to compute accurate rates
    std::unordered_map<std::string, std::pair<size_t, std::chrono::steady_clock::time_point>> prev_vals;
    // additional windowed samples for multi-second averaging (e.g., 5s)
    const std::chrono::seconds window_sec(5);
    std::unordered_map<std::string, std::pair<size_t, std::chrono::steady_clock::time_point>> prev_vals_window;
    // EMA state for sustained/rolling-rate smoothing
    std::unordered_map<std::string, double> prev_rate_ema;
    std::unordered_map<std::string, std::chrono::steady_clock::time_point> prev_rate_ema_tp;
    // Print an initial status immediately so users see something during init
    // Returns the computed (rec_rate, ids_rate) so the loop can adapt the interval.
    auto build_and_emit = [&](bool flush_file) -> std::pair<double,double> {
        std::ostringstream ss;
        size_t stage_val = SIZE_MAX;
        // short key mapping for a compact, professional line
        auto short_key = [](const std::string &k)->std::string {
            if (k == "produced_chunks") return "prod";
            if (k == "processed_chunks") return "proc";
            if (k == "processed_records") return "rec";
            if (k == "parsed_ids") return "ids";
            if (k == "winners_found") return "win";
            if (k == "queue_size") return "q";
            if (k == "stage") return "stage";
            return k;
        };
        // Build a compact metrics section using the counters provided in g_counters
        std::vector<std::pair<std::string,size_t>> ordered_vals;
        ordered_vals.reserve(g_counters.size());
        for (const auto &p : g_counters) {
            size_t v = p.second ? p.second->load() : 0;
            ordered_vals.emplace_back(p.first, v);
            if (p.first == "stage") stage_val = v;
        }
        // Stage label
        std::string stage_label = "";
        if (stage_val != SIZE_MAX) {
            if (!g_stage_names.empty() && stage_val < g_stage_names.size()) stage_label = g_stage_names[stage_val];
            else stage_label = std::to_string(stage_val);
        }
    // dynamic line: show stage label but omit a timestamp to keep the live
    // progress line compact and avoid producing duplicate timestamped lines.
    ss << "[" << stage_label << "] ";
        bool first = true;
        // Display rounding for high-volume counters to avoid rapid visual churn.
        // parsed_ids is updated very frequently; round its visible value to
        // the nearest million and append '+' when the true value exceeds
        // the rounded display. This keeps the live line stable while rates
        // remain accurate.
        static const size_t IDS_ROUND = 1000000u;
        for (const auto &kv : ordered_vals) {
            const std::string &k = kv.first;
            size_t v = kv.second;
            if (!first) ss << " "; first = false;
            ss << short_key(k) << "=";
            if (k == "parsed_ids" && v >= IDS_ROUND) {
                size_t rounded = (v / IDS_ROUND) * IDS_ROUND;
                ss << human_int(rounded);
                if (v > rounded) ss << "+";
            } else {
                // default: show full human-friendly integer with separators
                ss << human_int(v);
            }
        }
    // Build a lookup for rate calculations
    std::unordered_map<std::string, size_t> vals;
    for (const auto &kv : ordered_vals) vals[kv.first] = kv.second;

    // throughput: compute records/sec and ids/sec using a sliding multi-second
    // window (window_sec) so rates are averaged over a few seconds instead of
    // instantaneously. We still keep precise per-sample prev_vals but compute
    // reported rates from prev_vals_window when available.
    double rec_rate = 0.0;
    double ids_rate = 0.0;
    auto now_tp = std::chrono::steady_clock::now();
    auto itpr = vals.find("processed_records");
    if (itpr != vals.end()) {
        size_t curr = itpr->second;
        // capture small-sample previous value (if any) before we update
        size_t prev_small = 0;
        std::chrono::steady_clock::time_point prev_small_tp;
        auto itpv = prev_vals.find("processed_records");
        if (itpv != prev_vals.end()) {
            prev_small = itpv->second.first;
            prev_small_tp = itpv->second.second;
        }
        // update fine-grained sample
        prev_vals["processed_records"] = { curr, now_tp };

        // try to compute windowed rate
        auto itw = prev_vals_window.find("processed_records");
        if (itw != prev_vals_window.end()) {
            size_t prev = itw->second.first;
            auto prev_tp = itw->second.second;
            double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(now_tp - prev_tp).count();
            if (elapsed >= double(window_sec.count())) {
                if (curr >= prev) rec_rate = double(curr - prev) / elapsed;
                // advance the window sample
                itw->second = { curr, now_tp };
            } else {
                // fallback to short-interval rate if available
                if (prev_small_tp.time_since_epoch().count() != 0) {
                    double small_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(now_tp - prev_small_tp).count();
                    if (small_elapsed >= 1e-6 && curr >= prev_small) rec_rate = double(curr - prev_small) / small_elapsed;
                }
            }
        } else {
            // first time: seed the window sample
            prev_vals_window["processed_records"] = { curr, now_tp };
            // if we had a previous small sample, compute a short-term rate
            if (prev_small_tp.time_since_epoch().count() != 0) {
                double small_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(now_tp - prev_small_tp).count();
                if (small_elapsed >= 1e-6 && curr >= prev_small) rec_rate = double(curr - prev_small) / small_elapsed;
            }
        }
    }
    auto itpi = vals.find("parsed_ids");
    if (itpi != vals.end()) {
        size_t curr = itpi->second;
        // capture small-sample previous value before update
        size_t prev_small = 0;
        std::chrono::steady_clock::time_point prev_small_tp;
        auto itpv = prev_vals.find("parsed_ids");
        if (itpv != prev_vals.end()) {
            prev_small = itpv->second.first;
            prev_small_tp = itpv->second.second;
        }
        prev_vals["parsed_ids"] = { curr, now_tp };
        auto itw = prev_vals_window.find("parsed_ids");
        if (itw != prev_vals_window.end()) {
            size_t prev = itw->second.first;
            auto prev_tp = itw->second.second;
            double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(now_tp - prev_tp).count();
            if (elapsed >= double(window_sec.count())) {
                if (curr >= prev) ids_rate = double(curr - prev) / elapsed;
                itw->second = { curr, now_tp };
            } else {
                if (prev_small_tp.time_since_epoch().count() != 0) {
                    double small_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(now_tp - prev_small_tp).count();
                    if (small_elapsed >= 1e-6 && curr >= prev_small) ids_rate = double(curr - prev_small) / small_elapsed;
                }
            }
        } else {
            prev_vals_window["parsed_ids"] = { curr, now_tp };
            if (prev_small_tp.time_since_epoch().count() != 0) {
                double small_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(now_tp - prev_small_tp).count();
                if (small_elapsed >= 1e-6 && curr >= prev_small) ids_rate = double(curr - prev_small) / small_elapsed;
            }
        }
    }
    // At this point rec_rate and ids_rate hold the instantaneous/windowed rates
    // Compute an exponential moving average (EMA) for each rate so the
    // displayed value is a sustained rolling mean rather than a noisy instant.
    auto apply_ema = [&](const std::string &key, double inst_rate)->double {
        double tau = double(window_sec.count()); // seconds
        auto now = now_tp;
        double ema = inst_rate;
        auto it = prev_rate_ema.find(key);
        if (it != prev_rate_ema.end()) {
            double prev_ema = it->second;
            auto it_tp = prev_rate_ema_tp.find(key);
            std::chrono::duration<double> dt = now - it_tp->second;
            double dt_s = dt.count();
            if (dt_s <= 0.0) dt_s = 1e-6;
            double alpha = 1.0 - std::exp(-dt_s / tau);
            ema = prev_ema + alpha * (inst_rate - prev_ema);
        }
        prev_rate_ema[key] = ema;
        prev_rate_ema_tp[key] = now;
        return ema;
    };

    double inst_rec = rec_rate;
    double inst_ids = ids_rate;
    rec_rate = apply_ema("processed_records", inst_rec);
    ids_rate = apply_ema("parsed_ids", inst_ids);

    // Choose which rate to show depending on current stage to avoid showing
    // an irrelevant 0.0 r/s during indexing. If we're in the "index" stage
    // prefer showing ids/s; if in shard/merge prefer r/s; otherwise show both.
    // Reduce visual churn on high-volume rates by rounding the displayed
    // ids/s value to a magnitude-aware unit (e.g., 10k when rate >100k).
    auto format_coarse_rate = [&](double rate)->std::string {
        if (rate <= 0.0) return human_rate(0.0);
        double unit = 1.0;
        if (rate >= 100000.0) unit = 10000.0; // show in 10k steps
        else if (rate >= 10000.0) unit = 1000.0; // 1k steps
        else if (rate >= 1000.0) unit = 100.0; // 100 steps
        else unit = 1.0;
        double rounded = std::floor(rate / unit) * unit;
        std::string s = human_rate(rounded);
        if (rate > rounded) s += "+";
        return s;
    };

    if (!stage_label.empty() && stage_label == "index") {
        ss << " | " << format_coarse_rate(ids_rate) << " ids/s";
    } else if (!stage_label.empty() && (stage_label == "shard" || stage_label == "merge")) {
        ss << " | " << format_coarse_rate(rec_rate) << " r/s";
    } else {
        ss << " | " << format_coarse_rate(rec_rate) << " r/s";
        if (ids_rate > 0.0) ss << ", " << format_coarse_rate(ids_rate) << " ids/s";
    }
    std::string line = ss.str();
        if (interactive) {
            std::string out = "\r" + line;
            if (prev_len > line.size()) out += std::string(prev_len - line.size(), ' ');
            prev_len = line.size();
            std::lock_guard<std::mutex> lk(g_log_mutex);
            std::cerr << out << std::flush;
            if (flush_file && g_log_file.is_open()) g_log_file << line << std::endl;
        } else {
            if (flush_file) log_info(line);
        }
        return {rec_rate, ids_rate};
    };
    // Adaptive emission: wait for wakeups but only emit when we have useful
    // information to show. Emit when either the stage counter is non-zero or
    // any provided counter is non-zero (so we don't print an "init" line).
    int curr_interval = std::max(interval_seconds, (int)window_sec.count());
    std::pair<double,double> rates{0.0,0.0};

    while (g_report_running.load()) {
        {
            std::unique_lock<std::mutex> lk(g_report_mutex);
            g_report_cv.wait_for(lk, std::chrono::seconds(curr_interval));
        }
        if (!g_report_running.load()) break;
        bool has_metrics = false;
        for (const auto &p : g_counters) {
            if (p.second && p.second->load() > 0) { has_metrics = true; break; }
        }
        if (g_stage_counter && g_stage_counter->load() == 0 && !has_metrics) {
            // still initializing and nothing to show yet
            continue;
        }
        rates = build_and_emit(true);
        double rec_rate = rates.first;
        double ids_rate = rates.second;
        // Adaptive interval logic: when throughput is very high, reduce update frequency
        int next_interval = interval_seconds;
        double peak = std::max(rec_rate, ids_rate);
        if (peak > 500000.0) next_interval = std::min(60, interval_seconds * 6);
        else if (peak > 100000.0) next_interval = std::min(30, interval_seconds * 3);
        else if (peak > 20000.0) next_interval = std::min(10, interval_seconds * 2);
        else next_interval = interval_seconds;
        // never update more frequently than the averaging window
        curr_interval = std::max(next_interval, (int)window_sec.count());
    }
    // When stopping, if interactive, print newline to end status line
    if (interactive) {
        std::lock_guard<std::mutex> lk(g_log_mutex);
        std::cerr << std::endl;
        if (g_log_file.is_open()) g_log_file << std::endl;
    }
}

void start_progress_report(const std::vector<std::pair<std::string, const std::atomic<size_t>*>>& counters, int interval_seconds, const std::vector<std::string>& stage_names) {
    stop_progress_report();
    g_counters = counters;
    g_stage_names = stage_names;
    // locate stage counter pointer for stage-aware logging
    g_stage_counter = nullptr;
    for (const auto &p : g_counters) {
        if (p.first == "stage") { g_stage_counter = p.second; break; }
    }
    g_report_running.store(true);
    g_report_thread = std::thread(reporter_loop, interval_seconds);
}

void stop_progress_report() {
    if (g_report_running.load()) {
        g_report_running.store(false);
        g_report_cv.notify_all();
        if (g_report_thread.joinable()) g_report_thread.join();
    }
    g_stage_counter = nullptr;
}

void notify_progress_report() {
    if (g_report_running.load()) {
        // notify the reporter thread to wake and re-draw the line
        g_report_cv.notify_all();
    }
}
