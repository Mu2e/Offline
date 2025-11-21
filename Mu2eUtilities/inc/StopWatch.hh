// A simple class to track processing times and report them
// Tracks at the 0.01 microsecond level, so not suitable for below this level
//
// Original author: Michael MacKenzie, 2025

#ifndef __MU2E_STOPWATCH__
#define __MU2E_STOPWATCH__

// C++ includes
#include <unordered_map>
#include <string>
#include <chrono>
#include <limits>
#include <format>
#include <iostream>

// ROOT includes
#include "TString.h"

namespace mu2e {
  class StopWatch {
  public:

    // a time tracker instance
    struct Timer_t {
      std::string _name;
      std::chrono::steady_clock::time_point _last_time;
      double _duration;
      double _max_duration;
      double _min_duration;
      unsigned _count;
      bool _active;
      double _calibration; // time offset per time check, in 0.01 us

      // for storing time in internal units
      constexpr static double _ns_to_unit = 0.1; // store time in 0.01 us
      constexpr static double _us_to_unit = 1000.*_ns_to_unit;
      constexpr static double _unit_to_us = 1./_ns_to_unit/1.e3;
      constexpr static double _unit_to_ms = 1./_ns_to_unit/1.e6;
      constexpr static double _unit_to_s  = 1./_ns_to_unit/1.e9;

      Timer_t(std::string name = "default", double calib = 0.1) : _name(name) { Reset(); SetCalibration(calib); }

      // Capture the duration time, increment the counter
      double Increment() {
        // If not active, start the clock
        if(!_active) {
          SetTime();
          return 0.;
        }

        // If active, return the duration and reset the clock
        const auto time_now = std::chrono::steady_clock::now();
        const double curr_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(time_now-_last_time).count()*_ns_to_unit - _calibration;
        _duration += curr_duration;
        ++_count;
        _min_duration = std::min(_min_duration, curr_duration);
        _max_duration = std::max(_max_duration, curr_duration);
        _last_time = time_now;
        return curr_duration*_unit_to_us; // return the time spent in the last iteration, in micro seconds
      }

      // Set the expected timing impact, given in us
      void SetCalibration(const double calib) { _calibration = calib*_us_to_unit; }

      // Initialize the starting point time, set active
      void SetTime() { _last_time = std::chrono::steady_clock::now(); _active = true;}

      // Store the duration (if running) and stop the clock
      double StopTime() {
        const double time = Increment();
        _active = false;
        return time;
      }

      // Reset the time
      void Reset()   { _duration = 0.; _count = 0; _max_duration = 0.; _min_duration = std::numeric_limits<double>::max(); SetTime(); _active = false; }

      // Accessors
      std::string Name() const { return _name; }
      unsigned Count  () const { return _count; }
      double   Time   () const { return _duration*_unit_to_s ; } // in seconds
      double   AvgTime() const { return (_count > 0) ? _duration / _count * _unit_to_us : 0.; } // in us
      double   AvgRate() const { return (_duration > 0.) ? _count / Time() : 0.; } // in Hz
      double   MaxTime() const { return _max_duration*_unit_to_us; } // in us
      double   MinTime() const { return (_count > 0) ? _min_duration*_unit_to_us : 0.; } // in us

    };

    // Constructor
    StopWatch() : _calibration(0.1) {}

    // Retrieve a timing instance, create it if needed
    Timer_t& Timer(std::string name) {
      if(!_timers.count(name)) _timers[name] = Timer_t(name, _calibration);
      return _timers[name];
    }

    // Increment a timer, returning the duration of the last iteration
    double Increment(std::string name) {
      const bool first = !_timers.count(name);
      Timer_t& timer = Timer(name);
      double duration = 0.;
      if(first) timer.SetTime(); // on the first call, just initialize the timer
      else duration = timer.Increment();
      return duration;
    }
    void IncrementAll() { for(auto timer : _timers) timer.second.Increment(); }
    void SetTime(std::string name) {
      Timer_t& timer = Timer(name);
      timer.SetTime();
    }
    double StopTime(std::string name) {
      Timer_t& timer = Timer(name);
      return timer.StopTime();
    }
    void StopAll() { for(auto timer : _timers) timer.second.StopTime(); }

    void Calibrate(int iterations = 1e5) {
      _calibration = 0.; // turn off calibration
      for(int itest = 0; itest < iterations; ++itest) Increment("AAA-TimeCalib");
      Timer("AAA-TimeCalib").StopTime();
      _calibration = Timer("AAA-TimeCalib").AvgTime(); // update the calibration
    }
    double Calibration() { return _calibration; }

    // Print a timer info
    void PrintTimer(const Timer_t& timer, std::ostream& os = std::cout, bool print_header = true) const {
      if(print_header) os << std::format("{:>30s}: {:>14s} {:>14s} {:>15s} {:>10s} {:>10s} {:>12s} {:>15s}\n", "Timer", "Time (s)", "Avg time (ms)", "Avg rate (kHz)",
                                         "Min (ms)", "Max (ms)", "Count", "Overhead (ms)");
      os << std::format("{:>30s}: {:14.5g} {:14.5g} {:15.5g} {:10.3g} {:10.3g} {:12d} {:15.3g}\n",
                        timer.Name().c_str(),
                        timer.Time(), // s
                        timer.AvgTime() / 1000., // ms
                        timer.AvgRate() / 1000., // kHz
                        timer.MinTime() / 1000., // ms
                        timer.MaxTime() / 1000., // ms
                        timer.Count(),
                        timer.Count()*_calibration / 1000. // ms
                        );
    }
    void Print(std::ostream& os = std::cout) const {
      if(_timers.size() == 0) return;

      // Sort by total time
      std::vector<const Timer_t*> timers;
      for(auto timer : _timers) {
        if(timer.first == "AAA-TimeCalib") continue; // ignore the calibration
        timers.push_back(&(_timers.at(timer.first))); // get a stable pointer
      }
      std::sort(timers.begin(), timers.end(),
                [&](const Timer_t* t_1, const Timer_t* t_2) {
                  return t_1->Time() > t_2->Time();// sort in reverse
                });

      // Print out each timer
      constexpr int line_width = 130;
      os << std::format("{:-<{}s}\n", "", line_width);
      bool first = true;
      for(const auto* timer : timers) {
        PrintTimer(*timer, os, first);
        first = false;
      }
      os << std::format("{:-<{}s}\n", "", line_width);
    }


  private:
    // List of time instances by name
    std::unordered_map<std::string, Timer_t> _timers; // for tracking processing time, use unordered for O(1) lookup time
    double _calibration; // in microseconds
  };

  std::ostream& operator<<(std::ostream& os, const StopWatch& watch) {
    watch.Print(os);
    return os;
  }
}

#endif
