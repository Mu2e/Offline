#include "Offline/TimeoutService/inc/TimeoutWatchdog.hh"

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Persistency/Provenance/ScheduleContext.h"

#include <cstdio>
#include <sstream>

namespace mu2e {

  TimeoutWatchdog::TimeoutWatchdog(Parameters const& config,
                                   art::ActivityRegistry& registry)
    : eventTimeoutMs_(config().eventTimeoutMs())
    , moduleTimeoutMs_(config().moduleTimeoutMs())
    , debugLevel_(config().debugLevel())
  {
    states_.expand_to_num_schedules();

    if(config().registerPreEventCallback()) {
      registry.sPreProcessEvent.watch([this](art::Event const& e, art::ScheduleContext const sc) {
        this->startEvent(e, sc.id());
      });
    }
  }

  void TimeoutWatchdog::startEvent(art::Event const& e, art::ScheduleID const sid) {
    State& state = states_.at(sid);
    if (state.event != e.event() || state.subRun != e.subRun() || state.run != e.run()) {
      // New event: reset cancellation state and apply event-level deadline.
      state.run = e.run();
      state.subRun = e.subRun();
      state.event = e.event();
      state.moduleLabel.clear();
      state.stopSource = std::stop_source{};
      state.stopToken = state.stopSource.get_token();

      if (eventTimeoutMs_ > 0.0) {
        Clock::time_point const now = Clock::now();
        state.eventDeadline = now + std::chrono::duration_cast<Clock::duration>(
                                                                                std::chrono::duration<double, std::milli>(eventTimeoutMs_));
      } else {
        state.eventDeadline.reset();
      }

      state.moduleDeadline.reset();

      if (debugLevel_ > 1) {
        std::printf("[TimeoutWatchdog::%s] schedule %u Event %u:%u:%u eventTimeoutMs=%.3f\n",
                    __func__,
                    static_cast<unsigned>(sid.id()),
                    state.run,
                    state.subRun,
                    state.event,
                    eventTimeoutMs_);
      }
    }
  }

  double TimeoutWatchdog::moduleBudgetFor_(std::optional<double> allowedTimeMs) const {
    // Module-provided budget takes precedence over service default.
    if (allowedTimeMs) return *allowedTimeMs;
    return moduleTimeoutMs_;
  }

  void TimeoutWatchdog::startModule(art::ScheduleID const sid,
                                    std::string const& moduleLabel,
                                    std::optional<double> allowedTimeMs) {
    // A stop request is sticky for the rest of the event: if an earlier module
    // timed out, check() returns true for this module immediately.
    State& state = states_.at(sid);
    state.moduleLabel = moduleLabel;
    if (state.stopToken.stop_requested()) { //FIXME: Need to decide if a previous module was timed out but the event wasn't, do we allow next modules to continue anyway
      if (debugLevel_ > 0) {
        std::printf("[TimeoutWatchdog::%s] Event %u:%u:%u stop already requested when starting module=%s\n",
                    __func__,
                    state.run,
                    state.subRun,
                    state.event,
                    state.moduleLabel.c_str());
      }
      // add source reset here if we want to ignore previous module's stop
    }

    double const budget = moduleBudgetFor_(allowedTimeMs);
    if (budget > 0.0) {
      Clock::time_point const now = Clock::now();
      state.moduleDeadline = now + std::chrono::duration_cast<Clock::duration>(
                                                                               std::chrono::duration<double, std::milli>(budget));
    } else {
      state.moduleDeadline.reset();
    }

    if (debugLevel_ > 1) {
      std::printf("[TimeoutWatchdog::%s] module=%s moduleTimeoutMs=%.3f\n",
                  __func__,
                  state.moduleLabel.c_str(),
                  budget);
    }
  }

  void TimeoutWatchdog::endModule(art::ScheduleID const sid) {
    State& state = states_.at(sid);
    if (debugLevel_ > 1 && !state.moduleLabel.empty()) {
      std::printf("[TimeoutWatchdog::%s] module=%s\n",
                  __func__,
                  state.moduleLabel.c_str());
    }

    // Clear module-only state; event-level deadline remains active.
    state.moduleDeadline.reset();
    state.moduleLabel.clear();
  }

  std::optional<TimeoutWatchdog::Clock::time_point>
  TimeoutWatchdog::eventDeadline(art::ScheduleID const sid) const { return states_.at(sid).eventDeadline; }

  std::optional<TimeoutWatchdog::Clock::time_point>
  TimeoutWatchdog::moduleDeadline(art::ScheduleID const sid) const { return states_.at(sid).moduleDeadline; }

  bool TimeoutWatchdog::check(art::ScheduleID const sid) {
    State& state = states_.at(sid);
    // Once cancellation is requested, preserve sticky behavior for callers.
    if (state.stopToken.stop_requested()) {
      if (debugLevel_ > 0) {
        std::printf("[TimeoutWatchdog::%s] Event %u:%u:%u stop already requested module=%s\n",
                    __func__,
                    state.run,
                    state.subRun,
                    state.event,
                    state.moduleLabel.c_str());
      }
      return true;
    }

    Clock::time_point const now = Clock::now();

    bool const timedOutEvent = state.eventDeadline && (now > *state.eventDeadline);
    bool const timedOutModule = state.moduleDeadline && (now > *state.moduleDeadline);

    if (!timedOutEvent && !timedOutModule) {
      return false;
    }

    // Transition to cancelled state so subsequent checks can short-circuit quickly.
    state.stopSource.request_stop();

    if (debugLevel_ > 0) {
      std::printf("[TimeoutWatchdog::%s] Event %u:%u:%u timeout module=%s eventExceeded=%d moduleExceeded=%d\n",
                  __func__,
                  state.run,
                  state.subRun,
                  state.event,
                  state.moduleLabel.c_str(),
                  timedOutEvent ? 1 : 0,
                  timedOutModule ? 1 : 0);
    }

    return true;
  }

  std::stop_token TimeoutWatchdog::stopToken(art::ScheduleID const sid) const {
    return states_.at(sid).stopToken;
  }

} // namespace mu2e
