#include "Offline/TimeoutService/inc/TimeoutWatchdog.hh"

#include "art/Framework/Services/Registry/ActivityRegistry.h"

#include <cstdio>
#include <sstream>

namespace mu2e {

  TimeoutWatchdog::TimeoutWatchdog(Parameters const& config,
                                   art::ActivityRegistry& registry)
    : eventTimeoutMs_(config().eventTimeoutMs())
    , moduleTimeoutMs_(config().moduleTimeoutMs())
    , debugLevel_(config().debugLevel())
  {
    if(config().registerPreEventCallback()) {
      registry.sPreProcessEvent.watch([this](art::Event const& e, art::ScheduleContext) {
        this->startEvent(e);
      });
    }
  }

  void TimeoutWatchdog::startEvent(art::Event const& e) {
    if (tls_.event != e.event() || tls_.subRun != e.subRun() || tls_.run != e.run()) {
      // New event: reset cancellation state and apply event-level deadline.
      tls_.run = e.run();
      tls_.subRun = e.subRun();
      tls_.event = e.event();
      tls_.moduleLabel.clear();
      tls_.stopSource = std::stop_source{};
      tls_.stopToken = tls_.stopSource.get_token();

      if (eventTimeoutMs_ > 0.0) {
        Clock::time_point const now = Clock::now();
        tls_.eventDeadline = now + std::chrono::duration_cast<Clock::duration>(
                                                                               std::chrono::duration<double, std::milli>(eventTimeoutMs_));
      } else {
        tls_.eventDeadline.reset();
      }

      tls_.moduleDeadline.reset();

      if (debugLevel_ > 1) {
        std::printf("[TimeoutWatchdog::%s] Event %u:%u:%u eventTimeoutMs=%.3f\n",
                    __func__,
                    tls_.run,
                    tls_.subRun,
                    tls_.event,
                    eventTimeoutMs_);
      }
    }
  }

  double TimeoutWatchdog::moduleBudgetFor_(std::optional<double> allowedTimeMs) const {
    // Module-provided budget takes precedence over service default.
    if (allowedTimeMs) return *allowedTimeMs;
    return moduleTimeoutMs_;
  }

  void TimeoutWatchdog::startModule(std::string const& moduleLabel,
                                    std::optional<double> allowedTimeMs) {
    // Module scope starts a fresh stop token so checks are local to this invocation.
    tls_.moduleLabel = moduleLabel;
    if (tls_.stopToken.stop_requested()) { //FIXME: Need to decide if a previous module was timed out but the event wasn't, do we allow next modules to continue anyway
      if (debugLevel_ > 0) {
        std::printf("[TimeoutWatchdog::%s] Event %u:%u:%u stop already requested when starting module=%s\n",
                    __func__,
                    tls_.run,
                    tls_.subRun,
                    tls_.event,
                    tls_.moduleLabel.c_str());
      }
      // add source reset here if we want to ignore previous module's stop
    }

    double const budget = moduleBudgetFor_(allowedTimeMs);
    if (budget > 0.0) {
      Clock::time_point const now = Clock::now();
      tls_.moduleDeadline = now + std::chrono::duration_cast<Clock::duration>(
                                                                              std::chrono::duration<double, std::milli>(budget));
    } else {
      tls_.moduleDeadline.reset();
    }

    if (debugLevel_ > 1) {
      std::printf("[TimeoutWatchdog::%s] module=%s moduleTimeoutMs=%.3f\n",
                  __func__,
                  tls_.moduleLabel.c_str(),
                  budget);
    }
  }

  void TimeoutWatchdog::endModule() {
    if (debugLevel_ > 1 && !tls_.moduleLabel.empty()) {
      std::printf("[TimeoutWatchdog::%s] module=%s\n",
                  __func__,
                  tls_.moduleLabel.c_str());
    }

    // Clear module-only state; event-level deadline remains active.
    tls_.moduleDeadline.reset();
    tls_.moduleLabel.clear();
  }

  std::optional<TimeoutWatchdog::Clock::time_point>
  TimeoutWatchdog::eventDeadline() const { return tls_.eventDeadline; }

  std::optional<TimeoutWatchdog::Clock::time_point>
  TimeoutWatchdog::moduleDeadline() const { return tls_.moduleDeadline; }

  bool TimeoutWatchdog::check() {
    // Once cancellation is requested, preserve sticky behavior for callers.
    if (tls_.stopToken.stop_requested()) {
      if (debugLevel_ > 0) {
        std::printf("[TimeoutWatchdog::%s] Event %u:%u:%u stop already requested module=%s\n",
                    __func__,
                    tls_.run,
                    tls_.subRun,
                    tls_.event,
                    tls_.moduleLabel.c_str());
      }
      return true;
    }

    Clock::time_point const now = Clock::now();

    bool const timedOutEvent = tls_.eventDeadline && (now > *tls_.eventDeadline);
    bool const timedOutModule = tls_.moduleDeadline && (now > *tls_.moduleDeadline);

    if (!timedOutEvent && !timedOutModule) {
      return false;
    }

    // Transition to cancelled state so subsequent checks can short-circuit quickly.
    tls_.stopSource.request_stop();

    if (debugLevel_ > 0) {
      std::printf("[TimeoutWatchdog::%s] Event %u:%u:%u timeout module=%s eventExceeded=%d moduleExceeded=%d\n",
                  __func__,
                  tls_.run,
                  tls_.subRun,
                  tls_.event,
                  tls_.moduleLabel.c_str(),
                  timedOutEvent ? 1 : 0,
                  timedOutModule ? 1 : 0);
    }

    return true;
  }

  std::stop_token TimeoutWatchdog::stopToken() const {
    return tls_.stopToken;
  }

} // namespace mu2e
