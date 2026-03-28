#include "Offline/TimeoutService/inc/TimeoutWatchdog.hh"

#include <cstdio>
#include <sstream>

namespace mu2e {

// Thread-local state keeps timeout bookkeeping independent per worker thread.
thread_local TimeoutWatchdogService::State TimeoutWatchdogService::tls_{};

TimeoutWatchdogService::TimeoutWatchdogService(Parameters const& config,
                                               art::ActivityRegistry&)
  : eventTimeoutMs_(config().eventTimeoutMs())
  , moduleTimeoutMs_(config().moduleTimeoutMs())
  , debugLevel_(config().debugLevel())
{}

static std::string eventIdString_(art::Event const& e) {
  std::ostringstream os;
  os << e.id().run() << ":" << e.id().subRun() << ":" << e.id().event();
  return os.str();
}

void TimeoutWatchdogService::startEvent(art::Event const& e) {
  std::string const id = eventIdString_(e);
  if (tls_.eventId != id) {
    // New event: reset cancellation state and apply event-level deadline.
    tls_.eventId = id;
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
      std::printf("[TimeoutWatchdogService::%s] event=%s eventTimeoutMs=%.3f\n",
                  __func__,
                  tls_.eventId.c_str(),
                  eventTimeoutMs_);
    }
  }
}

double TimeoutWatchdogService::moduleBudgetFor_(std::optional<double> allowedTimeMs) const {
  // Module-provided budget takes precedence over service default.
  if (allowedTimeMs) return *allowedTimeMs;
  return moduleTimeoutMs_;
}

void TimeoutWatchdogService::startModule(std::string const& moduleLabel,
                                         std::optional<double> allowedTimeMs) {
  // Module scope starts a fresh stop token so checks are local to this invocation.
  tls_.moduleLabel = moduleLabel;
  tls_.stopSource = std::stop_source{};
  tls_.stopToken = tls_.stopSource.get_token();

  double const budget = moduleBudgetFor_(allowedTimeMs);
  if (budget > 0.0) {
    Clock::time_point const now = Clock::now();
    tls_.moduleDeadline = now + std::chrono::duration_cast<Clock::duration>(
      std::chrono::duration<double, std::milli>(budget));
  } else {
    tls_.moduleDeadline.reset();
  }

  if (debugLevel_ > 1) {
    std::printf("[TimeoutWatchdogService::%s] module=%s moduleTimeoutMs=%.3f\n",
                __func__,
                tls_.moduleLabel.c_str(),
                budget);
  }
}

void TimeoutWatchdogService::endModule() {
  if (debugLevel_ > 2 && !tls_.moduleLabel.empty()) {
    std::printf("[TimeoutWatchdogService::%s] module=%s\n",
                __func__,
                tls_.moduleLabel.c_str());
  }

  // Clear module-only state; event-level deadline remains active.
  tls_.moduleDeadline.reset();
  tls_.moduleLabel.clear();
}

std::optional<TimeoutWatchdogService::Clock::time_point>
TimeoutWatchdogService::eventDeadline() const { return tls_.eventDeadline; }

std::optional<TimeoutWatchdogService::Clock::time_point>
TimeoutWatchdogService::moduleDeadline() const { return tls_.moduleDeadline; }

bool TimeoutWatchdogService::check() const {
  // Once cancellation is requested, preserve sticky behavior for callers.
  if (tls_.stopToken.stop_requested()) {
    if (debugLevel_ > 0) {
      std::printf("[TimeoutWatchdogService::%s] stop already requested event=%s module=%s\n",
                  __func__,
                  tls_.eventId.c_str(),
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
    std::printf("[TimeoutWatchdogService::%s] timeout event=%s module=%s eventExceeded=%d moduleExceeded=%d\n",
                __func__,
                tls_.eventId.c_str(),
                tls_.moduleLabel.c_str(),
                timedOutEvent ? 1 : 0,
                timedOutModule ? 1 : 0);
  }

  return true;
}

std::stop_token TimeoutWatchdogService::stopToken() const {
  return tls_.stopToken;
}

} // namespace mu2e
