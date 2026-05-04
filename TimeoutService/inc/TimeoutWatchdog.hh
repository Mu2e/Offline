//
// An art service to provide cooperative timeout checks for event/module processing.
//

#ifndef TimeoutService_TimeoutWatchdog_hh
#define TimeoutService_TimeoutWatchdog_hh

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceTable.h"
#include "fhiclcpp/types/Atom.h"

#include <chrono>
#include <optional>
#include <stop_token>
#include <string>

namespace art {
  class ActivityRegistry;
}

namespace mu2e {

class TimeoutWatchdog {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    // Time budgets in milliseconds. Zero or negative disables the corresponding timeout.
    fhicl::Atom<double> eventTimeoutMs{
      Name("eventTimeoutMs"),
      Comment("Event timeout budget in milliseconds; <= 0 disables event-level timeout"),
      0.0};
    fhicl::Atom<double> moduleTimeoutMs{
      Name("moduleTimeoutMs"),
      Comment("Default module timeout budget in milliseconds; <= 0 disables module-level timeout"),
      0.0};
    fhicl::Atom<bool> registerPreEventCallback{
      Name("registerPreEventCallback"),
      Comment("Whether to register the pre-event callback to start the event timer. If false, the service user must call startEvent() directly."),
      true};
    fhicl::Atom<int> debugLevel{
      Name("debugLevel"),
      Comment("Service debug verbosity: 0=off, 1=timeouts, 2=module/event boundaries"),
      0};

  };

  using Parameters = art::ServiceTable<Config>;

  TimeoutWatchdog(Parameters const& config, art::ActivityRegistry&);
  TimeoutWatchdog(TimeoutWatchdog const&) = delete;
  TimeoutWatchdog& operator=(TimeoutWatchdog const&) = delete;
  TimeoutWatchdog(TimeoutWatchdog&&) = delete;
  TimeoutWatchdog& operator=(TimeoutWatchdog&&) = delete;

  // --- called by modules (cooperative) ---
  void startEvent(art::Event const& e);
  void startModule(std::string const& moduleLabel,
                   std::optional<double> allowedTimeMs = std::nullopt);
  void endModule();

  // Returns true once the current event/module has exceeded a configured budget.
  bool check();

  // Stop token that can be used in downstream cooperative cancellation points.
  std::stop_token stopToken() const;

  // helpers
  std::optional<std::chrono::steady_clock::time_point> eventDeadline() const;
  std::optional<std::chrono::steady_clock::time_point> moduleDeadline() const;

  // RAII guard declaration (defined below)
  class ModuleGuard;

private:
  using Clock = std::chrono::steady_clock;

  struct State {
    // Identifiers for diagnostics
    art::RunNumber_t run; // event ID
    art::SubRunNumber_t subRun;
    art::EventNumber_t event;
    std::string moduleLabel; // current module

    // Deadlines (if enabled)
    std::optional<Clock::time_point> eventDeadline;
    std::optional<Clock::time_point> moduleDeadline;

    // Cancellation state for cooperative checks.
    std::stop_source stopSource;
    std::stop_token stopToken;

    State()
      : run(0), subRun(0), event(0)
      , moduleLabel()
      , eventDeadline()
      , moduleDeadline()
      , stopSource()
      , stopToken(stopSource.get_token()) {}
  };

  double moduleBudgetFor_(std::optional<double> allowedTimeMs) const;

  State tls_;

  double eventTimeoutMs_;
  double moduleTimeoutMs_;
  int debugLevel_;
};

class TimeoutWatchdog::ModuleGuard {
public:
  ModuleGuard(TimeoutWatchdog& svc,
              art::Event const& e,
              std::string const& moduleLabel,
              std::optional<double> allowedTimeMs = std::nullopt)
    : svc_{svc}
  {
    svc_.startEvent(e); // Only necessary if the service is configured with registerPreEventCallback=false
    svc_.startModule(moduleLabel, allowedTimeMs);
  }

  ~ModuleGuard() noexcept {
    // Never throw from destructor
    svc_.endModule();
  }

  bool check() const { return svc_.check(); }

  std::stop_token stopToken() const { return svc_.stopToken(); }

private:
  TimeoutWatchdog& svc_;
};

} // namespace mu2e

DECLARE_ART_SERVICE(mu2e::TimeoutWatchdog, SHARED)

#endif /* TimeoutService_TimeoutWatchdog_hh */
