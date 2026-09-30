//
// An art service to provide cooperative timeout checks for event/module processing.
//
// Deadlines and cancellation state are kept per art schedule: every call names the
// schedule it acts on, so concurrent schedules never touch the same state.
//

#ifndef TimeoutService_TimeoutWatchdog_hh
#define TimeoutService_TimeoutWatchdog_hh

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceTable.h"
#include "art/Utilities/PerScheduleContainer.h"
#include "art/Utilities/ScheduleID.h"
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
  // sid is the schedule processing the event: scheduleID() in a legacy module,
  // ProcessingFrame::scheduleID() in a shared one.
  void startEvent(art::Event const& e, art::ScheduleID sid);
  void startModule(art::ScheduleID sid,
                   std::string const& moduleLabel,
                   std::optional<double> allowedTimeMs = std::nullopt);
  void endModule(art::ScheduleID sid);

  // Returns true once the current event/module has exceeded a configured budget.
  bool check(art::ScheduleID sid);

  // Stop token that can be used in downstream cooperative cancellation points.
  std::stop_token stopToken(art::ScheduleID sid) const;

  // helpers
  std::optional<std::chrono::steady_clock::time_point> eventDeadline(art::ScheduleID sid) const;
  std::optional<std::chrono::steady_clock::time_point> moduleDeadline(art::ScheduleID sid) const;

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

  // One State per schedule, sized at construction and never resized, so each
  // schedule reads and writes only its own element.
  art::PerScheduleContainer<State> states_;

  double eventTimeoutMs_;
  double moduleTimeoutMs_;
  int debugLevel_;
};

class TimeoutWatchdog::ModuleGuard {
public:
  ModuleGuard(TimeoutWatchdog& svc,
              art::Event const& e,
              art::ScheduleID sid,
              std::string const& moduleLabel,
              std::optional<double> allowedTimeMs = std::nullopt)
    : svc_{svc}
    , sid_{sid}
  {
    svc_.startEvent(e, sid_); // Only necessary if the service is configured with registerPreEventCallback=false
    svc_.startModule(sid_, moduleLabel, allowedTimeMs);
  }

  ~ModuleGuard() noexcept {
    // Never throw from destructor
    svc_.endModule(sid_);
  }

  bool check() const { return svc_.check(sid_); }

  std::stop_token stopToken() const { return svc_.stopToken(sid_); }

private:
  TimeoutWatchdog& svc_;
  art::ScheduleID sid_;
};

} // namespace mu2e

DECLARE_ART_SERVICE(mu2e::TimeoutWatchdog, SHARED)

#endif /* TimeoutService_TimeoutWatchdog_hh */
