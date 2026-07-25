#include "Offline/DbService/inc/RunInfo.hh"

#include "Offline/GeneralUtilities/inc/TimeUtility.hh"

#include "cetlib_except/exception.h"

#include <ctime>

using namespace mu2e;

//**************************************************
int RunInfo::isDone(int timeout_seconds) const {
  // Check if we have the required data
  if (_transitions.empty() || _subruns.empty()) {
    throw cet::exception("RUNINFO_INCOMPLETE")
        << "RunInfo::isDone() cannot determine status: "
        << "transitions or subruns list is empty\n";
  }

  // Get the last transition
  const RunTransition& lastTransition = _transitions.back();
  int lastTransitionType = lastTransition.typeId();

  // Check for user stop (transitions 0, 1, 6, 7)
  if (lastTransitionType == 0 || lastTransitionType == 1 ||
      lastTransitionType == 6 || lastTransitionType == 7) {
    return 1;  // user stop
  }

  // Check for error stop (transition 2)
  if (lastTransitionType == 2) {
    return 2;  // error stop
  }

  // Check for timeout (probable crash)
  // Get the last transition time
  std::string lastTransTime = lastTransition.transitionTime();
  std::time_t lastTransTimeT;
  int rc = TimeUtility::parseTimeTZ(lastTransTime, lastTransTimeT);
  if (rc != 0) {
    throw cet::exception("RUNINFO_BAD_TIME")
        << "RunInfo::isDone() cannot parse transition time: " << lastTransTime
        << "\n";
  }

  // Get the last subrun stop time
  const RunSubRun& lastSubrun = _subruns.back();
  std::time_t lastSubrunTimeT = lastSubrun.stopTimeUnix();

  // Use the more recent of the two times
  std::time_t lastActivityTime = std::max(lastTransTimeT, lastSubrunTimeT);

  // Get current time
  std::time_t currentTime = std::time(nullptr);

  // Check if timeout exceeded
  if ((currentTime - lastActivityTime) > timeout_seconds) {
    return 3;  // probably crash end
  }

  // Still running
  return 0;
}
