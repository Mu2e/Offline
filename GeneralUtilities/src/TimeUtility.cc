#include "Offline/GeneralUtilities/inc/TimeUtility.hh"
#include <iomanip>

using namespace mu2e;

int TimeUtility::parseTimeTZ(std::string const& stime, std::time_t& ttime) {

  ttime = 0;

  // must at least be a date
  if (stime.size() < 10) {
    return 1;
  }

  std::string ftime(24, ' ');
  std::string zero = "2000-01-01 00:00:00";
  // copy data and time part
  for (size_t ipt = 0; ipt < 19; ipt++) {
    if (ipt != 10) {
      if (ipt < stime.size()) {
        ftime[ipt] = stime[ipt];
      } else {
        ftime[ipt] = zero[ipt];
      }
    }
  }

  // skip ns if present
  size_t ipt = 19;
  if (stime[19] == '.') {
    // advance past ns
    while (stime[ipt] != '-' && stime[ipt] != '+' && ipt < stime.size())
      ipt++;
  }

  // sign and hour digits of time zone offset
  int tzi = 0;
  if (stime.size() >= ipt + 3) {
    std::string tz = stime.substr(ipt, 3);
    tzi = std::stoi(tz);
  }

  std::tm tt{}; // must initialize to zero
  std::istringstream ss(ftime);
  ss >> std::get_time(&tt, "%Y-%m-%d %H:%M:%S");
  if (ss.fail()) {
    return 1;
  }
  // these don't usually get set, but make sure they aren't used
  tt.tm_isdst = 0;
  tt.tm_gmtoff = 0;
  // here is the main point - interpret date-time as UTC,
  // convert to epoch time, then correct for explicit tz offset
  ttime = timegm(&tt) - tzi * 3600;

  return 0;
}

std::string TimeUtility::reformat1(const std::string& stime) {
    std::string result = stime;

    size_t dotPos = result.find('.');

    if (dotPos != std::string::npos) { // Found a decimal point
        size_t endPos = dotPos + 1;

        while (endPos < result.length() && std::isdigit(result[endPos])) {
            endPos++;
        }

        result.erase(dotPos, endPos - dotPos); // Remove the decimal part
    }
    result[10] = 'T';
    return result;
}
