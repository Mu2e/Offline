#ifndef DbService_RunSelect_hh
#define DbService_RunSelect_hh

// holds selection criteria for runs for runTool

#include <string>

namespace mu2e {

class RunSelect {
 public:

  // run:
  //     runnumber, like "110000", means just this one run
  //     start-stop (inclusive) like "110000-120000",
  //     empty means all
  // last: return only last N runs, 0 means all
  // type, like "4" measn only runs of the type (see flags), empty is all
  // time, TIME1 or TIME1/TIME2  for runs since a time, or in the range
  // days only runs take in the last N days, 0 means all
  RunSelect(const std::string& run, int last, const std::string& type,
            const std::string& time, int days) :
      _run(run),
      _last(last), _type(type), _time(time), _days(days) {}

  const std::string& run() const { return _run; }
  int last() const { return _last; }
  const std::string& type() const { return _type; }
  const std::string& time() const { return _time; }
  int days() const { return _days; }

 private:
  std::string _run;
  int _last;
  std::string _type;
  std::string _time;
  int _days;
};

}  // namespace mu2e

#endif
