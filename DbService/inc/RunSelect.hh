#ifndef DbService_RunSelect_hh
#define DbService_RunSelect_hh

// holds selection criteria for runs for runTool

#include <string>

namespace mu2e {

class RunSelect {
 public:

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
