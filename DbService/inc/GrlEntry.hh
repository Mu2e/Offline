#ifndef DbService_GrlEntry_hh
#define DbService_GrlEntry_hh

// holds a goodrun entry - a name of the bit word, the values of the bit word
// and the interval of validity

#include "Offline/DbTables/inc/DbIoV.hh"
#include <string>
#include <vector>

namespace mu2e {

class GrlEntry {
 public:
  typedef std::vector<GrlEntry> EntryVec;

  GrlEntry(std::string const& name, DbIoV const& iov, int value = 0,
           int eid = -1, int retired = 0, std::string const& create_user = "",
           std::string const& create_time = "") :
      _name(name),
      _iov(iov), _value(value), _eid(eid), _retired(retired),
      _create_user(create_user), _create_time(create_time) {}

  const std::string& name() const { return _name; }
  const DbIoV& iov() const { return _iov; }
  int value() const { return _value; }
  int eid() const { return _eid; }
  int retired() const { return _retired; }
  const std::string& createUser() const { return _create_user; }
  const std::string& createTime() const { return _create_time; }

 private:
  std::string _name;
  DbIoV _iov;
  int _value;
  int _eid;
  int _retired;
  std::string _create_user;
  std::string _create_time;
};

}  // namespace mu2e

#endif
