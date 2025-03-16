#ifndef DbService_GrlHeader_hh
#define DbService_GrlHeader_hh

// holds a defintion of a goodrun list - name, etc
// locked==0 is unlocked, locked >0 is locked

#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

class GrlHeader {
 public:
  typedef std::vector<GrlHeader> HeaderVec;
  GrlHeader(const std::string& name, int locked = 0, int lid = -1,
            std::string const& create_time = "",
            std::string const& create_user = "") :
      _name(name),
      _locked(locked), _lid(lid), _create_time(create_time),
      _create_user(create_user) {}
  const std::string& name() const { return _name; }
  int locked() const { return _locked; }
  int lid() const { return _lid; }
  const std::string& createTime() const { return _create_time; }
  const std::string& createUser() const { return _create_user; }
  std::string formatted() const;

 private:
  std::string _name;
  int _locked;
  int _lid;
  std::string _create_time;
  std::string _create_user;
};

}  // namespace mu2e

#endif
