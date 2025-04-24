#ifndef DbService_GrlWord_hh
#define DbService_GrlWord_hh

// definition of a word that can have bits set for a GRL entry
// it is just the name of the word

#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

class GrlWord {
 public:
  typedef std::vector<GrlWord> WordVec;

  GrlWord(std::string const& name, std::string const& description = "",
          std::string const& create_time = "",
          std::string const& create_user = "") :
      _name(name),
      _description(description), _create_time(create_time),
      _create_user(create_user) {}

  std::string const& name() const { return _name; }
  std::string const& description() const { return _description; }
  std::string const& createTime() const { return _create_time; }
  std::string const& createUser() const { return _create_user; }

 private:
  std::string _name;
  std::string _description;
  std::string _create_time;
  std::string _create_user;
};

}  // namespace mu2e

#endif
