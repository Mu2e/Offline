#ifndef DbService_GrlBit_hh
#define DbService_GrlBit_hh

// definition of a bits that can be set within a given word

#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

class GrlBit {
 public:
  typedef std::vector<GrlBit> BitVec;

  GrlBit(std::string const& name, std::string const& bitname, int bitnumber,
         std::string const& description = "",
         std::string const& create_time = "",
         std::string const& create_user = "") :
      _name(name),
      _bitname(bitname), _bitnumber(bitnumber), _description(description),
      _create_time(create_time), _create_user(create_user) {}

  std::string const& name() const { return _name; }
  std::string const& bitname() const { return _bitname; }
  int bitnumber() const { return _bitnumber; }
  std::string const& description() const { return _description; }
  std::string const& createTime() const { return _create_time; }
  std::string const& createUser() const { return _create_user; }

 private:
  std::string _name;
  std::string _bitname;
  int _bitnumber;
  std::string _description;
  std::string _create_time;
  std::string _create_user;
};

}  // namespace mu2e

#endif
