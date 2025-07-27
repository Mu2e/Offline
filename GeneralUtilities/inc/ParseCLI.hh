#ifndef GeneralUtilities_ParseCLI_hh
#define GeneralUtilities_ParseCLI_hh

//
// a class to define and interpret command line arguments,
//   can work with or without a subcommand structure
//

#include "Offline/GeneralUtilities/inc/StringVec.hh"
#include <iostream>
#include <string>

namespace mu2e {

class ParseCLI {

public:
  // autohelp means respond to "-h" or "--help" with appropriate output from declarations
  // verbose is for debugging
  ParseCLI(const std::string& helpstr = "", bool autohelp = true, bool verbose = false) :
      _helpstr(helpstr), _autohelp(autohelp), _verbose(verbose) {}

  // required for each subcommand, including the blank global subcommand ""
  int addSubcommand(const std::string& subcommand, const std::string& helpstr = "");
  // subcommand or "" for global switches
  // name = for referring to the switch, always referenced together with subcommand
  // sname = short name, like "v" for CLI "-v"
  // lname = long name, like "verbose" for CLI "--verbose"
  // args = true if this switch should have a following argument
  // helpstr = help string for this switch
  // defaultstr = default switch argument
  // repeated = true if swtich can be repeated to increase count or make a list,
  //       if false, last instance will be used
  // required = parse will fail if not present and no default provided
  int addSwitch(const std::string& subcommand, const std::string& name, const std::string& sname,
                const std::string& lname, bool args = false, const std::string& helpstr = "",
                const std::string& defaultstr = "", bool repeated = false,
                bool required = false);

  // do not remove the command name in the first arg
  int setArgs(int argc, char** argv);

  const std::string& subcommand() const { return _subcommand; }
  const StringVec& positionals() const { return _positionals; }

  int getCount(const std::string& subcommand, const std::string& name) const;
  std::string getString(const std::string& subcommand, const std::string& name) const;
  StringVec getStrings(const std::string& subcommand, const std::string& name) const;
  int getInt(const std::string& subcommand, const std::string& name) const;
  std::vector<int> getInts(const std::string& subcommand, const std::string& name) const;
  float getFloat(const std::string& subcommand, const std::string& name) const;
  std::vector<float> getFloats(const std::string& subcommand, const std::string& name) const;

  // useful to see the result of defining and parsing the args
  void print();

private:
  static constexpr size_t BADIND{9999};
  size_t findItem(const std::string& subcommand, const std::string& name) const;
  void autohelp() const;

  std::string _helpstr;
  bool _autohelp;
  bool _verbose;

  struct subitem {
    std::string subcommand;
    std::string helpstr;
  };

  struct item {
    std::string subcommand;
    std::string name;
    std::string sname;
    std::string lname;
    bool args;
    std::string helpstr;
    std::string defaultstr;
    bool repeated;
    bool required;
    StringVec values;
  };
  std::vector<subitem> _subs;
  std::vector<item> _items;
  std::string _command;
  std::string _subcommand;
  StringVec _positionals;
};

} // namespace mu2e

#endif
