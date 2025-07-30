#include "Offline/GeneralUtilities/inc/ParseCLI.hh"
#include <span>

using namespace mu2e;

/*********************************************************/
int ParseCLI::addSubcommand(const std::string& subcommand, const std::string& helpstr) {
  subitem ss;
  ss.subcommand = subcommand;
  ss.helpstr = helpstr;
  for (const auto& si : _subs) {
    if (ss.subcommand == si.subcommand) {
      std::cout << "Error: subcommand " << ss.subcommand << " is repeated\n";
      return 1;
    }
  }
  _subs.emplace_back(ss);
  return 0;
}

/*********************************************************/
int ParseCLI::addSwitch(const std::string& subcommand, const std::string& name,
                        const std::string& sname, const std::string& lname, bool args,
                        const std::string& helpstr, const std::string& defaultstr, bool repeated,
                        bool required) {
  item ii;
  ii.subcommand = subcommand;
  ii.name = name;
  ii.sname = sname;
  ii.lname = lname;
  ii.args = args;
  ii.helpstr = helpstr;
  ii.defaultstr = defaultstr;
  ii.repeated = repeated;
  ii.required = required;
  for (auto& jj : _items) {
    if (ii.subcommand == jj.subcommand) {
      if (ii.name == jj.name) {
        std::cout << "Error: subcommand " << ii.subcommand << " has repeated switch " << ii.name
                  << "\n";
        return 1;
      }
    }
  }
  _items.push_back(ii);
  return 0;
}

/*********************************************************/
void ParseCLI::print() {
  std::cout << "current subcommand : " << _subcommand << "\n";
  if (_subs.size() > 0) {
    std::cout << "subcommands : \n";
    for (const auto& x : _subs) {
      if(x.subcommand.empty()) {
        std::cout << "  " << "(blank)" << "\n";
      } else {
        std::cout << "  " << x.subcommand << "\n";
      }
    }
  }
  if (_positionals.size() > 0) {
    std::cout << "positionals : \n";
    for (const auto& x : _positionals) {
      std::cout << "  " << x << "\n";
    }
  }
  for (const auto& ii : _items) {
    std::cout << "name         : " << ii.name << "\n";
    std::cout << "  subcommand : " << ii.subcommand << "\n";
    std::cout << "  sname      : " << ii.sname << "\n";
    std::cout << "  lname      : " << ii.lname << "\n";
    std::cout << "  args       : " << (ii.args ? "t" : "f") << "\n";
    std::cout << "  help       : " << ii.helpstr << "\n";
    std::cout << "  default    : " << ii.lname << "\n";
    std::cout << "  repeated   : " << (ii.repeated ? "t" : "f") << "\n";
    std::cout << "  required   : " << (ii.required ? "t" : "f") << "\n";
    std::cout << "  values     : \n";
    for (const auto& x : ii.values) {
      std::cout << "    " << x << "\n";
    }
  }
}

/*********************************************************/
int ParseCLI::setArgs(int argc, char** argv) {
  StringVec a0, a1, aa;

  if (_autohelp) { // create help switches
    for (const auto& ss : _subs)
      addSwitch(ss.subcommand, "help", "h", "help", false, "print help","");
  }

  std::span<char*> args_span(argv, static_cast<size_t>(argc));

  _command = args_span[0];
  size_t last_slash_pos = _command.rfind('/');

  if (last_slash_pos != std::string::npos) _command = _command.substr(last_slash_pos + 1);

  for (size_t i = 1; i < args_span.size(); ++i) {
    a0.emplace_back(args_span[i]); // Safely access via span::at()
  }

  // loop over a0 and expand args as needed into words in a1
  // "-ab=value" into "-a -b value" and --key=value into "--key value"
  for (const auto& a : a0) {
    if (a[0] == '-') {
      // split -k=value or --key=value into two words
      size_t ie = a.find('=');
      std::string p1, p2;
      if (ie != std::string::npos) {
        p1 = a.substr(0, ie);
        p2 = a.substr(ie + 1);
      } else {
        p1 = a;
      }
      if (p1.size() > 1 && p1[1] != '-') { // if -k, not --key
        // splits -abc into -a -b -c
        for (size_t i = 1; i < p1.size(); i++) {
          std::string ss = "-" + std::string(1, p1[i]);
          a1.emplace_back(ss); // save -a
        }
      } else {
        a1.emplace_back(p1); // save --key
      }
      if (p2.size() > 0)
        a1.emplace_back(p2); // save value from -k=value
    } else {
      // this is word that doesn't begin with "-", save it
      a1.emplace_back(a);
    }
  }

  // now all words are in order in a1
  if (_verbose) {
    std::cout << "start parsed words:\n";
    for (const auto& a : a1) {
      std::cout << "  " << a << "\n";
    }
    std::cout << "end parsed words\n";
  }

  // now compare a1[] words against the defined switches
  // put results in _items.values
  size_t i = 0;
  while (i < a1.size()) {
    std::string a = a1[i];
    if (_verbose)
      std::cout << "processing : " << a << "\n";
    if (a[0] != '-') {
      // if first non-switch then it is the subcommand, else a positional arg
      if (_subcommand.empty()) {
        if (_verbose)
          std::cout << "found subcommand : " << a << "\n";

        _subcommand = a;

        // check for valid subcommand
        bool found = false;
        for (const auto& i : _items) {
          if (_subcommand == i.subcommand)
            found = true;
        }
        if (!found) {
          std::cout << "Error: subcommand " << _subcommand << " was not declared\n";
          return 1;
        }

      } else {
        if (_verbose)
          std::cout << "adding postional : " << a << "\n";
        _positionals.emplace_back(a);
      }
    } else {                             // word starts with dash
      std::string sname, lname;          // for current switch
      if (a.size() > 1 && a[1] == '-') { // words starts with two dashes
        lname = a.substr(2);
      } else { // one letter switch
        sname = a.substr(1);
      }

      // find matching declared switch
      size_t iit{BADIND}, j{0}; // iit will be index in _items
      while (iit == BADIND && j < _items.size()) {
        const item& itemj = _items[j];
        bool cmatch = (_subcommand == itemj.subcommand);
        bool nmatch =
            ((!sname.empty() && sname == itemj.sname) || (!lname.empty() && lname == itemj.lname));
        if (cmatch && nmatch) {
          if (_verbose)
            std::cout << "arg \"" << a << "\" matches to " << itemj.subcommand << "," << itemj.name
                      << "\n";
          iit = j;
        }
        j++;
      } // end iit while loop

      if (iit == BADIND) {
        std::cout << "Error: switch " << a << " does not match any declared switch\n";
        return 1;
      } else {
        item& vitem = _items[iit]; // valid item
        if (vitem.args) {          // if arg expected
          if (i >= (a1.size() - 1)) {
            std::cout << "Error: switch " << a << " does not have an argument\n";
            return 1;
          }
          if (a1[i + 1][0] == '-') {
            std::cout << "Error: switch " << a << " does not have an argument\n";
            return 1;
          }
          if (_verbose)
            std::cout << "switch \"" << a << "\" has arg \"" << a1[i + 1] << "\"\n";
          vitem.values.emplace_back(a1[i + 1]);
          i++;   // extra increment to go past arg
        } else { // does not expect an arg
          vitem.values.emplace_back("FOUND");
        }
      }
    }    // endif arg starts with dash
    i++; // increment to next word in a1 args
  }      // end i while

  for (const auto& i : _items) {
    if (i.required) {
      if (i.values.size() == 0 && i.defaultstr.empty()) {
        std::cout << "Error: required subcommand " << i.subcommand << " name " << i.name
                  << " switch was required but not present\n";
      }
    }
  }

  if (_verbose)
    print();

  autohelp();

  return 0;
}

/*********************************************************/
void ParseCLI::autohelp() const {
  if (!_autohelp)
    return;
  bool need = false;
  for (const auto& ii : _items) {
    if (ii.name == "help" && getCount(ii.subcommand, ii.name) > 0)
      need = true;
  }
  if (!need)
    return;

  std::cout << "\n";
  if (_subcommand.empty()) {
    std::cout << "   " << _command << " [GLOBAL OPTIONS] [SUBCOMMAND] [SUBCOMMAND OPTIONS]"
              << "\n\n";
    std::cout << "   " << _helpstr << "\n\n";
  } else {
    std::cout << _command << " " << _subcommand << " [SUBCOMMAND OPTIONS]"
              << "\n";
    for (const auto& ss : _subs) {
      if (ss.subcommand == _subcommand) {
        std::cout << "\n"
                  << "   " << ss.helpstr << "\n\n";
      }
    }
  }

  for (const auto& ii : _items) {
    if (ii.subcommand == _subcommand) {
      std::cout << "     ";
      if (!ii.sname.empty())
        std::cout << "-" << ii.sname;
      if (!ii.sname.empty() && !ii.lname.empty())
        std::cout << ",";
      if (!ii.sname.empty())
        std::cout << "--" << ii.lname << "  ";
      std::cout << ii.helpstr;
      StringVec notes;
      if(ii.repeated) notes.emplace_back("repeated");
      if(ii.required) notes.emplace_back("required");
      if(!ii.defaultstr.empty()) notes.emplace_back("default="+ii.defaultstr);
      if(notes.size()>0) {
        std::cout << " (";
        for(size_t i=0; i<notes.size(); i++) {
          if(i!=0) std::cout << ",";
          std::cout << notes[i];
        }
        std::cout << ")";
      }
      std::cout << "\n";
    }
  }

  if (_subcommand.empty()) {
    std::cout << "\n   subcommands:\n";
    for (const auto& ss : _subs) {
      if (!ss.subcommand.empty()) {
        std::cout << "   " << ss.subcommand << " - " << ss.helpstr << "\n";
      }
    }
  }
  std::cout << "\n";

  return;
}

/*********************************************************/
size_t ParseCLI::findItem(const std::string& subcommand, const std::string& name) const {
  for (size_t i = 0; i < _items.size(); i++) {
    const item& itemi = _items[i];
    if (subcommand == itemi.subcommand && name == itemi.name)
      return i;
  }
  std::cout << "Error: can't find subcommand=" << subcommand << " and name " << name << "\n";
  return BADIND;
}

/*********************************************************/
int ParseCLI::getCount(const std::string& subcommand, const std::string& name) const {
  size_t i = findItem(subcommand, name);
  if (i == BADIND)
    return 0;
  return int(_items[i].values.size());
}

/*********************************************************/
std::string ParseCLI::getString(const std::string& subcommand, const std::string& name) const {
  size_t i = findItem(subcommand, name);
  if (i == BADIND)
    return std::string();
  const item& itemi = _items[i];
  if (itemi.values.size() == 0) {
    return itemi.defaultstr;
  }
  if(itemi.repeated) {
    return itemi.values[0];
  } else {
    return itemi.values.back();
  }
}

/*********************************************************/
StringVec ParseCLI::getStrings(const std::string& subcommand, const std::string& name) const {
  size_t i = findItem(subcommand, name);
  if (i == BADIND)
    return StringVec();
  return _items[i].values;
}

/*********************************************************/
int ParseCLI::getInt(const std::string& subcommand, const std::string& name) const {
  std::string str = getString(subcommand,name);
  return std::stoi(str);
}

/*********************************************************/
std::vector<int> ParseCLI::getInts(const std::string& subcommand, const std::string& name) const {
  std::vector<int> intv;
  StringVec strs = getStrings(subcommand,name);
  for(const auto& istr : strs) intv.emplace_back(std::stoi(istr));
  return intv;
}

/*********************************************************/
float ParseCLI::getFloat(const std::string& subcommand, const std::string& name) const {
  std::string fstr = getString(subcommand,name);
  return std::stof(fstr);
}

/*********************************************************/
std::vector<float> ParseCLI::getFloats(const std::string& subcommand, const std::string& name) const {
  std::vector<float> floatv;
  StringVec strs = getStrings(subcommand,name);
  for(const auto& istr : strs) floatv.emplace_back(std::stof(istr));
  return floatv;
}
