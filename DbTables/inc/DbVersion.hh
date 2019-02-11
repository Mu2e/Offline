#ifndef DbTables_DbVersion_hh
#define DbTables_DbVersion_hh

#include <string>
#include <boost/algorithm/string.hpp>

namespace mu2e {

  class DbVersion {
  public:
    DbVersion():_purpose(std::string("none")),_major(-1),
	      _minor(-1),_extension(-1) {}
    DbVersion(std::string const& purpose, int major, 
	    int minor, int extension):_purpose(purpose),
	     _major(major),_minor(minor),_extension(extension){}
    DbVersion(std::string const& purpose, std::string version):
    _purpose(purpose),_major(-1),_minor(-1),_extension(-1){
      if( version[0]=='v' || version[0]=='V' ) 
	version = version.substr(1,version.size()-1);
      std::vector<std::string> cont;
      boost::split(cont, version, boost::is_any_of("_-/: "));
      for(auto& s: cont) boost::algorithm::trim(s);

      if(cont.size()>0) {
	if(cont[0]=="*" || cont[0]=="") {
	  _major = -1;
	} else {
	  _major = std::stoi(cont[0]);
	}
      }
      if(cont.size()>1) {
	if(cont[1]=="*" || cont[1]=="") {
	  _minor = -1;
	} else {
	  _minor = std::stoi(cont[1]);
	}
      }
      if(cont.size()>2) {
	if(cont[2]=="*" || cont[2]=="") {
	  _extension = -1;
	} else {
	  _extension = std::stoi(cont[2]);
	}
      }
    }
    std::string const& purpose() const {return _purpose;}
    int major() const { return _major;}
    int minor() const { return _minor;}
    int extension() const { return _extension;}
    std::string to_string(std::string const& sep = std::string("/")) const {
      return _purpose+" "
	+std::to_string(_major)+sep
	+std::to_string(_minor)+sep
	+std::to_string(_extension);}
  private:
    std::string _purpose;
    int _major;
    int _minor;
    int _extension;
  };
}
#endif
