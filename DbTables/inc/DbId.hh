#ifndef DbTables_DbId_hh
#define DbTables_DbId_hh

#include <string>

namespace mu2e {

  class DbId {
  public:
    DbId() { setDb("mu2e_conditions_prd"); }
    DbId(std::string const& name) { setDb(name.c_str()); }
    DbId(const char* name) { setDb(name); }
    const std::string& name() const { return _name; }
    const std::string& host() const { return _host; }
    const std::string& port() const { return _port; }
    const std::string& url() const { return _url; }
    const std::string& urlNoCache() const { return _urlNoCache; }
    void setDb(std::string const& name);
  private:
    std::string _name;
    std::string _host;
    std::string _port;
    std::string _url;
    std::string _urlNoCache;
  };
}
#endif
