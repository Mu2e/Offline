#ifndef DbTables_DbId_hh
#define DbTables_DbId_hh

#include <string>

namespace mu2e {

  class DbId {
  public:
    DbId() {}
    DbId(std::string const& name, std::string const& host, std::string const& port,
	 std::string const& url, std::string const& urlNoCache):
      _name(name),_host(host),_port(port),_url(url),_urlNoCache(urlNoCache) {}
    const std::string& name() const { return _name; }
    const std::string& host() const { return _host; }
    const std::string& port() const { return _port; }
    const std::string& url() const { return _url; }
    const std::string& urlNoCache() const { return _urlNoCache; }
  private:
    std::string _name;
    std::string _host;
    std::string _port;
    std::string _url;
    std::string _urlNoCache;
  };
}
#endif
