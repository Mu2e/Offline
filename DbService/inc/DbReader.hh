#ifndef DbService_DbReader_hh
#define DbService_DbReader_hh

#include <string>
#include <chrono>
#include <curl/curl.h>
#include "DbTables/inc/DbId.hh"
#include "DbTables/inc/DbValCache.hh"


namespace mu2e {
  class DbReader {
  public:
    DbReader(const DbId& id = DbId());
    ~DbReader();
    const DbId& id() const { return _id; }

    struct DbQueryResult {
      std::string reply;
    };

    int query(std::string& csv, const std::string& select, 
	      const std::string& table, const std::string& where="",
	      const std::string& order="");

    int fillTableByCid(DbTable::ptr_t ptr, int cid);
    int fillValTables(DbValCache& vcache);

    std::string& lastError() { return _lastError; }
    double lastTime() { return _lastTime.count()*1.0e-6; } // seconds
    double totalTime() { return _totalTime.count()*1.0e-6; } // seconds

    void setDbId(DbId id) { _id = id; }
    void setTimeout(float timeout=3600) { _timeout = timeout; }
    void setAbortOnFail(bool aof=true) { _abortOnFail = aof; }
    void setUseCache(bool uc=true) { _useCache = uc; }
    void setVerbose(int verbose) { _verbose=verbose; }
    void setTimeVerbose(int timeVerbose) { _timeVerbose=timeVerbose; }

  private:

    DbId _id;
    CURL *_curl_handle;
    CURLcode _curl_ret;
    DbQueryResult _result;
    float _timeout; // in s
    std::chrono::microseconds _lastTime;
    std::chrono::microseconds _totalTime;
    std::string _lastError;
    bool _removeHeader;
    bool _abortOnFail;
    bool _useCache; // default = true, false = do not use web cache
    int _verbose;
    int _timeVerbose;
  };
}
#endif
