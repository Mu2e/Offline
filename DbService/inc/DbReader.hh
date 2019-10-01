#ifndef DbService_DbReader_hh
#define DbService_DbReader_hh
//
// use curl to get data from the condtions database web interface
// 
// To use:
// create with default (production database) or specify another database,
// then call query method.  Everything goes to the web cache url except
// val.* IOV tables
//
// this code will retry up to the timeout, then abort
// if you want to handle the failure, set setAbortOnFail(false)
//
#include <string>
#include <chrono>
#include <curl/curl.h>
#include "DbTables/inc/DbId.hh"
#include "DbTables/inc/DbValCache.hh"


namespace mu2e {
  class DbReader {
  public:

    struct QueryForm {
      std::string csv; // the answer
      std::string select;
      std::string table;
      std::string where;
      std::string order;
    };

    DbReader(const DbId& id = DbId());
    ~DbReader();
    const DbId& id() const { return _id; }

    // run query url, answer returned in csv
    int query(std::string& csv, const std::string& select, 
	      const std::string& table, const std::string& where="",
	      const std::string& order="");
    // an alternate form
    int query(QueryForm& qf);

    // this keeps the socket open between url's, so it is more efficient
    int multiQuery(std::vector<QueryForm>& qfv);

    int fillTableByCid(DbTable::ptr_t ptr, int cid);
    int fillValTables(DbValCache& vcache);

    std::string& lastError() { return _lastError; }
    double lastTime() { return _lastTime.count()*1.0e-6; } // seconds
    double totalTime() { return _totalTime.count()*1.0e-6; } // seconds

    void setDbId(DbId id) { _id = id; }
    // time to keep retrying to read data from the web server
    void setTimeout(float timeout=3600) { _timeout = timeout; }
    // stop on errors, set false to handle errors in caller
    void setAbortOnFail(bool aof=true) { _abortOnFail = aof; }
    // use cached queries where possible (non-val tables)
    void setUseCache(bool uc=true) { _useCache = uc; }
    // if 0, skip cache, go to DB. If non-zero, use cache, but
    // renew cache every lifetime (integer seconds)
    void setCacheLifetime(int clt=0) { _cacheLifetime = clt; }
    void setVerbose(int verbose) { _verbose = verbose; }
    void setTimeVerbose(int timeVerbose) { _timeVerbose = timeVerbose; }

  private:

    // for internal use with curl
    struct DbQueryResult {
      std::string reply;
    };

    // these pieces are put together in the query routines
    int openHandle();
    int closeHandle();
    int queryCore(std::string& csv, const std::string& select, 
	      const std::string& table, const std::string& where="",
	      const std::string& order="");

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
    int _cacheLifetime;
    int _verbose;
    int _timeVerbose;
  };
}
#endif
