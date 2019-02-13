#ifndef DbService_DbCurl_hh
#define DbService_DbCurl_hh

//
// the curl library requires the response to be put in a function
// with a particular signature and behavior.
// This is used by DbTables/inc/DbReader.hh for curl
//

namespace mu2e {
  static size_t curlWriteCallback(void *contents, size_t size, size_t nmemb, void *userp)  {
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
  }
}

#endif
