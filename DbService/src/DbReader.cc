#include <unistd.h>
#include <iostream>
#include <iomanip>
#include "cetlib_except/exception.h"
#include "DbService/inc/DbReader.hh"
#include "DbService/inc/DbCurl.hh"
#include "DbTables/inc/DbUtil.hh"


using namespace std;

mu2e::DbReader::DbReader():_curl_handle(nullptr),
		_timeout(3600),_totalTime(0),_removeHeader(true),
	        _abortOnFail(true),_useCache(true),_cacheLifetime(0),
					 _verbose(0),_timeVerbose(0),_saveCsv(true) {

  // allocates memory for curl
  curl_global_init(CURL_GLOBAL_ALL);
}

mu2e::DbReader::~DbReader() {
  // free memory
  curl_global_cleanup();
}

int mu2e::DbReader::query(std::string& csv, 
			  const std::string& select, 
			  const std::string& table, 
			  const std::string& where,
			  const std::string& order) {

  int rc;

  // reserve resources, alloc memory
  rc = openHandle();
  if (rc!=0) return rc;

  // issue the http call
  rc = queryCore(csv,select,table,where,order);

  // always close handles after single queries
  // otherwise there can be too many open sockets
  // in grid jobs
  closeHandle();

  return rc;
}

int mu2e::DbReader::query(QueryForm& qf) {
  return query(qf.csv,qf.select,qf.table,qf.where,qf.order);
}

int mu2e::DbReader::multiQuery(std::vector<QueryForm>& qfv) {

  int rc = openHandle();
  if (rc!=0) return rc;

  auto iter = qfv.begin();
  while(rc==0 && iter!=qfv.end()) {
    rc = queryCore(iter->csv,iter->select,iter->table,
		   iter->where,iter->order);
    iter++;
  }

  closeHandle();

  return rc;

}



int mu2e::DbReader::fillTableByCid(DbTable::ptr_t ptr, int cid) {
  std::string csv;
  std::string where="cid:eq:"+std::to_string(cid);
  int rc = query(csv,ptr->query(),ptr->dbname(),where);
  if(rc!=0) return rc;
  ptr->fill(csv,_saveCsv);
  return 0;
}


int mu2e::DbReader::fillValTables(DbValCache& vcache) {
  std::string csv;
  int rc;

  auto start_time = std::chrono::high_resolution_clock::now();

  // do all reads at once, for efficiency
  std::vector<QueryForm> qfv(11);

  // DbTables to fill
  ValTables tables;
  ValCalibrations calibrations;
  ValIovs iovs;
  ValGroups groups;
  ValGroupLists grouplists;
  ValPurposes purposes;
  ValLists lists;
  ValTableLists tablelists;
  ValVersions versions;
  ValExtensions extensions;
  ValExtensionLists extensionlists;

  // load up the queries
  qfv[0].select = tables.query();
  qfv[0].table = tables.dbname();
  qfv[1].select = calibrations.query();
  qfv[1].table = calibrations.dbname();
  qfv[2].select = iovs.query();
  qfv[2].table = iovs.dbname();
  qfv[3].select = groups.query();
  qfv[3].table = groups.dbname();
  qfv[4].select = grouplists.query();
  qfv[4].table = grouplists.dbname();
  qfv[5].select = purposes.query();
  qfv[5].table = purposes.dbname();
  qfv[6].select = lists.query();
  qfv[6].table = lists.dbname();
  qfv[7].select = tablelists.query();
  qfv[7].table = tablelists.dbname();
  qfv[8].select = versions.query();
  qfv[8].table = versions.dbname();
  qfv[9].select = extensions.query();
  qfv[9].table = extensions.dbname();
  qfv[10].select = extensionlists.query();
  qfv[10].table = extensionlists.dbname();

  rc = multiQuery(qfv);
  if(rc!=0) return rc;

  tables.fill(qfv[0].csv,_saveCsv);
  vcache.setValTables(tables);
  
  calibrations.fill(qfv[1].csv,_saveCsv);
  vcache.setValCalibrations(calibrations);

  iovs.fill(qfv[2].csv,_saveCsv);
  vcache.setValIovs(iovs);

  groups.fill(qfv[3].csv,_saveCsv);
  vcache.setValGroups(groups);

  grouplists.fill(qfv[4].csv,_saveCsv);
  vcache.setValGroupLists(grouplists);

  purposes.fill(qfv[5].csv,_saveCsv);
  vcache.setValPurposes(purposes);

  lists.fill(qfv[6].csv,_saveCsv);
  vcache.setValLists(lists);

  tablelists.fill(qfv[7].csv,_saveCsv);
  vcache.setValTableLists(tablelists);

  versions.fill(qfv[8].csv,_saveCsv);
  vcache.setValVersions(versions);

  extensions.fill(qfv[9].csv,_saveCsv);
  vcache.setValExtensions(extensions);

  extensionlists.fill(qfv[10].csv,_saveCsv);
  vcache.setValExtensionLists(extensionlists);

  auto end_time = std::chrono::high_resolution_clock::now();
  _lastTime = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
  _totalTime += _lastTime;

  if(_timeVerbose>0) {
    std::cout<<"DbReader::fillValCache took " <<
      std::setprecision(6) << _lastTime.count()*1.0e-6 <<" s" << std::endl;
  }
  if(_verbose>3) {
    std::cout<<"DbReader::fillValCache results " << std::endl;
    vcache.print();
  }


  return 0;

}

int mu2e::DbReader::openHandle() {

  if(_id.name().empty()) {
      throw cet::exception("DBREADER_DBID NOT_SET") 
	<< "DbReader found the DbId was not set\n";
  }

  _curl_handle = curl_easy_init();

  if (!_curl_handle) {
    if(_abortOnFail) {
      throw cet::exception("DBREADER_FAILED_CURL_INIT") << 
	"DbReader failed to initialize curl" <<"\n";
    } else {
      if(_verbose>0) {
	std::cout << "DbReader failed to initialize curl" << std::endl;
      }
      return 1;
    }
  }

  // send all data to this function
  curl_easy_setopt(_curl_handle, CURLOPT_WRITEFUNCTION, mu2e::curlWriteCallback);
  // we pass our 'response' struct to the callback function
  curl_easy_setopt(_curl_handle, CURLOPT_WRITEDATA, (void *)&_result);
  // some servers don't like requests that are made without a user-agent 
  // field, so we provide one
  curl_easy_setopt(_curl_handle, CURLOPT_USERAGENT, "libcurl-agent/1.0");
  // Enable redirection
  curl_easy_setopt(_curl_handle, CURLOPT_FOLLOWLOCATION, 1);
  // do not require certificate verification
  curl_easy_setopt(_curl_handle, CURLOPT_SSL_VERIFYPEER, 0);

  return 0;

}

int mu2e::DbReader::queryCore(std::string& csv, 
			      const std::string& select, 
			      const std::string& table, 
			      const std::string& where,
			      const std::string& order) {

  std::string url;
  bool qcache = true;
  bool qval = table.substr(0,3)=="val";
  if( qval ) {
    // val tables describe IoV, so need to be up-to-date
    if(_cacheLifetime<=0) qcache = false;
  } else {
    if( !_useCache ) qcache = false;
  }
  if(qcache) {
    // result is cached on the web server
    // calibration table entries don't change
    url = _id.url();
  } else {
    // result will come directly from the db
    url = _id.urlNoCache();
  }
  url.append("t=");
  url.append(table);
  url.append("&c=");
  url.append(select);
  if(!where.empty()) {
    url.append("&w=");
    url.append(where);
  }
  if(!order.empty()) {
    url.append("&o=");
    url.append(order);
  }

  if(_verbose>3) {
    std::string time = DbUtil::timeString();
    std::cout << "DbReader " << time 
	      <<"  url="<<url << std::endl;
  }

  string urlfinal = url;
  if(qval && _cacheLifetime>0) {
    // add a field, rounded to the nearest n second, so that
    // cache has an effective lifetime by changing every n sec
    time_t tt;
    time(&tt);
    ostringstream ss;
    ss << tt/_cacheLifetime;
    urlfinal.append("&z=");
    urlfinal.append(ss.str());
  }
  // put the url in the handle
  curl_easy_setopt(_curl_handle, CURLOPT_URL, urlfinal.c_str());
  // return an error when the http header returns an error (like bad gateway)
  curl_easy_setopt(_curl_handle, CURLOPT_FAILONERROR, 1L);

  auto start_time = std::chrono::high_resolution_clock::now();
  _curl_ret = CURLE_RECV_ERROR;

  // I am forced to use random numbers for time to retry - initialize
  srandom( start_time.time_since_epoch().count()%RAND_MAX );

  int itry = 0;
  int nsec = 0;
  double sleep_time = 0.0; // in seconds
  while(_curl_ret != CURLE_OK && nsec < _timeout) {
    if(itry>0) {
      double st = ((double)random()/(double)RAND_MAX) * (1 << itry);
      sleep( st );
      sleep_time += st;
    }
    _result.reply.clear();
    if(_verbose>5) std::cout << "DbReader start perform with itry="
			     <<itry << std::endl;
    // the actual GET over the network
    _curl_ret = curl_easy_perform(_curl_handle);
    // save error state
    _lastError = curl_easy_strerror(_curl_ret);

    int http_code = 0;
    curl_easy_getinfo (_curl_handle, CURLINFO_RESPONSE_CODE, &http_code);

    if (_curl_ret == CURLE_OK && http_code == 200) {
      if(_verbose>5) {
	std::cout << "DbReader curlOutput=";
	std::size_t iend = _result.reply.size();
	if(_verbose>9 || iend<400) {
	  std::cout << _result.reply << std::endl;
	} else {
	  std::cout << _result.reply.substr(0,200) << std::endl;
	  std::cout << " ... " << std::endl;
	  std::cout << _result.reply.substr(iend-200,200) << std::endl;
	}
	std::cout << "DbReader return code=" << _lastError 
		  << "   http code: " << http_code << std::endl;
      }
    } else {
      if(_verbose>0) std::cout << "DbReader try: " << itry
			       << "   lastError: " << _lastError 
			       << "   http code: " << http_code << std::endl;
    }

    auto try_time = std::chrono::high_resolution_clock::now();
    nsec = ( std::chrono::duration_cast<std::chrono::seconds>(try_time - start_time) ).count();

    itry++;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  _lastTime = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
  _totalTime += _lastTime;

  if(_timeVerbose>3 || _verbose>3) {
    std::string time = DbUtil::timeString();

    std::cout<<"DbReader "<< time << " "
	     << std::setprecision(6) << _lastTime.count()*1.0e-6 <<" s " 
	     << "to read " << table 
	     << " on try " << itry << " with " 
	     << std::setprecision(6) << sleep_time << " s sleeping"
	     << std::endl;
  }

  if (_curl_ret != CURLE_OK) {
    if (_abortOnFail) {
      throw cet::exception("DBREADER_FAILED_CURL") << 
	"DbReader failed to read database " << _id.name() 
		   << ", last error: "<<_lastError <<"\n";
    }
    return 1;
  }



  // move the pointer to string memory from the result to user's csv
  csv = std::move(_result.reply);
  _result.reply.clear();

  // queries return with a line of column titles, remove it
  if(_removeHeader) {
    auto firstN = csv.find("\n");
    if(firstN != std::string::npos) {
      csv = csv.substr(firstN+1, csv.size()-firstN-1);
    }
  }

  return 0;
}

int mu2e::DbReader::closeHandle() {

  // close the socket and release other dynamic resources
  curl_easy_cleanup(_curl_handle);
  // cleanup frees the memory held by the handle
  _curl_handle = nullptr;

  return 0;
}
