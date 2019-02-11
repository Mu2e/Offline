#include <unistd.h>
#include <iostream>
#include <iomanip>
#include "cetlib_except/exception.h"
#include "DbService/inc/DbReader.hh"
#include "DbService/inc/DbCurl.hh"

mu2e::DbReader::DbReader(const DbId& id):_id(id),_maxTries(10),
		_totalTime(0),_removeHeader(true),_abortOnFail(true),
		_useCache(true),_verbose(0),_timeVerbose(0) {

  curl_global_init(CURL_GLOBAL_ALL);
  _curl_handle = curl_easy_init();

  if (!_curl_handle) {
    throw cet::exception("DBREADER_FAILED_CURL_INIT") << 
      "DbReader failed to initialize curl" <<"\n";
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
}

mu2e::DbReader::~DbReader() {
    curl_global_cleanup();
}

int mu2e::DbReader::query(std::string& csv, 
			  const std::string& select, 
			  const std::string& table, 
			  const std::string& where,
			  const std::string& order) {

  // result is cached on the web server
  // calibration table entries don't change
  std::string url = _id.url();
  if( !_useCache || table.substr(0,3)=="val" ) {
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

  //http://dbdata0vm.fnal.gov:9091/QE/mu2e/dev/app/SQ/query?t=testt&c=*&dbname=mu2e_conditions_dev


  if(_verbose>3) std::cout << "DbReader url="<<url << std::endl;
  curl_easy_setopt(_curl_handle, CURLOPT_URL, url.c_str());
  // return an error when the http header returns an error (like bad gateway)
  curl_easy_setopt(_curl_handle, CURLOPT_FAILONERROR, 1L);

  uint itry = 0;
  auto start_time = std::chrono::high_resolution_clock::now();
  _curl_ret = CURLE_RECV_ERROR;

  while(_curl_ret != CURLE_OK && itry < _maxTries) {
    if(itry>0) sleep(itry);
    _result.reply.clear();
    if(_verbose>5) std::cout << "DbReader start perform with itry="<<itry << std::endl;
    _curl_ret = curl_easy_perform(_curl_handle);
    _lastError = curl_easy_strerror(_curl_ret);
    if (_curl_ret == CURLE_OK) {
      if(_verbose>5) {
	std::cout << "DbReader curlOutput=" << _result.reply << std::endl;
	std::cout << "DbReader return code=" << _lastError << std::endl;
      }
    } else {
      if(_verbose>0) std::cout << "DbReader lastError=" << _lastError << std::endl;
    }
    //curl_easy_cleanup(_curl_handle);
    itry++;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  _lastTime = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
  _totalTime += _lastTime;

  if(_timeVerbose>3) {
    std::cout<<"DbReader::query took " <<
      std::setprecision(6) << _lastTime.count()*1.0e-6 <<" s " 
	     << "to read " << table << std::endl;
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

int mu2e::DbReader::fillTableByCid(DbTable::ptr_t ptr, int cid) {
  std::string csv;
  std::string where="cid:eq:"+std::to_string(cid);
  int rc = query(csv,ptr->query(),ptr->dbname(),where);
  if(rc!=0) return rc;
  ptr->fill(csv);
  return 0;
}


int mu2e::DbReader::fillValTables(DbValCache& vcache) {
  std::string csv;
  int rc;

  auto start_time = std::chrono::high_resolution_clock::now();

  ValTables tables;
  rc = query(csv,tables.query(),tables.dbname());
  if(rc!=0) return rc;
  tables.fill(csv);
  vcache.setValTables(tables);
  
  ValCalibrations calibrations;
  rc = query(csv,calibrations.query(),calibrations.dbname());
  if(rc!=0) return rc;
  calibrations.fill(csv);
  vcache.setValCalibrations(calibrations);

  ValIovs iovs;
  rc = query(csv,iovs.query(),iovs.dbname());
  if(rc!=0) return rc;
  iovs.fill(csv);
  vcache.setValIovs(iovs);

  ValGroups groups;
  rc = query(csv,groups.query(),groups.dbname());
  if(rc!=0) return rc;
  groups.fill(csv);
  vcache.setValGroups(groups);

  ValGroupLists grouplists;
  rc = query(csv,grouplists.query(),grouplists.dbname());
  if(rc!=0) return rc;
  grouplists.fill(csv);
  vcache.setValGroupLists(grouplists);

  ValPurposes purposes;
  rc = query(csv,purposes.query(),purposes.dbname());
  if(rc!=0) return rc;
  purposes.fill(csv);
  vcache.setValPurposes(purposes);

  ValLists lists;
  rc = query(csv,lists.query(),lists.dbname());
  if(rc!=0) return rc;
  lists.fill(csv);
  vcache.setValLists(lists);

  ValTableLists tablelists;
  rc = query(csv,tablelists.query(),tablelists.dbname());
  if(rc!=0) return rc;
  tablelists.fill(csv);
  vcache.setValTableLists(tablelists);

  ValVersions versions;
  rc = query(csv,versions.query(),versions.dbname());
  if(rc!=0) return rc;
  versions.fill(csv);
  vcache.setValVersions(versions);

  ValExtensions extensions;
  rc = query(csv,extensions.query(),extensions.dbname());
  if(rc!=0) return rc;
  extensions.fill(csv);
  vcache.setValExtensions(extensions);

  ValExtensionLists extensionlists;
  rc = query(csv,extensionlists.query(),extensionlists.dbname());
  if(rc!=0) return rc;
  extensionlists.fill(csv);
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

