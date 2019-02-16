#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include "cetlib_except/exception.h"
#include "DbTables/inc/DbIoV.hh"

void mu2e::DbIoV::subtract(DbIoV const& iov, uint32_t run, uint32_t subrun) {

  // check for no overlap
  if(iov.endRun()<startRun() || 
     (iov.endRun()==startRun() && iov.endSubrun()<startSubrun()) ) return;
  if(iov.startRun()>endRun() || 
     (iov.startRun()==endRun() && iov.startSubrun()>endSubrun()) ) return;

  //check if iov start/end are inside this interval
  bool sin = inInterval(iov.startRun(),iov.startSubrun());
  bool ein = inInterval(iov.endRun(),iov.endSubrun());

  int method = 0;
  if(! sin && ! ein) { // subtract all
    method = 0;
  } else if(sin && ! ein) { //subtract middle to end
    method = 1;
  } else if(ein && ! sin) { // subtract start to middle
    method = 2;
  } else { 
    // subtract a piece in the middle, keep the early or 
    // late fragment based on the run/subrun
    method = 1;
    if(run > iov.startRun() || 
       (run==iov.startRun() && subrun>=iov.startSubrun() ) ) method = 2;
  }

  if(method==0) { // subtract the whole interval
    set(0,0,0,0);
  } else if(method==1) { // subtract middle to end
    uint32_t er = iov.startRun();
    uint32_t es = iov.startSubrun();
    if(es==0) {
      er--;
      es = maxsr;
    } else {
      es--;
    }
    set(startRun(),startSubrun(),er,es);
    return;
  } else if(method==2) { // subtract start to middle
    uint32_t sr = iov.endRun();
    uint32_t ss = iov.endSubrun();
    if(ss==maxsr) {
      sr++;
      ss = 0;
    } else {
      ss++;
    }
    set(sr,ss,endRun(),endSubrun());
  }

  return;

}


void mu2e::DbIoV::overlap(DbIoV const& iov) {

  // check for no overlap
  bool null = false;
  if(iov.endRun()<startRun() || 
     (iov.endRun()==startRun() && iov.endSubrun()<startSubrun()) ) null = true;
  if(iov.startRun()>endRun() || 
     (iov.startRun()==endRun() && iov.startSubrun()>endSubrun()) ) null = true;
  if(null) {
    set(0,0,0,0);
    return;
  }

  //check if iov start/end are inside this interval
  bool sin = inInterval(iov.startRun(),iov.startSubrun());
  bool ein = inInterval(iov.endRun(),iov.endSubrun());

  if(sin && ein) { // overlap is iov
    set(iov.startRun(),iov.startSubrun(),iov.endRun(),iov.endSubrun());
  } else if(sin && ! ein) { //overlap middle to end
    set(iov.startRun(),iov.startSubrun(),endRun(),endSubrun());
  } else if(ein && ! sin) { // overlap start to middle
    set(startRun(),startSubrun(),iov.endRun(),iov.endSubrun());
  }
  // last case (!sin && !ein) is no change
  return;
}


std::string mu2e::DbIoV::simpleString() const {
  std::ostringstream ss;
  ss << startRun() << " " << startSubrun()  << " "
     << endRun() << " " << endSubrun();
  return ss.str();
}

std::string mu2e::DbIoV::to_string(bool compress) const {
  std::ostringstream ss;
  int w = ( compress ? 0 : 6 );
  ss << std::setw(w) << startRun() << ":";
  ss << std::setw(w) << startSubrun() << "-";
  ss << std::setw(w) << endRun() << ":";
  ss << std::setw(w) << endSubrun();
  return ss.str();
}

void mu2e::DbIoV::setByString(std::string iovstr) {
  boost::trim(iovstr); // remove leading, trailing whitespace
  boost::to_upper(iovstr);
  if(iovstr.empty()||iovstr=="MAX"||iovstr=="ALL") {
    setMax();
    return;
  }
  if(iovstr=="EMPTY") {
    set(0,0,0,0);
    return;
  }
  std::vector<std::string> words;
  std::string start,end;
  if(iovstr.find('-')!=std::string::npos) { // has a dash
    boost::split(words,iovstr, boost::is_any_of("-"), 
		 boost::token_compress_off);
    start = words[0];
    end = words[1];
  } else {
    // this allows strings with spaces like
    // "1 2 3 4" -> start = 1:2, end=3:4
    boost::split(words,iovstr, boost::is_any_of(" \t"), 
		 boost::token_compress_on);
    if(words.size()==4) {
      start = words[0]+":"+words[1];
      end = words[2]+":"+words[3];
    } else if(words.size()==2) {
      start = words[0];
      end = words[1];
    } else if(words.size()==1) {
      start = words[0];
      end = words[0];
    } else {
      throw cet::exception("DBIOV_BAD_INIT_STRING") 
	<< "DbIoV::setByString cannot interpret string: " << iovstr << "\n";
    }
  }

  // require a string for start and end
  if(start.empty()||end.empty()) {
    throw cet::exception("DBIOV_EMPTY_INIT_STRING") 
      << "DbIoV::setByString found start or end point was blank: " << iovstr << "\n";
  }

  boost::split(words,start, boost::is_any_of(":"),boost::token_compress_on);
  uint32_t startr = 0, startsr = 0;
  if(words.size()>=1) {
    boost::to_upper(words[0]);
    if(words[0]!="MIN") {
      startr  = std::stoi(words[0]);
    }
  }
  if(words.size()>=2) {
    boost::to_upper(words[1]);
    if(words[1]!="MIN") {
      startsr = std::stoi(words[1]);
    }
  }
  boost::split(words,end, boost::is_any_of(":"),boost::token_compress_on);
  uint32_t endr = maxr, endsr = maxsr;
  if(words.size()>=1) {
    boost::to_upper(words[0]);
    if(words[0]!="MAX") {
      endr  = std::stoi(words[0]);
    }
  }
  if(words.size()>=2) {
    boost::to_upper(words[1]);
    if(words[1]!="MAX") {
      endsr = std::stoi(words[1]);
    }
  }

  if( endr<startr || (endr==startr && endsr<startsr) ) {
    throw cet::exception("DBIOV_NEGATIVE_INIT_STRING") 
      << "DbIoV::setByString found end point was before start point: " << iovstr << "\n";
  }

  set(startr,startsr,endr,endsr);

}
