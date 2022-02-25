
#include "Offline/DbTables/inc/DbSet.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <iostream>

namespace mu2e {

  //***********************************************

void DbSet::add(int tid, int cid, DbIoV const& iov) { 
  _emap[tid].emplace_back(cid, iov); 
}

  //***********************************************

DbSet::EIoV DbSet::find(int tid, uint32_t run, uint32_t subrun) const {

  auto iter = _emap.find(tid);
  if (iter == _emap.end()) {
    return EIoV();
  }
  
  for (auto const& ee : iter->second) {
    if (ee.iov().inInterval(run, subrun)) {
      return ee;
    }
  }

  // if no return above, then find a nearby entry, if requested
  if(_nearestMatch) {
    EIoV br;
    int brr=br.iov().maxRun(), bsr=br.iov().maxSubrun();
    for(auto const& r : iter->second) { // find closest entry
      int dr = run - r.iov().endRun();
      int ds = subrun + r.iov().maxSubrun() - r.iov().endSubrun();
      if(dr==0) ds = subrun - r.iov().endSubrun();
      if( dr>0 && ds>0 && (dr<brr || ( dr==brr && ds<bsr) ) ) {
        brr = dr;
        bsr = ds;
        br = r;
      }
    } // loop over Rows for table type
    if(br.cid()>0) { // then something was found
      return br;
    }
  } // try nearest match

  return EIoV();

}

  //***********************************************

void DbSet::clear() {
  _emap.clear();
}

  //***********************************************

void DbSet::print() const {
  std::cout << "  tid       valid range        cid" << std::endl;
  for (auto const& p : _emap) {
    int tid = p.first;
    for (auto r : p.second) {
      std::cout << std::setw(5) << tid << std::setw(20) << r.iov() << std::setw(8) << std::right
                << r.cid() << std::endl;
    }
  }
}

  //***********************************************

void DbSet::printShort() const {
  std::cout << "  tid     N IoV" << std::endl;
  for (auto const& p : _emap) {
    std::cout << std::setw(5) << p.first 
              << std::setw(8) << p.second.size() << std::endl;
  }
}


} // namespace mu2e
