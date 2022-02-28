
#include "Offline/DbTables/inc/DbCache.hh"
#include <iostream>
#include <iomanip>

void mu2e::DbCache::add(int cid, mu2e::DbTable::cptr_t const& ptr) { 

  _tables[cid] = ptr;

  _cidFifo.push(cid);
  _nAdded++;
  _size += ptr->size();
  if(_size>_hwmSize) _hwmSize = _size;
  if(_nAdded%_purgeInterval==0) purge();

}

mu2e::DbTable::cptr_t mu2e::DbCache::get(int cid) {
  auto it = _tables.find(cid);
  if(it != _tables.end()) {
    return it->second;
  } else {
    return mu2e::DbTable::cptr_t(nullptr);
  }
}

void mu2e::DbCache::purge() {

  if(_size < _limitSize) return;

  _nPurges++;
  while(_size>_purgeEnd*_limitSize) {
    int cid = _cidFifo.front();
    _cidFifo.pop();
    auto it = _tables.find(cid);
    if(it != _tables.end()) {
      _size =- it->second->size();
      _nPurged++;
      _tables.erase(it);
    }
  }

  return;
}


void mu2e::DbCache::print() const {
  for(auto t: _tables) {
    std::cout << std::setw(6) << t.first << " " 
	      << std::setw(15) << t.second->name() << std::endl;
  }
}

void mu2e::DbCache::printStats() const {

  std::cout << "    cache nTable : " << _nAdded << "\n";
  std::cout << "    cache HWM    : " << _hwmSize << " b\n";
  std::cout << "    cache nPurges : " << _nPurges << "\n";
  std::cout << "    cache nPurged : " << _nPurged << "\n";

}
