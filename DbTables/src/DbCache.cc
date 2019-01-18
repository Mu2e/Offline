#include <algorithm>
#include "DbTables/inc/DbCache.hh"

void mu2e::DbCache::add(int cid, mu2e::DbTable::cptr_t const& ptr) { 
  _tables[cid] = ptr;
}

mu2e::DbTable::cptr_t mu2e::DbCache::get(int cid) {
  auto it = _tables.find(cid);
  if(it != _tables.end()) {
    return it->second;
  } else {
    return mu2e::DbTable::cptr_t(nullptr);
  }
}

int mu2e::DbCache::purge(const size_t target) {
  // loop over tables, remove oldest
  size_t current = size();
  if(current<1.1*target) return 0;

  return 0;
}

size_t mu2e::DbCache::size() {
  size_t s=0;
  for(auto const& t : _tables) s += t.second->size();
  return s;
}


void mu2e::DbCache::print() const {
  for(auto t: _tables) {
    std::cout << std::setw(6) << t.first << " " 
	      << std::setw(15) << t.second->name() << std::endl;
  }
}
