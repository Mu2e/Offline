#include "Offline/DbService/inc/OpenHandle.hh"

using namespace mu2e;

//**************************************************
template <typename T>
OpenHandle<T>::OpenHandle(const std::string& database) :
  _tool(database), _name(T::cxname), _cid(0) {}

//**************************************************
template <typename T>
int OpenHandle<T>::update(uint32_t run, uint32_t subrun) {
  int rc = 0;
  if (_iov.inInterval(run, subrun)) return rc;

  std::string csv;
  _tool.table(_name, run, subrun, csv, _cid, _iov);
  _table = std::make_shared<T>();
  _table->fill(csv);

  return 0;
}

//**************************************************
template <typename T>
const T& OpenHandle<T>::table() {
  return *_table;
}
