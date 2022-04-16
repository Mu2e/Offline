#ifndef DbTables_DbSet_hh
#define DbTables_DbSet_hh

// A calibration set as it is stored and maintained in the
// database includes pointers from a purpose, version, extension,
// and groups down to the list of IoVs.  This container
// holds just the list of IoVs which is the result of flattening the
// heirarchical structure.  For faster and simpler lookup and manipulation

#include "Offline/DbTables/inc/DbIoV.hh"
#include <map>
#include <vector>

namespace mu2e {

class DbSet {
 public:
  // IoV extended with its related cid
  class EIoV {
   public:
    EIoV() : _cid(-1) {}
    EIoV(int cid, DbIoV const& iov) : _cid(cid), _iov(iov) {}
    int cid() const { return _cid; }
    DbIoV const& iov() const { return _iov; }

   private:
    int _cid;
    DbIoV _iov;
  };

  typedef std::map<int, std::vector<EIoV>> EIoVMap;

  DbSet() : _nearestMatch(false) {}

  void setNearestMatch(bool b) { _nearestMatch = b; }
  void add(int tid, int cid, DbIoV const& iov);
  // return appropriate EIoV (IoV,cid) for given tid and run,subrun
  EIoV find(int tid, uint32_t run, uint32_t subrun) const;
  EIoVMap const& emap() const { return _emap; }

  void clear();
  void print() const;
  void printShort() const;

 private:
  // for each tid key, hold a vector of (IoV,cid)
  std::map<int, std::vector<EIoV>> _emap;
  // whether or not to accept nearest matches in cid find
  bool _nearestMatch;
};

}  // namespace mu2e

#endif
