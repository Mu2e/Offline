//
//  Define the CRV test module geometry for KinKal
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKalGeom_CRV_hh
#define KinKalGeom_CRV_hh
#include "KinKal/Geometry/Rectangle.hh"
#include <vector>
#include <memory>
namespace mu2e {
  namespace KKGeom {
    class CRV {
    public:
      using RecPtr = std::shared_ptr<KinKal::Rectangle>;
      using RecPtrVec = std::vector<RecPtr>;
      // default constructor (empty)
      CRV();
      // construct with rectangles representing the midplane of CRV sectors (shields)
      CRV(RecPtrVec const& sectors,std::vector<std::string> snames) :
        sectors_(sectors),snames_(snames) {}
      // accessors
      auto const& sectors() const { return sectors_; }
      auto const& sector(size_t isector) const { return sectors_.at(isector); }
      // access sectors can by name: this might be inefficient. return value must
      // be tested against the sector vector (could be null)
      RecPtrVec::const_iterator sector(std::string const& sector_name) const {
        auto retval = sectors_.end();
        auto ifnd = std::find(snames_.begin(),snames_.end(),sector_name);
        if(ifnd != snames_.end()){
          auto index = std::distance(snames_.begin(),ifnd);
          retval = sectors_.begin() + index;
        }
        return retval;
      }

      private:
      RecPtrVec sectors_;
      std::vector<std::string> snames_;
    };
  }
}

#endif
