#ifndef CaloConditions_CaloDAQConditions_hh
#define CaloConditions_CaloDAQConditions_hh

//
// CaloDAQConditions stores channel maps vs detector map
//

// C++ includes
#include <vector>

// Mu2e includes
#include "CalorimeterGeom/inc/CaloInfo.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  class CaloDAQConditions : virtual public ProditionsEntity {
  public:

    typedef std::shared_ptr<CaloDAQConditions> ptr_t;
    typedef std::shared_ptr<const CaloDAQConditions> cptr_t;

    //CaloDAQConditions():_name("CaloDAQConditions") {}
    constexpr static const char* cxname = {"CaloDAQConditions"};
    
    // construct with constants, then some values are computed and filled below
    CaloDAQConditions(std::vector<uint16_t> DIRAC2CaloMap, std::vector<uint16_t> Calo2DIRACMap) :
       ProditionsEntity(cxname),
       // _name("CaloDAQConditions"),
      _DIRAC2CaloMap(DIRAC2CaloMap), _Calo2DIRACMap(Calo2DIRACMap){}

    virtual ~CaloDAQConditions() {}
    
    //ora ..    std::string const& name() const { return _name; }
    void print(std::ostream& os) const;

    uint16_t packetIdTocaloRoId(uint16_t packetId) const;
    uint16_t caloRoIdToPacketId(uint16_t caloRoId) const;
    //
    // all of these must be called to fill this object ...
    //
    // From Dirac POinter to Calo Offline roID
    void setDIRAC2CaloMap(std::vector<uint16_t> DIRAC2CaloMap) { _DIRAC2CaloMap = DIRAC2CaloMap; }
    // From Calo Offline roID to Dirac Pointer
    void setCalo2DIRACMap(std::vector<uint16_t> Calo2DIRACMap) { _Calo2DIRACMap = Calo2DIRACMap; }

  private:
    
    //ora ... std::string _name;

    std::vector<uint16_t> _DIRAC2CaloMap;
    std::vector<uint16_t> _Calo2DIRACMap;

  };
  
}

#endif

