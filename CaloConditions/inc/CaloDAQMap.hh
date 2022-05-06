#ifndef CaloConditions_CaloDAQMap_hh
#define CaloConditions_CaloDAQMap_hh

//
// CaloDAQMap stores channel maps vs detector map
//

// C++ includes
#include <vector>

// Mu2e includes
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  class CaloDAQMap : virtual public ProditionsEntity {
  public:

    typedef std::shared_ptr<CaloDAQMap> ptr_t;
    typedef std::shared_ptr<const CaloDAQMap> cptr_t;

    //CaloDAQMap():_name("CaloDAQMap") {}
    constexpr static const char* cxname = {"CaloDAQMap"};

    // construct with constants, then some values are computed and filled below
    CaloDAQMap(std::vector<uint16_t> DIRAC2CaloMap, std::vector<uint16_t> Calo2DIRACMap) :
       ProditionsEntity(cxname),
       // _name("CaloDAQMap"),
      _DIRAC2CaloMap(DIRAC2CaloMap), _Calo2DIRACMap(Calo2DIRACMap){}

    virtual ~CaloDAQMap() {}

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

