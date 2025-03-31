#ifndef CosmicRayShieldGeom_CRSScintillatorShield_hh
#define CosmicRayShieldGeom_CRSScintillatorShield_hh
//
// Representation of one ScintillatorShield in CosmicRayShield.
//
//
// Original author KLG based on Rob Kutschke's Device
//

// C++ includes
#include <vector>

// Mu2e includes
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorShieldId.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"


namespace mu2e
{

  class CRSScintillatorShield
  {

    friend class CosmicRayShieldMaker;

    //disable default constructor
    CRSScintillatorShield();

    public:

    CRSScintillatorShield(CRSScintillatorShieldId const & id,
                          std::string const & name,
                          const std::shared_ptr<CRSScintillatorBarDetail> barDetails,
                          const std::string &absorberMaterialName, const std::string &aluminumSheetMaterialName, const std::string &FEBMaterialName,
                          CRSScintillatorShieldId precedingSector, int sectorType, int countersPerModule);

    // Accept the compiler generated destructor, copy constructor and assignment operators

    // Accessors
    const CRSScintillatorShieldId id() const { return _id;}

    const std::string& getName() const {return _name;}

    int nModules() const {return int(_modules.size());}

    const std::vector<CRSScintillatorModule>& getCRSScintillatorModules() const
    {
      return _modules;
    }

    const CRSScintillatorModule& getModule( int n) const
    {
      return _modules.at(n);
    }

    const CRSScintillatorModule& getModule( const CRSScintillatorModuleId& moduleid ) const
    {
      return _modules.at(moduleid.getModuleNumber());
    }

    const CRSScintillatorLayer& getLayer( const CRSScintillatorLayerId& lid ) const
    {
      return _modules.at(lid.getModuleNumber()).getLayer(lid);
    }

    CRSScintillatorBar const & getBar( const CRSScintillatorBarId& bid ) const
    {
      return _modules.at(bid.getModuleNumber()).getBar(bid);
    }

    const CRSScintillatorBarDetail& getCRSScintillatorBarDetail() const
    {
      return *_barDetails;
    }

    const std::string &getAbsorberMaterialName() const
    {
      return _absorberMaterialName;
    }

    const std::string &getAluminumSheetMaterialName() const
    {
      return _absorberMaterialName;
    }

    const std::string &getFEBMaterialName() const
    {
      return _FEBMaterialName;
    }

    // Formatted string embedding the id of the shield.
    std::string name( std::string const & base ) const;

    // On readback from persistency, recursively recompute mutable members.
    //    void fillPointers( const CosmicRayShield& cosmicRayShield ) const;

    private:

    CRSScintillatorShieldId _id;

    std::string _name;

    std::vector<CRSScintillatorModule> _modules;

    // Detailed info about scintillator bars assuming they are the same for all bars within a shield
    std::shared_ptr<CRSScintillatorBarDetail> _barDetails;

    std::string _absorberMaterialName;
    std::string _aluminumSheetMaterialName;
    std::string _FEBMaterialName;


    // Information needed for the coincidence finder
    public:
    const CRSScintillatorShieldId getPrecedingSector() const {return _precedingSector;}
    const int getSectorType() const {return _sectorType;}
    const int getCountersPerModule() const {return _countersPerModule;}

    private:
    CRSScintillatorShieldId _precedingSector; //the sector id which precedes this CRV sector
    int _sectorType; //e.g. R=1, L=2, T=3, ...
    int _countersPerModule;  //per layer

  };

} //namespace mu2e

#endif /* CosmicRayShieldGeom_CRSScintillatorShield_hh */
