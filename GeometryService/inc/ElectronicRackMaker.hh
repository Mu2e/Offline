#ifndef GeometryServices_ElectronicRackMaker_hh
#define GeometryServices_ElectronicRackMaker_hh
//
// Class to construct and return ElectronicRack
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ElectronicRack;
  class SimpleConfig;

  class ElectronicRackMaker {
  public:

    static std::unique_ptr<ElectronicRack>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* GeometryServices_ElectronicRackMaker_hh */
