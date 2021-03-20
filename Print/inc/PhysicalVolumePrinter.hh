//
//  Utility class to print PhysicalVolume
// 
#ifndef Print_inc_PhysicalVolumePrinter_hh
#define Print_inc_PhysicalVolumePrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class PhysicalVolumePrinter : public ProductPrinter {
  public:

    PhysicalVolumePrinter() { }
    PhysicalVolumePrinter(const Config& conf):ProductPrinter(conf) { }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override {};
    void PrintSubRun(art::SubRun const& subrun,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<PhysicalVolumeInfoMultiCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<PhysicalVolumeInfoMultiCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const PhysicalVolumeInfoMultiCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const PhysicalVolumeInfoSingleStage& obj,
	       int ind = -1, std::ostream& os = std::cout);
    void Print(const std::pair<cet::map_vector_key, mu2e::PhysicalVolumeInfo>& obj, 
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);
    void PrintPVListHeader(std::ostream& os = std::cout);

  };

}
#endif
