
#include "Print/inc/PhysicalVolumePrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::PhysicalVolumePrinter::PrintSubRun(art::SubRun const& subrun,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<PhysicalVolumeInfoMultiCollection> > vah = subrun.getMany<PhysicalVolumeInfoMultiCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = subrun.getValidHandle<PhysicalVolumeInfoMultiCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::PhysicalVolumePrinter::Print(
	     const art::Handle<PhysicalVolumeInfoMultiCollection>& handle,
	     std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::PhysicalVolumePrinter::Print(
             const art::ValidHandle<PhysicalVolumeInfoMultiCollection>& handle,
	     std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle,os);
}

void 
mu2e::PhysicalVolumePrinter::Print(
  const PhysicalVolumeInfoMultiCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "PhysicalVolumeInfoMultiCollection has " << coll.size() << " stages\n";
  if(verbose()==1) PrintListHeader(os);
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++, os);
}

void 
mu2e::PhysicalVolumePrinter::Print(
	     const PhysicalVolumeInfoSingleStage& obj, 
	     int ind, std::ostream& os) {
  if(verbose()<1) return;
  const PhysicalVolumeInfoSingleStage& ss = obj;
  int n = ss.size();
  if(verbose()>=1) {
    if(verbose()>1) PrintListHeader(os);
    std::cout << std::setw(3) << ind << std::setw(8) << n << std::endl;
  }
  if(verbose()>1) {
    PrintPVListHeader(os);
    int i=0;
    for(auto& pv: ss) {
      Print(pv,i++,os);
    }
  }

}


void 
mu2e::PhysicalVolumePrinter::Print(const 
	   std::pair<cet::map_vector_key, mu2e::PhysicalVolumeInfo>& obj, 
	   int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  int i = obj.first.asInt();
  const PhysicalVolumeInfo& pv = obj.second;
  os 
    << " " << std::setw(5) << i << std::setw(5) << pv.copyNo() 
    << "  "
    << " " << std::setw(43) << pv.name() 
    << " " << std::setw(23) << pv.materialName() 
    << std::endl;

}

void 
mu2e::PhysicalVolumePrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::PhysicalVolumePrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind    count\n";

}

void 
mu2e::PhysicalVolumePrinter::PrintPVListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << " ind    id copyNo                                     name                    material\n";
}

