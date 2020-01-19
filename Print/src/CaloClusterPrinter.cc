
#include "Print/inc/CaloClusterPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::CaloClusterPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CaloClusterCollection> > vah;
    event.getManyByType(vah);
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CaloClusterCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::CaloClusterPrinter::Print(const art::Handle<CaloClusterCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloClusterPrinter::Print(const art::ValidHandle<CaloClusterCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloClusterPrinter::Print(const CaloClusterCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CaloClusterCollection has " << coll.size() << " clusters\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::CaloClusterPrinter::Print(const art::Ptr<CaloCluster>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::CaloClusterPrinter::Print(const mu2e::CaloCluster& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  if( obj.energyDep() < _eCut ) return;

  CLHEP::Hep3Vector const&  cog = obj.cog3Vector();

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  if(verbose()==1) {
    os 
      << " " << std::setw(5) << obj.diskId()
      << " " << std::setw(5) << obj.size()
      << " " 
      << " " << std::setw(8) << std::setprecision(1) << obj.time()
      << " " << std::setw(8) << std::setprecision(1) << obj.energyDep()
      << "  " 
      << " " << std::setw(8) << std::setprecision(1) << cog.x()
      << " " << std::setw(8) << std::setprecision(1) << cog.y()
      << " " << std::setw(8) << std::setprecision(1) << cog.z()
      << std::endl;
  } else if(verbose()==2) {
     os 
       << "  secId: " << std::setw(8) << obj.diskId()
       << "  size: " << std::setw(4) << obj.size()
       << "  time: " << std::setw(8) << std::setprecision(1) << obj.time()
       << "  eDep: " << std::setw(8) << std::setprecision(1) << obj.energyDep() << "\n";

     os
       << "  e1: " << std::setw(8) << std::setprecision(1) << obj.e1()
       << "  e9: " << std::setw(8) << std::setprecision(1) << obj.e9()
       << "  e25: " << std::setw(8) << std::setprecision(1) << obj.e25()
       << "  2ndMom: " << std::setw(8) << std::setprecision(1) << obj.secondMoment()
       << "  angle: " << std::setw(8) << std::setprecision(1) << obj.angle() << "\n";
     os
       << "  COG: " << std::setw(8) << std::setprecision(1) << cog.x()
       << " " << std::setw(8) << std::setprecision(1) << cog.y()
       << " " << std::setw(8) << std::setprecision(1) << cog.z() 
       << "    split: " << std::setw(4) << obj.isSplit()
       << "\n";
     os 
       << " CrystalHits:";
       for(auto& ic: obj.caloCrystalHitsPtrVector()) {
	 os << " " << (*ic).id();
       }
      os << "\n";

  }
}

void 
mu2e::CaloClusterPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::CaloClusterPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind  secId   size    time    energy            position\n";

}

