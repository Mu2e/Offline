
#include "Offline/Print/inc/CaloClusterMCPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void
mu2e::CaloClusterMCPrinter::Print(art::Event const& event,
                                std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CaloClusterMCCollection> > vah = event.getMany<CaloClusterMCCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CaloClusterMCCollection>(tag);
      Print(ih);
    }
  }
}

void
mu2e::CaloClusterMCPrinter::Print(const art::Handle<CaloClusterMCCollection>& handle,
                                std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void
mu2e::CaloClusterMCPrinter::Print(const art::ValidHandle<CaloClusterMCCollection>& handle,
                                std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void
mu2e::CaloClusterMCPrinter::Print(const CaloClusterMCCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CaloClusterMCCollection has " << coll.size() << " clusters\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void
mu2e::CaloClusterMCPrinter::Print(const art::Ptr<CaloClusterMC>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void
mu2e::CaloClusterMCPrinter::Print(const mu2e::CaloClusterMC& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  if( obj.totalEnergyDep() < _eCut ) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  if(verbose()==1) {
    os
      << " " << std::setw(5) << obj.caloHitMCs().size()
      << " " << std::setw(5) << obj.energyDeposits().size()
      << " "
      << " " << std::setw(8) << std::setprecision(1) << obj.totalEnergyDep()
      << " " << std::setw(8) << std::setprecision(1) << obj.totalEnergyDepG4()
      << std::endl;
  } else if(verbose()==2) {
     os
       << "  caloHitMCs:";
     for(auto const& h : obj.caloHitMCs() )  os << " " << h.key();
     os << std::endl;

     os << "     SimPart  rel   time    eDep    eDepG4    mom" << std::endl;
     for(auto const& cemc: obj.energyDeposits()) {
       os
         << "     " << std::setw(5) << cemc.sim().key()
         << " " << std::setw(5) << cemc.rel().relationship()
         << " " << std::setw(8) << std::setprecision(1) << cemc.time()
         << " " << std::setw(8) << std::setprecision(1) << cemc.energyDep()
         << " " << std::setw(8) << std::setprecision(1) << cemc.energyDepG4()
         << " " << std::setw(8) << std::setprecision(1) << cemc.momentumIn()
         << std::endl;
     }   // loop over eDep MC
  } // verbose == 2
}

void
mu2e::CaloClusterMCPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void
mu2e::CaloClusterMCPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind   nHits   nDeps   eDep eDepG4\n";

}

