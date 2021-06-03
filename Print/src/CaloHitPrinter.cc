
#include "Print/inc/CaloHitPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::CaloHitPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CaloHitCollection> > vah = event.getMany<CaloHitCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CaloHitCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::CaloHitPrinter::Print(const art::Handle<CaloHitCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloHitPrinter::Print(const art::ValidHandle<CaloHitCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloHitPrinter::Print(const CaloHitCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CaloHitCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::CaloHitPrinter::Print(const art::Ptr<CaloHit>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::CaloHitPrinter::Print(const mu2e::CaloHit& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  if( obj.energyDep() < _eCut ) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  os 
    << " " << std::setw(5) << obj.crystalID()
    << " " 
    << " " << std::setw(8) << std::setprecision(1) << obj.time()
    << " " << std::setw(8) << std::setprecision(1) << obj.energyDep()
    << " " << std::setw(8) << std::setprecision(1) << obj.energyDepTot()
    << " " << std::setw(5) << std::setprecision(1) << obj.nSiPMs()
    << std::endl;

  if(verbose()<2) return;

  os << "     SiPM   time   energy  chi2" << std::endl;
  for(auto const& artrdigi: obj.recoCaloDigis()) {
    if (artrdigi.isAvailable()) {
      auto const& rdigi = *artrdigi;
      os 
	<< "     " << std::setw(5) << rdigi.SiPMID() 
	<< " " 
	<< " " << std::setw(8) << std::setprecision(1) << rdigi.time()
	<< " " << std::setw(8) << std::setprecision(1) << rdigi.energyDep()
	<< " " << std::setw(8) << std::setprecision(1) << rdigi.chi2()
	<< std::endl;
    } else {
      os <<"     CaloRecoDigi Ptr !isAvailable" << std::endl;
    } //endif
  } // loop on recoCaloDigi

}

void 
mu2e::CaloHitPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::CaloHitPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind     id     time     energy    eTot   nSiPMI\n";

}

