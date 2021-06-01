
#include "Print/inc/TrkCaloIntersectPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::TrkCaloIntersectPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<TrkCaloIntersectCollection> > vah = event.getMany<TrkCaloIntersectCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<TrkCaloIntersectCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::TrkCaloIntersectPrinter::Print(const art::Handle<TrkCaloIntersectCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::TrkCaloIntersectPrinter::Print(const art::ValidHandle<TrkCaloIntersectCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::TrkCaloIntersectPrinter::Print(const TrkCaloIntersectCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "TrkCaloIntersectCollection has " << coll.size() << " intersections\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::TrkCaloIntersectPrinter::Print(const art::Ptr<TrkCaloIntersect>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::TrkCaloIntersectPrinter::Print(const mu2e::TrkCaloIntersect& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  KalRepPtr const&  tptr = obj.trk();
  KalRepPtr::key_type tkey = 0;
  if(tptr) tkey = tptr.key();

  os 
    << " " << std::setw(5) << obj.diskId()
    << " " << std::setw(6) << tkey
    << " " << std::setw(7) << obj.trkId()
    << "  "
    << " " << std::setw(8) << std::setprecision(1) << obj.pathLengthEntrance()
    << " " << std::setw(8) << std::setprecision(1) 
          << obj.pathLenghtEntranceErr()
    << "    " << std::setw(8) << std::setprecision(1) << obj.pathLengthExit()
    << std::endl;

}

void 
mu2e::TrkCaloIntersectPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::TrkCaloIntersectPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind   secId   trkKey  trkId  path_ent path_ent_err path_exit\n";

}

