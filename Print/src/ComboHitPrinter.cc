
#include "Print/inc/ComboHitPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::ComboHitPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<ComboHitCollection> > vah;
    event.getManyByType(vah);
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<ComboHitCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::ComboHitPrinter::Print(const art::Handle<ComboHitCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::ComboHitPrinter::Print(const art::ValidHandle<ComboHitCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::ComboHitPrinter::Print(const ComboHitCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "ComboHitCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::ComboHitPrinter::Print(const art::Ptr<ComboHit>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::ComboHitPrinter::Print(const mu2e::ComboHit& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  if(verbose()==1) {
    os 
      << " " << std::setw(5) << obj.nCombo()
      << " " << std::setw(5) << obj.nStrawHits()
      << " " 
      << " " << std::setw(8) << std::setprecision(3) << obj.pos().x()
      << " " << std::setw(8) << std::setprecision(3) << obj.pos().y()
      << " " << std::setw(9) << std::setprecision(3) << obj.pos().z()
      << " " << std::setw(7) << std::setprecision(1) << obj.time()
      << " " << std::setw(8) << std::setprecision(5) << obj.energyDep()
      << " " << std::setw(8) << std::setprecision(4) << obj.qual()
      << std::endl;
  } else if(verbose()==2) {

    os
      << "  StrawId: " << std::setw(5) << obj.strawId().asUint16()
      << "   StrawHitFlag: ";
    for(auto sn: obj.flag().bitNames()) { 
      if(obj.flag().hasAnyProperty(StrawHitFlag(sn.first))) 
	os << " " << sn.first;
    }
    os << "\n";
    os
      << "   nCombo: " << std::setw(2) << obj.nCombo()
      << " nStraw: " << std::setw(2) << obj.nStrawHits()
      << " time: " << std::setw(7) << std::setprecision(1) << obj.time()
      << " E: " << std::setw(8) << std::setprecision(5) << obj.energyDep()
      << " qual: " << std::setw(7) << std::setprecision(4) << obj.qual()
      << std::endl;
    os
      << "   wireRes: " << std::setw(8) << std::setprecision(3) << obj.wireRes()
      << " transRes: " << std::setw(8) << std::setprecision(3) << obj.transRes()
      << " wireDist: " << std::setw(8) << std::setprecision(3) << obj.wireDist()
      << "\n";
    os
      << "   pos: " << std::setw(8) << std::setprecision(3) << obj.pos().x()
      << " " << std::setw(8) << std::setprecision(3) << obj.pos().y()
      << " " << std::setw(9) << std::setprecision(3) << obj.pos().z()
      << "  dir: " << std::setw(8) << std::setprecision(3) << obj.wdir().x()
      << " " << std::setw(8) << std::setprecision(3) << obj.wdir().y()
      << " " << std::setw(9) << std::setprecision(3) << obj.wdir().z()
      << "\n";
    os 
      << "   indexArray:";
    for(auto ii: obj.indexArray()) os << " " << ii;
    os << std::endl;
  }

}

void 
mu2e::ComboHitPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::ComboHitPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os <<  "ind  nCombo nStraw   x        y         z        t        E       qual\n";

}

