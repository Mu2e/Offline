#include "Offline/Print/inc/StrawGasStepPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/ConditionsService/inc/AcceleratorParams.hh"
#include <iomanip>
#include <string>

void mu2e::StrawGasStepPrinter::Print(art::Event const& event,
                                      std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<StrawGasStepCollection> > vah =
        event.getMany<StrawGasStepCollection>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<StrawGasStepCollection>(tag);
      Print(ih);
    }
  }
}

void mu2e::StrawGasStepPrinter::Print(
    const art::Handle<StrawGasStepCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::StrawGasStepPrinter::Print(
    const art::ValidHandle<StrawGasStepCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::StrawGasStepPrinter::Print(const StrawGasStepCollection& coll,
                                      std::ostream& os) {
  if (verbose() < 0 ) return;
  os << "StrawGasStepCollection has " << coll.size() << " steps\n";
  if (verbose() >= 1) PrintListHeader();
  int i = 0;
  for (const auto& obj : coll) Print(obj, i++);
}

void mu2e::StrawGasStepPrinter::Print(const art::Ptr<StrawGasStep>& obj,
                                      int ind, std::ostream& os) {
  if (verbose() < 1) return;
  Print(*obj, ind);
}

void mu2e::StrawGasStepPrinter::Print(const mu2e::StrawGasStep& obj, int ind,
                                      std::ostream& os) {
  if (verbose() < 1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if (ind >= 0) os << std::setw(4) << ind;

  long unsigned int pkey = 0;
  auto const& pptr = obj.simParticle();
  if (pptr) pkey = pptr->id().asUint();
  ConditionsHandle<AcceleratorParams> accPar("ignored");
  double mbtime = accPar->deBuncherPeriod;

  os << ", " << std::setw(5) << pkey << ", " << std::setw(8)
     << obj.strawId().asUint16() << ", " << std::setw(5) << std::setprecision(5)
     << obj.ionizingEdep() << ", " << std::setw(5) << std::setprecision(3)
     << obj.stepLength() << ", " << std::setw(5) << std::setprecision(2)
     << obj.time() << ", " << std::setw(7) << std::setprecision(2)
     << fmod(obj.time(),mbtime) << ", " << std::setw(5) << std::setprecision(2)
     << obj.momentum().R() << ", " << std::setw(5) << std::setprecision(2)
     // << obj.momentum().z() << ", " << std::setw(5) << std::setprecision(2)
     << obj.position() // << ", " << std::setw(5) << std::setprecision(2)
     // << obj.startPosition().x() << ", " << std::setw(5) << std::setprecision(2)
     // << obj.startPosition().y() << ", " << std::setw(5) << std::setprecision(2)
     // << obj.startPosition().y() << ", " << std::setw(5) << std::setprecision(2)
     <<std::endl;
}

void mu2e::StrawGasStepPrinter::PrintHeader(const std::string& tag,
                                            std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}

void mu2e::StrawGasStepPrinter::PrintListHeader(std::ostream& os) {
  if (verbose() < 1) return;
  os << " ind SimPart StrawID eDep    length    time   time%mbtime   position \n";
  // os << " ind SimPart StrawID eDep    length    time   time%mbtime    ptot   pz   position   x   y   z \n";
}
