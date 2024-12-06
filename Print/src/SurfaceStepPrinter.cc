#include "Offline/Print/inc/SurfaceStepPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include <iomanip>
#include <string>

void mu2e::SurfaceStepPrinter::Print(art::Event const& event,
                                      std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<SurfaceStepCollection> > vah =
        event.getMany<SurfaceStepCollection>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<SurfaceStepCollection>(tag);
      Print(ih);
    }
  }
}

void mu2e::SurfaceStepPrinter::Print(
    const art::Handle<SurfaceStepCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::SurfaceStepPrinter::Print(
    const art::ValidHandle<SurfaceStepCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::SurfaceStepPrinter::Print(const SurfaceStepCollection& coll,
                                      std::ostream& os) {
  if (verbose() < 0 ) return;
  os << "SurfaceStepCollection has " << coll.size() << " steps\n";
  if (verbose() >= 1) PrintListHeader();
  int i = 0;
  for (const auto& obj : coll) Print(obj, i++);
}

void mu2e::SurfaceStepPrinter::Print(const art::Ptr<SurfaceStep>& obj,
                                      int ind, std::ostream& os) {
  if (verbose() < 1) return;
  Print(*obj, ind);
}

void mu2e::SurfaceStepPrinter::Print(const mu2e::SurfaceStep& obj, int ind,
                                      std::ostream& os) {
  if (verbose() < 1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if (ind >= 0) os << std::setw(4) << ind;

  long unsigned int pkey = 0;
  auto const& pptr = obj.simParticle();
  if (pptr) pkey = pptr->id().asUint();
  double mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();

  os << ", " << std::setw(5) << pkey << ", "
     << std::setw(8) << obj.surfaceId() << "  "
     << std::setw(5) << std::setprecision(5)     << obj.energyDeposit() << "  "
     << std::setw(5) << std::setprecision(3)     << obj.pathLength() << " "
     << std::setw(5) << std::setprecision(2)     << obj.time() << " "
     << std::setw(7) << std::setprecision(2)     << fmod(obj.time(),mbtime) << " "
     << std::setw(5)  << std::setprecision(1) << obj.startPosition().x() << ","
     << std::setw(5)  << std::setprecision(1) << obj.startPosition().y() << ","
     << std::setw(5)  << std::setprecision(1) << obj.startPosition().z() << " "
     << std::setw(5)  << std::setprecision(1) << obj.momentum().x() << ","
     << std::setw(5)  << std::setprecision(1) << obj.momentum().y() << ","
     << std::setw(5)  << std::setprecision(1) << obj.momentum().z()
     << std::endl;
}

void mu2e::SurfaceStepPrinter::PrintHeader(const std::string& tag,
                                            std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}

void mu2e::SurfaceStepPrinter::PrintListHeader(std::ostream& os) {
  if (verbose() < 1) return;
  os << " ind      SimPart        SurfaceId    eDep     length    time time%mb pos  X    Y    Z   mom  X     Y     Z\n";
}
