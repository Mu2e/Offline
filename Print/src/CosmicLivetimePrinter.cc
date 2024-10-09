#include "Offline/Print/inc/CosmicLivetimePrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::CosmicLivetimePrinter::PrintSubRun(art::SubRun const& subrun,
                                              std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<CosmicLivetime> > vah =
        subrun.getMany<CosmicLivetime>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    unsigned ncosmic(0);
    for (const auto& tag : tags()) {
      auto ih = subrun.getValidHandle<CosmicLivetime>(tag);
      Print(ih);
      ++ncosmic;
      nPrimaries_ += ih->primaries();
      livetime_ += ih->liveTime();
    }
    if(ncosmic> 1) std::cout << ">1 Cosmic Livetime!" << std::endl;
  }
}

void mu2e::CosmicLivetimePrinter::Print(art::Event const& event,
                                        std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<CosmicLivetime> > vah =
        event.getMany<CosmicLivetime>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<CosmicLivetime>(tag);
      Print(ih);
    }
  }
}

void mu2e::CosmicLivetimePrinter::Print(
    const art::Handle<CosmicLivetime>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::CosmicLivetimePrinter::Print(
    const art::ValidHandle<CosmicLivetime>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::CosmicLivetimePrinter::Print(const mu2e::CosmicLivetime& obj,
                                        int ind, std::ostream& os) {
  ++nsub_;
  nPrimaries_ += obj.primaries();
  livetime_ += obj.liveTime();
  if (verbose() < 1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);

  os << "  N Primaries:" << obj.primaries() << "  Area:" << obj.area()
     << "  LowE:" << obj.lowE() << "  HighE:" << obj.highE()
     << "  FluxConstant:" << obj.fluxConstant()
     << "  Livetime:" << obj.liveTime() << std::endl;
}

void mu2e::CosmicLivetimePrinter::PrintHeader(const std::string& tag,
                                              std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}

void mu2e::CosmicLivetimePrinter::PrintEndJob(std::ostream& os) {
  os << "Processed " << nsub_ << " Subruns for a total of " << std::setprecision(0) << nPrimaries_ << " primaries and " << livetime_ << " seconds total livetime" << std::endl;
}
