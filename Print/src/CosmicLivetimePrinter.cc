#include "Offline/Print/inc/CosmicLivetimePrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include "canvas/Persistency/Common/Sampled.h"
#include "canvas/Persistency/Provenance/SampledInfo.h"
#include <iomanip>
#include <string>

void mu2e::CosmicLivetimePrinter::Print(art::Event const& event, std::ostream& os) {}

void mu2e::CosmicLivetimePrinter::PrintSubRun(art::SubRun const& subrun,
                                              std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<CosmicLivetime> > vcl =
        subrun.getMany<CosmicLivetime>();
    for (auto const& cl : vcl) Print(cl);
    // also look for sampled instances
    auto vscl = subrun.getMany<art::Sampled<CosmicLivetime>>();
    if(vscl.size() > 0){
      for (auto const& scl : vscl){
        std::cout << "SampledCosmicLivetime with tag " << scl->originalInputTag() << std::endl;
        auto sinfomh = subrun.getHandle<art::SampledSubRunInfo>("SamplingInput");
        if(sinfomh.isValid()){
          auto const& sinfom = *sinfomh;
          for(auto sinfoit = sinfom.begin(); sinfoit != sinfom.end(); ++sinfoit) {
            if(sinfoit->first.find("Cosmic") != std::string::npos){
              std::cout << "With SampledSubRunInfo entry for dataset " << sinfoit->first << " Has the following CosmicLivetimes: " << std::endl;
              for(auto const& sr : sinfoit->second.ids){
                std::cout << sr << " : ";
                auto sclp = scl->get(sinfoit->first,sr);
                if(!sclp.empty()) Print(*sclp);
              }
            }
          }
        }
      }
    }

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
  if(nsub_ > 0){
    os << "Processed " << nsub_ << " Subruns for a total of " << std::setprecision(0) << nPrimaries_ << " primaries and " << livetime_ << " seconds total livetime" << std::endl;
  }
}
