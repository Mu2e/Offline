#include "Offline/Print/inc/PrescaleFilterFractionPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include "canvas/Persistency/Common/Sampled.h"
#include "canvas/Persistency/Provenance/SampledInfo.h"
#include <iomanip>
#include <string>

void mu2e::PrescaleFilterFractionPrinter::Print(art::Event const& event, std::ostream& os) {}

void mu2e::PrescaleFilterFractionPrinter::PrintSubRun(art::SubRun const& subrun, std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<PrescaleFilterFraction> > ffl = subrun.getMany<PrescaleFilterFraction>();
    for (auto const& cl : ffl) Print(cl);
    // also look for sampled instances
    auto ffs = subrun.getMany<art::Sampled<PrescaleFilterFraction>>();
    if(ffs.size() > 0){
      for (auto const& ff : ffs){
        std::cout << "PrescaleFilterFraction with tag " << ff->originalInputTag() << std::endl;
        auto sinfomh = subrun.getHandle<art::SampledSubRunInfo>("SamplingInput");
        if(sinfomh.isValid()){
          auto const& sinfom = *sinfomh;
          for(auto sinfoit = sinfom.begin(); sinfoit != sinfom.end(); ++sinfoit) {
            if(sinfoit->first.find("Cosmic") != std::string::npos){
              std::cout << "With SubRunInfo entry for dataset " << sinfoit->first << " Has the following PrescaleFilterFractions: " << std::endl;
              for(auto const& sr : sinfoit->second.ids){
                std::cout << sr << " : ";
                auto ffp = ff->get(sinfoit->first,sr);
                if(!ffp.empty()) Print(*ffp);
              }
            }
          }
        }
      }
    }

  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ff = subrun.getValidHandle<PrescaleFilterFraction>(tag);
      Print(ff);
    }
  }
}

void mu2e::PrescaleFilterFractionPrinter::Print(
    const art::Handle<PrescaleFilterFraction>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::PrescaleFilterFractionPrinter::Print(
    const art::ValidHandle<PrescaleFilterFraction>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::PrescaleFilterFractionPrinter::Print(const mu2e::PrescaleFilterFraction& obj,
                                        int ind, std::ostream& os) {
  os << std::setiosflags(std::ios::fixed | std::ios::right);
  os << " Fraction passing filter " << obj.filterFraction()  <<  " N Seen " << obj.nSeen()
    << " with prescale fraction " << obj.prescaleFraction();
  if(obj.chained())os << " chained";
  os << std::endl;
}

void mu2e::PrescaleFilterFractionPrinter::PrintHeader(const std::string& tag,
                                              std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}
