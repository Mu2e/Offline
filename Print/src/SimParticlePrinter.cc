
#include "Offline/Print/inc/SimParticlePrinter.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::SimParticlePrinter::Print(art::Event const& event,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<SimParticleCollection> > vah =
        event.getMany<SimParticleCollection>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<SimParticleCollection>(tag);
      Print(ih);
    }
  }
}

void mu2e::SimParticlePrinter::Print(
    const art::Handle<SimParticleCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::SimParticlePrinter::Print(
    const art::ValidHandle<SimParticleCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::SimParticlePrinter::Print(const SimParticleCollection& coll,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  os << "SimParticleCollection has " << coll.size() << " particles\n";
  if (verbose() == 1) PrintListHeader();
  int i = 0;
  for (const auto& obj : coll) Print(obj.second, i++, obj.first.asUint());
}

void mu2e::SimParticlePrinter::Print(const art::Ptr<SimParticle>& obj, int ind,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  Print(*obj, ind, obj.key());
}

void mu2e::SimParticlePrinter::Print(const mu2e::SimParticle& obj, int ind,
                                     std::size_t key, std::ostream& os) {
  if (verbose() < 1) return;

  if (obj.startMomentum().vect().mag() < _pCut) return;
  if ((abs(obj.pdgId()) == PDGCode::e_minus || obj.pdgId() == PDGCode::gamma) &&
      obj.startMomentum().vect().mag() < _emPCut)
    return;
  if (_primaryOnly && (!obj.isPrimary())) return;

  art::Ptr<SimParticle> const& pptr = obj.parent();
  int pkey = -1;
  if (pptr) pkey = int(pptr.key());

  art::Ptr<GenParticle> const& gptr = obj.genParticle();
  std::string gid("none");
  if (gptr) gid = gptr->generatorId().name();

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if (ind >= 0) os << std::setw(4) << ind;

  if (verbose() == 1) {
    os << " " << std::setw(7)  << key
       << " " << std::setw(7)  << pkey
       << " " << std::setw(11) << obj.pdgId()
       << " " << std::setw(8)  << std::setprecision(1) << obj.startPosition().x()
       << " " << std::setw(8)  << std::setprecision(1) << obj.startPosition().y()
       << " " << std::setw(8)  << std::setprecision(1) << obj.startPosition().z()
       << " " << std::setw(9)  << std::setprecision(1) << obj.startMomentum().vect().mag()
       << " " << std::setw(9)  << std::setprecision(1) << obj.startMomentum().vect().perp()
       << " " << std::setw(8)  << std::setprecision(1) << obj.startGlobalTime()
       << " " << std::setw(8)  << std::setprecision(1) << obj.startProperTime()
       << "   "
       << " " << std::setw(8)  << std::setprecision(1) << obj.endPosition().x()
       << " " << std::setw(8)  << std::setprecision(1) << obj.endPosition().y()
       << " " << std::setw(8)  << std::setprecision(1) << obj.endPosition().z()
       << " " << std::setw(9)  << std::setprecision(1) << obj.endMomentum().vect().mag()
       << " " << std::setw(9)  << std::setprecision(1) << obj.endMomentum().vect().perp()
       << " " << std::setw(8)  << std::setprecision(1) << obj.endGlobalTime()
       << " " << std::setw(8)  << std::setprecision(1) << obj.endProperTime()
       << " " << std::setw(6)  << obj.endVolumeIndex()
       << "  "
       << " " << std::setiosflags(std::ios::left) << obj.stoppingCode().name()
       << std::endl;
  } else if (verbose() == 2) {
    os << "  id: " << std::setw(8) << key << " pdgId: " << std::setw(4)
       << obj.pdgId() << " parentKey: " << std::setw(8) << pkey
       << " genId: " << std::setiosflags(std::ios::left) << gid << "\n";
    os << "  Start  gtime: " << std::setw(8) << std::setprecision(1)
       << obj.startGlobalTime() << "  ptime: " << std::setw(8)
       << std::setprecision(1) << obj.startProperTime()
       << "  volIdx: " << std::setw(4) << obj.startVolumeIndex()
       << "  g4status: " << std::setw(4) << obj.startG4Status() << "\n";
    os << "  Start  creationCode: " << std::setiosflags(std::ios::left)
       << obj.creationCode().name()
       << "    realCreationCode: " << std::setiosflags(std::ios::left)
       << obj.originParticle().creationCode().name() << "\n";
    os << "  Start pos: " << std::setw(8) << std::setprecision(1)
       << obj.startPosition().x() << " " << std::setw(8) << std::setprecision(1)
       << obj.startPosition().y() << " " << std::setw(8) << std::setprecision(1)
       << obj.startPosition().z() << "  mom: " << std::setw(8)
       << std::setprecision(1) << obj.startMomentum().x() << " " << std::setw(8)
       << std::setprecision(1) << obj.startMomentum().y() << " " << std::setw(8)
       << std::setprecision(1) << obj.startMomentum().z() << "\n";
    os << "  End  gtime: " << std::setw(8) << std::setprecision(1)
       << obj.endGlobalTime() << "  ptime: " << std::setw(8)
       << std::setprecision(1) << obj.endProperTime()
       << "  volIdx: " << std::setw(4) << obj.endVolumeIndex()
       << "  g4status: " << std::setw(4) << obj.endG4Status() << "\n";
    os << "  End  stopCode: " << std::setiosflags(std::ios::left)
       << obj.stoppingCode().name() << "  stepKE: " << std::setw(8)
       << std::setprecision(1) << obj.preLastStepKineticEnergy()
       << "  nStep: " << std::setw(4) << obj.nSteps() << "\n";
    os << "  End pos: " << std::setw(8) << std::setprecision(1)
       << obj.endPosition().x() << " " << std::setw(8) << std::setprecision(1)
       << obj.endPosition().y() << " " << std::setw(8) << std::setprecision(1)
       << obj.endPosition().z() << "  mom: " << std::setw(8)
       << std::setprecision(1) << obj.endMomentum().x() << " " << std::setw(8)
       << std::setprecision(1) << obj.endMomentum().y() << " " << std::setw(8)
       << std::setprecision(1) << obj.endMomentum().z() << "\n";

    if (!obj.daughterIds().empty()) {
      os << "  daughters: ";
      for (auto& d : obj.daughterIds()) os << " " << d.asUint();
      os << std::endl;
    }

  }  // end if verbose
}

void mu2e::SimParticlePrinter::PrintHeader(const std::string& tag,
                                           std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}

void mu2e::SimParticlePrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind      key    parent    pdgId        Start  Position            P           pT    Time    ProperTime     End Position               P           pT    Time   ProperTime vol   process\n";
}
