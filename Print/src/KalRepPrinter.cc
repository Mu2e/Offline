
#include "Print/inc/KalRepPrinter.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "RecoDataProducts/inc/TrackSummary.hh"
#include "art/Framework/Principal/Provenance.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "Mu2eUtilities/inc/toHepPoint.hh"

#include <string>

void 
mu2e::KalRepPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<KalRepPtrCollection> > vah = event.getMany<KalRepPtrCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<KalRepPtrCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::KalRepPrinter::Print(const art::Handle<KalRepPtrCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::KalRepPrinter::Print(const art::ValidHandle<KalRepPtrCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::KalRepPrinter::Print(const KalRepPtrCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "KalRepPtrCollection has " << coll.size() << " tracks\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::KalRepPrinter::Print(const art::Ptr<KalRep>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::KalRepPrinter::Print(const KalRep& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  //reuse code from TrackSummaryMaker
  const KalRep* krep = &obj;

  GeomHandle<DetectorSystem> det;
  GeomHandle<VirtualDetector> vdg;
  if(!krep->fitCurrent()) {
    os << " KalRep::fitCurrent() is false " << std::endl;
    return;
  }
  TrackSummary sum(krep->fitStatus().success(),
		   krep->charge(), krep->nActive(),
		   krep->nDof(), krep->chisq(),
		   krep->t0().t0(), krep->t0().t0Err(),
		   krep->flt0());

  CLHEP::Hep3Vector entpos = 
    det->toDetector(vdg->getGlobal(VirtualDetectorId::TT_FrontPA));
  double zent = entpos.z();
  double firsthitfltlen = krep->firstHit()->kalHit()->hit()->fltLen();
  double lasthitfltlen = krep->lastHit()->kalHit()->hit()->fltLen();
  double entlen = std::min(firsthitfltlen,lasthitfltlen);
  if(!TrkHelixUtils::findZFltlen(krep->traj(),zent,entlen,0.1)) {
    throw cet::exception("RUNTIME")<<"Error from findZFltlen()\n";
  }

  double loclen(0.0);
  TrackSummary::HelixParams helix(*krep->localTrajectory(entlen,loclen));
  TrackSummary::TrackStateAtPoint st(helix,
				     fromHepPoint(krep->position(entlen)),
				     krep->momentum(entlen),
				     krep->momentumErr(entlen).covMatrix(),
				     krep->arrivalTime(entlen),
				     entlen
				     );
  
  sum.addState(st);

  _tsprinter.setVerbose(verbose());
  _tsprinter.Print(sum);
  
}

void 
mu2e::KalRepPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::KalRepPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  _tsprinter.PrintListHeader(os);
}


