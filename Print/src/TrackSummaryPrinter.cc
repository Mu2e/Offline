
#include "Print/inc/TrackSummaryPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>

void 
mu2e::TrackSummaryPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<TrackSummaryCollection> > vah = event.getMany<TrackSummaryCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<TrackSummaryCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::TrackSummaryPrinter::Print(const art::Handle<TrackSummaryCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::TrackSummaryPrinter::Print(const art::ValidHandle<TrackSummaryCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::TrackSummaryPrinter::Print(const TrackSummaryCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "TrackSummaryCollection has " << coll.size() << " tracks\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::TrackSummaryPrinter::Print(const art::Ptr<TrackSummary>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::TrackSummaryPrinter::Print(const mu2e::TrackSummary& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  if( obj.states().size()==0 ) {
    os 
      << " a track with fitstatus " << obj.fitstatus()
      << " has no trajectory"
      << std::endl;
    return;
  }

  // only look at the first state, usually the tracker entrance point
  const mu2e::TrackSummary::TrackStateAtPoint& state = obj.states()[0];
  const mu2e::TrackSummary::HelixParams& helix = state.helix();

  if(verbose()==1) {
    os 
      << " " << std::setw(5) << obj.fitstatus()
      << " " << std::setw(8) << std::setprecision(3) << obj.fitcon()
      << " " << std::setw(5) << obj.charge()
      << " " << std::setw(8) << std::setprecision(1) << state.momentum().mag()
      << " " << std::setw(8) << std::setprecision(3) << state.momentumError()
      << " " << std::setw(8) << std::setprecision(4) << helix.tanDip()
      << " " << std::setw(7) << std::setprecision(1) << helix.d0()
      << " " << std::setw(7) << std::setprecision(1) << helix.dOut()
      << " " << std::setw(7) << std::setprecision(1) << obj.t0()
      << " " << std::setw(7) << std::setprecision(4) << obj.nactive()
      << std::endl;
  } else if(verbose()>=2) {

    os 
      << "  "
      << " fitStatus: " << std::setw(3) << obj.fitstatus()
      << "  charge: " << std::setw(3) << obj.charge()
      << "  nactive: " << std::setw(3) << obj.nactive()
      << "  chi2: " << std::setw(8) << std::setprecision(2) << obj.chi2()
      << "  ndof: " << std::setw(3) << obj.ndof()
      << "\n";
    
    os
      << "  "
      << "  fitcon: " << std::setw(8) << std::setprecision(3) << obj.fitcon()
      << "  t0: " << std::setw(8) << std::setprecision(1) << obj.t0()
      << "  t0Err: " << std::setw(8) << std::setprecision(3) << obj.t0Err()
      << "  flt0: " << std::setw(8) << std::setprecision(1) << obj.flt0()
      << "\n";

    //const std::vector<TrackStateAtPoint>& 
    int i=0;
    for(const auto& st: obj.states()) {
      const auto& hh = st.helix();
      os
	<< "  State " << i++ << "\n"
	<< "  "
	<< "  d0: "     << std::setw(7) << std::setprecision(1) << hh.d0()
	<< "  phi0: "   << std::setw(7) << std::setprecision(3) << hh.phi0()
	<< "  omega: "  << std::setw(7) << std::setprecision(5) << hh.omega()
	<< "  z0: "     << std::setw(7) << std::setprecision(1) << hh.z0()
	<< "  tanDip: " << std::setw(7) << std::setprecision(3) << hh.tanDip()
	<< "\n";
      os
	<< "  "
	<< "  dOut: "     << std::setw(8) << std::setprecision(1) << hh.dOut()
	<< "  radius: "   << std::setw(8) << std::setprecision(1) << hh.radius()
	<< "  wavelength: "  << std::setw(8) << std::setprecision(1) 
	<< hh.wavelength()
	<< "\n";
      os
	<< "  "
	<< "  pos: " << std::setw(8) << std::setprecision(1) 
	   << st.position().x()
	<< " " << std::setw(8) << std::setprecision(1) << st.position().y()
	<< " " << std::setw(8) << std::setprecision(1) << st.position().z()
	<< "  mom: " << std::setw(8) << std::setprecision(3) 
	   << st.momentum().x()
	<< " " << std::setw(8) << std::setprecision(3) << st.momentum().y()
	<< " " << std::setw(8) << std::setprecision(3) << st.momentum().z()
	<< "\n";
      os
	<< "  momErr: "   << std::setw(8) << std::setprecision(3) 
           << st.momentumError()
	<< "  arrtime: "   << std::setw(8) << std::setprecision(1) 
           << st.arrivalTime()
	<< "  fltTime: "   << std::setw(8) << std::setprecision(1) 
           << st.flightLength()
	<< "  costh: "  << std::setw(8) << std::setprecision(4) << st.costh()
	<< "\n";
      if(verbose()>=3) {
	os
	  << "     d0: " << std::setw(8) << std::setprecision(1) << hh.d0() 
	    << " +/- "  << std::setw(8) << std::setprecision(1)
	    << sqrt(hh.covariance()(1,1)) << "\n";
	os
	  << "   phi0: " << std::setw(8) << std::setprecision(4) << hh.phi0() 
	    << " +/- "  << std::setw(8) << std::setprecision(4)
	    << sqrt(hh.covariance()(2,2)) << "\n";
	os
	  << "  omega: " << std::setw(8) << std::setprecision(6) << hh.omega() 
	    << " +/- "  << std::setw(8) << std::setprecision(6)
	  << sqrt(hh.covariance()(3,3)) << "\n";
	os
	  << "     z0: " << std::setw(8) << std::setprecision(1) << hh.z0() 
	    << " +/- "  << std::setw(8) << std::setprecision(1)
	    << sqrt(hh.covariance()(4,4)) << "\n";
	os
	  << " tanDip: " << std::setw(8) << std::setprecision(4) << hh.tanDip() 
	    << " +/- "  << std::setw(8) << std::setprecision(4)
	    << sqrt(hh.covariance()(5,5)) << "\n";
      }
      if(verbose()>=4) {
	os << "  Helix covariance:\n";
	PrintMatrix(hh.covariance(),os,0);

	os << "  Helix correlations:\n";
	PrintMatrix(hh.covariance(),os,1);
      }
      
    } // end loop over states

  /*
    << "  "
    << " " << std::setw(8) << std::setprecision(1) << obj.momentum().z()
    << " " << std::setw(7) << std::setprecision(1) << obj.properTime()
    << "   "
    << " " << std::setiosflags(std::ios::left) << obj.generatorId().name()
    << std::endl;

 const double maxd = th._d0 + 2./th._om;
        TrackSummary sum(krep->fitStatus().success(),
                         krep->charge(), krep->nActive(),
                         krep->nDof(), krep->chisq(),
                         krep->t0().t0(), krep->t0().t0Err(),
                         krep->flt0());

        // The following code is based on KalFitMC.cc
        CLHEP::Hep3Vector entpos = det->toDetector(vdg->getGlobal(VirtualDetectorId::TT_FrontPA));
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

    */

  }
}

void 
mu2e::TrackSummaryPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::TrackSummaryPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind  status   fitcon  chrg     p      pErr   tanDip     d0     dMax     t0    nActive\n";
}

