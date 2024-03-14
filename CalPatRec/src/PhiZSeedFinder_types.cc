///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalPatRec/inc/ChannelID.hh"
#include "Offline/CalPatRec/inc/PhiZSeedFinder_types.hh"

using CalPatRec::ChannelID;

namespace mu2e {
  namespace PhiZSeedFinderTypes {

    float stationZ   [kNStations];
//-----------------------------------------------------------------------------
    Data_t::Data_t() {
      for (int is=0; is<kNStations; is++) {
        //        fListOfSeeds    [is].reserve(100);
      }
    }
//-----------------------------------------------------------------------------
// only fListOfSeeds contains new'ed memory,
// the ManagedList destructor should handle that
//-----------------------------------------------------------------------------
    Data_t::~Data_t() {
    }

//-----------------------------------------------------------------------------
// do not deallocate memory used by fListOfSeeds, re-use it
//-----------------------------------------------------------------------------
    void Data_t::InitEvent(const art::Event* Evt, int DebugLevel) {
      event      = Evt;
      debugLevel = DebugLevel;

      for (int is=0; is<kNStations; is++) {
//-----------------------------------------------------------------------------
// re-initialize faces ... no timing bins
//-----------------------------------------------------------------------------
        for (int face=0; face<kNFaces; face++) {
          FaceZ_t* fz = &fFaceData[is][face];
          fz->fHitData.clear() ;
        }
      }
//-----------------------------------------------------------------------------
// in case of a vector, 'clear()' erases it
//-----------------------------------------------------------------------------
//      fListOfDeltaCandidates.clear();
//-----------------------------------------------------------------------------
// proton candidates will be reinitialized one by one, as needed
// just set the number of used ones to zero...
//-----------------------------------------------------------------------------
//      fListOfProtonCandidates.clear();
    }
//-----------------------------------------------------------------------------
// called just once
//-----------------------------------------------------------------------------
    void Data_t::InitGeometry() {
      mu2e::GeomHandle<mu2e::Tracker> tH;
      tracker     = tH.get();

      // mu2e::GeomHandle<mu2e::DiskCalorimeter> cH;
      // calorimeter = cH.get();

//-----------------------------------------------------------------------------
// define station Z coordinates and calculate the time-of-flight between
// the station and each calorimeter disk for a typical mu-->e conversion electron
//-----------------------------------------------------------------------------
      for (unsigned ipl=0; ipl<tracker->nPlanes(); ipl += 2) {
        const Plane* p1 = &tracker->getPlane(ipl);
        const Plane* p2 = &tracker->getPlane(ipl+1);
        int ist         = ipl/2;
        stationZ[ist]   = (p1->origin().z()+p2->origin().z())/2;
      }
//-----------------------------------------------------------------------------
// per-panel constants
//-----------------------------------------------------------------------------
      ChannelID cx, co;

      int npl = tracker->nPlanes();
      for (int ipl=0; ipl< npl; ipl++) {
        const Plane* pln = &tracker->getPlane(ipl);
        int  ist = ipl/2;

        for (unsigned ipn=0; ipn<pln->nPanels(); ipn++) {
          const Panel* panel = &pln->getPanel(ipn);

          int face;
          if (panel->id().getPanel() % 2 == 0) face = 0;
          else                                 face = 1;

          cx.Station   = ist;
          cx.Plane     = ipl % 2;
          cx.Face      = face;
          cx.Panel     = ipn;

          ChannelID::orderID (&cx, &co);

          FaceZ_t*  fz = &fFaceData[co.Station][co.Face];
          Pzz_t*    pz = fz->Panel(co.Panel);
          pz->fID      = 3*co.Face+co.Panel;

          pz->wx  = panel->straw0Direction().x();
          pz->wy  = panel->straw0Direction().y();
//-----------------------------------------------------------------------------
// can't simply define nx,ny via the wx,wy - there is a couple of sign inversions:
// - plane=0,face=1 is flipped wrt plane=0,face=0
// - plane=1 is flipped wrt plane-0
//-----------------------------------------------------------------------------
          float phi = panel->straw0MidPoint().phi();
          pz->nx    = cos(phi);
          pz->ny    = sin(phi);
          pz->z     = (panel->getStraw(0).getMidPoint().z()+panel->getStraw(1).getMidPoint().z())/2.;

          if (co.Panel == 0) fz->z = pz->z;
        }
        stationUsed[ist] = 1;
      }
//-----------------------------------------------------------------------------
// precalculate overlaps between panels - 12x12 matrix,
// use first aand second stations for that
// a panel covers 120 deg, two panels overlap if n1*n2 > -sin(30 deg) = -0.5
//-----------------------------------------------------------------------------
      for (int ich=0; ich<2; ich++) {
        for (int f1=0; f1<3; f1++) {
          FaceZ_t* fz1 = &fFaceData[ich][f1];
          for (int ip1=0; ip1<3; ip1++) {
            Pzz_t* pz1 = fz1->Panel(ip1);
            for (int f2=f1+1; f2<4; f2++) {
              FaceZ_t* fz2 = &fFaceData[ich][f2];
              for (int ip2=0; ip2<3; ip2++) {
                Pzz_t* pz2 = fz2->Panel(ip2);
                float n1n2 = pz1->nx*pz2->nx+pz1->ny*pz2->ny;
                if (n1n2 >= -0.5-1.e-6) {
                  panelOverlap[ich][pz1->fID][pz2->fID] = 1;
                  panelOverlap[ich][pz2->fID][pz1->fID] = 1;
                }
                else {
                  panelOverlap[ich][pz1->fID][pz2->fID] = 0;
                  panelOverlap[ich][pz2->fID][pz1->fID] = 0;
                }
              }
            }
          }
        }
      }
    }

  }
}
