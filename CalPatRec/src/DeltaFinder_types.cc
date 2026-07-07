//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
#include "Offline/GeometryService/inc/GeomHandle.hh"
// #include "Offline/GeometryService/inc/DetectorSystem.hh"

#include "Offline/CalPatRec/inc/DeltaFinder_types.hh"

using CLHEP::Hep3Vector;

namespace mu2e {

  namespace DeltaFinderTypes {

    float stationZ   [kNStations];

//-----------------------------------------------------------------------------
    FaceZ_t::FaceZ_t() {
      memset(fFirst ,0,kMaxNTimeBins*sizeof(int));
      memset(fLast  ,0,kMaxNTimeBins*sizeof(int));
      // memset(fPFirst,0,kMaxNTimeBins*sizeof(int));
      // memset(fPLast ,0,kMaxNTimeBins*sizeof(int));
      fHitData.reserve(10);
      fProtonHitData.reserve(10);
    }

//-----------------------------------------------------------------------------
    Data_t::Data_t() {
      for (int is=0; is<kNStations; is++) {
        fListOfSeeds    [is].reserve(100);
      }
      for(int i = 0; i < 2; ++i) {
        for(int j = 0; j < 12; ++j) {
          for(int k = 0; k < 12; ++k) panelOverlap[i][j][k] = 0;
        }
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

      // Check if it's an On-Spill or Off-Spill event, to determine the time bins
      if(ewm && ewm->spillType() == EventWindowMarker::onspill) {
        _nTimeBins = int(2000.f / timeBin);
        if(_nTimeBins > kMaxNTimeBins) _nTimeBins = kMaxNTimeBins;
      } else { // default to the maximum
        _nTimeBins = kMaxNTimeBins;
      }

      for (int is=0; is<kNStations; is++) {
        fListOfSeeds       [is].clear();
        fListOfProtonSeeds [is].clear();
        fListOfComptonSeeds[is].clear();
//-----------------------------------------------------------------------------
// re-initialize faces
//-----------------------------------------------------------------------------
        if(doTiming > 2) watch->SetTime("FaceZ_t::clear");
        for (int face=0; face<kNFaces; face++) {
          FaceZ_t* fz = &fFaceData[is][face];
          fz->fHitData.clear() ;
          fz->fProtonHitData.clear() ;
          memset(fz->fFirst ,0xff,_nTimeBins*sizeof(int));
          memset(fz->fLast  ,0xff,_nTimeBins*sizeof(int));
          // memset(fz->fPFirst,0xff,_nTimeBins*sizeof(int));
          // memset(fz->fPLast ,0xff,_nTimeBins*sizeof(int));
        }
        if(doTiming > 2) watch->StopTime("FaceZ_t::clear");
      }
//-----------------------------------------------------------------------------
// in case of a vector, 'clear()' erases it
//-----------------------------------------------------------------------------
      fListOfDeltaCandidates.clear();
//-----------------------------------------------------------------------------
// proton candidates will be reinitialized one by one, as needed
// just set the number of used ones to zero...
//-----------------------------------------------------------------------------
      fListOfProtonCandidates.clear();
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
//-----------------------------------------------------------------------------
    int Data_t::nSeedsTot() {
      int n(0);
      for (int i=0; i<kNStations; i++) n += fListOfSeeds[i].N();
      return n;
    }

//-----------------------------------------------------------------------------
    void Data_t::printHitData(HitData_t* Hd, const char* Option) {
      const mu2e::ComboHit* ch0 = &(*chcol)[0];

      const mu2e::ComboHit* ch = Hd->fHit;

      int index = ch-ch0;

      printf("index   sid  Stn:Pln:Pnl:Str seed delta  Time     TCorr     eDep    chi2min     sigw\n");
      printf("------------------------------------------\n");
      printf("%5i %5i ",
             index,ch->strawId().asUint16());

      printf(" %3i %3i %3i %3i",
             ch->strawId().station(),
             ch->strawId().plane(),
             ch->strawId().panel(),
             ch->strawId().straw());

      int seed_index = (Hd->fSeed) ? Hd->fSeed->Index() : -1;

      printf("%5i %5i",seed_index, Hd->fDeltaIndex);
      printf("  %8.3f %8.3f %8.5f",ch->time(),ch->correctedTime(),ch->energyDep());
      printf(" %8.2f %8.2f\n",Hd->fChi2Min,ch->posRes(ComboHit::wire));
    }

//-----------------------------------------------------------------------------
    void Data_t::printDeltaSeed(DeltaSeed* Seed, const char* Option) {
      printf("---------------------------------------------");
      printf("-------------------------------------------------------------------------------------\n");
      printf("index good type delta  SHID  SHID  SHID  SHID");
      printf("  chi2all/N chi2par/N chi2perp/N   chi21   chi22 mintime  maxtime     X        Y         Z   nch nsh\n");
      printf("---------------------------------------------");
      printf("-------------------------------------------------------------------------------------\n");

      printf("%5i  %03i %4i %5i",Seed->Index(),Seed->fGood,Seed->fType,Seed->fDeltaIndex);
//-----------------------------------------------------------------------------
// print hit ID's in each face
//-----------------------------------------------------------------------------
      for (int face=0; face<kNFaces; face++) {
        const HitData_t* hd = Seed->HitData(face);
        if (hd == nullptr) printf(" %5i",-1);
        else {
          const ComboHit* hit = hd->fHit;
          printf(" %5i",hit->strawId().asUint16());
        }
      }

      printf(" %8.2f   %7.2f %7.2f %7.2f %7.2f",Seed->Chi2TotN(),Seed->Chi2ParN(),Seed->Chi2PerpN(),Seed->fChi21,Seed->fChi22);
      printf("%8.1f %8.1f",Seed->MinHitTime(),Seed->MaxHitTime());
      printf(" %8.3f %8.3f %9.3f",Seed->CofM.x(),Seed->CofM.y(),Seed->CofM.z());
      // printf("%4i",Seed->fNFacesWithHits);
      printf("%4i",Seed->nHits());
      printf("%4i",Seed->nStrawHits());
      printf("\n");
    }

//-----------------------------------------------------------------------------
    void Data_t::printDeltaCandidate(DeltaCandidate* Delta, const char* Option) {

      printf("------------------------------------------------------------------------------------------------------\n");
      printf("      i    nh  ns s1  s2     X       Y         Z          chi21   chi22   htmin   htmax   t0     \n");
      printf("------------------------------------------------------------------------------------------------------\n");
      printf(":dc:%05i %3i",Delta->Index(),Delta->fNHits);
      printf(" %3i",Delta->fNSeeds);
      printf(" %2i  %2i %7.2f %7.2f %9.2f",Delta->fFirstStation,Delta->fLastStation,
             Delta->CofM.x(),Delta->CofM.y(),Delta->CofM.z());
      printf("\n");
      printf("------------------------------------------------------------------------------------------------------\n");

      for (int is=Delta->fFirstStation;is<=Delta->fLastStation; is++) {
        DeltaSeed* ds = Delta->Seed(is);
        if (ds != NULL) {

          int face0 = ds->SFace(0);
          int face1 = ds->SFace(1);

          const HitData_t* hd0 = ds->HitData(face0);
          const HitData_t* hd1 = (face1 >= 0) ? ds->HitData(face1) : nullptr;

          printf("          %3i    %3i:%03i",ds->fNHits,is,ds->Index());
          printf(" %7.2f %7.2f %9.2f",ds->CofM.x(),ds->CofM.y(),ds->CofM.z());
          float chi22 = (hd1) ? hd1->fChi2Min : -1;
          printf(" %7.1f %7.1f",hd0->fChi2Min, chi22);
          printf(" %7.1f %7.1f",ds->MinHitTime(),ds->MaxHitTime());
          printf(" %7.1f ",Delta->T0(is));

          printf("(");
          for (int face=0; face<kNFaces; face++) {
            const HitData_t* hd = ds->HitData(face);
            if (hd == nullptr) printf(" %5i",-1);
            else {
              const ComboHit* hit = hd->fHit;
              printf(" %5i",hit->strawId().asUint16());
            }
            if (face != kNFaces-1) printf(",");
          }

          printf(")\n");
        }
      }
    }

//-----------------------------------------------------------------------------
// used in DeltaFinderAna only
//-----------------------------------------------------------------------------
    int findIntersection(const HitData_t* Hd1, const HitData_t* Hd2, Intersection_t* Result) {
      double x1, y1, x2, y2, nx1, ny1, nx2, ny2;

      const XYZVectorF& p1 = Hd1->fHit->pos();
      x1 =  p1.x();
      y1 =  p1.y();

      const XYZVectorF& p2 = Hd2->fHit->pos();
      x2 =  p2.x();
      y2 =  p2.y();

      XYZVectorF wdir1 = Hd1->fHit->uDir();
      nx1 = wdir1.x();
      ny1 = wdir1.y();

      XYZVectorF wdir2 = Hd2->fHit->uDir();
      nx2 = wdir2.x();
      ny2 = wdir2.y();

      double n1n2  = nx1*nx2+ny1*ny2;
      double r12n1 = (x1-x2)*nx1+(y1-y2)*ny1;
      double r12n2 = (x1-x2)*nx2+(y1-y2)*ny2;
//-----------------------------------------------------------------------------
// intersection point, in 2D two lines always intersect
// wd1, wd2 - distances to the hits, sign convention: delta=hit-intersection
//-----------------------------------------------------------------------------
      double t1   = (r12n2*n1n2-r12n1)/(1-n1n2*n1n2);
      Result->wd1 = -t1;

      Result->x = x1+nx1*t1;
      Result->y = y1+ny1*t1;
      Result->z = (p1.z()+p2.z())/2;

      double t2   = (r12n2-n1n2*r12n1)/(1-n1n2*n1n2);
      Result->wd2 = -t2;

      return 0;
    }

    // Utility function to convert a float to a sortable unsigned 32-bit integer key.
    uint32_t floatToSortableInt(float f) {
      uint32_t u;
      // Type-pun safe conversion using memcpy (avoids strict aliasing violations)
      memcpy(&u, &f, sizeof(float));

      // If the number is negative, flip all bits *except* the sign bit logic
      // The most significant bit indicates the sign (1 for negative)
      if (u & 0x80000000) {
        // For negative numbers, flip all *other* bits so they sort in reverse order
        // within the negative range and before positive numbers in a standard
        // unsigned integer sort.
        u = ~u;
      } else {
        // For positive numbers, simply flip the sign bit itself to move them
        // to the upper half of the unsigned integer range.
        u |= 0x80000000;
      }
      return u;
    }

    // Function to perform a single pass of Counting Sort on a specific byte
    void countingSortPass(std::vector<const ComboHit*>& input,
                          std::vector<const ComboHit*>& output,
                          const int byte_shift) {

      constexpr int max_8_bit = 256; // uint8_t(-1) is 255
      uint32_t count[max_8_bit];
      memset(count, 0, sizeof(count));
      const size_t n = input.size();
      if(n == 0) return;

      // Count the occurences of each byte value
      for (size_t i = 0; i < n; i++) {
        const float f = input[i]->time();
        const uint32_t u = floatToSortableInt(f);
        const uint8_t byte = (u >> byte_shift) & 0xFF;
        count[byte]++;
      }

      // Re-map count[byte] to last location of that byte section
      for (int i = 1; i < max_8_bit; i++) {
        count[i] += count[i - 1];
      }

      // Fill the output array by iterating backwards in the location map
      for (size_t i = n - 1; ; i--) {
        const float f = input[i]->time();
        const uint32_t u = floatToSortableInt(f);
        const uint8_t byte = (u >> byte_shift) & 0xFF;
        output[count[byte] - 1] = input[i]; // put the value in its mapped location by byte
        count[byte]--; // next byte occurrence will go 1 before this one
        if(i == 0) break; // do here by hand to avoid size_t(-1)
      }
    }
  }
}
