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
    Data_t::Data_t() {
      for (int i=0; i<kNStations; i++) {
        listOfSeeds[i] = new TClonesArray("mu2e::DeltaSeed",50);
        listOfSeeds[i]->SetOwner(kTRUE);
      }

      for (int is=0; is<kNStations; is++) {
        for (int face=0; face<kNFaces; face++) {
          for (int ip=0; ip<kNPanelsPerFace; ip++) {
            PanelZ_t* pz = &oTracker[is][face][ip];
            pz->fHitData = new std::vector<HitData_t> ;
          }
        }
      }
    }

//-----------------------------------------------------------------------------
    Data_t::~Data_t() {
      for (int i=0; i<kNStations; i++) {
        delete listOfSeeds[i];
      }

      for (int is=0; is<kNStations; is++) {
        for (int face=0; face<kNFaces; face++) {
          for (int ip=0; ip<kNPanelsPerFace; ip++) {
            PanelZ_t* pz = &oTracker[is][face][ip];
            delete pz->fHitData;
          }
        }
      }
    }

//-----------------------------------------------------------------------------
    void Data_t::orderID(ChannelID* X, ChannelID* O) {
      if (X->Panel % 2 == 0) X->Face = 0;
      else                   X->Face = 1; // define original face

      O->Station = X->Station;            // stations already ordered
      O->Plane   = X->Plane;              // planes already ordered, but not necessary for ordered construct

      if (X->Station % 2 == 0) {
        if (X->Plane == 0) O->Face = 1 - X->Face;
        else               O->Face = X->Face + 2;
      }
      else {
        if (X->Plane == 0) O->Face = X->Face;
        else               O->Face = 3 - X->Face; // order face
      }

      O->Panel = int(X->Panel/2);                // order panel

      int n = X->Station + X->Plane + X->Face;   // pattern has no intrinsic meaning, just works
      if (n % 2 == 0) O->Layer = 1 - X->Layer;
      else            O->Layer = X->Layer;       // order layer
    }

//-----------------------------------------------------------------------------
  void Data_t::deOrderID(ChannelID* X, ChannelID* O) {

    X->Station = O->Station;

    X->Plane   = O->Plane;

    if(O->Station % 2 ==  0) {
      if(O->Plane == 0) X->Face = 1 - O->Face;
      else X->Face = O->Face - 2;
    }
    else {
      if(O->Plane == 0) X->Face = O->Face;
      else X->Face = 3 - O->Face;
    }

    if(X->Face == 0) X->Panel = O->Panel * 2;
    else X->Panel = 1 + (O->Panel * 2);

    int n = X->Station + X->Plane + X->Face;
    if(n % 2 == 0) X->Layer = 1 - O->Layer;
    else X->Layer = O->Layer;
  }

//-----------------------------------------------------------------------------
    void Data_t::InitEvent(const art::Event* Evt, int DebugLevel) {
      event      = Evt;
      fNSeeds    = 0;
      debugLevel = DebugLevel;

      for (int is=0; is<kNStations; is++) {
        nseeds_per_station[is] = 0;
        listOfSeeds[is]->Clear();
      }

      listOfDeltaCandidates.clear();
    }

//-----------------------------------------------------------------------------
    void Data_t::InitGeometry() {
      mu2e::GeomHandle<mu2e::Tracker> tH;
      tracker     = tH.get();

      mu2e::GeomHandle<mu2e::DiskCalorimeter> cH;
      calorimeter = cH.get();

      ChannelID cx, co;
      int       nDisks    = calorimeter->nDisk();

      double    disk_z[2] = {0, 0};                            // in the tracker frame

      for (int i=0; i<nDisks; ++i){
        Hep3Vector gpos = calorimeter->disk(i).geomInfo().origin();
        Hep3Vector tpos = calorimeter->geomUtil().mu2eToTracker(gpos);
        disk_z[i] = tpos.z();
      }
//-----------------------------------------------------------------------------
// define station Z coordinates and calculate the time-of-flight between
// the station and each calorimeter disk for a typical mu-->e conversion electron
//-----------------------------------------------------------------------------
      for (unsigned ipl=0; ipl<tracker->nPlanes(); ipl += 2) {
        const Plane* p1 = &tracker->getPlane(ipl);
        const Plane* p2 = &tracker->getPlane(ipl+1);
        int ist         = ipl/2;
        stationZ[ist]   = (p1->origin().z()+p2->origin().z())/2;

        for (int iDisk=0; iDisk<nDisks; ++iDisk){
          stationToCaloTOF[iDisk][ist] = (disk_z[iDisk] - stationZ[ist])/sin(meanPitchAngle)/CLHEP::c_light;
        }
      }

      float     z_tracker_center(0.);
      // int       nPlanesPerStation(2);

//-----------------------------------------------------------------------------
// per-panel TOF corrections
//-----------------------------------------------------------------------------
      int npl = tracker->nPlanes();
      for (int ipl=0; ipl< npl; ipl++) {
        const Plane* pln = &tracker->getPlane(ipl);
        int  ist = ipl/2;

        for (unsigned ipn=0; ipn<pln->nPanels(); ipn++) {
          const Panel* panel = &pln->getPanel(ipn);

          int face;
          if (panel->id().getPanel() % 2 == 0) face = 0;
          else                                 face = 1;

          for (unsigned il=0; il<panel->nLayers(); ++il) {
            cx.Station   = ist;
            cx.Plane     = ipl % 2;
            cx.Face      = face;
            cx.Panel     = ipn;
            cx.Layer     = il;
            orderID (&cx, &co);
            int os       = co.Station;
            int of       = co.Face;
            int op       = co.Panel;
            PanelZ_t* pz = &oTracker[os][of][op];
            pz->fPanel   = panel;
//-----------------------------------------------------------------------------
// panel caches phi of its center and the z, phi runs from 0 to 2*pi
//-----------------------------------------------------------------------------
            pz->wx  = panel->straw0Direction().x();
            pz->wy  = panel->straw0Direction().y();
            pz->phi = panel->straw0MidPoint().phi();
            pz->z   = (panel->getStraw(0).getMidPoint().z()+panel->getStraw(1).getMidPoint().z())/2.;

            int  uniqueFaceId = ipl*mu2e::StrawId::_nfaces + of;
            faceTOF[uniqueFaceId] = (z_tracker_center - pz->z)/sin(meanPitchAngle)/CLHEP::c_light;
          }
        }
        stationUsed[ist] = 1;
      }
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
      printf("  chi2all/N chi2perp/N   chi21   chi22 mintime  maxtime     X        Y         Z   nfwh nch nsh\n");
      printf("---------------------------------------------");
      printf("-------------------------------------------------------------------------------------\n");

      printf("%5i  %03i %4i %5i",Seed->Index(),Seed->fGood,Seed->fType,Seed->fDeltaIndex);
//-----------------------------------------------------------------------------
// print hit ID's in each face
//-----------------------------------------------------------------------------
      for (int face=0; face<kNFaces; face++) {
        const HitData_t* hd = Seed->hitData[face];
        if (hd == nullptr) printf(" %5i",-1);
        else {
          const ComboHit* hit = hd->fHit;
          printf(" %5i",hit->strawId().asUint16());
        }
      }

      printf(" %8.2f   %7.2f %7.2f %7.2f",Seed->Chi2AllN(),Seed->Chi2PerpN(),Seed->fChi21,Seed->fChi22);
      printf("%8.1f %8.1f",Seed->MinHitTime(),Seed->MaxHitTime());
      printf(" %8.3f %8.3f %9.3f",Seed->CofM.x(),Seed->CofM.y(),Seed->CofM.z());
      printf("%4i",Seed->fNFacesWithHits);
      printf("%4i",Seed->NHits());
      printf("%4i",Seed->NStrawHits());
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
        DeltaSeed* ds = Delta->seed[is];
        if (ds != NULL) {

          int face0 = ds->SFace(0);
          int face1 = ds->SFace(1);

          const HitData_t* hd0 = ds->HitData(face0);
          const HitData_t* hd1 = (face1 >= 0) ? ds->HitData(face1) : nullptr;

          printf("          %3i  %3i    %3i:%03i",ds->fNHits,ds->fNHitsCE,is,ds->Index());
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
    int findIntersection(const HitData_t* Hd1, const HitData_t* Hd2, Intersection_t* Result) {
      double x1, y1, x2, y2, nx1, ny1, nx2, ny2;

      const XYZVectorF& p1 = Hd1->fHit->pos();
      x1 =  p1.x();
      y1 =  p1.y();

      const XYZVectorF& p2 = Hd2->fHit->pos();
      x2 =  p2.x();
      y2 =  p2.y();

      const XYZVectorF& wdir1 = Hd1->fHit->wdir();
      nx1 = wdir1.x();
      ny1 = wdir1.y();

      const XYZVectorF& wdir2 = Hd2->fHit->wdir();
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

//-----------------------------------------------------------------------------
// testOrderID & testdeOrderID not used in module, only were used to make sure OrderID and deOrderID worked as intended
//-----------------------------------------------------------------------------
    void Data_t::testdeOrderID() {

      ChannelID x, o;

      for (int s=0; s<2; ++s) {
        for (int f=0; f<4; ++f) {
          for (int pa=0; pa<3; ++pa) {
            for (int l=0; l<2; ++l) {

              o.Station          = s;
              o.Face             = f;
              if (f < 2) o.Plane = 0;
              else       o.Plane = 1;
              o.Panel            = pa;
              o.Layer            = l;

              deOrderID(&x, &o);

              printf(" testdeOrderID: Initial(station = %i, plane = %i, face = %i, panel = %i, layer = %i)",
                     x.Station,x.Plane,x.Face,x.Panel,x.Layer);
              printf("  Ordered(station = %i, plane = %i, face = %i, panel = %i, layer = %i)\n",
                     o.Station,o.Plane,o.Face,o.Panel,o.Layer);
            }
          }
        }
      }
    }

//-----------------------------------------------------------------------------
    void Data_t::testOrderID() {

      ChannelID x, o;

      for (int s=0; s<2; ++s) {
        for (int pl=0; pl<2; ++pl) {
          for (int pa=0; pa<6; ++pa) {
            for (int l=0; l<2; ++l) {
              x.Station = s;
              x.Plane   = pl;
              x.Panel   = pa;
              x.Layer   = l;
              orderID(&x, &o);
              printf(" testOrderID: Initial(station = %i, plane = %i, face = %i, panel = %i, layer = %i)",
                     x.Station, x.Plane, x.Face, x.Panel, x.Layer);
              printf("  Ordered(station = %i, plane = %i, face = %i, panel = %i, layer = %i)\n",
                     o.Station, o.Plane, o.Face, o.Panel, o.Layer);
            }
          }
        }
      }
    }
  }
}
