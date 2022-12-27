//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////


#include "Offline/CalPatRec/inc/DeltaFinder_types.hh"

namespace mu2e {
  namespace DeltaFinderTypes {

    DeltaCandidate::DeltaCandidate() {
      fIndex  = -1;
      for(int s=0; s<kNStations; ++s) {
        dxy    [s] = -1;
        seed   [s] = NULL;
        fT0Min [s] = -1.e10;
        fT0Max [s] =  1.e10;
      }
      fFirstStation = 999;
      fLastStation  =  -1;
      fMcPart       = NULL;
      fNHits        = 0;
      fNStrawHits   = 0;
      fNHitsCE      = 0;
      n_seeds       = 0;
    }

    DeltaCandidate::DeltaCandidate(int Index, DeltaSeed* Seed, int Station) {
      fIndex  = Index;
      for(int s=0; s<kNStations; ++s) {
        dxy    [s] = -1;
        seed   [s] = NULL;
        fT0Min [s] = -1.e10;
        fT0Max [s] =  1.e10;
      }
      fFirstStation = 999;
      fLastStation  =  -1;
      fMcPart       = NULL;
      fNHits        = 0;
      fNStrawHits   = 0;
      fNHitsCE      = 0;
      n_seeds       = 0;

      if (Seed) AddSeed(Seed,Station);
    }

//-----------------------------------------------------------------------------
    void DeltaCandidate::AddSeed(DeltaSeed* Seed, int Station) {
      seed[Station]         = Seed;
      n_seeds              += 1;

      if (fFirstStation > Station) fFirstStation = Station;
      if (fLastStation  < Station) fLastStation  = Station;
//-----------------------------------------------------------------------------
// ** FIXME recalculate the center of gravity - dont' need to be exact here
// seeds with NHits=1 don't have their center of gravity defined
//-----------------------------------------------------------------------------
      if (Seed->fNFacesWithHits > 1) {
        float x = (CofM.x()*fNHits+Seed->CofM.x()*Seed->NHits())/(fNHits+Seed->NHits());
        float y = (CofM.y()*fNHits+Seed->CofM.y()*Seed->NHits())/(fNHits+Seed->NHits());
        CofM.SetX(x);
        CofM.SetY(y);
      }
      fNHits               += Seed->NHits();
      fNStrawHits          += Seed->NStrawHits();
//-----------------------------------------------------------------------------
// update t0min and t0max
// FIXME - need more reasonable limits
//-----------------------------------------------------------------------------
      float t0min = Seed->T0Min();
      float t0max = Seed->T0Max();
      float t0    = (t0min+t0max)/2;
      float dt    = 50-(t0max-t0min)/2;
      fT0Min[Station] = t0-dt;
      fT0Max[Station] = t0+dt;
//-----------------------------------------------------------------------------
// finally, set seed DeltaIndex
//-----------------------------------------------------------------------------
      Seed->fDeltaIndex = fIndex;
    }

//-----------------------------------------------------------------------------
    void DeltaCandidate::MergeDeltaCandidate(DeltaCandidate* Delta) {
      int is1 = Delta->FirstStation();
      int is2 = Delta->LastStation ();
      for (int is=is1; is<=is2; is++) {
        if (Delta->seed[is] == nullptr)                               continue;
        if (seed[is] != nullptr) {
          printf("ERROR in DeltaCandidate::MergeDeltaCandidate: ");
          printf("merged DC also has a segment in station %i\n",is);
                                                                      continue;
        }
        seed[is] = Delta->seed[is];
        seed[is]->fDeltaIndex = fIndex;
//-----------------------------------------------------------------------------
// increment hit count only if a seed has been addded
//-----------------------------------------------------------------------------
        n_seeds              += 1;
        fNHits               += Delta->seed[is]->NHits();
        fNStrawHits          += Delta->seed[is]->NStrawHits();

        if (fFirstStation > is) fFirstStation = is;
        if (fLastStation  < is) fLastStation  = is;
      }

      float x = (CofM.x()*fNHits+Delta->CofM.x()*Delta->NHits())/(fNHits+Delta->NHits());
      float y = (CofM.y()*fNHits+Delta->CofM.y()*Delta->NHits())/(fNHits+Delta->NHits());
      CofM.SetX(x);
      CofM.SetY(y);
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
  }
}
