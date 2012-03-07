//
// Calorimeter's data cluster container
//
// $Id: CaloCluster.hh,v 1.2 2012/03/07 18:00:38 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/03/07 18:00:38 $
//
// Original author G. Pezzullo & G. Tassielli
//


#ifndef RecoDataProducts_CaloCluster_hh
#define RecoDataProducts_CaloCluster_hh

// Mu2e includes:
#include "art/Persistency/Common/Ptr.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"


// C++ includes
#include <vector>

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"


namespace mu2e {


typedef art::Ptr< CaloCrystalHit>                  CaloCrystalHitPtr;
typedef std::vector<CaloCrystalHitPtr>      CaloCrystalHitPtrVector;


struct CaloCluster{

private:


public:
        int                                            vaneId;    //number of vane
        size_t                                    clusterSize;
        float                                            time;    //(ns)
        //float                                           _dt;    //(ns)
        int                                            cogRow;    //running from 0, 1, ... (nCrystalR - 1)
        int                                         cogColumn;    //running from 0, 1, ... (nCrystalZ - 1)
        float                                       energyDep;    //(MeV)
        CLHEP::Hep3Vector                          cog3Vector;    //center of gravity of the cluster in the mu2e frame
        CLHEP::Hep3Vector                     cog3VectorError;    //
        CaloCrystalHitPtrVector      caloCrystalHitsPtrVector;

        // public:

        CaloCluster():
                vaneId(0),
                clusterSize(0),
                time(0.),
                cogRow(0),
                cogColumn(0),
                //_dt(0.),
                energyDep(0.){
        }

        CaloCluster(int                                    iVane):
                vaneId(iVane),
                clusterSize(0),
                time(0.),
                cogRow(0),
                cogColumn(0),
                //_dt(0.),
                energyDep(0.){
        }

        CaloCluster(int                                    iVane,
                        float                                   time,
                        //float                                     dt,
                        float                                 energy,
                        CaloCrystalHitPtrVector              caloCrystalHits):
                                vaneId(iVane),
                                time(time),
                                //_dt(dt),
                                energyDep(energy){

                for(CaloCrystalHitPtrVector::iterator it = caloCrystalHits.begin(); it!= caloCrystalHits.end(); ++it){
                        caloCrystalHitsPtrVector.push_back(*it);
                }
                clusterSize=caloCrystalHits.size();
                cogRow=0;
                cogColumn=0;
        }

        void AddHit (CaloCrystalHitPtr &a){
                caloCrystalHitsPtrVector.push_back(a);
                time*=(float)clusterSize;
                ++clusterSize;
                time += a->time();
                time /= (float)clusterSize;
                energyDep += a->energyDep();
        }

};

} // namespace mu2e

#endif /* RecoDataProducts_CaloCluster_hh */
