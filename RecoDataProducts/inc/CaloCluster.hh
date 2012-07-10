//
// Calorimeter's data cluster container
//
// $Id: CaloCluster.hh,v 1.4 2012/07/10 00:02:20 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:20 $
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
        int                                            _vaneId;    //number of vane
        float                                            _time;    //(ns)
        int                                            _cogRow;    //running from 0, 1, ... (nCrystalR - 1)
        int                                         _cogColumn;    //running from 0, 1, ... (nCrystalZ - 1)
        float                                       _energyDep;    //(MeV)
        CLHEP::Hep3Vector                          _cog3Vector;    //center of gravity of the cluster in the localVaneFrame
        CLHEP::Hep3Vector                     _cog3VectorError;
        CaloCrystalHitPtrVector      _caloCrystalHitsPtrVector;

public:
        CaloCluster():
                _vaneId(0),
                _time(0.),
                _cogRow(0),
                _cogColumn(0),
                _energyDep(0.){
        }

        CaloCluster(int                                    iVane):
                _vaneId(iVane),
                _time(0.),
                _cogRow(0),
                _cogColumn(0),
                _energyDep(0.){
        }

        CaloCluster(int                                    iVane,
                        float                                   time,
                        float                                 energy,
                        CaloCrystalHitPtrVector              caloCrystalHits):
                                _vaneId(iVane),
                                _time(time),
                                _energyDep(energy),
                                _caloCrystalHitsPtrVector(caloCrystalHits){}

        void print( std::ostream& ost = std::cout, bool doEndl = true ) const;


        void AddHit (CaloCrystalHitPtr &a);

        //Accessors
        inline int                                                   vaneId() const{return _vaneId;}//index of the vane
        inline float                                                   time() const{return _time;}//(ns)
        float                                                       timeErr() const;//ns
        float                                             timeFasterCrystal() const;//ns
        float                                          timeFasterCrystalErr() const;//ns
        inline int                                                   cogRow() const{return _cogRow;}//running from 0, 1, ... (nCrystalR - 1)
        inline int                                                cogColumn() const{return _cogColumn;}//running from 0, 1, ... (nCrystalZ - 1)
        inline float                                              energyDep() const{return _energyDep;}//(MeV)
        float                                                  energyDepErr() const;
        inline CLHEP::Hep3Vector                                 cog3Vector() const{return _cog3Vector;}
        inline CLHEP::Hep3Vector                            cog3VectorError() const{return _cog3VectorError;}
        inline CaloCrystalHitPtrVector const&      caloCrystalHitsPtrVector() const{return _caloCrystalHitsPtrVector;}
        inline size_t                                                  size() const{return _caloCrystalHitsPtrVector.size();}
        float                                                     showerDir() const;
        float                                                  errShowerDir() const;
        int                                                           wSize() const;
        int                                                           vSize() const;
        int                                              cryEnergydepMaxRow() const;
        int                                           cryEnergydepMaxColumn() const;

        //Setting parameters
        void SetVaneId ( int vane) ;
        void SetTime (float time) ;
        void SetCogRow ( int row) ;
        void SetCogColumn ( int column) ;
        void SetEnergyDep ( float energyDep) ;
        void SetCog3Vector ( CLHEP::Hep3Vector cog3Vector) ;
        void SetCog3VectorError ( CLHEP::Hep3Vector cog3VectorErr) ;

};

} // namespace mu2e

#endif /* RecoDataProducts_CaloCluster_hh */
