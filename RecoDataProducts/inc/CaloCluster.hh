//
// Calorimeter's data cluster container
//
// $Id: CaloCluster.hh,v 1.3 2012/03/19 19:35:42 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/03/19 19:35:42 $
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
        //size_t                                    clusterSize;
        float                                            _time;    //(ns)
        //float                                           _dt;    //(ns)
        int                                            _cogRow;    //running from 0, 1, ... (nCrystalR - 1)
        int                                         _cogColumn;    //running from 0, 1, ... (nCrystalZ - 1)
        float                                       _energyDep;    //(MeV)
        CLHEP::Hep3Vector                          _cog3Vector;    //center of gravity of the cluster in the mu2e frame
        CLHEP::Hep3Vector                     _cog3VectorError;    //
        CaloCrystalHitPtrVector      _caloCrystalHitsPtrVector;
        float                                       _showerDir; //angular coefficient which identifies the direction of the shower in the local vane reference.
                                                                //the linear equqation is Column = showerDir*Row + offSet. the offset colud be found using the
                                                                //information that the line must pass through the cog coordinates
        float                                    _errShowerDir;


public:
        CaloCluster():
                _vaneId(0),
                //clusterSize(0),
                _time(0.),
                _cogRow(0),
                _cogColumn(0),
                //_dt(0.),
                _energyDep(0.),
                _showerDir(0.0),
                _errShowerDir(0.0){
        }

        CaloCluster(int                                    iVane):
                _vaneId(iVane),
                //rootclusterSize(0),
                _time(0.),
                _cogRow(0),
                _cogColumn(0),
                //_dt(0.),
                _energyDep(0.),
                _showerDir(0.0),
                _errShowerDir(0.0){
        }

        CaloCluster(int                                    iVane,
                        float                                   time,
                        //float                                     dt,
                        float                                 energy,
                        CaloCrystalHitPtrVector              caloCrystalHits):
                                _vaneId(iVane),
                                _time(time),
                                //_dt(dt),
                                _energyDep(energy),
                                _caloCrystalHitsPtrVector(caloCrystalHits){

                //                for(CaloCrystalHitPtrVector::iterator it = caloCrystalHits.begin(); it!= caloCrystalHits.end(); ++it){
                //                        caloCrystalHitsPtrVector.push_back(*it);
                //                }
                //clusterSize=caloCrystalHits.size();
//                _cogRow=0;
//                _cogColumn=0;
        }

        void print( std::ostream& ost = std::cout, bool doEndl = true ) const;


        void AddHit (CaloCrystalHitPtr &a) ;
//        {
//                _caloCrystalHitsPtrVector.push_back(a);
//                //time*=(float)clusterSize;
//                _time*=(float)caloCrystalHitsPtrVector.size();
//                //++clusterSize;
//                _time += a->time();
//                //time /= (float)clusterSize;
//                _time /= (float)caloCrystalHitsPtrVector.size();
//                _energyDep += a->energyDep();
//        }
        //Accessors
        inline int                                            vaneId() const{return _vaneId;}    //number of vane
        //size_t                                    clusterSize;
        inline float                                            time() const{return _time;}    //(ns)
        //float                                           _dt;    //(ns)
        inline int                                            cogRow() const{return _cogRow;}    //running from 0, 1, ... (nCrystalR - 1)
        inline int                                         cogColumn() const{return _cogColumn;}    //running from 0, 1, ... (nCrystalZ - 1)
        inline float                                       energyDep() const{return _energyDep;}    //(MeV)
        inline CLHEP::Hep3Vector                          cog3Vector() const{return _cog3Vector;}    //center of gravity of the cluster in the mu2e frame
        inline CLHEP::Hep3Vector                     cog3VectorError() const{return _cog3VectorError;}    //
        inline CaloCrystalHitPtrVector const&      caloCrystalHitsPtrVector() const{return _caloCrystalHitsPtrVector;}
        inline size_t                                           size() const{return _caloCrystalHitsPtrVector.size();}
        inline float                                       showerDir() const{return _showerDir;}    //(MeV)
        inline float                                    errShowerDir() const{return _errShowerDir;}    //(MeV)

        //Setting parameters
        void SetVaneId ( int vane) ;
        void SetTime (float time) ;
        void SetCogRow ( int row) ;
        void SetCogColumn ( int column) ;
        void SetEnergyDep ( float energyDep) ;
        void SetCog3Vector ( CLHEP::Hep3Vector cog3Vector) ;
        void SetCog3VectorError ( CLHEP::Hep3Vector cog3VectorErr) ;
        void SetShowerDir(float dir);
        void SetErrShowerDir(float errDir);
        //void AddHit (CaloCrystalHitPtr &a, float t, )


};

} // namespace mu2e

#endif /* RecoDataProducts_CaloCluster_hh */
