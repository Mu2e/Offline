#ifndef RecoDataProducts_CaloCluster_hh
#define RecoDataProducts_CaloCluster_hh
//
// Calorimeter's data cluster container
//
// $Id: CaloCluster.hh,v 1.6 2013/03/08 01:22:32 echenard Exp $
// $Author: echenard $
// $Date: 2013/03/08 01:22:32 $
//
// Original author G. Pezzullo & G. Tassielli
//

// Mu2e includes:
#include "art/Persistency/Common/Ptr.h"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"

// C++ includes
#include <vector>
#include <ostream>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {



   class CaloCluster {
       

       public:
    
            typedef art::Ptr< CaloCrystalHit>          CaloCrystalHitPtr;
            typedef std::vector<CaloCrystalHitPtr>     CaloCrystalHitPtrVector;

	    CaloCluster() : 
	       _sectionId(-1),_time(-1),_energyDep(-1),_e1(-1),_e9(-1),_e25(-1),_secondMoment(-1),_angle(0),_cog3Vector(CLHEP::Hep3Vector(0,0,0)),
	       _caloCrystalHitsPtrVector(),_isSplit(false)
	    {}

	    CaloCluster(int iSection, double time, double energy, CaloCrystalHitPtrVector caloCrystalHits, bool isSplit) : 
	       _sectionId(iSection),_time(time),_energyDep(energy),_e1(-1),_e9(-1),_e25(-1),_secondMoment(-1),_angle(0),
	       _cog3Vector(CLHEP::Hep3Vector(0,0,0)),_caloCrystalHitsPtrVector(caloCrystalHits),_isSplit(isSplit)
	    {}

	    void print(std::ostream& ost = std::cout) const;


	    //Accessors
	    int                                           sectionId() const{return _sectionId;}       
	    int                                                size() const{return _caloCrystalHitsPtrVector.size();}
	    double                                             time() const{return _time;}            
	    double                                        energyDep() const{return _energyDep;}       
	    double                                               e1() const{return _e1;}       
	    double                                               e9() const{return _e9;}       
	    double                                              e25() const{return _e25;}       
	    double                                     secondMoment() const{return _secondMoment;}       
	    double                                            angle() const{return _angle;}       
	    CLHEP::Hep3Vector const&                     cog3Vector() const{return _cog3Vector;}      
	    CaloCrystalHitPtrVector const& caloCrystalHitsPtrVector() const{return _caloCrystalHitsPtrVector;}
	    bool                                            isSplit() const{return _isSplit;} 


	    //Setting parameters
	    void cog3Vector(CLHEP::Hep3Vector cog3Vector)      {_cog3Vector = cog3Vector;}
	    void energyRing(double e1, double e9, double e25)  {_e1 = e1;_e9 = e9;_e25 = e25;}
	    void secondMoment(double second)                   {_secondMoment = second;}
	    void angle(double angle)                           {_angle = angle;}
         


	 private:
	 
	    int                      _sectionId;   
	    double                   _time;       
	    double                   _energyDep;
	    double                   _e1;
	    double                   _e9;
	    double                   _e25;	      
	    double                   _secondMoment;	      
	    double                   _angle;	      
	    CLHEP::Hep3Vector        _cog3Vector; 
	    CaloCrystalHitPtrVector  _caloCrystalHitsPtrVector;
	    bool                     _isSplit;    

   };

} 

#endif /* RecoDataProducts_CaloCluster_hh */
