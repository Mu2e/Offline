 //
 // $Id: CaloClusterTools.cc,v 1.3 2013/03/05 20:33:25 aluca Exp $
 // $Author: aluca $
 // $Date: 2013/03/05 20:33:25 $
 //
 // Original author G. Pezzullo, A. Luca' & G. Tassielli
 //

 // C++ includes
 #include <ostream>

 // Framework includes.
 #include "cetlib/exception.h"

 // Mu2e includes
 #include "CaloCluster/inc/CaloClusterTools.hh"
 #include "CalorimeterGeom/inc/VaneCalorimeter.hh"
 #include "GeometryService/inc/GeomHandle.hh"
 #include <cmath>

// Other includes
#include "cetlib/pow.h"

 using namespace std;

 namespace mu2e {
   double timeScale(double energy){
     double res = 0.0;
     double timeScale = 1.2;//[ns]
     res = timeScale/std::sqrt(energy);
     return res;
   }
   void CaloClusterTools::print( std::ostream& ost, bool doEndl ) const{
     _cluster.print(ost, false);
     ost << " timeErr"         << timeErr() << "\n"
	 << " energyDepErr: "  << energyDepErr();
     if ( doEndl ){
       ost << endl;
    }
   }
   
   CaloClusterTools::CaloClusterTools(CaloCluster const &clu):
     _cluster(clu),
     _calorimeter(*GeomHandle<VaneCalorimeter>()){
   }
   
   double CaloClusterTools::timeFasterCrystal() const{
     
     double time = 0.0;
     double tmp =0.0;
     
     if(_cluster.caloCrystalHitsPtrVector().size() == 0) return time;
     
     time = _cluster.caloCrystalHitsPtrVector().at(0)->time();
     for( size_t itCD=0; itCD<_cluster.caloCrystalHitsPtrVector().size(); ++itCD){
       tmp =_cluster.caloCrystalHitsPtrVector().at(itCD)->time();
       if(tmp<time){
	 time = tmp;
       }
     }
     
     return time;
   }
   
   double CaloClusterTools::timeFasterCrystalErr() const{
     double time = 0.0
       , tmp =0.0
       , energy = 0.0;
     
     if(_cluster.caloCrystalHitsPtrVector().size() == 0) return time;
     
     time = _cluster.caloCrystalHitsPtrVector().at(0)->time();
     for( size_t itCD=0; itCD<_cluster.caloCrystalHitsPtrVector().size(); ++itCD){
       tmp =_cluster.caloCrystalHitsPtrVector().at(itCD)->time();
       if(tmp<time){
	 time = tmp;
	 energy = _cluster.caloCrystalHitsPtrVector().at(itCD)->energyDep();
       }
     }
     time = timeScale(energy);
     return time;
   }
   
   double   CaloClusterTools::energyDepErr( ) const{
     
     double stocasticTerm = 0.014;
     double argo =0.0;
     argo = _cluster.energyDep()*0.001;//conversion in GeV
     argo = stocasticTerm/pow(argo, 0.25);
     argo *= _cluster.energyDep();
     
     return argo;
   }
   
   double CaloClusterTools::timeErr( ) const{
     double  tmpEq = 0.;
     
     double weight = 0.0;
     double tmp =0.0;
     
     for( size_t itCD=0; itCD<_cluster.caloCrystalHitsPtrVector().size(); ++itCD){
       tmp = (timeScale(_cluster.caloCrystalHitsPtrVector().at(itCD)->energyDep() ) )*( _cluster.caloCrystalHitsPtrVector().at(itCD)->energyDep() );
       weight += std::pow(tmp, 2);
       
       tmpEq += std::pow(_cluster.caloCrystalHitsPtrVector().at(itCD)->energyDep(), 2.);
     }
     
     weight /= tmpEq;
     weight = std::sqrt(weight);weight = std::sqrt(weight);
     
     return weight;
   }
   
   
   
   double CaloClusterTools::showerDir( ) const{
     //Get handle to calorimeter
     GeomHandle<VaneCalorimeter> cg;
     
     double R =0.;
     double Z =0., V= 0.0, W=0.0, VW=0.0, VQ = 0.0, RQ=0.0, RZ=0.0;
     
     for( size_t i=0; i< _cluster.caloCrystalHitsPtrVector().size(); ++i){
       std::vector<art::Ptr<CaloHit> > const& ROIds = (*_cluster.caloCrystalHitsPtrVector().at(i)).readouts();
       CaloHit const& thehit = *ROIds.at(0);
       
       //Get Z and R from readout
       double tZ = cg->crystalZByRO(thehit.id());
       double tR = cg->crystalRByRO(thehit.id());
       
       V += tR;
       W += tZ;
       VW += tR*tZ;
       VQ += pow(tR,2);
       
       //Multiply Z and R for the crystal's energy, which is the weight we use in that algorithm
       tZ *=(_cluster.caloCrystalHitsPtrVector().at(i))->energyDep();
       tR*=(_cluster.caloCrystalHitsPtrVector().at(i))->energyDep();
       
       Z += tZ;
       RQ += pow(tR, 2);
       R += tR;
       RZ += tZ*tR;
     }
     
     int size = _cluster.caloCrystalHitsPtrVector().size();
     double m = V*W - size*VW;
     m /= (pow(V,2) - size*VQ);

     double D = RQ - pow(R,2);
     double errM = (1- pow(m,2)*D);
     if(errM>=0.0){
       errM = sqrt(errM);
	 }else{
       errM = 1e-4;
     }
     
     
     return m;
   }
   
   double CaloClusterTools::errShowerDir( ) const{
     GeomHandle<VaneCalorimeter> cg;

     double R =0.;
     double Z =0., V= 0.0, W=0.0, VW=0.0, VQ = 0.0, RQ=0.0, RZ=0.0;
     
     for( size_t i=0; i< _cluster.caloCrystalHitsPtrVector().size(); ++i){
       std::vector<art::Ptr<CaloHit> > const& ROIds = (*(_cluster.caloCrystalHitsPtrVector().at(i)) ).readouts();
       CaloHit const& thehit = *ROIds.at(0);
       
       //Get Z and R from readout
       double tZ = cg->crystalZByRO(thehit.id());
       double tR = cg->crystalRByRO(thehit.id());
       
       V += tR;
       W += tZ;
       VW += tR*tZ;
       VQ += pow(tR,2);
       
       //Multiply Z and R for the crystal's energy, which is the weight we use in that algorithm
       tZ *=(_cluster.caloCrystalHitsPtrVector().at(i))->energyDep();
       tR*=(_cluster.caloCrystalHitsPtrVector().at(i))->energyDep();
       
       Z += tZ;
       RQ += pow(tR, 2);
       R += tR;
       RZ += tZ*tR;
       
       //Calculate the sum of the square of the weight for derive the error of the weighted mean
     }
     
     int size = _cluster.caloCrystalHitsPtrVector().size();
     double m = V*W - size*VW;
     m /= (pow(V,2) - size*VQ);
     
     double D = RQ - pow(R,2);
     double errM = (1- pow(m,2)*D);
     if(errM>=0.0){
       errM = sqrt(errM);
     }else{
       errM = 1e-4;
     }
     
     return errM;
   }
   
   int CaloClusterTools::wSize( ) const{
     GeomHandle<VaneCalorimeter> cg;
	 int tmpWmin = cg->nCrystalZ(), tmpWmax = 0;
	 int res = -1;
	 for(size_t i=0; i<_cluster.caloCrystalHitsPtrVector().size(); ++i){
	   
	   CaloCrystalHit const & hit = *(_cluster.caloCrystalHitsPtrVector().at(i));
	   
	   std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
	   if(ROIds.size()<1 ) continue;
	   
	   CaloHit const & thehit = *ROIds.at(0);
	   
	   
	   if(cg->crystalZByRO(thehit.id()) > tmpWmax){
	     tmpWmax = cg->crystalZByRO(thehit.id());
	   }
	   if(cg->crystalZByRO(thehit.id()) < tmpWmin){
	     tmpWmin = cg->crystalZByRO(thehit.id());
	   }
	 }//end loop on _caloCrystalHitsPtr
	 
	 res = tmpWmax - tmpWmin + 1;
	 return res;
	 
 }
   
   int CaloClusterTools::vSize( ) const{
     GeomHandle<VaneCalorimeter> cg;
     int tmpVmin = cg->nCrystalR(), tmpVmax = 0;
     int res = -1;
     for(size_t i=0; i<_cluster.caloCrystalHitsPtrVector().size(); ++i){
       
       CaloCrystalHit const & hit = *(_cluster.caloCrystalHitsPtrVector().at(i));
       
       std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
		 if(ROIds.size()<1 ) continue;
		 
		 CaloHit const & thehit = *ROIds.at(0);
		 
		 
		 if(cg->crystalRByRO(thehit.id()) > tmpVmax){
		   tmpVmax = cg->crystalRByRO(thehit.id());
		 }
		 if(cg->crystalRByRO(thehit.id()) < tmpVmin){
		   tmpVmin = cg->crystalRByRO(thehit.id());
		 }
     }//end loop on _caloCrystalHitsPtr
     
     res = tmpVmax - tmpVmin + 1;
     return res;
     
   }

   
   
   int CaloClusterTools::cryEnergydepMaxRow( ) const{
     GeomHandle<VaneCalorimeter> cg;
	 double tmpMaxE = 0.0;
	 int row = 0;
	 for(size_t i=0; i<_cluster.caloCrystalHitsPtrVector().size(); ++i){
	   
	   CaloCrystalHit const & hit = *(_cluster.caloCrystalHitsPtrVector().at(i));
	   if(hit.energyDep() > tmpMaxE){
	     tmpMaxE = hit.energyDep();
	     std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
	     if(ROIds.size()<1 ) continue;
	     
	     CaloHit const & thehit = *ROIds.at(0);
	     
	     
	     row    =  cg->crystalRByRO(thehit.id());
	     
	     
	     
	     
	   }
	 }//end loop on _caloCrystalHitsPtr
	 return row;
   }
   
   int CaloClusterTools::cryEnergydepMaxColumn( ) const{
	 GeomHandle<VaneCalorimeter> cg;
	 double tmpMaxE = 0.0;
	 int column = 0;
	 for(size_t i=0; i<_cluster.caloCrystalHitsPtrVector().size(); ++i){
	   
	   CaloCrystalHit const& hit = *(_cluster.caloCrystalHitsPtrVector().at(i));
	   
	   if(hit.energyDep() > tmpMaxE){
	     tmpMaxE = hit.energyDep();
	     std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
	     if(ROIds.size()<1 ) continue;
	     
	     
	     CaloHit const & thehit = *ROIds.at(0);
	     
	     column   =  cg->crystalZByRO(thehit.id());
	     
	     
	   }
	 }//end loop on _caloCrystalHitsPtr
	 return column;
   }
   
 } // namespace mu2e
