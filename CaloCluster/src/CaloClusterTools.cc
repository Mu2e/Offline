//
// $Id: CaloClusterTools.cc,v 1.2 2012/09/08 02:24:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/09/08 02:24:25 $
//
// Original author G. Pezzullo, A. Luca' & G. Tassielli
//

// C++ includes
#include <ostream>

// Mu2e includes
#include "CaloCluster/inc/CaloClusterTools.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include <cmath>

// Other includes
#include "cetlib/pow.h"

using namespace std;

namespace mu2e {

  namespace {
    double timeScale(double energy){
      double res = 0.0;
      double timeScale = 1.2;//[ns]
      res = timeScale/std::sqrt(energy);
      return res;
    }
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
      , tmp = 0.0
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

  double   CaloClusterTools::energyDepErr() const{

    double stocasticTerm = 0.014;
    double argo =0.0;
    argo = _cluster.energyDep()*0.001;//conversion in GeV
    argo = stocasticTerm/pow(argo, 0.25);
    argo *= _cluster.energyDep();

    return argo;
  }

  double CaloClusterTools::timeErr() const{
    double  tmpEq = 0.;

    double weight = 0.0;
    double tmp =0.0;

    for( size_t itCD=0; itCD<_cluster.caloCrystalHitsPtrVector().size(); ++itCD){
      tmp = (timeScale(_cluster.caloCrystalHitsPtrVector().at(itCD)->energyDep() ) )*( _cluster.caloCrystalHitsPtrVector().at(itCD)->energyDep() );
      weight += cet::square(tmp);

      tmpEq += cet::square(_cluster.caloCrystalHitsPtrVector().at(itCD)->energyDep());
    }

    weight /= tmpEq;
    weight = std::sqrt(weight);weight = std::sqrt(weight);

    return weight;
  }


  double CaloClusterTools::showerDir() const{

    double R =0.;
    double Z =0., V= 0.0, W=0.0, VW=0.0, VQ = 0.0, RQ=0.0, RZ=0.0;

    for( size_t i=0; i< _cluster.caloCrystalHitsPtrVector().size(); ++i){
      std::vector<art::Ptr<CaloHit> > const& ROIds = (*_cluster.caloCrystalHitsPtrVector().at(i)).readouts();
      CaloHit const& thehit = *ROIds.at(0);

      //Get Z and R from readout
      double tZ = _calorimeter.getCrystalZByRO(thehit.id());
      double tR = _calorimeter.getCrystalRByRO(thehit.id());

      V += tR;
      W += tZ;
      VW += tR*tZ;
      VQ += cet::square(tR);

      //Multiply Z and R for the crystal's energy, which is the weight we use in that algorithm
      tZ *=(_cluster.caloCrystalHitsPtrVector().at(i))->energyDep();
      tR*=(_cluster.caloCrystalHitsPtrVector().at(i))->energyDep();

      Z += tZ;
      RQ += cet::square(tR);
      R += tR;
      RZ += tZ*tR;
    }

    int size = _cluster.caloCrystalHitsPtrVector().size();
    double m = V*W - size*VW;
    m /= (cet::square(V) - size*VQ);

    double D = RQ - cet::square(R);
    double errM = (1- cet::square(m)*D);
    if(errM>=0.0){
      errM = sqrt(errM);
    }else{
      errM = 1e-4;
    }


    cout << "Shower dir: " << m << endl;
    return m;
  }

  double CaloClusterTools::errShowerDir() const{

    double R =0.;
    double Z =0., V= 0.0, W=0.0, VW=0.0, VQ = 0.0, RQ=0.0, RZ=0.0;

    for( size_t i=0; i< _cluster.caloCrystalHitsPtrVector().size(); ++i){
      std::vector<art::Ptr<CaloHit> > const& ROIds = (*(_cluster.caloCrystalHitsPtrVector().at(i)) ).readouts();
      CaloHit const& thehit = *ROIds.at(0);

      //Get Z and R from readout
      double tZ = _calorimeter.getCrystalZByRO(thehit.id());
      double tR = _calorimeter.getCrystalRByRO(thehit.id());

      V += tR;
      W += tZ;
      VW += tR*tZ;
      VQ += cet::square(tR);

      //Multiply Z and R for the crystal's energy, which is the weight we use in that algorithm
      tZ *=(_cluster.caloCrystalHitsPtrVector().at(i))->energyDep();
      tR*=(_cluster.caloCrystalHitsPtrVector().at(i))->energyDep();

      Z += tZ;
      RQ += cet::square(tR);
      R += tR;
      RZ += tZ*tR;

      //Calculate the sum of the square of the weight for derive the error of the weighted mean
    }

    int size = _cluster.caloCrystalHitsPtrVector().size();
    double m = V*W - size*VW;
    m /= (cet::square(V) - size*VQ);

    double D = RQ - cet::square(R);
    double errM = (1- cet::square(m)*D);
    if(errM>=0.0){
      errM = sqrt(errM);
    }else{
      errM = 1e-4;
    }

    return errM;
  }

  int CaloClusterTools::wSize() const{

    int tmpWmin = _calorimeter.nCrystalZ(), tmpWmax = 0;
    int res = -1;
    for(size_t i=0; i<_cluster.caloCrystalHitsPtrVector().size(); ++i){

      CaloCrystalHit const & hit = *(_cluster.caloCrystalHitsPtrVector().at(i));

      std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
      if(ROIds.size()<1 ) continue;

      CaloHit const & thehit = *ROIds.at(0);


      if(_calorimeter.getCrystalZByRO(thehit.id()) > tmpWmax){
        tmpWmax = _calorimeter.getCrystalZByRO(thehit.id());
      }
      if(_calorimeter.getCrystalZByRO(thehit.id()) < tmpWmin){
        tmpWmin = _calorimeter.getCrystalZByRO(thehit.id());
      }
    }//end loop on _caloCrystalHitsPtr

    res = tmpWmax - tmpWmin + 1;
    return res;

  }

  int CaloClusterTools::vSize() const{

    int tmpVmin = _calorimeter.nCrystalR(), tmpVmax = 0;
    int res = -1;
    for(size_t i=0; i<_cluster.caloCrystalHitsPtrVector().size(); ++i){

      CaloCrystalHit const & hit = *(_cluster.caloCrystalHitsPtrVector().at(i));

      std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
      if(ROIds.size()<1 ) continue;

      CaloHit const & thehit = *ROIds.at(0);


      if(_calorimeter.getCrystalRByRO(thehit.id()) > tmpVmax){
        tmpVmax = _calorimeter.getCrystalRByRO(thehit.id());
      }
      if(_calorimeter.getCrystalRByRO(thehit.id()) < tmpVmin){
        tmpVmin = _calorimeter.getCrystalRByRO(thehit.id());
      }
    }//end loop on _caloCrystalHitsPtr

    res = tmpVmax - tmpVmin + 1;
    return res;

  }



  int CaloClusterTools::cryEnergydepMaxRow() const{

    double tmpMaxE = 0.0;
    int row = 0;
    for(size_t i=0; i<_cluster.caloCrystalHitsPtrVector().size(); ++i){

      CaloCrystalHit const & hit = *(_cluster.caloCrystalHitsPtrVector().at(i));
      if(hit.energyDep() > tmpMaxE){
        tmpMaxE = hit.energyDep();
        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
        if(ROIds.size()<1 ) continue;

        CaloHit const & thehit = *ROIds.at(0);


        row    =  _calorimeter.getCrystalRByRO(thehit.id());




      }
    }//end loop on _caloCrystalHitsPtr
    return row;
  }

  int CaloClusterTools::cryEnergydepMaxColumn() const{

    double tmpMaxE = 0.0;
    int column = 0;
    for(size_t i=0; i<_cluster.caloCrystalHitsPtrVector().size(); ++i){

      CaloCrystalHit const& hit = *(_cluster.caloCrystalHitsPtrVector().at(i))
        ;
      if(hit.energyDep() > tmpMaxE){
        tmpMaxE = hit.energyDep();
        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
        if(ROIds.size()<1 ) continue;


        CaloHit const & thehit = *ROIds.at(0);

        column   =  _calorimeter.getCrystalZByRO(thehit.id());


      }
    }//end loop on _caloCrystalHitsPtr
    return column;
  }

} // namespace mu2e
