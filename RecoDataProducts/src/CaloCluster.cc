//
// $Id: CaloCluster.cc,v 1.2 2012/07/10 00:02:20 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:20 $
//
// Original author G. Pezzullo
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib/exception.h"

// Mu2e includes
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include <cmath>

using namespace std;

namespace mu2e {

float timeScale(float energy){
        float res = 0.0;
        float timeScale = 1.2;//[ns]
        res = timeScale/std::sqrt(energy);
        return res;
}

float CaloCluster::timeFasterCrystal() const{

        float time = 0.0;
        float tmp =0.0;

        if(_caloCrystalHitsPtrVector.size() == 0) return time;

        time = _caloCrystalHitsPtrVector.at(0)->time();
        for( size_t itCD=0; itCD<_caloCrystalHitsPtrVector.size(); ++itCD){
                tmp =_caloCrystalHitsPtrVector.at(itCD)->time();
                if(tmp<time){
                        time = tmp;
                }
        }

        return time;
}

float CaloCluster::timeFasterCrystalErr() const{
        float time = 0.0
                        , tmp =0.0
                        , energy = 0.0;

        if(_caloCrystalHitsPtrVector.size() == 0) return time;

        time = _caloCrystalHitsPtrVector.at(0)->time();
        for( size_t itCD=0; itCD<_caloCrystalHitsPtrVector.size(); ++itCD){
                tmp =_caloCrystalHitsPtrVector.at(itCD)->time();
                if(tmp<time){
                        time = tmp;
                        energy = _caloCrystalHitsPtrVector.at(itCD)->energyDep();
                }
        }
        time = timeScale(energy);
        return time;
}

float   CaloCluster::energyDepErr() const{

        float stocasticTerm = 0.014;
        float argo =0.0;
        argo = _energyDep/1000.0;//conversion in GeV
        argo = stocasticTerm/pow(argo, 0.25);
        argo *= _energyDep;

        return argo;
}

float CaloCluster::timeErr() const{
        double  tmpEq = 0.;

        float weight = 0.0;
        float tmp =0.0;

        for( size_t itCD=0; itCD<_caloCrystalHitsPtrVector.size(); ++itCD){
                tmp = (timeScale(_caloCrystalHitsPtrVector.at(itCD)->energyDep() ) )*( _caloCrystalHitsPtrVector.at(itCD)->energyDep() );
                weight += std::pow(tmp, 2);

                tmpEq += std::pow(_caloCrystalHitsPtrVector.at(itCD)->energyDep(), 2.);
        }

        weight /= tmpEq;
        weight = std::sqrt(weight);weight = std::sqrt(weight);

        return weight;
}



float CaloCluster::showerDir() const{
        //Get handle to calorimeter
        GeomHandle<Calorimeter> cg;

        float R =0.;
        float Z =0., V= 0.0, W=0.0, VW=0.0, VQ = 0.0, RQ=0.0, RZ=0.0;

        for( size_t i=0; i< _caloCrystalHitsPtrVector.size(); ++i){
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*_caloCrystalHitsPtrVector.at(i)).readouts();
                CaloHit const& thehit = *ROIds.at(0);

                //Get Z and R from readout
                double tZ = cg->getCrystalZByRO(thehit.id());
                double tR = cg->getCrystalRByRO(thehit.id());

                V += tR;
                W += tZ;
                VW += tR*tZ;
                VQ += pow(tR,2);

                //Multiply Z and R for the crystal's energy, which is the weight we use in that algorithm
                tZ *=(*_caloCrystalHitsPtrVector.at(i)).energyDep();
                tR*=(*_caloCrystalHitsPtrVector.at(i)).energyDep();

                Z += tZ;
                RQ += pow(tR, 2);
                R += tR;
                RZ += tZ*tR;
        }

        int size = _caloCrystalHitsPtrVector.size();
        float m = V*W - size*VW;
        m /= (pow(V,2) - size*VQ);

        float D = RQ - pow(R,2);
        float errM = (1- pow(m,2)*D);
        if(errM>=0.0){
                errM = sqrt(errM);
        }else{
                errM = 1e-4;
        }


        return m;
}

float CaloCluster::errShowerDir() const{
        GeomHandle<Calorimeter> cg;

        float R =0.;
        float Z =0., V= 0.0, W=0.0, VW=0.0, VQ = 0.0, RQ=0.0, RZ=0.0;

        for( size_t i=0; i< _caloCrystalHitsPtrVector.size(); ++i){
                std::vector<art::Ptr<CaloHit> > const& ROIds = (*_caloCrystalHitsPtrVector.at(i)).readouts();
                CaloHit const& thehit = *ROIds.at(0);

                //Get Z and R from readout
                double tZ = cg->getCrystalZByRO(thehit.id());
                double tR = cg->getCrystalRByRO(thehit.id());

                V += tR;
                W += tZ;
                VW += tR*tZ;
                VQ += pow(tR,2);

                //Multiply Z and R for the crystal's energy, which is the weight we use in that algorithm
                tZ *=(*_caloCrystalHitsPtrVector.at(i)).energyDep();
                tR*=(*_caloCrystalHitsPtrVector.at(i)).energyDep();

                Z += tZ;
                RQ += pow(tR, 2);
                R += tR;
                RZ += tZ*tR;

                //Calculate the sum of the square of the weight for derive the error of the weighted mean
        }

        int size = _caloCrystalHitsPtrVector.size();
        float m = V*W - size*VW;
        m /= (pow(V,2) - size*VQ);

        float D = RQ - pow(R,2);
        float errM = (1- pow(m,2)*D);
        if(errM>=0.0){
                errM = sqrt(errM);
        }else{
                errM = 1e-4;
        }

        return errM;
}

int CaloCluster::wSize() const{
        GeomHandle<Calorimeter> cg;
        int tmpWmin = cg->nCrystalZ(), tmpWmax = 0;
        int res = -1;
        for(size_t i=0; i<_caloCrystalHitsPtrVector.size(); ++i){

                CaloCrystalHit const & hit = (*_caloCrystalHitsPtrVector.at(i));

                std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                if(ROIds.size()<1 ) continue;

                CaloHit const & thehit = *ROIds.at(0);


                if(cg->getCrystalZByRO(thehit.id()) > tmpWmax){
                        tmpWmax = cg->getCrystalZByRO(thehit.id());
                }
                if(cg->getCrystalZByRO(thehit.id()) < tmpWmin){
                        tmpWmin = cg->getCrystalZByRO(thehit.id());
                }
        }//end loop on _caloCrystalHitsPtr

        res = tmpWmax - tmpWmin + 1;
        return res;

}

int CaloCluster::vSize() const{
        GeomHandle<Calorimeter> cg;
        int tmpVmin = cg->nCrystalR(), tmpVmax = 0;
        int res = -1;
        for(size_t i=0; i<_caloCrystalHitsPtrVector.size(); ++i){

                CaloCrystalHit const & hit = (*_caloCrystalHitsPtrVector.at(i));

                std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                if(ROIds.size()<1 ) continue;

                CaloHit const & thehit = *ROIds.at(0);


                if(cg->getCrystalRByRO(thehit.id()) > tmpVmax){
                        tmpVmax = cg->getCrystalRByRO(thehit.id());
                }
                if(cg->getCrystalRByRO(thehit.id()) < tmpVmin){
                        tmpVmin = cg->getCrystalRByRO(thehit.id());
                }
        }//end loop on _caloCrystalHitsPtr

        res = tmpVmax - tmpVmin + 1;
        return res;

}



int CaloCluster::cryEnergydepMaxRow() const{
        GeomHandle<Calorimeter> cg;
        double tmpMaxE = 0.0;
        int row = 0;
        for(size_t i=0; i<_caloCrystalHitsPtrVector.size(); ++i){

                CaloCrystalHit const & hit = (*_caloCrystalHitsPtrVector.at(i));
                if(hit.energyDep() > tmpMaxE){
                        tmpMaxE = hit.energyDep();
                        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                        if(ROIds.size()<1 ) continue;

                        CaloHit const & thehit = *ROIds.at(0);


                        row    =  cg->getCrystalRByRO(thehit.id());




                }
        }//end loop on _caloCrystalHitsPtr
        return row;
}

int CaloCluster::cryEnergydepMaxColumn() const{
        GeomHandle<Calorimeter> cg;
        double tmpMaxE = 0.0;
        int column = 0;
        for(size_t i=0; i<_caloCrystalHitsPtrVector.size(); ++i){

                CaloCrystalHit const& hit = (*_caloCrystalHitsPtrVector.at(i));
                if(hit.energyDep() > tmpMaxE){
                        tmpMaxE = hit.energyDep();
                        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                        if(ROIds.size()<1 ) continue;


                        CaloHit const & thehit = *ROIds.at(0);

                        column   =  cg->getCrystalZByRO(thehit.id());


                }
        }//end loop on _caloCrystalHitsPtr
        return column;
}

void CaloCluster::SetVaneId (int vane) {
        _vaneId = vane;
}

void CaloCluster::SetTime (float time) {
        _time = time;
}

void CaloCluster::SetCogRow ( int row) {
        _cogRow = row;
}

void CaloCluster::SetCogColumn ( int column) {
        _cogColumn = column;
}

void CaloCluster::SetEnergyDep ( float energyDep) {
        _energyDep = energyDep;
}

void CaloCluster::SetCog3Vector ( CLHEP::Hep3Vector cog3Vector) {
        _cog3Vector = cog3Vector;
}

void CaloCluster::SetCog3VectorError ( CLHEP::Hep3Vector cog3VectorErr) {
        _cog3VectorError = cog3VectorErr;
}

void CaloCluster::AddHit (CaloCrystalHitPtr &a) {

        _caloCrystalHitsPtrVector.push_back(a);

        _time *=_energyDep;
        _time += (a->time())*(a->energyDep());

        _energyDep += a->energyDep();
        _time /= _energyDep;
}

// Print the information found in this hit.
void CaloCluster::print( ostream& ost, bool doEndl ) const {

        ost << "CaloCluster :   "
                        << " vane: "          << _vaneId
                        << " time: "          << _time
                        << " timeErr"         << timeErr()<<endl
                        << " energyDep: "     << _energyDep
                        << " energyDepErr: "     << energyDepErr()<<endl
                        << " COGrow: "        << _cogRow
                        << " COGcolumn: "     << _cogColumn<<endl
                        << " COG3Vector.u: "    << _cog3Vector.x()
                        << " COG3Vector.v: "    << _cog3Vector.y()
                        << " COG3Vector.w: "    << _cog3Vector.z()<<endl
                        << " COG3VectorError.u: "    << _cog3VectorError.x()
                        << " COG3VectorError.v: "    << _cog3VectorError.y()
                        << " COG3VectorError.w: "    << _cog3VectorError.z()
                        << " size: "          << _caloCrystalHitsPtrVector.size();

        if ( doEndl ){
                ost << endl;
        }

}

} // namespace mu2e
