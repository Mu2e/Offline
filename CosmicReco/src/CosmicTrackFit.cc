// Author: S Middleton
// Date: March 2019
// Purpose: Holds functions for the fitting of Cosmic Tracks in tracker

// Mu2e Cosmics:
#include "Offline/CosmicReco/inc/CosmicTrackFit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"

// art
#include "canvas/Persistency/Common/Ptr.h"

//Mu2e General:
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

//Fitting
#include "Offline/Mu2eUtilities/inc/ParametricFit.hh"
#include "Offline/Mu2eUtilities/inc/BuildLinearFitMatrixSums.hh"
#include "Offline/CosmicReco/inc/MinuitDriftFitter.hh"
#include "Offline/CosmicReco/inc/DriftFitUtils.hh"
#include "Offline/DataProducts/inc/GenVector.hh"

//ROOT:
#include "TMatrixD.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/DistFunc.h"

// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

//String:
#include <string>

using namespace std;
using namespace boost::accumulators;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;
using namespace ROOT::Math;

std::vector<double> TrackerConstraints(1500, 1500);

struct ycomp {
    bool operator()(XYZVectorF const& p1, XYZVectorF p2) { return p1.y() > p2.y(); }
  };


namespace mu2e
{
    CosmicTrackFit::CosmicTrackFit(const Config& conf) :
        _Npara (conf.Npara()),
        _diag (conf.diag()),
        _debug  (conf.debug()),
        _dontuseflag (conf.dontuseflag()),
        _minnsh   (conf.minnsh()),
        _minnch  (conf.minnch()),
        _n_outliers (conf.n_outliers()),
        _maxniter (conf.maxniter()),
        _max_seed_chi2 (conf.max_seed_chi2()),
        _max_chi2_change (conf.max_chi2_change()),
        _max_position_deviation (conf.max_position_deviation()),
        _maxHitDOCA (conf.maxHitDOCA()),
        _maxLogL (conf.maxLogL()),
        _gaussTres (conf.gaussTres()),
        _maxTres (conf.maxTres()),
        _maxd (conf.maxd()),
        _maxpull (conf.maxpull()),
        _useTSeedDirection (conf.UseTSeedDirection())
        {}

    /* ---------------Initialize Fit----------------//
    //----------------------------------------------*/
    bool CosmicTrackFit::initCosmicTrack(const char* title, CosmicTrackSeed& tseed, ComboHitCollection &combohits) {

        bool is_ok(false);
        RunFitChi2(title, tseed, combohits);
        is_ok = tseed._status.hasAllProperties(TrkFitFlag::helixOK);
        return is_ok;
  }

  std::vector<XYZVectorF> SortPoints(std::vector<XYZVectorF> pointY){
        std::vector<XYZVectorF> sortedPoints;
          std::sort(pointY.begin(), pointY.end(),ycomp());
          for (unsigned i=0; i<pointY.size(); ++i) {
                      sortedPoints.push_back(pointY[i]);
              }
              return sortedPoints;
   }

  /*-------------Line Direction-------------------------//
    Range of methods for finding track directions - some are redundent-check!
  //----------------------------------------------*/

  XYZVectorF CosmicTrackFit::InitLineDirection(const ComboHit *ch0, const ComboHit *chN) {
      float tx = chN->pos().x() - ch0->pos().x();
      float ty = chN->pos().y() - ch0->pos().y();
      float tz = chN->pos().z() - ch0->pos().z();
      XYZVectorF track(tx,ty,tz);

      return track.Unit();
    }


    XYZVectorF CosmicTrackFit::LineDirection(double a1, double b1, const ComboHit *ch0, const ComboHit *chN, XYZVectorF ZPrime) {
      XYZVectorF track(a1,b1,1);
      return track.Unit();
    }

    XYZVectorF CosmicTrackFit::ConvertPointToDetFrame(XYZVectorF vec){
        Hep3Vector vec1(vec.x(),vec.y(),vec.z());
        GeomHandle<DetectorSystem> det;
        Hep3Vector vec2 = det->toDetector(vec1);
        XYZVectorF XYZ(vec2.x(), vec2.y(), vec2.z());
        return XYZ;

    }

    void CosmicTrackFit::BeginFit(const char* title, CosmicTrackSeed &tseed, art::Event const& event, ComboHitCollection const& chcol, std::vector<StrawHitIndex> &panelHitIdxs){
      ComboHitCollection combohits;
      combohits.setParent(chcol.parent());
      for (size_t i=0;i<panelHitIdxs.size();i++){
        combohits.push_back(chcol[panelHitIdxs[i]]);
      }
      tseed._status.clear(TrkFitFlag::helixOK);
      bool init(false);
      if (!tseed._status.hasAllProperties(TrkFitFlag::circleInit)) {
        init = true;
        if (initCosmicTrack( title, tseed, combohits)){
          tseed._status.merge(TrkFitFlag::circleInit);
        }
        else {
          FillTrackHitCollections(tseed, event, chcol, panelHitIdxs);
          return;
        }
      }
      if (!init)RunFitChi2(title, tseed, combohits);

      FillTrackHitCollections(tseed, event, chcol, panelHitIdxs);

      return;
    }

    void CosmicTrackFit::FillTrackHitCollections(CosmicTrackSeed &tseed, art::Event const& event, ComboHitCollection const& chcol, std::vector<StrawHitIndex> &panelHitIdxs){
      ComboHitCollection::CHCIter chids;
      chcol.fillComboHits( panelHitIdxs, chids);
      for (auto const& it : chids){
        tseed._straw_chits.push_back(*it);
      }
      for(size_t ich= 0; ich<tseed._straw_chits.size(); ich++){
        std::vector<StrawHitIndex> shitids;
        tseed._straw_chits.fillStrawHitIndices(ich, shitids);

        /*
        for(auto const& ids : shitids){
          size_t    istraw   = (ids);
          TrkStrawHitSeed tshs;
          tshs._index  = istraw;
          tshs._t0 = tseed._t0;
          tseed._trkstrawhits.push_back(tshs);
        }
        */
      }

    }

   void CosmicTrackFit::RunFitChi2(const char* title, CosmicTrackSeed& tseed, ComboHitCollection &combohits) {
        CosmicTrack* track = &tseed._track;
        tseed._status.merge(TrkFitFlag::helixOK);
        tseed._status.merge(TrkFitFlag::helixConverged);
        FitAll(title, tseed, combohits, track);

    }



    void CosmicTrackFit::FitAll(const char* title, CosmicTrackSeed &tseed, ComboHitCollection &combohits, CosmicTrack* cosmictrack){

     ::BuildLinearFitMatrixSums S;
     ComboHit   *hitP1(0), *hitP2(0);
     size_t nHits (combohits.size());
     int DOF = (nHits);// - (_Npara);
     const ComboHit* ch0 = &combohits[0];
     const ComboHit* chN = &combohits[combohits.size()-1];

     //Step 1: Get Initial Estimate of track direction
     XYZVectorF ZPrime;
     if (_useTSeedDirection){
       ZPrime = tseed._track.FitEquation.Dir;
     }else{
       ZPrime = InitLineDirection(ch0, chN);
     }
     std::vector<XYZVectorF> AxesList = ParametricFit::GetAxes(ZPrime);
     TrackAxes InitAxes = ParametricFit::GetTrackAxes(ZPrime);

    //Step 2: Loop over hits and get track parameters based on above estimated track direction
    for (size_t f1=0; f1<nHits; ++f1){
      hitP1 = &combohits[f1];
      if (!use_hit(*hitP1) && hitP1->nStrawHits() < _minnsh)  continue;
      XYZVectorF point(hitP1->pos().x(),hitP1->pos().y(),hitP1->pos().z());
      std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP1, InitAxes._XDoublePrime, InitAxes._YDoublePrime);
      S.addPoint(point, InitAxes._XDoublePrime, InitAxes._YDoublePrime, InitAxes._ZPrime, ErrorsXY[0], ErrorsXY[1]);
     }

      //Step 3: Get the first estmiate of track parameters, get a updated track direction vector from these parameters
      double a0 = S.GetAlphaX()[0][0];
      double a1 = S.GetAlphaX()[1][0];
      double b0 = S.GetAlphaY()[0][0];
      double b1 = S.GetAlphaY()[1][0];
      XYZVectorF Direction(a1,b1,1);
      XYZVectorF UpdatedTrackDirection =Direction.Unit();

      //Step 4: Update axes and store them as initial axes for plotting
      ZPrime = UpdatedTrackDirection;
      TrackAxes Axes = ParametricFit::GetTrackAxes(ZPrime);
      TrackParams InitParams(a0,a1,b0,b1);

      cosmictrack->SetInitParams(InitParams);
      cosmictrack->SetInitCoordSystem(Axes);
      cosmictrack->SetFitCoordSystem(Axes);
      //Step 5: Loop for initial diagnostics
      if(_debug>0){
        cosmictrack->set_initchisq_dofY(S.GetChi2Y()/abs(DOF));
        cosmictrack->set_initchisq_dofX(S.GetChi2X()/abs(DOF));
        cosmictrack->set_initchisq_dof(S.GetTotalChi2()/abs(DOF));
      }

      //Step 6: Begin iteration for finding the best track fit possible.
      unsigned niter(0);
      ::BuildLinearFitMatrixSums S_niteration;
      bool converged = false;
      CosmicTrack* BestTrack = cosmictrack;
      double chi2_best_track = 10000000;//chosen arbitary high number

      //Step 7: find best prev chi2
      while(niter < _maxniter && converged==false){
              niter +=1;
         double previous_chi2, changed_chi2;
              if (niter == 1) {
                       previous_chi2 = S.GetTotalChi2()/abs(DOF);
                       chi2_best_track = previous_chi2;
       }
              else {
                       previous_chi2 = S_niteration.GetTotalChi2()/abs(DOF);
              }
             //Step 8: Remove previous sums and use previously stored axes from updated last iteration:
             S_niteration.clear();
             Axes = ParametricFit::GetTrackAxes(cosmictrack->FitParams.Direction());
             cosmictrack->set_niter(niter );
             for (size_t f4=0; f4 < nHits; ++f4){

                   hitP2 = &combohits[f4];
                    if (((!use_hit(*hitP2) ) && (hitP2->nStrawHits() < _minnsh) )) continue;
              XYZVectorF point(hitP2->pos().x(),hitP2->pos().y(),hitP2->pos().z());
              std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP2, Axes._XDoublePrime, Axes._YDoublePrime);
              S_niteration.addPoint(point, Axes._XDoublePrime, Axes._YDoublePrime, Axes._ZPrime, ErrorsXY[0],ErrorsXY[1]);
              }
          //Get new parameters
               a0 = S_niteration.GetAlphaX()[0][0];
          a1 = S_niteration.GetAlphaX()[1][0];
          b0 = S_niteration.GetAlphaY()[0][0];
          b1 = S_niteration.GetAlphaY()[1][0];

              //Step 9: Get new direction:
           //XYZVectorF DirectionSecond(a1,b1,1);
           //XYZVectorF UpdatedTrackDirectionSecond = DirectionSecond.Unit();

            //Step 10 - Get New Axes:
            //Axes = ParametricFit::GetTrackAxes(UpdatedTrackDirectionSecond);
           cosmictrack->SetFitCoordSystem(Axes);
           //Step 11: Update Parameters:
           TrackParams FitParams(a0,a1,b0,b1);
           cosmictrack->SetFitParams(FitParams);

           //Step 12: Update Chi2:
           if(_debug > 0){
                      cosmictrack->Diag.FinalChiTot = (S_niteration.GetTotalChi2()/abs(DOF));
                   cosmictrack->Diag.FinalChiY = S_niteration.GetChi2Y()/abs(DOF);
                   cosmictrack->Diag.FinalChiX = S_niteration.GetChi2X()/abs(DOF);
           }

           if(abs((a0) - (InitParams.A0)) > _max_position_deviation or abs((b0) - (InitParams.B0)) > _max_position_deviation){continue;}

           //Step13: Check if Chi2 is improved:, If so call this "BestTrack"
           float updated_chi2 = S_niteration.GetTotalChi2()/abs(DOF);
           changed_chi2 = chi2_best_track - updated_chi2;

           if (updated_chi2 < chi2_best_track  ){

                //BestTrack->clear_diag();
                TrackParams FitParams(a0,a1,b0,b1);
                   BestTrack->SetFitParams(FitParams);
                BestTrack->SetFitCoordSystem(Axes);

                      XYZVectorF EndTrackPosition(BestTrack->FitParams.A0,BestTrack->FitParams.B0,0);
                      TrackEquation EndTrack(EndTrackPosition, BestTrack->FitCoordSystem._ZPrime);
                      BestTrack->SetFitEquation(EndTrack);

                   //Step 14: Save this Chi2:
                if(_debug > 0){
                        BestTrack->Diag.FinalChiTot = (S_niteration.GetTotalChi2()/abs(DOF));
                           BestTrack->Diag.FinalChiY = S_niteration.GetChi2Y()/abs(DOF);
                        BestTrack->Diag.FinalChiX = S_niteration.GetChi2X()/abs(DOF);

                }
                chi2_best_track = updated_chi2;
                cosmictrack=BestTrack;
              }

         if( abs(changed_chi2) < _max_chi2_change ){
                  converged = true;
                  cosmictrack->converged = true;
               }

           if(niter == _maxniter && converged ==false ){
                     tseed._status.clear(TrkFitFlag::helixOK);
                     tseed._status.clear(TrkFitFlag::helixConverged);
                     continue;
                 }

   }
   cosmictrack=BestTrack;
   if(cosmictrack->converged and  _diag > 0){
         for (size_t f5=0; f5<nHits; ++f5){
                     hitP2 = &combohits[f5];
                if (((!use_hit(*hitP2) ) && (hitP2->nStrawHits() < _minnsh) )) continue;
                XYZVectorF point(hitP2->pos().x(),hitP2->pos().y(),hitP2->pos().z());
                XYZVectorF point_prime(point.Dot(BestTrack->FitCoordSystem._XDoublePrime), point.Dot(BestTrack->FitCoordSystem._YDoublePrime), point.Dot(BestTrack->FitCoordSystem._ZPrime));

                std::vector<double> ErrorsXY = ParametricFit::GetErrors(hitP1, BestTrack->FitCoordSystem._XDoublePrime, BestTrack->FitCoordSystem._YDoublePrime);

                      float newRx = ParametricFit::GetResidualX(BestTrack->FitParams.A0, BestTrack->FitParams.A1, point_prime);
                float newRy = ParametricFit::GetResidualY(BestTrack->FitParams.B0, BestTrack->FitParams.B1, point_prime);

                if(abs(newRx)/ErrorsXY[0] > _maxpull or abs(newRy)/ErrorsXY[1] > _maxpull ){ cosmictrack->converged = false;}

                BestTrack->SetCovarience(S_niteration.GetCovX()[0][0],S_niteration.GetCovX()[0][1], S_niteration.GetCovX()[1][0], S_niteration.GetCovX()[1][1], S_niteration.GetCovY()[0][0], S_niteration.GetCovY()[0][1], S_niteration.GetCovY()[1][0], S_niteration.GetCovY()[1][1]);

                      cosmictrack = BestTrack;

              }
    }
     S.clear();

     if(cosmictrack->converged == true){

          ConvertFitToDetectorFrame(cosmictrack->FitCoordSystem, cosmictrack->FitParams.Position(), cosmictrack->FitParams.Direction(), cosmictrack, true, false);

     }

   }

/*------------Translate fit back into XYZ and/or the detector frame------//
Using matrices to ctransform from local to global coordinates
//-----------------------------------------------------------------------*/
void CosmicTrackFit::ConvertFitToDetectorFrame(TrackAxes axes, XYZVectorF Position, XYZVectorF Direction, CosmicTrack* cosmictrack, bool seed, bool det){
        TMatrixD A(3,3);
        A[0][0] = axes._XDoublePrime.X();
        A[0][1] = axes._YDoublePrime.X();
        A[0][2] = axes._ZPrime.X();

        A[1][0] = axes._XDoublePrime.Y();
        A[1][1] = axes._YDoublePrime.Y();
        A[1][2] = axes._ZPrime.Y();

        A[2][0] = axes._XDoublePrime.Z();
        A[2][1] = axes._YDoublePrime.Z();
        A[2][2] = axes._ZPrime.Z();

        TMatrixD P(3,1);
        P[0][0] = Position.X();
        P[1][0] = Position.Y();
        P[2][0] = Position.Z();

        TMatrixD D(3,1);
        D[0][0] = Direction.X();
        D[1][0] = Direction.Y();
        D[2][0] = Direction.Z();

        TMatrixD PXYZ(A*P);
        TMatrixD DXYZ(A*D);

        DXYZ[0][0] = DXYZ[0][0]/ DXYZ[2][0];
        DXYZ[1][0] = DXYZ[1][0]/ DXYZ[2][0];
        DXYZ[2][0] = DXYZ[2][0]/ DXYZ[2][0];

        PXYZ[0][0] = PXYZ[0][0] - DXYZ[0][0]*PXYZ[2][0]/ DXYZ[2][0];
        PXYZ[1][0] = PXYZ[1][0]- DXYZ[1][0]*PXYZ[2][0]/ DXYZ[2][0];
        PXYZ[2][0] = PXYZ[2][0]- DXYZ[2][0]*PXYZ[2][0]/ DXYZ[2][0];
        /*
        //TMatrixD sigmaPos(3,1);
        TMatrixD Cov(3,3);
        TMatrixD At(A);
        At.T();
        for(unsigned i =0; i< 3; i++){
                for(unsigned j =0 ; j <3; j++){

                    if(i==0 and j==0) {Cov[i][j] =  cosmictrack->FitParams.Covarience.sigA0; }
                    if(i==1 and j ==1) { Cov[i][j] = cosmictrack->FitParams.Covarience.sigA1;}
                    else {Cov[i][j] =0;}
                }
        }
        TMatrixD CovXYZ(A*Cov*At);
        for(unsigned i=0; i<3 ; i++){
                for(unsigned j=0; j<3; j++){
                        cout<<"i= "<<i<<" j= "<<j<<" "<<Cov[i][j]<<"=======> "<<CovXYZ[i][j]<<endl;

                }

        }
        */
        TMatrixD sigmaPos(3,1);
        sigmaPos[0][0] = cosmictrack->FitParams.Covarience.sigA0;
        sigmaPos[1][0] = cosmictrack->FitParams.Covarience.sigB0;
        sigmaPos[2][0] = 0;
        TMatrixD sigmaPosXYZ(A*sigmaPos);

        TMatrixD sigmaDir(3,1);
        sigmaDir[0][0] = cosmictrack->FitParams.Covarience.sigA1;
        sigmaDir[1][0] = cosmictrack->FitParams.Covarience.sigB1;
        sigmaDir[2][0] = 0;
        TMatrixD sigmaDirXYZ(A*sigmaDir);

        XYZVectorF Pos(PXYZ[0][0], PXYZ[1][0], PXYZ[2][0]);
        XYZVectorF Dir(DXYZ[0][0], DXYZ[1][0] , DXYZ[2][0]);

        if(seed == true){//is this the local coords at start?
                if (det == true){ // is this detector frame?
                        XYZVectorF PosInDet = ConvertPointToDetFrame(Pos);
                        TrackEquation XYZTrack(PosInDet, Dir);
                        cosmictrack->SetFitEquation(XYZTrack);


                } else {
                        TrackEquation XYZTrack(Pos, Dir);
                        cosmictrack->SetFitEquation(XYZTrack);

                }}
        else{
                TrackEquation XYZTrack(Pos, Dir);
                cosmictrack->SetMinuitEquation(XYZTrack);
        }

    }


    bool CosmicTrackFit::goodTrack(CosmicTrack& track)
    {
        if(track.Diag.FinalChiTot < _max_seed_chi2) return true;
        else return false;
    }

    bool CosmicTrackFit::use_hit(const ComboHit& thit) const
    {
        return (!thit._flag.hasAnyProperty(_dontuseflag));
    }

    bool CosmicTrackFit::use_track(double track_length) const
    {
        return (track_length > _maxd) ? false : true ;
    }

    void CosmicTrackFit::DriftFit(CosmicTrackSeed& tseed, StrawResponse const& _srep ){

      FitResult endresult = MinuitDriftFitter::DoFit(_diag, tseed, _srep, _tracker, _maxHitDOCA, _minnch, _maxLogL, _gaussTres, _maxTres);

      tseed._track.MinuitParams.A0 =  endresult.bestfit[0];//a0
      tseed._track.MinuitParams.A1 =  endresult.bestfit[1];//a1
      tseed._track.MinuitParams.B0 =  endresult.bestfit[2];//b0
      tseed._track.MinuitParams.B1 =  endresult.bestfit[3];//b1
      tseed._track.MinuitParams.T0 =  endresult.bestfit[4];//t0

      tseed._track.MinuitParams.deltaA0 =  endresult.bestfiterrors[0];//erra0
      tseed._track.MinuitParams.deltaA1 =  endresult.bestfiterrors[1];//erra1
      tseed._track.MinuitParams.deltaB0 =  endresult.bestfiterrors[2];//errb0
      tseed._track.MinuitParams.deltaB1 =  endresult.bestfiterrors[3];//errb1
      tseed._track.MinuitParams.deltaT0 =  endresult.bestfiterrors[4];//errt0

      if(endresult.bestfitcov.size() !=0 ){
        TrackCov Cov(endresult.bestfitcov[0], 0., 0., endresult.bestfitcov[1], endresult.bestfitcov[2],0.,0., endresult.bestfitcov[3]);
        tseed._track.MinuitParams.Covarience = Cov;
      }
      if(endresult.NLL !=0){ tseed._track.minuit_converged = true; }

      XYZVectorF X(1,0,0);
      XYZVectorF Y(0,1,0);
      XYZVectorF Z(0,0,1);

      TrackAxes XYZ(X,Y,Z);
      tseed._track.MinuitCoordSystem = XYZ;
      tseed._track.MinuitEquation.Pos = XYZVectorF(tseed._track.MinuitParams.A0,tseed._track.MinuitParams.B0,0);
      tseed._track.MinuitEquation.Dir = XYZVectorF(tseed._track.MinuitParams.A1,tseed._track.MinuitParams.B1,1);

      unsigned int n_outliers = 0;
      if(endresult.FullFitEndTimeResiduals.size() >0){
        for(unsigned i = 0; i< endresult.FullFitEndTimeResiduals.size()-1; i++){
          if( endresult.FullFitEndTimeResiduals[i] > _maxTres or isnan(endresult.FullFitEndTimeResiduals[i])==true){
            n_outliers +=1;
            tseed._straw_chits[i]._flag.merge(StrawHitFlag::outlier);

          }
        }

      }
      if( n_outliers  > _n_outliers) {
        tseed._track.minuit_converged = false;
      }

    }

}//end namespace
