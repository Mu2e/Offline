//Author: S Middleton
//Date: Oct 2019
//Purpose: An improved analyzer for Cosmics

#include <iostream>
#include <string>
#include <cmath>

// Cosmic Tracks:
#include "Offline/CosmicReco/inc/CosmicTrackFit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/CosmicReco/inc/CosmicTrackMCInfo.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/CosmicReco/inc/ComboHitInfoMC.hh"

//Mu2e Data Prods:
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/DataProducts/inc/GenVector.hh"

//Utilities
#include "Offline/CosmicReco/inc/DriftFitUtils.hh"
#include "Offline/Mu2eUtilities/inc/ParametricFit.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/Mu2eUtilities/inc/BuildLinearFitMatrixSums.hh"

// Mu2e diagnostics
#include "Offline/GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

// ROOT incldues
#include "TStyle.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TTree.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TProfile.h"

//Geom
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"

using namespace std;

namespace mu2e
{
        class CosmicTrackDetails : public art::EDAnalyzer {
        public:
        struct Config{
              using Name=fhicl::Name;
              using Comment=fhicl::Comment;
              fhicl::Atom<int> diag{Name("diagLevel"), Comment("set to 1 for info"),1};
              fhicl::Atom<bool> mcdiag{Name("mcdiag"), Comment("set on for MC info"),true};
              fhicl::Atom<art::InputTag> chtag{Name("ComboHitCollection"),Comment("tag for combo hit collection")};
              fhicl::Atom<art::InputTag> tctag{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
              fhicl::Atom<art::InputTag> costag{Name("CosmicTrackSeedCollection"),Comment("tag for cosmci track seed collection")};
              fhicl::Atom<art::InputTag> mcdigistag{Name("StrawDigiMCCollection"),Comment("StrawDigi collection tag"),"makeSD"};
            };
        typedef art::EDAnalyzer::Table<Config> Parameters;

        explicit CosmicTrackDetails(const Parameters& conf);
        virtual ~CosmicTrackDetails();
        virtual void beginJob() override;
        virtual void beginRun(const art::Run& r) override;
        virtual void analyze(const art::Event& e) override;
        virtual void endJob() override;
    private:

        Config _conf;

        int  _diag;
        bool _mcdiag;
        std::ofstream outputfile;
        art::InputTag   _chtag; //ComboHits
        art::InputTag   _tctag; //Timeclusters
        art::InputTag   _costag; //Straight tracks
        art::InputTag   _mcdigistag; //MC Digis

        const ComboHitCollection* _chcol;
        const CosmicTrackSeedCollection* _coscol;
        const TimeClusterCollection* _tccol;
        const StrawDigiMCCollection* _mcdigis;
        CosmicTrackMCInfo trueinfo;

        TTree* _cosmic_tree;

        Float_t _mc_phi_angle;
        Float_t _reco_phi_angle;
        Float_t _mc_theta_angle;
        Float_t _reco_theta_angle;

        Float_t _MinuitA0;
        Float_t _MinuitA1;
        Float_t _MinuitB1;
        Float_t _MinuitB0;

        Float_t _FitCovA0;
        Float_t _FitCovA1;
        Float_t _FitCovB0;
        Float_t _FitCovB1;
        Float_t _FitCovA0A1;
        Float_t _FitCovB0B1;

        Float_t _ErrorA0;
        Float_t _ErrorA1;
        Float_t _ErrorB1;
        Float_t _ErrorB0;

        Float_t _TrueA0;
        Float_t _TrueA1;
        Float_t _TrueB1;
        Float_t _TrueB0;

        Float_t _PhiSIM;
        Float_t _ThetaSIM;
        Float_t _MomentumSIM;

        float _FitDOCAs[8129];
        float _TrueDOCAs[8129];
        float _FitTOCAs[8129];
        float _RecoAmbig[8129];
        float _TrueTimeResiduals[8129];
        float _PullsX[8129];
        float _PullsY[8129];
        Float_t _NLL;
        float _TrueAmbig[8129];
        float _hit_time[8129];
        float _hit_drift_time[8129];
        float _RecoResiduals[8129];
        float _TrueResiduals[8129];
        float _AmbigRatio[8129];
        float _DriftDistance[8129];
        float _TrueDriftDistance[8129];

        Int_t _evt;

        Int_t _nsh, _nch; // # associated straw hits / event
        Int_t _ntc; // # clusters/event
        int _nhits[8129]; // # hits used
        Int_t _n_panels; // # panels
        Int_t _n_stations; // # stations
        Int_t _n_planes; // # stations

        int _nused;
        Float_t _cluster_time;

        ProditionsHandle<Tracker> _alignedTracker_h;
        ProditionsHandle<StrawResponse> _strawResponse_h;

        Bool_t _StraightTrackInit, _StraightTrackConverged, _StraightTrackOK, _hitsOK;
        Int_t _strawid;
        vector<ComboHitInfoMC> _chinfomc;

        CosmicTrackMCInfo FitMC(const StrawDigiMCCollection*& _mcdigis);
        CosmicTrackMCInfo FillDriftMC(ComboHit const& chi, double reco_ambig, CosmicTrackMCInfo info, double t0, const Tracker* tracker);
              bool findData(const art::Event& evt);

    };

    CosmicTrackDetails::CosmicTrackDetails(const Parameters& conf) :
        art::EDAnalyzer(conf),
        _diag (conf().diag()),
        _mcdiag (conf().mcdiag()),
        _chtag (conf().chtag()),
        _tctag (conf().tctag()),
        _costag (conf().costag()),
        _mcdigistag (conf().mcdigistag())
        {}

    CosmicTrackDetails::~CosmicTrackDetails(){}

    void CosmicTrackDetails::beginJob() {

        if(_diag > 0){
                art::ServiceHandle<art::TFileService> tfs;
                _cosmic_tree=tfs->make<TTree>("cosmic_tree"," Diagnostics for Cosmic Track Fitting");
                _cosmic_tree->Branch("nused",  &_nused ,   "nused/I");
                _cosmic_tree->Branch("evt",&_evt,"evt/I");
                _cosmic_tree->Branch("nhits",&_nhits,"nhits[nused]/I");
                _cosmic_tree->Branch("StrawHitsInEvent", &_nsh, "StrawHitsInEvent/I");
                _cosmic_tree->Branch("ComboHitsInEvent", &_nch, "ComboHitsInEvent/I");
                _cosmic_tree->Branch("PanelsCrossedInEvent", &_n_panels, "PanelsCrossedInEvent/I");
                _cosmic_tree->Branch("PlanesCrossedInEvent", &_n_planes, "PlanesCrossedInEvent/I");
                _cosmic_tree->Branch("StatonsCrossedInEvent", &_n_stations, "StationsCrossedInEvent/I");
                _cosmic_tree->Branch("TimeClustersInEvent", &_ntc, "TimeClusterInEvent/I");
                _cosmic_tree->Branch("hit_time", &_hit_time, "hit_time[nused]/F");
                _cosmic_tree->Branch("hit_drit_time", &_hit_drift_time, "hit_drift_time[nused]/F");
                _cosmic_tree->Branch("hitsOK",&_hitsOK,"hitsOK/B");
                _cosmic_tree->Branch("StraightTrackInit",&_StraightTrackInit,"StraightTrackInit/B");
                _cosmic_tree->Branch("StraightTrackOK",&_StraightTrackOK,"StraightTrackOK/B");
                _cosmic_tree->Branch("StraightTrackConverged",&_StraightTrackConverged,"StraightTrackConverged/B");

                _cosmic_tree->Branch("MinuitA0",&_MinuitA0,"MinuitA0/F");
                _cosmic_tree->Branch("MinuitA1",&_MinuitA1,"MinuitA1/F");
                _cosmic_tree->Branch("MinuitB0",&_MinuitA0,"MinuitB0/F");
                _cosmic_tree->Branch("MinuitB1",&_MinuitA1,"MinuitB1/F");


                _cosmic_tree->Branch("ErrorA0",&_ErrorA0,"ErrorA0/F");
                _cosmic_tree->Branch("ErrorA1",&_ErrorA1,"ErrorA1/F");
                _cosmic_tree->Branch("ErrorB0",&_ErrorA0,"ErrorB0/F");
                _cosmic_tree->Branch("ErrorB1",&_ErrorA1,"ErrorB1/F");

                _cosmic_tree->Branch("FitCovA0",&_FitCovA0,"FitCovA0/F");
                _cosmic_tree->Branch("FitCovA1",&_FitCovA1,"FitCovA1/F");
                _cosmic_tree->Branch("FitCovB0",&_FitCovA0,"FitCovB0/F");
                _cosmic_tree->Branch("FitCovB1",&_FitCovB1,"FitCovB1/F");
                _cosmic_tree->Branch("FitCovA0A1",&_FitCovA0A1,"FitCovA0A1/F");
                _cosmic_tree->Branch("FitCovB0B1",&_FitCovB0B1,"FitCovB0B1/F");

                _cosmic_tree->Branch("RecoPhi",&_reco_phi_angle, "RecoPhi/F");
                _cosmic_tree->Branch("RecoTheta",&_reco_theta_angle, "RecoTheta/F");

                _cosmic_tree->Branch("FitDOCAs",&_FitDOCAs,"FitDOCAs[nused]/F");
                _cosmic_tree->Branch("RecoAmbig",&_RecoAmbig,"RecoAmbig[nused]/F");
                _cosmic_tree->Branch("FitTOCAs",&_FitTOCAs,"FitTOCAs[nused]/F");
                _cosmic_tree->Branch("PullsX",&_PullsX,"PullsX[nused]/F");
                _cosmic_tree->Branch("PullsY",&_PullsY,"PullsY[nused]/F");
                _cosmic_tree->Branch("RecoResiduals",&_RecoResiduals,"RecoResiduals[nused]/F");

                _cosmic_tree->Branch("DriftDistance",&_DriftDistance,"DriftDistance[nused]/F");

                if(_mcdiag ){
                        _cosmic_tree->Branch("AmbigRatio",&_AmbigRatio,"AmbigRatio[nused]/F");
                        _cosmic_tree->Branch("TrueA0",&_TrueA0,"TrueA0/F");
                        _cosmic_tree->Branch("TrueA1",&_TrueA1,"TrueA1/F");
                        _cosmic_tree->Branch("TrueB0",&_TrueA0,"TrueB0/F");
                        _cosmic_tree->Branch("TrueB1",&_TrueA1,"TrueB1/F");

                        _cosmic_tree->Branch("TruePhi",&_mc_phi_angle, "TruePhi/F");
                        _cosmic_tree->Branch("TrueTheta",&_mc_theta_angle, "TrueTheta/F");
                        _cosmic_tree->Branch("TrueDOCAs",&_TrueDOCAs,"TrueDOCAs[nused]/F");
                        _cosmic_tree->Branch("TrueTimeResiduals",&_TrueTimeResiduals,"TrueTimeResiduals[nused]/F");
                        _cosmic_tree->Branch("TrueAmbig",&_TrueAmbig,"TrueAmbig[nused]/F");
                        _cosmic_tree->Branch("TrueResiduals",&_TrueResiduals,"TrueResiduals[nused]/F");
                        _cosmic_tree->Branch("TrueDriftDistance",&_TrueDriftDistance,"TrueDriftDistance[nused]/F");
                        _cosmic_tree->Branch("Momentum", &_MomentumSIM, "Momentum/F");
                }

        }
      }

    void CosmicTrackDetails::beginRun(const art::Run& run){}

    void CosmicTrackDetails::analyze(const art::Event& event) {

        const Tracker *tracker = _alignedTracker_h.getPtr(event.id()).get();
        StrawResponse const& srep = _strawResponse_h.get(event.id());

        _evt = event.id().event();

        if(!findData(event))
                throw cet::exception("RECO")<<"No Time Clusters in event"<< endl;

        _ntc = _tccol->size();
        _nch = _chcol->size();

        for(size_t itc=0; itc<_tccol->size();++itc){
                TimeCluster tc = (*_tccol)[itc];
                _cluster_time =  tc._t0._t0;

        }

        for(size_t ist = 0;ist < _coscol->size(); ++ist){

                CosmicTrackSeed sts =(*_coscol)[ist];
                CosmicTrack st = sts._track;
                double t0 = sts._t0.t0();
                TrkFitFlag const& status = sts._status;

                if (!status.hasAllProperties(TrkFitFlag::helixOK) ){ continue; }
                if(st.converged == false or st.minuit_converged  == false) { continue; }

                std::vector<int> panels, planes, stations;

                _reco_phi_angle=acos(st.FitEquation.Dir.x()/st.FitEquation.Dir.Mag2());
                _reco_theta_angle=acos(st.FitEquation.Dir.y()/sqrt(st.FitEquation.Dir.Mag2()));
                _MinuitA0=(st.MinuitParams.A0);
                _MinuitA1=(st.MinuitParams.A1);
                _MinuitB1=(st.MinuitParams.B1);
                _MinuitB0=(st.MinuitParams.B0);

                _ErrorA0=(st.MinuitParams.deltaA0);
                _ErrorA1=(st.MinuitParams.deltaA1);
                _ErrorB1=(st.MinuitParams.deltaB1);
                _ErrorB0=(st.MinuitParams.deltaB0);

                _FitCovA0  = st.MinuitParams.Covarience.sigA0;
                _FitCovA1 = st.MinuitParams.Covarience.sigA1;
                _FitCovB0 = st.MinuitParams.Covarience.sigB0;
                _FitCovB1 = st.MinuitParams.Covarience.sigB1;
                _FitCovA0A1 = st.MinuitParams.Covarience.sigA0A1;
                _FitCovB0B1 = st.MinuitParams.Covarience.sigB0B1;

               if(_mcdiag){

                        trueinfo = FitMC(_mcdigis);

                        _mc_phi_angle=(trueinfo.TruePhi);
                        _mc_theta_angle=(trueinfo.TrueTheta);

                        _TrueA1=(trueinfo.TrueFitEquation.Dir.X());
                        _TrueB1=(trueinfo.TrueFitEquation.Dir.Y());
                        _TrueA0=(trueinfo.TrueFitEquation.Pos.X());
                        _TrueB0=(trueinfo.TrueFitEquation.Pos.Y());
                }
                for(size_t i=0; i<sts._straw_chits.size();i++){
                        ComboHit const& chit = sts._straw_chits[i];
                        _nused = i;

                        _nhits[_nused] = chit.nStrawHits();
                        panels.push_back(chit.strawId().panel());
                        planes.push_back(chit.strawId().plane());
                        stations.push_back(chit.strawId().station());
                        _hit_time[_nused] = chit.time();
                        _hit_drift_time[_nused] = chit.driftTime();

                        if(_mcdiag){
                                trueinfo = FitMC(_mcdigis);
                                trueinfo = FillDriftMC(chit, _RecoAmbig[_nused], trueinfo, t0, tracker);
                                    double truedoca = DriftFitUtils::GetTestDOCA(chit, trueinfo.TrueFitEquation.Pos.X(), trueinfo.TrueFitEquation.Dir.X(), trueinfo.TrueFitEquation.Pos.Y(),trueinfo.TrueFitEquation.Dir.Y(),  tracker);

                                double fitdoca = DriftFitUtils::GetTestDOCA(chit, st.MinuitParams.A0,st.MinuitParams.A1, st.MinuitParams.B0, st.MinuitParams.B1,  tracker);

                                double recoambig = DriftFitUtils::GetAmbig(chit, st.MinuitParams.A0,st.MinuitParams.A1, st.MinuitParams.B0, st.MinuitParams.B1, tracker);

                                double res = DriftFitUtils::GetRPerp(srep, chit, st.MinuitParams.A0,st.MinuitParams.A1, st.MinuitParams.B0, st.MinuitParams.B1, tracker);

                                double trueambig = DriftFitUtils::GetAmbig(chit, trueinfo.TrueFitEquation.Pos.X(), trueinfo.TrueFitEquation.Dir.X(), trueinfo.TrueFitEquation.Pos.Y(),trueinfo.TrueFitEquation.Dir.Y(), tracker);

                                double trueres = DriftFitUtils::GetRPerp(srep,chit, trueinfo.TrueFitEquation.Pos.X(), trueinfo.TrueFitEquation.Dir.X(), trueinfo.TrueFitEquation.Pos.Y(),trueinfo.TrueFitEquation.Dir.Y(),  tracker);

                                double truedriftdistance = DriftFitUtils::GetDriftDistance(srep,chit, trueinfo.TrueFitEquation.Pos.X(), trueinfo.TrueFitEquation.Dir.X(), trueinfo.TrueFitEquation.Pos.Y(),trueinfo.TrueFitEquation.Dir.Y(),  tracker);

                                double driftdistance = DriftFitUtils::GetDriftDistance(srep, chit, st.MinuitParams.A0,st.MinuitParams.A1, st.MinuitParams.B0, st.MinuitParams.B1, tracker);

                                if(fitdoca < 2.5 and truedoca < 2.5 and abs(trueinfo.TrueFitEquation.Pos.X() ) < 5000 and abs(trueinfo.TrueFitEquation.Pos.Y()) < 5000 and abs(trueinfo.TrueFitEquation.Dir.X())< 5 and abs(trueinfo.TrueFitEquation.Dir.Y())< 5){
                                        _RecoAmbig[_nused] = recoambig;
                                        _FitDOCAs[_nused] = recoambig*fitdoca;
                                        _RecoResiduals[_nused] = res ;
                                        _FitTOCAs[_nused] = fitdoca / 0.0625;
                                        _TrueDOCAs[_nused] = trueambig*truedoca;
                                        _TrueTimeResiduals[_nused] = truedoca / 0.0625;
                                        _TrueAmbig[_nused] = trueambig;
                                        _TrueResiduals[_nused] = trueres ;
                                        _MomentumSIM = trueinfo.TrueMomentum;
                                        _AmbigRatio[_nused] = trueambig/recoambig;
                                        _TrueDriftDistance[_nused] = truedriftdistance;
                                        _DriftDistance[_nused] = driftdistance;
                           }
                }
                }
                _hitsOK = status.hasAllProperties(TrkFitFlag::hitsOK);
                if(status.hasAllProperties(TrkFitFlag::Straight)){
                        _StraightTrackOK = status.hasAllProperties(TrkFitFlag::helixOK);
                        _StraightTrackConverged = status.hasAllProperties(TrkFitFlag::helixConverged);
                        _StraightTrackInit = status.hasAllProperties(TrkFitFlag::circleInit);
                }

                _n_panels = std::set<float>( panels.begin(), panels.end() ).size();
                _n_planes = std::set<float>( planes.begin(), planes.end() ).size();
                _n_stations = std::set<float>( stations.begin(), stations.end() ).size();
                if(st.minuit_converged == true){
                         _cosmic_tree->Fill();
               }
        }

}


    void CosmicTrackDetails::endJob() {}

    CosmicTrackMCInfo CosmicTrackDetails::FitMC(const StrawDigiMCCollection*& _mcdigis){

        ::BuildLinearFitMatrixSums S;
        CosmicTrackMCInfo TrackTrueInfo;

        StrawDigiMC hitP1;
        StrawDigiMC first = (*_mcdigis)[0];

        auto const& spmcp0= first.earlyStrawGasStep();
        XYZVectorF pos0(spmcp0->position().x(), spmcp0->position().y(), spmcp0->position().z());
        XYZVectorF dir0(spmcp0->momentum().x(), spmcp0->momentum().y(), spmcp0->momentum().z());

        for(size_t ich = 0;ich < _mcdigis->size(); ++ich){
            hitP1 = (*_mcdigis)[ich];

            auto const& spmcp= hitP1.earlyStrawGasStep();
            XYZVectorF posN(spmcp->position().x(), spmcp->position().y(), spmcp->position().z());

            XYZVectorF ZPrime = (spmcp->momentum().Unit());

            TrackAxes TrueAxes = ParametricFit::GetTrackAxes(ZPrime);
            TrackTrueInfo.TrueTrackCoordSystem = (TrueAxes);

            XYZVectorF point(posN.x(), posN.y(), posN.z());
            XYZVectorF X(1,0,0);
            XYZVectorF Y(0,1,0);
            XYZVectorF Z(0,0,1);
            S.addPoint( point, X,Y,Z, 1,1);

        }

        TrackParams RawTrueParams(S.GetAlphaX()[0][0], S.GetAlphaX()[1][0], S.GetAlphaY()[0][0], S.GetAlphaY()[1][0]);

        XYZVectorF TruePos(S.GetAlphaX()[0][0], S.GetAlphaY()[0][0], 0);

        XYZVectorF TrueDir(S.GetAlphaX()[1][0], S.GetAlphaY()[1][0], 1);
        TrueDir = TrueDir.Unit();
        TrueDir = TrueDir/TrueDir.Z();

        pos0.SetX(pos0.X()-(dir0.X()*pos0.Z()/dir0.Z()));
        pos0.SetY(pos0.Y()-(dir0.Y()*pos0.Z()/dir0.Z()));
        pos0.SetZ(pos0.Z()-(dir0.Z()*pos0.Z()/dir0.Z()));
        dir0 = dir0/dir0.Z();

        TrackEquation TrueTrack(pos0, dir0);

        TrackTrueInfo.TrueFitEquation = (TrueTrack);
        TrackTrueInfo.TruePhi =(atan(TrueDir.y()/TrueDir.x()));
        TrackTrueInfo.TrueTheta = (acos(TrueDir.x()/sqrt(TrueDir.Mag2())));


        int n{-1};
        for(size_t ich = 0;ich < _mcdigis->size(); ++ich){
                StrawDigiMC digimc= (*_mcdigis)[ich];
                auto const& sim_cal = digimc.strawGasStep(mu2e::StrawEnd::cal)->simParticle();
                auto const& sim_hv  = digimc.strawGasStep(mu2e::StrawEnd::cal)->simParticle();
                if ( sim_cal == sim_hv ){

                  if (n == 0){
                        TrackTrueInfo.TrueMomentum = sqrt(digimc.strawGasStep(mu2e::StrawEnd::cal)->momentum().mag2());;
                            TrackTrueInfo.TrueThetaSIM = acos(digimc.strawGasStep(mu2e::StrawEnd::cal)->momentum().z()/TrackTrueInfo.TrueMomentum );
                            TrackTrueInfo.TruePhiSIM = atan(digimc.strawGasStep(mu2e::StrawEnd::cal)->momentum().y()/digimc.strawGasStep(mu2e::StrawEnd::cal)->momentum().x());
                 }
                 n++;
      }
    }

    return TrackTrueInfo;
   }

    CosmicTrackMCInfo CosmicTrackDetails::FillDriftMC(ComboHit const& chit, double RecoAmbig, CosmicTrackMCInfo info, double t0,  const Tracker* tracker){
        double true_doca = DriftFitUtils::GetTestDOCA(chit, info.TrueFitEquation.Pos.X(), info.TrueFitEquation.Dir.X(), info.TrueFitEquation.Pos.Y(),info.TrueFitEquation.Dir.Y(),  tracker);
        double trueambig = DriftFitUtils::GetAmbig(chit, info.TrueFitEquation.Pos.X(), info.TrueFitEquation.Dir.X(), info.TrueFitEquation.Pos.Y(),info.TrueFitEquation.Dir.Y(),  tracker);
        info.Ambig.push_back(trueambig);
        info.TrueDOCA.push_back(true_doca);

        return info;
    }

    bool CosmicTrackDetails::findData(const art::Event& evt){
        _chcol = 0;
        _tccol = 0;
        _coscol = 0;
        auto chH = evt.getValidHandle<ComboHitCollection>(_chtag);
        _chcol = chH.product();
        auto tcH = evt.getValidHandle<TimeClusterCollection>(_tctag);
        _tccol =tcH.product();
        auto stH = evt.getValidHandle<CosmicTrackSeedCollection>(_costag);
        _coscol =stH.product();
        if(_mcdiag){
                _mcdigis=0;
                auto mcdH = evt.getValidHandle<StrawDigiMCCollection>(_mcdigistag);
                _mcdigis = mcdH.product();
        }
        return _chcol != 0 && _tccol!=0 && _coscol !=0 && (_mcdigis != 0 || !_mcdiag);
       }

}

using mu2e::CosmicTrackDetails;
DEFINE_ART_MODULE(CosmicTrackDetails)
