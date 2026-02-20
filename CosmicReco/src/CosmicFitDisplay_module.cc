//TGeo:
#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TGeoMatrix.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

//Mu2e Data Prods:
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
// Mu2e Utilities
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/Mu2eUtilities/inc/ParametricFit.hh"

//Mu2e Tracker Geom:
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"

// Mu2e diagnostics
#include "Offline/GeneralUtilities/inc/ParameterSetHelpers.hh"

//Cosmics:
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

// ROOT incldues
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TH2D.h"
#include "TF1.h"
#include "TH3D.h"
#include "Rtypes.h"
#include "TApplication.h"
#include "TArc.h"
#include "TTUBE.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TPolyMarker.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TText.h"
#include "TRotMatrix.h"


using namespace std;

namespace mu2e
{
  class CosmicFitDisplay : public art::EDAnalyzer {
    public:
        struct Config{
              using Name=fhicl::Name;
              using Comment=fhicl::Comment;
              fhicl::Atom<bool> mcdiag{Name("mcdiag"), Comment("set on for MC info"),true};
              fhicl::Atom<art::InputTag> chtag{Name("ComboHitCollection"),Comment("tag for combo hit collection")};
              fhicl::Atom<art::InputTag> tctag{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
              fhicl::Atom<art::InputTag> sttag{Name("CosmicTrackSeedCollection"),Comment("tag for cosmci track seed collection")};
              fhicl::Atom<bool> doDisplay{Name("doDisplay"),Comment("use display"), false};
              fhicl::Atom<bool> clickToAdvance{Name("clickToAdvance"),Comment("next event"), false};
         };
        typedef art::EDAnalyzer::Table<Config> Parameters;

              explicit CosmicFitDisplay(const Parameters& conf);

              virtual ~CosmicFitDisplay();
              virtual void beginJob() override;
              virtual void analyze(const art::Event& e) override;
    private:
        Config _conf;
              bool _mcdiag;
              Int_t _evt;

        // The module label of this instance of this module.
        std::string moduleLabel_;
        const ComboHitCollection* _chcol;
        const CosmicTrackSeedCollection* _coscol;
        //For Event Displays:
        TApplication* application_;
        TDirectory*   directory_ = nullptr;
        TTree* _cosmic_fit;
        TCanvas*      canvas_ = nullptr;
        TH2D* _display = nullptr;
        TNtuple* _ntTrack = nullptr;
        TNtuple* _ntHit = nullptr;

        art::InputTag   _chtag;//combo
        art::InputTag   _tctag;//timeclusters
        art::InputTag   _sttag;//Straight tracks
        bool doDisplay_;
        bool clickToAdvance_;

        Int_t _strawid;
        void plot2dDriftCompare(const art::Event& event);
        void plot2d(const art::Event& evt);
        void plotTrackerElements(const art::Event& event);
        void plot3dPrimes(const art::Event& evt);
        void plot3dXYZ(const art::Event& evt);
        void Geom3D(const art::Event& evt);
        std::vector<double> GetMaxAndMin(std::vector<double> myvector);
        bool findData(const art::Event& evt);
};

    CosmicFitDisplay::CosmicFitDisplay(const Parameters& conf) :
        art::EDAnalyzer(conf),
        _mcdiag (conf().mcdiag()),
        _chtag (conf().chtag()),
        _tctag (conf().tctag()),
        _sttag (conf().sttag()),
        doDisplay_ (conf().doDisplay()),
        clickToAdvance_ (conf().clickToAdvance())
        {}

    CosmicFitDisplay::~CosmicFitDisplay(){}

    void CosmicFitDisplay::beginJob() {
      // create diagnostics if requested...
      if ( !doDisplay_ ) return;
      art::ServiceHandle<art::TFileService> tfs;
      directory_ = gDirectory;
      // If needed, create the ROOT interactive environment. See note 1.
      if ( !gApplication ){
              int    tmp_argc(0);
              char** tmp_argv(0);
              application_ = new TApplication( "noapplication", &tmp_argc, tmp_argv );
      }
      // Create a canvas with a guaranteed unique name; the module label is unique within a job.
      TString name  = "canvas_"     + moduleLabel_;
      TString title = "Canvas for " + moduleLabel_;
      int window_size_x(1300);
      int window_size_y(600);
      canvas_ = tfs->make<TCanvas>(name,title,window_size_x,window_size_y);
      canvas_->Divide(2,2);

      }
      void CosmicFitDisplay::analyze(const art::Event& event) {
      //Call one or more of the macros from below here, for example:
        plot2dDriftCompare(event);
      }

/*
==============================NOTE:============================================//
Below here are a series of macros -  they are not glamorous but they produce useful debugging plots
================================================================================//
*/

      void CosmicFitDisplay::plot2dDriftCompare(const art::Event& event){
        _evt = event.id().event();
        findData(event);

        std::vector<double>  xseed, yseed, zseed, xdrift, ydrift, zdrift, a0seed, a1seed, b0seed, b1seed, a0drift, a1drift, b0drift, b1drift;
        std::vector<XYZVectorF> xprimes_seed, yprimes_seed, zprimes_seed, xprimes_drift,  yprimes_drift,  zprimes_drift;

        //find time clusters:
            unsigned  _ncosmics = _coscol->size();
        unsigned _nch = _chcol->size();
        //loop over tracks:

        for(size_t i =0; i < _ncosmics; i++){
                CosmicTrackSeed track =(*_coscol)[i];
                if(track._track.converged == false){continue;}
                if(track._track.minuit_converged == false){continue;}

                xprimes_seed.push_back(track._track.FitCoordSystem._XDoublePrime);
                yprimes_seed.push_back(track._track.FitCoordSystem._YDoublePrime);
                zprimes_seed.push_back(track._track.FitCoordSystem._ZPrime);

                xprimes_drift.push_back(track._track.MinuitCoordSystem._XDoublePrime);
                yprimes_drift.push_back(track._track.MinuitCoordSystem._YDoublePrime);
                zprimes_drift.push_back(track._track.MinuitCoordSystem._ZPrime);

                if(isnan(track._track.FitParams.A0) == true && isnan(track._track.FitParams.A1) == true && isnan(track._track.FitParams.B0) == true && isnan(track._track.FitParams.B1) == true) continue;

                        a0seed.push_back(track._track.FitParams.A0);
                        a1seed.push_back(track._track.FitParams.A1);
                        b0seed.push_back(track._track.FitParams.B0);
                        b1seed.push_back(track._track.FitParams.B1);
                        a0drift.push_back(track._track.MinuitParams.A0);
                        a1drift.push_back(track._track.MinuitParams.A1);
                        b0drift.push_back(track._track.MinuitParams.B0);
                        b1drift.push_back(track._track.MinuitParams.B1);

                        for(size_t i =0; i < _nch; i++){
                                ComboHit const& chit =(*_chcol)[i];
                                xseed.push_back(chit.pos().Dot(xprimes_seed[0]));
                                yseed.push_back(chit.pos().Dot(yprimes_seed[0]));
                                zseed.push_back(chit.pos().Dot(zprimes_seed[0]));
                                xdrift.push_back(chit.pos().Dot(xprimes_drift[0]));
                                ydrift.push_back(chit.pos().Dot(yprimes_drift[0]));
                                zdrift.push_back(chit.pos().Dot(zprimes_drift[0]));


                        }//end hit loop

                      double minz_seed = *std::min_element(zseed.begin(), zseed.end());
                      double maxz_seed = *std::max_element(zseed.begin(), zseed.end());
                      double minx_seed = *std::min_element(xseed.begin(), xseed.end());
                      double maxx_seed = *std::max_element(xseed.begin(), xseed.end());
                      double miny_seed = *std::min_element(yseed.begin(), yseed.end());
                      double maxy_seed = *std::max_element(yseed.begin(), yseed.end());

                      double minz_drift = *std::min_element(zdrift.begin(), zdrift.end());
                      double maxz_drift = *std::max_element(zdrift.begin(), zdrift.end());
                      double minx_drift = *std::min_element(xdrift.begin(), xdrift.end());
                      double maxx_drift = *std::max_element(xdrift.begin(), xdrift.end());
                      double miny_drift = *std::min_element(ydrift.begin(), ydrift.end());
                      double maxy_drift = *std::max_element(ydrift.begin(), ydrift.end());

                      GeomHandle<Tracker> th;
                      const Tracker* tracker = th.get();
                      TubsParams envelope(tracker->g4Tracker()->getInnerTrackerEnvelopeParams());
                      //double zlimit{envelope.zHalfLength()};
                      if (doDisplay_) {

                              std::cout << "Run: " << event.id().run()
                           << "  Subrun: " << event.id().subRun()
                           << "  Event: " << event.id().event()<<std::endl;
                              TLine  major_error_line, minor_error_line, out_line, fit_to_trackxprime, fit_to_trackyprime;
                              TPolyMarker poly;
                              TText  text;
                              TEllipse strawXsec;
                              TArc arc;
                              canvas_->SetTitle("foo title");
                              auto pad = canvas_->cd(1);
                              pad->Clear();
                              canvas_->SetTitle("bar title");

                              auto xzplot = pad->DrawFrame(minz_drift-100,minx_drift-100, maxz_drift+100, maxx_drift+150);
                              xzplot->GetYaxis()->SetTitleOffset(1.25);
                              xzplot->SetTitle( "Drift Fit X'' vs Z'; Z'(mm);X''(mm)");
                                    //arc.SetFillStyle(0);
                                    //arc.DrawArc(0.,0., envelope.outerRadius());
                                    //arc.DrawArc(0.,0., envelope.innerRadius());
                              fit_to_trackxprime.SetLineColor(kYellow);
                              fit_to_trackxprime.SetLineColor(kGreen);

                              poly.SetMarkerSize(1);
                              poly.SetMarkerStyle(4);

                              int ihit =0;
                              for(size_t i =0; i < _nch; i++){

                                        ComboHit const& chit =(*_chcol)[i];
                                              ihit +=1;
                                              if (ihit == 5) continue;
                                              if (ihit == 10) {
                                                      ihit = ihit+1;
                                              }
                                              if (ihit == 13) {
                                                      ihit = 29;
                                              }
                                        auto const& p = chit.pos();
                                        auto w = chit.uDir();
                                        auto const& s = chit.wireRes();
                                        auto const& t = chit.transRes();
                                        double x0prime{(p.Dot(xprimes_drift[0]))} ;
                                        poly.SetMarkerColor(ihit);
                                        major_error_line.SetLineColor(ihit);
                                        minor_error_line.SetLineColor(ihit);
                                        double z0prime{(p.Dot(zprimes_drift[0]))};
                                        poly.DrawPolyMarker( 1, &z0prime, &x0prime );
                                        XYZVectorF major = (s*w);
                                        XYZVectorF minor = GenVector::ZDir().Cross(w) * t;
                                        double major2 = (s*w).Mag2();
                                        double minor2 = (GenVector::ZDir().Cross(w) * t).Mag2();

                                        double major_x1 = p.Dot(xprimes_drift[0])+sqrt(major2)*(major.Unit()).Dot(xprimes_drift[0]);
                                        double major_x2 = p.Dot(xprimes_drift[0])-sqrt(major2)*(major.Unit()).Dot(xprimes_drift[0]);
                                        double major_z1 = p.Dot(zprimes_drift[0])+sqrt(major2)*(major.Unit()).Dot(zprimes_drift[0]);
                                        double major_z2 = p.Dot(zprimes_drift[0])-sqrt(major2)*(major.Unit()).Dot(zprimes_drift[0]);
                                        double minor_x1 = p.Dot(xprimes_drift[0])+sqrt(minor2)*(minor.Unit()).Dot(xprimes_drift[0]);
                                        double minor_x2 = p.Dot(xprimes_drift[0])-sqrt(minor2)*(minor.Unit()).Dot(xprimes_drift[0]);
                                        double minor_z1 = p.Dot(zprimes_drift[0])+sqrt(minor2)*(minor.Unit()).Dot(zprimes_drift[0]);
                                        double minor_z2 = p.Dot(zprimes_drift[0])-sqrt(minor2)*(minor.Unit()).Dot(zprimes_drift[0]);

                                        major_error_line.DrawLine( major_z1, major_x1, major_z2, major_x2);
                                        minor_error_line.DrawLine( minor_z1, minor_x1, minor_z2, minor_x2);

                                }//end hit loop

                               if(a1drift.size() > 0){

                                TF1 *trackline_xprime = new TF1("line", "[0]+[1]*x", minz_drift,maxz_drift);
                                trackline_xprime->SetParameter(0, a0drift[0]);
                                trackline_xprime->SetParameter(1, a1drift[0]);

                                trackline_xprime->SetLineColor(6);
                                trackline_xprime->Draw("same");

                                }

                              pad = canvas_->cd(2);
                              pad->Clear();
                              auto yzplot = pad->DrawFrame(minz_drift-100,miny_drift-100, maxz_drift+100, maxy_drift+150);
                              yzplot->GetYaxis()->SetTitleOffset(1.25);
                              yzplot->SetTitle( "Drift Fit Y'' vs Z'; Z'(mm);Y''(mm)");
                              ihit = 0;
                                    for(size_t i =0; i < _nch; i++){
                                ComboHit const& chit =(*_chcol)[i];
                                              ihit+=1;
                                              if (ihit == 5) continue;
                                              if (ihit == 10) {
                                                      ihit = ihit+1;
                                              }
                                              if (ihit == 13) {
                                                      ihit = 29;
                                              }
                                        auto const& p = chit.pos();
                                        auto w = chit.uDir();
                                        auto const& s = chit.wireRes();
                                        double y0prime{(p.Dot(yprimes_drift[0]))} ;
                                        double z0prime{(p.Dot(zprimes_drift[0]))};
                                        poly.SetMarkerColor(ihit);
                                        major_error_line.SetLineColor(ihit);
                                        poly.DrawPolyMarker( 1, &z0prime, &y0prime );
                                        double y1 = p.Dot(yprimes_drift[0])+s*w.Dot(yprimes_drift[0]);
                                        double y2 = p.Dot(yprimes_drift[0])-s*w.Dot(yprimes_drift[0]);
                                        double z1 = p.Dot(zprimes_drift[0])+s*w.Dot(zprimes_drift[0]);
                                        double z2 = p.Dot(zprimes_drift[0])-s*w.Dot(zprimes_drift[0]);
                                        major_error_line.DrawLine( z1, y1, z2, y2);

                                      }
                                      if(b1drift.size() > 0){

                                        TF1 *trackline_yprime = new TF1("line", "[0]+[1]*x",minz_drift,maxz_drift);
                                        trackline_yprime->SetParameter(0, b0drift[0]);
                                        trackline_yprime->SetParameter(1, b1drift[0]);
                                        trackline_yprime->SetLineColor(6);
                                        trackline_yprime->Draw("same");

                                          }
                                          pad = canvas_->cd(3);
                                          pad->Clear();
                                    xzplot = pad->DrawFrame(minz_seed-100,minx_seed, maxz_seed+100, maxx_seed+150);
                                          xzplot->GetYaxis()->SetTitleOffset(1.25);
                                    xzplot->SetTitle( "Seed Fit X'' vs Z'; Z'(mm);Y''(mm)");
                                          ihit =0;
                                    for(size_t i =0; i < _nch; i++){

                                        ComboHit const& chit =(*_chcol)[i];
                                              ihit +=1;
                                              if (ihit == 5) continue;
                                              if (ihit == 10) {
                                                      ihit = ihit+1;
                                              }
                                              if (ihit == 13) {
                                                      ihit = 29;
                                              }
                                        auto const& p = chit.pos();
                                        auto w = chit.uDir();
                                        auto const& s = chit.wireRes();
                                        auto const& t = chit.transRes();
                                        double x0prime{(p.Dot(xprimes_seed[0]))} ;
                                        poly.SetMarkerColor(ihit);
                                        major_error_line.SetLineColor(ihit);
                                        minor_error_line.SetLineColor(ihit);
                                        double z0prime{(p.Dot(zprimes_seed[0]))};
                                        poly.DrawPolyMarker( 1, &z0prime, &x0prime );
                                        XYZVectorF major = (s*w);
                                        XYZVectorF minor = GenVector::ZDir().Cross(w) * t;
                                        double major2 = (s*w).Mag2();
                                        double minor2 = (GenVector::ZDir().Cross(w) * t).Mag2();

                                        double major_x1 = p.Dot(xprimes_seed[0])+sqrt(major2)*(major.Unit()).Dot(xprimes_seed[0]);
                                        double major_x2 = p.Dot(xprimes_seed[0])-sqrt(major2)*(major.Unit()).Dot(xprimes_seed[0]);
                                        double major_z1 = p.Dot(zprimes_seed[0])+sqrt(major2)*(major.Unit()).Dot(zprimes_seed[0]);
                                        double major_z2 = p.Dot(zprimes_seed[0])-sqrt(major2)*(major.Unit()).Dot(zprimes_seed[0]);
                                        double minor_x1 = p.Dot(xprimes_seed[0])+sqrt(minor2)*(minor.Unit()).Dot(xprimes_seed[0]);
                                        double minor_x2 = p.Dot(xprimes_seed[0])-sqrt(minor2)*(minor.Unit()).Dot(xprimes_seed[0]);
                                        double minor_z1 = p.Dot(zprimes_seed[0])+sqrt(minor2)*(minor.Unit()).Dot(zprimes_seed[0]);
                                        double minor_z2 = p.Dot(zprimes_seed[0])-sqrt(minor2)*(minor.Unit()).Dot(zprimes_seed[0]);

                                        major_error_line.DrawLine( major_z1, major_x1, major_z2, major_x2);
                                        minor_error_line.DrawLine( minor_z1, minor_x1, minor_z2, minor_x2);

                                }//end hit loop

                               if(a1seed.size() > 0){

                                TF1 *trackline_xprime = new TF1("line", "[0]+[1]*x", minz_seed,maxz_seed);
                                trackline_xprime->SetParameter(0, a0seed[0]);
                                trackline_xprime->SetParameter(1, a1seed[0]);

                                trackline_xprime->SetLineColor(6);
                                trackline_xprime->Draw("same");

                                }
                              pad = canvas_->cd(4);
                              pad->Clear();
                              yzplot = pad->DrawFrame(minz_seed-100,miny_seed-100, maxz_seed+100, maxy_seed+150);
                              yzplot->GetYaxis()->SetTitleOffset(1.25);
                              yzplot->SetTitle( "Seed Fit Y'' vs Z'; Z'(mm);Y''(mm)");
                              ihit = 0;
                                    for(size_t i =0; i < _nch; i++){
                                ComboHit const& chit =(*_chcol)[i];
                                              ihit+=1;
                                              if (ihit == 5) continue;
                                              if (ihit == 10) {
                                                      ihit = ihit+1;
                                              }
                                              if (ihit == 13) {
                                                      ihit = 29;
                                              }
                                        auto const& p = chit.pos();
                                        auto w = chit.uDir();
                                        auto const& s = chit.wireRes();
                                        double y0prime{(p.Dot(yprimes_seed[0]))} ;
                                        double z0prime{(p.Dot(zprimes_seed[0]))};
                                        poly.SetMarkerColor(ihit);
                                        major_error_line.SetLineColor(ihit);
                                        poly.DrawPolyMarker( 1, &z0prime, &y0prime );
                                        double y1 = p.Dot(yprimes_seed[0])+s*w.Dot(yprimes_seed[0]);
                                        double y2 = p.Dot(yprimes_seed[0])-s*w.Dot(yprimes_seed[0]);
                                        double z1 = p.Dot(zprimes_seed[0])+s*w.Dot(zprimes_seed[0]);
                                        double z2 = p.Dot(zprimes_seed[0])-s*w.Dot(zprimes_seed[0]);
                                        major_error_line.DrawLine( z1, y1, z2, y2);
                                      }
                                      if(b1seed.size() > 0){
                                        TF1 *trackline_yprime = new TF1("line", "[0]+[1]*x",minz_seed, maxz_seed);
                                        trackline_yprime->SetParameter(0, b0seed[0]);
                                        trackline_yprime->SetParameter(1, b1seed[0]);
                                        trackline_yprime->SetLineColor(6);
                                        trackline_yprime->Draw("same");
                                          }
                                          ostringstream title;
                                      title << "Run: " << event.id().run()
                                      << "  Subrun: " << event.id().subRun()
                                      << "  Event: " << event.id().event()<<".root";

                                      text.SetTextAlign(11);
                                      text.DrawTextNDC( 0., 0.01, title.str().c_str() );

                                      canvas_->Modified();
                                      canvas_->Update();
                                      canvas_->SaveAs(title.str().c_str());
                                      if ( clickToAdvance_ ){
                                        cerr << "Double click in the Canvas " << moduleLabel_ << " to continue:" ;
                                        gPad->WaitPrimitive();
                                            } else{
                                        char junk;
                                        cerr << "Enter any character to continue: ";
                                        cin >> junk;
                                      }
                                              cerr << endl;

                }//end cosmic loop
                }//display
                }//end drift plotting



      void CosmicFitDisplay::plot2d(const art::Event& event){
        _evt = event.id().event();
        findData(event);

        std::vector<double> initpullsx, initpullsy, pullsx, pullsy, x, y, z, rawx, rawy, rawz, xinit, yinit, zinit, out_x, out_y, out_z, a0, a1, b0, b1,a0init, a1init, b0init, b1init, chi_dof_XDoublePrimeZPrime, chi_dof_YDoublePrimeZPrime, initchi_dof_YDoublePrimeZPrime, initchi_dof_XDoublePrimeZPrime;

        std::vector<XYZVectorF> xprimes, yprimes, zprimes, xprimesinit,  yprimesinit,  zprimesinit;

         //find time clusters:
            unsigned  _ncosmics = _coscol->size();
        unsigned _nch = _chcol->size();
        //loop over tracks:
        for(size_t i =0; i < _ncosmics; i++){

                CosmicTrackSeed track =(*_coscol)[i];
                if(track._track.converged == false){continue;}
                xprimes.push_back(track._track.FitCoordSystem._XDoublePrime);
                yprimes.push_back(track._track.FitCoordSystem._YDoublePrime);
                zprimes.push_back(track._track.FitCoordSystem._ZPrime);

                xprimesinit.push_back(track._track.InitCoordSystem._XDoublePrime);
                yprimesinit.push_back(track._track.InitCoordSystem._YDoublePrime);
                zprimesinit.push_back(track._track.InitCoordSystem._ZPrime);

                if(isnan(track._track.FitParams.A0) == true && isnan(track._track.FitParams.A1) == true && isnan(track._track.FitParams.B0) == true && isnan(track._track.FitParams.B1) == true) continue;

                        a0.push_back(track._track.FitParams.A0);
                        a1.push_back(track._track.FitParams.A1);
                        b0.push_back(track._track.FitParams.B0);
                        b1.push_back(track._track.FitParams.B1);
                        a0init.push_back(track._track.InitParams.A0);
                        a1init.push_back(track._track.InitParams.A1);
                        b0init.push_back(track._track.InitParams.B0);
                        b1init.push_back(track._track.InitParams.B1);

                        chi_dof_XDoublePrimeZPrime.push_back(track._track.Diag.FinalChiX);
                        chi_dof_YDoublePrimeZPrime.push_back(track._track.Diag.FinalChiY);

                        initchi_dof_XDoublePrimeZPrime.push_back(track._track.Diag.InitialChiX);
                        initchi_dof_YDoublePrimeZPrime.push_back(track._track.Diag.InitialChiY);



        }
      if(xprimes.size() >0){
        // loop over combo hits

        for(size_t i =0; i < _nch; i++){
                ComboHit const& chit =(*_chcol)[i];
                x.push_back(chit.pos().Dot(xprimes[0]));
                y.push_back(chit.pos().Dot(yprimes[0]));
                z.push_back(chit.pos().Dot(zprimes[0]));
                xinit.push_back(chit.pos().Dot(xprimesinit[0]));
                yinit.push_back(chit.pos().Dot(yprimesinit[0]));
                zinit.push_back(chit.pos().Dot(zprimesinit[0]));
                rawx.push_back(chit.pos().x());
                rawy.push_back(chit.pos().y());
                rawz.push_back(chit.pos().z());

        }


        double minz = *std::min_element(z.begin(), z.end());
        double maxz = *std::max_element(z.begin(), z.end());

        double minx = *std::min_element(x.begin(), x.end());
        double maxx = *std::max_element(x.begin(), x.end());
        double miny = *std::min_element(y.begin(), y.end());
        double maxy = *std::max_element(y.begin(), y.end());

        double minzinit = *std::min_element(zinit.begin(), zinit.end());
        double maxzinit = *std::max_element(zinit.begin(), zinit.end());

        double minxinit = *std::min_element(xinit.begin(), xinit.end());
        double maxxinit = *std::max_element(xinit.begin(), xinit.end());
        double minyinit = *std::min_element(yinit.begin(), yinit.end());
        double maxyinit = *std::max_element(yinit.begin(), yinit.end());

        double minrawz = *std::min_element(rawz.begin(), rawz.end());
        double maxrawz = *std::max_element(rawz.begin(), rawz.end());

        double minrawx = *std::min_element(rawx.begin(), rawx.end());
        double maxrawx = *std::max_element(rawx.begin(), rawx.end());
        double minrawy = *std::min_element(rawy.begin(), rawy.end());
        double maxrawy = *std::max_element(rawy.begin(), rawy.end());

        GeomHandle<Tracker> tracker;
        GeomHandle<Straw> straw;
        GeomHandle<Plane> plane;
        GeomHandle<Panel> panel;
        // Annulus of a cylinder that bounds the tracker/straw info:
        TubsParams envelope(tracker->g4Tracker()->getInnerTrackerEnvelopeParams());

        if (doDisplay_) {

              std::cout << "Run: " << event.id().run()
           << "  Subrun: " << event.id().subRun()
           << "  Event: " << event.id().event()<<std::endl;
              TLine  major_error_line, minor_error_line, out_line, fit_to_trackxprime, fit_to_trackyprime;
              TArc   arc;
              TPolyMarker poly, out;
              TBox   box;
              TText  text;

              canvas_->SetTitle("foo title");
              auto pad = canvas_->cd(1);
              pad->Clear();
              canvas_->SetTitle("bar title");

              auto xzplot = pad->DrawFrame(minz-100,minx-100, maxz+100, maxx+150);
              xzplot->GetYaxis()->SetTitleOffset(1.25);
              xzplot->SetTitle( "Final Iteration Fit X'' vs Z'; Z'(mm);X''(mm)");
              arc.SetFillStyle(0);
              //arc.DrawArc(0.,0., envelope.outerRadius());
              //arc.DrawArc(0.,0., envelope.innerRadius());
              fit_to_trackxprime.SetLineColor(kYellow);
              fit_to_trackxprime.SetLineColor(kGreen);

              poly.SetMarkerSize(1);
              poly.SetMarkerStyle(4);

              int ihit =0;
              for(size_t i =0; i < _nch; i++){

                ComboHit const& chit =(*_chcol)[i];
                              ihit +=1;
                              if (ihit == 5) continue;
                              if (ihit == 10) {
                                      ihit = ihit+1;
                              }
                              if (ihit == 13) {
                                      ihit = 29;
                              }
                        auto const& p = chit.pos();
                        auto w = chit.uDir();
                        auto const& s = chit.wireRes();
                        auto const& t = chit.transRes();
                        double x0prime{(p.Dot(xprimes[0]))} ;
                        poly.SetMarkerColor(ihit);
                        major_error_line.SetLineColor(ihit);
                        minor_error_line.SetLineColor(ihit);
                        double z0prime{(p.Dot(zprimes[0]))};
                        poly.DrawPolyMarker( 1, &z0prime, &x0prime );
                        XYZVectorF major = (s*w);
                        XYZVectorF minor = GenVector::ZDir().Cross(w) * t;
                        double major2 = (s*w).Mag2();
                        double minor2 = (GenVector::ZDir().Cross(w) * t).Mag2();

                        double major_x1 = p.Dot(xprimes[0])+sqrt(major2)*(major.Unit()).Dot(xprimes[0]);
                        double major_x2 = p.Dot(xprimes[0])-sqrt(major2)*(major.Unit()).Dot(xprimes[0]);
                        double major_z1 = p.Dot(zprimes[0])+sqrt(major2)*(major.Unit()).Dot(zprimes[0]);
                        double major_z2 = p.Dot(zprimes[0])-sqrt(major2)*(major.Unit()).Dot(zprimes[0]);
                        double minor_x1 = p.Dot(xprimes[0])+sqrt(minor2)*(minor.Unit()).Dot(xprimes[0]);
                        double minor_x2 = p.Dot(xprimes[0])-sqrt(minor2)*(minor.Unit()).Dot(xprimes[0]);
                        double minor_z1 = p.Dot(zprimes[0])+sqrt(minor2)*(minor.Unit()).Dot(zprimes[0]);
                        double minor_z2 = p.Dot(zprimes[0])-sqrt(minor2)*(minor.Unit()).Dot(zprimes[0]);

                        major_error_line.DrawLine( major_z1, major_x1, major_z2, major_x2);
                        minor_error_line.DrawLine( minor_z1, minor_x1, minor_z2, minor_x2);
                        /*
                        TLatex latex;
                        stringstream pulls;
                        pulls<<pullsx[ihit-1];
                        pulls.precision(2);
                        const char* str_pulls = pulls.str().c_str();
                           latex.SetTextSize(0.05);
                           latex.SetTextColor(ihit);
                           latex.SetTextAlign(13);  //align at top

                           if(i%2 == 0){
                           latex.DrawLatex(z0prime-10, x0prime+75,str_pulls);
                           }
                           if(i%2 != 0){
                           latex.DrawLatex(z0prime-10, x0prime-75,str_pulls);
                           }
                        */
              }

              if(a1.size() > 0){

                TF1 *trackline_xprime = new TF1("line", "[0]+[1]*x", minz,maxz);
                trackline_xprime->SetParameter(0, a0[0]);
                trackline_xprime->SetParameter(1, a1[0]);

                trackline_xprime->SetLineColor(6);
                trackline_xprime->Draw("same");
                TLegend *leg = new TLegend(0.1,0.8,0.3,0.9);
                leg->AddEntry("#Chi^{2}/N = ", "#Chi^{2}/N =",  "");
                stringstream chi_info;
                chi_info<< chi_dof_XDoublePrimeZPrime[0];
                const char* str_chi_info = chi_info.str().c_str();
                leg->AddEntry("" ,str_chi_info,  "");
                leg->Draw("same");
              }
              pad = canvas_->cd(2);
              pad->Clear();
              auto yzplot = pad->DrawFrame(minz-100,miny-100, maxz+100, maxy+150);
              yzplot->GetYaxis()->SetTitleOffset(1.25);
              yzplot->SetTitle( "Final Iteration Fit Y'' vs Z'; Z'(mm);Y''(mm)");
              ihit = 0;
              for(size_t i =0; i < _nch; i++){
                ComboHit const& chit =(*_chcol)[i];
                              ihit+=1;
                              if (ihit == 5) continue;
                              if (ihit == 10) {
                                      ihit = ihit+1;
                              }
                              if (ihit == 13) {
                                      ihit = 29;
                              }
                        auto const& p = chit.pos();
                        auto w = chit.uDir();
                        auto const& s = chit.wireRes();
                        double y0prime{(p.Dot(yprimes[0]))} ;
                        double z0prime{(p.Dot(zprimes[0]))};
                        poly.SetMarkerColor(ihit);
                        major_error_line.SetLineColor(ihit);
                        poly.DrawPolyMarker( 1, &z0prime, &y0prime );
                        double y1 = p.Dot(yprimes[0])+s*w.Dot(yprimes[0]);
                        double y2 = p.Dot(yprimes[0])-s*w.Dot(yprimes[0]);
                        double z1 = p.Dot(zprimes[0])+s*w.Dot(zprimes[0]);
                        double z2 = p.Dot(zprimes[0])-s*w.Dot(zprimes[0]);
                        major_error_line.DrawLine( z1, y1, z2, y2);
                        /*
                        TLatex latex;
                        stringstream pulls;
                        pulls<<pullsy[ihit-1];
                        pulls.precision(2);
                        const char* str_pulls = pulls.str().c_str();
                           latex.SetTextSize(0.05);
                           latex.SetTextColor(ihit);
                           latex.SetTextAlign(13);  //align at top

                           if(i%2 == 0){
                           latex.DrawLatex(z0prime-10, y0prime+75,str_pulls);
                           }
                           if(i%2 != 0){
                           latex.DrawLatex(z0prime-10, y0prime-75,str_pulls);
                           }
                           */
              }
              if(b1.size() > 0){

                TF1 *trackline_yprime = new TF1("line", "[0]+[1]*x",minz,maxz);
                trackline_yprime->SetParameter(0, b0[0]);
                trackline_yprime->SetParameter(1, b1[0]);
                trackline_yprime->SetLineColor(6);
                trackline_yprime->Draw("same");
                TLegend *leg = new TLegend(0.1,0.8,0.2,0.9);
                leg->AddEntry("#Chi^{2}/N ", "#Chi^{2}/N",  "");
                stringstream chi_info;
                chi_info<< chi_dof_YDoublePrimeZPrime[0];
                const char* str_chi_info = chi_info.str().c_str();
                leg->AddEntry("" ,str_chi_info,  "");
                leg->Draw("same");

              }

              pad = canvas_->cd(3);
              pad->Clear();
              auto xzplot_init = pad->DrawFrame(minzinit-100,minxinit-100, maxzinit+100, maxxinit+150);
              xzplot_init->GetYaxis()->SetTitleOffset(1.25);
              xzplot_init->SetTitle( "Initial Fit (chi 2 min) X'' vs Z'; Z'(mm);X''(mm)");
              ihit = 0;
              for(size_t i =0; i < _nch; i++){
                ComboHit const& chit =(*_chcol)[i];
                              ihit+=1;
                              if (ihit == 5) continue;
                              if (ihit == 10) {
                                      ihit = ihit+1;
                              }
                              if (ihit == 13) {
                                      ihit = 29;
                              }
                        auto const& p = chit.pos();
                        auto w = chit.uDir();
                        auto const& s = chit.wireRes();
                        double x0primeinit{(p.Dot(xprimesinit[0]))} ;
                        double z0primeinit{(p.Dot(zprimesinit[0]))};
                        poly.SetMarkerColor(ihit);
                        major_error_line.SetLineColor(ihit);
                        poly.DrawPolyMarker( 1, &z0primeinit, &x0primeinit );
                        double x1 = p.Dot(xprimesinit[0])+s*w.Dot(xprimesinit[0]);
                        double x2 = p.Dot(xprimesinit[0])-s*w.Dot(xprimesinit[0]);
                        double z1 = p.Dot(zprimesinit[0])+s*w.Dot(zprimesinit[0]);
                        double z2 = p.Dot(zprimesinit[0])-s*w.Dot(zprimesinit[0]);
                        major_error_line.DrawLine( z1, x1, z2, x2);
                        /*
                        TLatex latex;
                        stringstream pulls;
                        pulls<<initpullsx[ihit-1];
                        pulls.precision(2);
                        const char* str_pulls = pulls.str().c_str();
                           latex.SetTextSize(0.05);
                           latex.SetTextColor(ihit);
                           latex.SetTextAlign(13);  //align at top

                           if(i%2 == 0){
                           latex.DrawLatex(z0primeinit-10, x0primeinit+75,str_pulls);
                           }
                           if(i%2 != 0){
                           latex.DrawLatex(z0primeinit-10, x0primeinit-75,str_pulls);
                           }
                           */
              }
              if(a1.size() > 0){

                TF1 *trackline_yprime = new TF1("line", "[0]+[1]*x",minzinit,maxzinit);
                trackline_yprime->SetParameter(0, a0init[0]);
                trackline_yprime->SetParameter(1, a1init[0]);
                trackline_yprime->SetLineColor(6);
                //trackline_yprime->Draw("same");
                TLegend *leg = new TLegend(0.1,0.8,0.2,0.9);
                leg->AddEntry("#Chi^{2}/N = ", "#Chi^{2}/N =",  "");
                stringstream chi_info;
                chi_info<< initchi_dof_XDoublePrimeZPrime[0];
                const char* str_chi_info = chi_info.str().c_str();
                leg->AddEntry("" ,str_chi_info,  "");
                leg->Draw("same");

              }

              pad = canvas_->cd(4);
              pad->Clear();
              auto yzplot_init = pad->DrawFrame(minzinit-100,minyinit-100, maxzinit+100, maxyinit+150);
              yzplot_init->GetYaxis()->SetTitleOffset(1.25);
              yzplot_init->SetTitle( "Initial (chi 2 min) Fit Y'' vs Z'; Z'(mm);Y''(mm)");
              ihit = 0;
              for(size_t i =0; i < _nch; i++){
                ComboHit const& chit =(*_chcol)[i];
                              ihit+=1;
                              if (ihit == 5) continue;
                              if (ihit == 10) {
                                      ihit = ihit+1;
                              }
                              if (ihit == 13) {
                                      ihit = 29;
                              }
                        auto const& p = chit.pos();
                        auto w = chit.uDir();
                        auto const& s = chit.wireRes();
                        double y0prime{(p.Dot(yprimesinit[0]))} ;
                        double z0prime{(p.Dot(zprimesinit[0]))};
                        poly.SetMarkerColor(ihit);
                        major_error_line.SetLineColor(ihit);
                        poly.DrawPolyMarker( 1, &z0prime, &y0prime );
                        double y1 = p.Dot(yprimesinit[0])+s*w.Dot(yprimesinit[0]);
                        double y2 = p.Dot(yprimesinit[0])-s*w.Dot(yprimesinit[0]);
                        double z1 = p.Dot(zprimesinit[0])+s*w.Dot(zprimesinit[0]);
                        double z2 = p.Dot(zprimesinit[0])-s*w.Dot(zprimesinit[0]);
                        major_error_line.DrawLine( z1, y1, z2, y2);
                        /*
                        TLatex latex;
                        stringstream pulls;
                        pulls<<initpullsy[ihit-1];
                        pulls.precision(2);
                        const char* str_pulls = pulls.str().c_str();
                           latex.SetTextSize(0.05);
                           latex.SetTextColor(ihit);
                           latex.SetTextAlign(13);  //align at top
                           if(i%2 == 0){
                           latex.DrawLatex(z0prime-10, y0prime+75,str_pulls);
                           }
                           if(i%2 != 0){
                           latex.DrawLatex(z0prime-10, y0prime-75,str_pulls);
                           }
                           */
              }
              if(b1.size() > 0){

                TF1 *trackline_yprime = new TF1("line", "[0]+[1]*x",minzinit,maxzinit);
                trackline_yprime->SetParameter(0, b0init[0]);
                trackline_yprime->SetParameter(1, b1init[0]);
                trackline_yprime->SetLineColor(6);
                //trackline_yprime->Draw("same");
                TLegend *leg = new TLegend(0.1,0.8,0.2,0.9);
                leg->AddEntry("#Chi^{2}/N = ", "#Chi^{2}/N =",  "");
                stringstream chi_info;
                chi_info<< initchi_dof_YDoublePrimeZPrime[0];
                const char* str_chi_info = chi_info.str().c_str();
                leg->AddEntry("" ,str_chi_info,  "");
                leg->Draw("same");
              }

              pad = canvas_->cd(5);
              pad->Clear();
              auto rawxzplot = pad->DrawFrame(minrawz-100,minrawx-100, maxrawz+100, maxrawx+100);
              rawxzplot->GetYaxis()->SetTitleOffset(1.25);
              rawxzplot->SetTitle( "Raw Hits X vs Z; Z(mm);X(mm)");
              ihit = 0;
              for(size_t i =0; i < _nch; i++){
                ComboHit const& chit =(*_chcol)[i];
                              ihit +=1;
                              if (ihit == 5) continue;
                              if (ihit == 10) {
                                      ihit = ihit+1;
                              }
                              if (ihit == 13) {
                                      ihit = 29;
                              }
                        auto const& p = chit.pos();
                        auto w = chit.uDir();
                        auto const& s = chit.wireRes();
                        double y0prime{p.x()} ;
                        double z0prime{p.z()};
                        poly.SetMarkerColor(ihit);
                        poly.DrawPolyMarker( 1, &z0prime, &y0prime );
                        double y1 = p.x()+s*w.x();
                        double y2 = p.x()-s*w.x();
                        double z1 = p.z()+s*w.z();
                        double z2 = p.z()-s*w.z();
                        major_error_line.SetLineColor(ihit);
                        major_error_line.DrawLine( z1, y1, z2, y2);
              }

              TF1 *trackline_x = new TF1("line", "[0]+[1]*x",  minrawz, maxrawz);
              trackline_x->SetParameter(0, a0[0]*xprimes[0].x());
              trackline_x->SetParameter(1, a1[0]);
              trackline_x->SetLineColor(6);
              //trackline_x->Draw("same");

              pad = canvas_->cd(6);
              pad->Clear();
              auto rawyzplot = pad->DrawFrame(minrawz-100,minrawy-100, maxrawz+100, maxrawy+100);//-plotLimits,-plotLimits,plotLimits,plotLimits);//pad->DrawFrame(minrawz-100, minrawy-100, maxrawz+100, maxrawy+100);
              rawyzplot->GetYaxis()->SetTitleOffset(1.25);
              rawyzplot->SetTitle( "Raw Hits Y vs Z; Z(mm);Y(mm)");
              ihit=0;
              for(size_t i =0; i < _nch; i++){
                ComboHit const& chit =(*_chcol)[i];
                              ihit +=1;
                              if (ihit == 5) continue;
                              if (ihit == 10) {
                                      ihit = ihit+1;
                              }
                        auto const& p = chit.pos();
                        auto w = chit.uDir();
                        auto const& s = chit.wireRes();
                        double y0prime{p.y()} ;
                        double z0prime{p.z()};
                        poly.SetMarkerColor(ihit);
                        poly.DrawPolyMarker( 1, &z0prime, &y0prime );
                        double y1 = p.y()+s*w.y();
                        double y2 = p.y()-s*w.y();
                        double z1 = p.z()+s*w.z();
                        double z2 = p.z()-s*w.z();
                        major_error_line.SetLineColor(ihit);
                        major_error_line.DrawLine( z1, y1, z2, y2);
              }

              TF1 *trackline_y = new TF1("line", "[0]+[1]*x", minrawz, maxrawz);
              trackline_y->SetParameter(0,b0[0]);
              trackline_y->SetParameter(1, b1[0]);
              trackline_y->SetLineColor(6);
              //trackline_y->Draw("same");


              ostringstream title;
              title << "Run: " << event.id().run()
              << "  Subrun: " << event.id().subRun()
              << "  Event: " << event.id().event()<<".root";

              text.SetTextAlign(11);
              text.DrawTextNDC( 0., 0.01, title.str().c_str() );

              canvas_->Modified();
              canvas_->Update();
              canvas_->SaveAs(title.str().c_str());
              if ( clickToAdvance_ ){
                cerr << "Double click in the Canvas " << moduleLabel_ << " to continue:" ;
                gPad->WaitPrimitive();
                    } else{
                char junk;
                cerr << "Enter any character to continue: ";
                cin >> junk;
              }
                      cerr << endl;
        }//End Diag
      }
      }//End 2d

      //Plots XZ or YZ of hits with straws:
      void CosmicFitDisplay::plotTrackerElements(const art::Event& event){
        _evt = event.id().event();
        findData(event);

         //find time clusters:
            unsigned  _ncosmics = _coscol->size();
        unsigned _nch = _chcol->size();
        //loop over tracks:
        std::vector<double> x,y,z;
        std::vector<XYZVectorF> xprimes, yprimes, zprimes;
        for(size_t i =0; i < _ncosmics; i++){

          CosmicTrackSeed track =(*_coscol)[i];
          if(track._track.converged == false){continue;}

          for(size_t i =0; i < _nch; i++){
                ComboHit const& chit =(*_chcol)[i];
                x.push_back(chit.pos().x());
                y.push_back(chit.pos().y());
                z.push_back(chit.pos().z());

            }
            double minz = *std::min_element(z.begin(), z.end());
                  double maxz = *std::max_element(z.begin(), z.end());
                  //double minx = *std::min_element(x.begin(), x.end());
            //double maxx = *std::max_element(x.begin(), x.end());
            double miny = *std::min_element(y.begin(), y.end());
            double maxy = *std::max_element(y.begin(), y.end());
                   GeomHandle<Tracker> th;
                   const Tracker* tracker = th.get();
            // Annulus of a cylinder that bounds the tracker/straw info:
            TubsParams envelope(tracker->g4Tracker()->getInnerTrackerEnvelopeParams());
            if (doDisplay_) {
                TPolyMarker poly, wirecentre;
                      TBox   box;
                TText  text;
                TLine error_line;
                TEllipse strawXsec;
                poly.SetMarkerStyle(2);
                //poly.SetMarkerSize(straw_radius);

                auto pad = canvas_;
                      pad->Clear();
                      canvas_->SetTitle("bar title");
                auto yzplot = pad->DrawFrame(minz-20, miny-20, maxz+20, maxy+20);
                      yzplot->GetYaxis()->SetTitleOffset(1.25);
                int ihit = 0;
                for(size_t i =0; i < _nch; i++){
                        ComboHit const& chit =(*_chcol)[i];
                              ihit +=1;
                              if (ihit == 5) continue;
                              if (ihit == 10) {
                                      ihit = ihit+1;
                              }

                        auto const& p = chit.pos();
                        auto w = chit.uDir();
                        auto const& s = chit.wireRes();
                        double y0prime{p.y()} ;
                        double z0prime{p.z()};
                        poly.SetMarkerColor(ihit);
                        poly.DrawPolyMarker( 1, &z0prime, &y0prime );

                        double y1 = p.y()+s*w.y();
                        double y2 = p.y()-s*w.y();
                        double z1 = p.z()+s*w.z();
                        double z2 = p.z()-s*w.z();
                        double const straw_radius = tracker->strawOuterRadius();
                        cout<<straw_radius<<endl;
                        Straw const&  straw = tracker->getStraw(chit.strawId());
                        auto spos = straw.getMidPoint();
                              auto sdir = straw.getDirection();
                              auto wpos = spos + chit.wireDist()*sdir;

                         wirecentre.SetMarkerColor(ihit);
                         double y0wire{wpos.y()} ;
                        double z0wire{wpos.z()};
                        wirecentre.DrawPolyMarker( 1,&z0wire , &y0wire );
                        error_line.SetLineColor(ihit);
                        error_line.DrawLine( z1, y1, z2, y2);
                        strawXsec.SetFillStyle(ihit);
                        strawXsec.DrawEllipse(wpos.z(), wpos.y(), straw_radius, straw_radius, 0 ,360, 0);
                   }
                   ostringstream title;
                      title << "Run: " << event.id().run()
                      << "  Subrun: " << event.id().subRun()
                      << "  Event: " << event.id().event()<<".root";

                      text.SetTextAlign(11);
                      text.DrawTextNDC( 0., 0.01, title.str().c_str() );
                      canvas_->SetTitle("foo title");

                      canvas_->Modified();
                      canvas_->Update();
                      canvas_->SaveAs(title.str().c_str());
                      if ( clickToAdvance_ ){
                        cerr << "Double click in the Canvas " << moduleLabel_ << " to continue:" ;
                        gPad->WaitPrimitive();
                            } else{
                        char junk;
                        cerr << "Enter any character to continue: ";
                        cin >> junk;
                      }
                              cerr << endl;
                        }//End Diag
                 }
      }//End 2d


      void CosmicFitDisplay::plot3dXYZ(const art::Event& event){

        _evt = event.id().event();  // add event id

        auto comboHits  = event.getValidHandle<ComboHitCollection>( _chtag );
        auto Tracks  = event.getValidHandle<CosmicTrackSeedCollection>( _sttag );
        GeomHandle<Tracker> tracker;

        // Annulus of a cylinder that bounds the tracker
        TubsParams envelope(tracker->g4Tracker()->getInnerTrackerEnvelopeParams());

        std::vector<double> x, y, z, a0, a1, b0, b1;
        std::vector<XYZVectorF> xprimes, yprimes, zprimes, initial_track_direction;
        x.reserve(comboHits->size());
        y.reserve(comboHits->size());
        z.reserve(comboHits->size());
        a1.reserve(Tracks->size());
        b1.reserve(Tracks->size());
        a0.reserve(Tracks->size());
        b0.reserve(Tracks->size());

        //Get track coordinate system and fit parameters
        for(auto const& track: *Tracks){
        xprimes.push_back(track._track.FitCoordSystem._XDoublePrime);
        yprimes.push_back(track._track.FitCoordSystem._YDoublePrime);
        zprimes.push_back(track._track.FitCoordSystem._ZPrime);
        initial_track_direction.push_back(track._track.InitParams.Direction() );
        a1.push_back(track._track.FitParams.A0);
        a0.push_back(track._track.FitParams.A1);
        b1.push_back(track._track.FitParams.B0);
        b0.push_back(track._track.FitParams.B1);
        }

        if(xprimes.size() >0){
        for(auto const& chit : *comboHits){

        x.push_back(chit.pos().x());
        y.push_back(chit.pos().y());
        z.push_back(chit.pos().z());
        }

        if (doDisplay_) {
              std::cout << "Run: " << event.id().run()
           << "  Subrun: " << event.id().subRun()
           << "  Event: " << event.id().event()<<std::endl;
              TTUBE tube;
              TText  text;
              TPolyMarker3D poly3D;

              // Draw the frame for the cylinders in plot:
              double plotLimits(1000.);
              double zlimit{envelope.zHalfLength()};
              //Create Pad:
              auto pad = canvas_;
              pad->Clear();
              canvas_->Draw();

              TH3D *xyzplot =new TH3D("XYZ","XYZ",1,-plotLimits,plotLimits,1,-plotLimits,plotLimits,1,-zlimit,zlimit);

              xyzplot->GetYaxis()->SetTitleOffset(1.5);
              xyzplot->GetXaxis()->SetTitleOffset(1.5);
                 xyzplot->GetZaxis()->SetTitleOffset(1.5);
              xyzplot->SetTitle( "Visualization in XYZ;X (mm);Y (mm); Z(mm)");
              xyzplot->SetMarkerStyle(2);
              xyzplot->SetMarkerSize(0.65);
              xyzplot->SetMarkerColor(kRed);
              xyzplot->SetStats(0);
              xyzplot->Draw();
              TTUBE *Tube = new TTUBE("","","", envelope.innerRadius(), envelope.outerRadius(), zlimit);
              Tube->SetFillColor(18);
              Tube->SetLineColor(18);
              Tube->Draw();
              poly3D.SetMarkerStyle(2);
              poly3D.SetMarkerSize(0.65);
              poly3D.SetMarkerColor(kRed);
              int index = 1;
              for ( auto const& chit : *comboHits ){
                        TPolyLine3D *errors = new TPolyLine3D;
                        errors->SetLineColor(kRed);
                        auto const& p = chit.pos();
                        auto w = chit.uDir();
                        auto const& s = chit.wireRes();
                        double x1 = p.x()+s*w.x();
                        double x2 = p.x()-s*w.x();
                        double z1 = p.z()+s*w.z();
                        double z2 = p.z()-s*w.z();
                        double y1 = p.y()+s*w.y();
                        double y2 = p.y()-s*w.y();
                        poly3D.SetPoint(index, p.x(), p.y(),p.z() );
                        errors->SetPoint(0, x1,y1,z1);
                        errors->SetNextPoint(x2,y2,z2);
                        index+=1;
                        errors->Draw("same");
                    }
              poly3D.Draw("same");
              if(a1.size() > 0 && b1.size() > 0){
                      TPolyLine3D *Track = new TPolyLine3D;
                      TPolyLine3D *InitTrack = new TPolyLine3D;
                Track->SetLineColor(kBlue);
                InitTrack->SetLineColor(kGreen);
                double tz1 = -1500;//z[0];
                double tz2 = 1500;//z[z.size()-1];

                double tx1 = a0[0] + a1[0]*tz1;
                double tx2 = a0[0] + a1[0]*tz2;
                double ty1 = b0[0] + b1[0]*tz1;//b1[0]*tz1+b0[0]*tx1+a0[0];
                double ty2 = b0[0] + b1[0]*tz2; //a1[0]*tz2+a0[0]*tx2+b0[0];
                Track->SetPoint(0, tx1, ty1, tz1);
                Track->SetNextPoint(tx2, ty2, tz2);
                Track->Draw("same");

              }

              canvas_->Update();

              //END 3D
              ostringstream title;


              title << "Run: " << event.id().run()
              << "  Subrun: " << event.id().subRun()
              << "  Event: " << event.id().event()<<".root";


              text.SetTextAlign(11);
              text.DrawTextNDC( 0., 0.01, title.str().c_str() );
              canvas_->Modified();

              canvas_->SaveAs(title.str().c_str());

              if ( clickToAdvance_ ){
                cerr << "Double click in the Canvas " << moduleLabel_ << " to continue:" ;
                gPad->WaitPrimitive();
                    } else{
                char junk;
                cerr << "Enter any character to continue: ";
                cin >> junk;
              }
                      cerr << endl;
            }

           }//end display
      }//end 3d


      void CosmicFitDisplay::plot3dPrimes(const art::Event& event){ //local coords

        _evt = event.id().event();  // add event id

        auto comboHits  = event.getValidHandle<ComboHitCollection>( _chtag );
        auto Tracks  = event.getValidHandle<CosmicTrackSeedCollection>( _sttag );
        GeomHandle<Tracker> tracker;

        // Annulus of a cylinder that bounds the tracker
        TubsParams envelope(tracker->g4Tracker()->getInnerTrackerEnvelopeParams());

        std::vector<double> x, y, z, a0, a1, b0, b1;
        std::vector<XYZVectorF> xprimes, yprimes, zprimes, initial_track_direction;
        x.reserve(comboHits->size());
        y.reserve(comboHits->size());
        z.reserve(comboHits->size());
        a1.reserve(Tracks->size());
        b1.reserve(Tracks->size());
        a0.reserve(Tracks->size());
        b0.reserve(Tracks->size());

        //Get track coordinate system and fit parameters
        for(auto const& track: *Tracks){
        xprimes.push_back(track._track.FitCoordSystem._XDoublePrime);
        yprimes.push_back(track._track.FitCoordSystem._YDoublePrime);
        zprimes.push_back(track._track.FitCoordSystem._ZPrime);
        initial_track_direction.push_back(track._track.InitParams.Direction() );
        a1.push_back(track._track.FitParams.A0);
        a0.push_back(track._track.FitParams.A1);
        b1.push_back(track._track.FitParams.B0);
        b0.push_back(track._track.FitParams.B1);
        }

        if(xprimes.size() >0){
        for(auto const& chit : *comboHits){

        x.push_back(chit.pos().Dot(xprimes[0]));
        y.push_back(chit.pos().Dot(yprimes[0]));
        z.push_back(chit.pos().Dot(zprimes[0]));
        }
       double minz = *std::min_element(z.begin(), z.end());
      double maxz = *std::max_element(z.begin(), z.end());
      double minx = *std::min_element(x.begin(), x.end());
      double maxx = *std::max_element(x.begin(), x.end());
      double miny = *std::min_element(y.begin(), y.end());
      double maxy = *std::max_element(y.begin(), y.end());
        if (doDisplay_) {
              std::cout << "Run: " << event.id().run()
           << "  Subrun: " << event.id().subRun()
           << "  Event: " << event.id().event()<<std::endl;
              TTUBE tube;
              TText  text;
              TPolyMarker3D poly3D;

              // Draw the frame for the cylinders in plot:
              //double plotLimits(1000.);
              //double zlimit{envelope.zHalfLength()};
              //Create Pad:
              auto pad = canvas_;
              pad->Clear();
              canvas_->Draw();

              TH3D *xyzplot =new TH3D("XYZ","XYZ",1,minz,maxz+500,1,-miny-500,maxy+500,1, minx-500, maxx+500);

              xyzplot->GetYaxis()->SetTitleOffset(1.5);
              xyzplot->GetXaxis()->SetTitleOffset(1.5);
                 xyzplot->GetZaxis()->SetTitleOffset(1.5);
              xyzplot->SetTitle( "Visualization in X''Y''Z';Z' (mm);Y'' (mm); X''(mm)");
              xyzplot->SetMarkerStyle(2);
              xyzplot->SetMarkerSize(0.65);
              xyzplot->SetMarkerColor(kRed);
              xyzplot->SetStats(0);
              xyzplot->Draw();

              poly3D.SetMarkerStyle(2);
              poly3D.SetMarkerSize(0.65);
              poly3D.SetMarkerColor(kRed);
              int index = 1;
              for ( auto const& chit : *comboHits ){
                        TPolyLine3D *errors = new TPolyLine3D;
                        errors->SetLineColor(kRed);
                        auto const& p = chit.pos();
                        double xdoubleprime{p.Dot(xprimes[0])};
                        double ydoubleprime{p.Dot(yprimes[0])};
                        double zprime{p.Dot(zprimes[0])};
                        auto w = chit.uDir();
                        auto const& s = chit.wireRes();
                        double x1 = p.Dot(xprimes[0])+(s*w).Dot(xprimes[0]);
                        double x2 = p.Dot(xprimes[0])-(s*w).Dot(xprimes[0]);
                        double z1 = p.Dot(zprimes[0])+(s*w).Dot(zprimes[0]);
                        double z2 = p.Dot(zprimes[0])-(s*w).Dot(zprimes[0]);
                        double y1 = p.Dot(yprimes[0])+(s*w).Dot(yprimes[0]);
                        double y2 = p.Dot(yprimes[0])-(s*w).Dot(yprimes[0]);
                        poly3D.SetPoint(index, zprime, ydoubleprime, xdoubleprime );
                        errors->SetPoint(0, z1,y1,x1);
                        errors->SetNextPoint(z2,y2,x2);
                        index+=1;
                        errors->Draw("same");
                    }
              poly3D.Draw("same");
              if(a1.size() > 0 && b1.size() > 0){
                      TPolyLine3D *Track = new TPolyLine3D;
                Track->SetLineColor(kBlue);

                double tz1 = minz;
                double tz2 = maxz;

                double tx1 = a0[0]+ a1[0]*tz1;
                double tx2 = a0[0] + a1[0]*tz2;
                double ty1 = b0[0] + b1[0]*tz1;
                double ty2 = b0[0] + b1[0]*tz2;

                Track->SetPoint(0, tz1, ty1, tx1);
                Track->SetNextPoint(tz2, ty2, tx2);
                Track->Draw("same");

              }
              canvas_->Update();
              //END 3D
              ostringstream title;


              title << "Run: " << event.id().run()
              << "  Subrun: " << event.id().subRun()
              << "  Event: " << event.id().event()<<".root";


              text.SetTextAlign(11);
              text.DrawTextNDC( 0., 0.01, title.str().c_str() );
              canvas_->Modified();

              canvas_->SaveAs(title.str().c_str());

              if ( clickToAdvance_ ){
                cerr << "Double click in the Canvas " << moduleLabel_ << " to continue:" ;
                gPad->WaitPrimitive();
                    } else{
                char junk;
                cerr << "Enter any character to continue: ";
                cin >> junk;
              }
                      cerr << endl;
            }
           }//end display
      }//end 3d

      //TODO :  below function is unfinished
      void CosmicFitDisplay::Geom3D(const art::Event& event){

          TGeoManager *geom = new TGeoManager("tracker_geom", "tracker_geom");
          //Materials/Media
          TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
          TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
          //Set Top Geom:
          TGeoVolume *top = geom->MakeBox("3DGeometry", Vacuum, 2000,2000,2000);
           geom->SetTopVolume(top);

           TGeoVolume *rootbox = geom->MakeBox("ROOT", Vacuum, 110., 50., 5.);
           rootbox->SetVisibility(kFALSE);
          //Rotate about z
          TGeoRotation *rot = new TGeoRotation("rot",0,90,0.);
          rot->RegisterYourself();

          GeomHandle<Tracker> tracker;

          // Annulus of a cylinder that bounds the tracker
          TubsParams envelope(tracker->g4Tracker()->getInnerTrackerEnvelopeParams());

          // Draw the frame for the cylinders in plot:
          double zlimit{envelope.zHalfLength()};
          double R_min{envelope.innerRadius()};
          double R_max{envelope.outerRadius()};


          TGeoVolume *tr = geom->MakeTube("tr", Vacuum, R_min, R_max, zlimit);
           tr->SetLineColor(kViolet);

_evt = event.id().event();  // add event id

            //get combo hits
           auto comboHits  = event.getValidHandle<ComboHitCollection>( _chtag );
           auto Tracks  = event.getValidHandle<CosmicTrackSeedCollection>( _sttag );



           // Create arrays for x,y,z:
           std::vector<double> x, y, z,a0,a1,b0,b1;
           x.reserve(comboHits->size());
           y.reserve(comboHits->size());
           z.reserve(comboHits->size());

          a1.reserve(Tracks->size());
          b1.reserve(Tracks->size());
          a0.reserve(Tracks->size());
          b0.reserve(Tracks->size());
          // loop over combo hits
          for(auto const& chit : *comboHits){
             x.push_back(chit.pos().x());
             y.push_back(chit.pos().y());
             z.push_back(chit.pos().z());
          }
          //loop over tracks:
          for(auto const& track: *Tracks){
             a1.push_back(track._track.FitParams.A1);
             a0.push_back(track._track.FitParams.A0);

             b1.push_back(track._track.FitParams.B1);
             b0.push_back(track._track.FitParams.B0);
          }

        //If Diag:
        if (doDisplay_) {
              std::cout << "Run: " << event.id().run()
           << "  Subrun: " << event.id().subRun()
           << "  Event: " << event.id().event()<<std::endl;


           //Add Node
           top->AddNode(tr,1, rot);


          //Close the geometry
          //geom->CloseGeometry();
          geom->SetVisLevel(4);
          top->Draw("ogle");
          if ( clickToAdvance_ ){
                cerr << "Double click in the Canvas " << moduleLabel_ << " to continue:" ;
                gPad->WaitPrimitive();
                    } else{
                char junk;
                cerr << "Enter any character to continue: ";
                cin >> junk;
              }
                      cerr << endl;

      }
   }

   std::vector<double> CosmicFitDisplay::GetMaxAndMin(std::vector<double> myvector){
        std::vector<double> MaxAndMin;
        double first = myvector[0];
            double smallest = first;
        double biggest = first;
            for (unsigned i=0; i<  myvector.size(); ++i){
        double element = myvector[i];
        if (element< smallest) {
            smallest = element;
         }
            }
        for (unsigned i=0; i<  myvector.size(); ++i){
        double element = myvector[i];
         if (element> biggest) {
            biggest = element;
         }
            }
    MaxAndMin.push_back(biggest);
    MaxAndMin.push_back(smallest);
    return MaxAndMin;
        }

bool CosmicFitDisplay::findData(const art::Event& evt){
        _chcol = 0;
        _coscol = 0;
        auto chH = evt.getValidHandle<ComboHitCollection>(_chtag);
        _chcol = chH.product();
        auto stH = evt.getValidHandle<CosmicTrackSeedCollection>(_sttag);
        _coscol =stH.product();

        return _chcol != 0 && _coscol !=0 ;
       }

}//End Namespace Mu2e

using mu2e::CosmicFitDisplay;
DEFINE_ART_MODULE(CosmicFitDisplay)



