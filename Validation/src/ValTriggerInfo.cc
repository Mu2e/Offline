
#include "Offline/Validation/inc/ValTriggerInfo.hh"


int mu2e::ValTriggerInfo::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hccs = tfs.make<TH1D>( "Ncss", "N caloClusters", 21, -0.5, 20.5);
  _htrk = tfs.make<TH1D>( "Ntrk", "N tracks", 21, -0.5, 20.5);
  _hhel = tfs.make<TH1D>( "Nhel", "N helixes", 21, -0.5, 20.5);
  _hhit = tfs.make<TH1D>( "Nhit", "N hitClusters", 21, -0.5, 20.5);
  _hcts = tfs.make<TH1D>( "Ncts", "N caloTrigSeeds", 21, -0.5, 20.5);
  _hcos = tfs.make<TH1D>( "Ncos", "N cosmics", 21, -0.5, 20.5);

  return 0;
}

int mu2e::ValTriggerInfo::fill(const mu2e::TriggerInfo & obj,
                                art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hccs->Fill(obj.caloClusters().size());
  _htrk->Fill(obj.tracks().size());
  _hhel->Fill(obj.helixes().size());
  _hhit->Fill(obj.hitClusters().size());
  _hcts->Fill(obj.caloTrigSeeds().size());
  _hcos->Fill(obj.cosmics().size());

  return 0;
}
