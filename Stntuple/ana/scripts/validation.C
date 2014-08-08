///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "Stntuple/ana/scripts/modules.hh"

def_name validation_001("validation");
def_name validation_002("wenu_filter");
def_name validation_003("ces_only");
def_name validation_003("debug_ntrk");
def_name val_stn("val_stn");
//-----------------------------------------------------------------------------
// standard validation studies
//-----------------------------------------------------------------------------
void validation() {
  stntuple::m_cal = (TCalAnaModule*) g.x->AddModule("TCalAnaModule",0);
  stntuple::m_ces = (TCesAnaModule*) g.x->AddModule("TCesAnaModule",0);
  stntuple::m_clu = (TClusterAnaModule*) g.x->AddModule("TClusterAnaModule",0);
  stntuple::m_cnv = (TConversionFilterModule*) g.x->AddModule("TConversionFilterModule",0);
  stntuple::m_emf = (TEmFilterModule*) g.x->AddModule("TEmFilterModule",0);
  stntuple::m_fwd = (TFwdDetAnaModule*) g.x->AddModule("TFwdDetAnaModule",0);
  stntuple::m_jet = (TJetAnaModule*) g.x->AddModule("TJetAnaModule",0);
  stntuple::m_lum = (TLumiMonModule*) g.x->AddModule("TLumiMonModule",0);
  stntuple::m_mc  = (TMcAnaModule*) g.x->AddModule("TMcAnaModule",0);
  stntuple::m_met = (TMetAnaModule*) g.x->AddModule("TMetAnaModule",0);
  stntuple::m_muo = (TMuoAnaModule*) g.x->AddModule("TMuoAnaModule",0);
  stntuple::m_pes = (TPesAnaModule*) g.x->AddModule("TPesAnaModule",0);
  stntuple::m_pho = (TPhotonAnaModule*) g.x->AddModule("TPhotonAnaModule",0);
  stntuple::m_svt = (TSvtAnaModule*) g.x->AddModule("TSvtAnaModule",0);
  stntuple::m_trk = (TTrackAnaModule*) g.x->AddModule("TTrackAnaModule",0);
  stntuple::m_trg = (TTrigAnaModule*) g.x->AddModule("TTrigAnaModule",0);
//-----------------------------------------------------------------------------
// WenuMonModule: require ELECTRON_CENTRAL_18 trigger path
//-----------------------------------------------------------------------------
  m_wen = (TWenuMonModule*) g.x->AddModule("TWenuMonModule",0);
  m_wen->AddL3TriggerName("ELECTRON_CENTRAL_18");
//-----------------------------------------------------------------------------
// JpsiMonModule: require event to pass one of the stream A dimuon trigger paths
//-----------------------------------------------------------------------------
  stntuple::m_jps = (TJpsiMonModule*) g.x->AddModule("TJpsiMonModule",0);
  stntuple::m_jps->AddL3TriggerName("JPSI_CMUP4_CMU1.5");
  stntuple::m_jps->AddL3TriggerName("JPSI_CMUCMU1.5");
  stntuple::m_jps->AddL3TriggerName("JPSI_CMUP4_CMX2");
  stntuple::m_jps->AddL3TriggerName("JPSI_CMU1.5_CMX2");
//---------------------------------------------------------------------------
// ValidatioModule: produce a set of histograms for doing the validation 
//---------------------------------------------------------------------------
  stntuple::m_val = (TValidationModule*) g.x->AddModule("TValidationModule",0);
}
//-----------------------------------------------------------------------------
// Wenu module only
//-----------------------------------------------------------------------------
void wenu_filter(int MinNTightElectrons) {
//-----------------------------------------------------------------------------
// WenuMonModule: require ELECTRON_CENTRAL_18 trigger path
//-----------------------------------------------------------------------------
  stntuple::m_wen = (TWenuMonModule*) g.x->AddModule("TWenuMonModule",1);
  stntuple::m_wen->AddL3TriggerName("ELECTRON_CENTRAL_18");
  stntuple::m_wen->SetMinNTightElectrons(MinNTightElectrons);
}


void ces_only(int PrintLevel) {
//-----------------------------------------------------------------------------
// WenuMonModule: require ELECTRON_CENTRAL_18 trigger path
//-----------------------------------------------------------------------------
  char s[100];
  if (PrintLevel > 0) {
    sprintf(s,"%i",PrintLevel);
    gSystem->Setenv("CesAna_PrintLevel",s);
  }
  stntuple::m_ces = (TCesAnaModule*) g.x->AddModule("TCesAnaModule",1);
}

void debug_ntrk(int Set, int NMin, int NMax) {
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  m_trk = (TTrackAnaModule*) g.x->AddModule("TTrackAnaModule",1);
  gSystem->Setenv("TrackAna_PrintLevel","101");
  stntuple::m_trk->SetSet_101(Set);
  stntuple::m_trk->SetMinNTracks_101(NMin);
  stntuple::m_trk->SetMaxNTracks_101(NMax);
}

void val_stn(int PdgCode = 11, int GeneratorCode = 2) {
//-----------------------------------------------------------------------------
// configure validation module
//-----------------------------------------------------------------------------
  m_val = (TValidationModule*) g.x->AddModule("TValidationModule",0);  
  m_val->SetPdgCode      (PdgCode);
  m_val->SetGeneratorCode(GeneratorCode);
}
