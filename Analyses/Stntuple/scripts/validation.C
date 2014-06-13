///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "global_vars.h"

def_name validation_001("validation");
def_name validation_002("wenu_filter");
def_name validation_003("ces_only");
def_name validation_003("debug_ntrk");
//-----------------------------------------------------------------------------
// standard validation studies
//-----------------------------------------------------------------------------
void validation() {
  g.x->AddModule(m_cal,0);
  g.x->AddModule(m_ces,0);
  g.x->AddModule(m_clu,0);
  g.x->AddModule(m_cnv,0);
  g.x->AddModule(m_emf,0);
  g.x->AddModule(m_fwd,0);
  g.x->AddModule(m_jet,0);
  g.x->AddModule(m_lum,0);
  g.x->AddModule(m_mc ,0);
  g.x->AddModule(m_met,0);
  g.x->AddModule(m_muo,0);
  g.x->AddModule(m_pes,0);
  g.x->AddModule(m_pho,0);
  g.x->AddModule(m_svt,0);
  g.x->AddModule(m_trk,0);
  g.x->AddModule(m_trg,0);
//-----------------------------------------------------------------------------
// WenuMonModule: require ELECTRON_CENTRAL_18 trigger path
//-----------------------------------------------------------------------------
  g.x->AddModule(m_wen,0);
  m_wen->AddL3TriggerName("ELECTRON_CENTRAL_18");
//-----------------------------------------------------------------------------
// JpsiMonModule: require event to pass one of the stream A dimuon trigger paths
//-----------------------------------------------------------------------------
  g.x->AddModule(m_jps,0);
  m_jps->AddL3TriggerName("JPSI_CMUP4_CMU1.5");
  m_jps->AddL3TriggerName("JPSI_CMUCMU1.5");
  m_jps->AddL3TriggerName("JPSI_CMUP4_CMX2");
  m_jps->AddL3TriggerName("JPSI_CMU1.5_CMX2");
}
//-----------------------------------------------------------------------------
// Wenu module only
//-----------------------------------------------------------------------------
void wenu_filter(int MinNTightElectrons) {
//-----------------------------------------------------------------------------
// WenuMonModule: require ELECTRON_CENTRAL_18 trigger path
//-----------------------------------------------------------------------------
  g.x->AddModule(m_wen,1);
  m_wen->AddL3TriggerName("ELECTRON_CENTRAL_18");
  m_wen->SetMinNTightElectrons(MinNTightElectrons);
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
  g.x->AddModule(m_ces,1);
}

void debug_ntrk(int Set, int NMin, int NMax) {
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  g.x->AddModule(m_trk,1);
  gSystem->Setenv("TrackAna_PrintLevel","101");
  m_trk->SetSet_101(Set);
  m_trk->SetMinNTracks_101(NMin);
  m_trk->SetMaxNTracks_101(NMax);
}
