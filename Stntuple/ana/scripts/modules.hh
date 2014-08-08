#ifndef stntuple_ana_scripts_modules_hh
#define stntuple_ana_scripts_modules_hh

class TClcAnaModule;
class TDebugModule;
class TLumiMonModule;
class TCalAnaModule;
class TCesAnaModule;
class TClusterAnaModule;
class TConversionFilterModule;
class TEmFilterModule;
class TFwdDetAnaModule;
class TJetAnaModule;
class TMcAnaModule;
class TMetAnaModule;
class TMuoAnaModule;
class TPesAnaModule;
class TPhotonAnaModule;
class TSvtAnaModule;
class TTrackAnaModule;
class TTrigAnaModule;
class TWenuMonModule;
class TJpsiMonModule;
class TValidationModule;

namespace stntuple {
  TClcAnaModule*           m_clc   = NULL;
  TCalAnaModule*           m_cal   = NULL;
  TCesAnaModule*           m_ces   = NULL;
  TClusterAnaModule*       m_clu   = NULL;
  TConversionFilterModule* m_cnv   = NULL;
  TDebugModule*            m_dbg   = NULL;
  TDFCModule*              m_dfc   = NULL;
  TEmFilterModule*         m_emf   = NULL;
  TFwdDetAnaModule*        m_fwd   = NULL;
  TJetAnaModule*           m_jet   = NULL;
  TLumiMonModule*          m_lum   = NULL;
  TMcAnaModule*            m_mc    = NULL;
  TMetAnaModule*           m_met   = NULL;   
  TMuoAnaModule*           m_muo   = NULL;
  TPesAnaModule*           m_pes   = NULL;
  TPhotonAnaModule*        m_pho   = NULL;
  TSvtAnaModule*           m_svt   = NULL;
  TTrackAnaModule*         m_trk   = NULL;
  TTrigAnaModule*          m_trg   = NULL;
  
  TWenuMonModule*          m_wen   = NULL;
  TJpsiMonModule*          m_jps   = NULL;

};


TValidationModule*       m_val   = NULL;
#endif
