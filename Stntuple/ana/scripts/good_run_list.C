//-----------------------------------------------------------------------------
//  calculate lumi for the ETF definition of a good run list
//  results:  lumi sum = 85630.539062(all)   50814.324219(good runs)
//  Grl = "ETF" or "WZ_PRD"
//-----------------------------------------------------------------------------

void good_run_list(TString Grl="ETF", int RMin=138425, int RMax=186598, int PrintAll = 0) {

  static TStnGoodRunList* grl = NULL;

  int *     grl_function(int);

  if (! grl) {
    grl = new TStnGoodRunList(Grl.Data());
    grl->Init();
  }
  else if (Grl != grl->GetName()) {
    delete grl;
    grl = new TStnGoodRunList(Grl.Data());
    grl->Init();
  }

  int good_run;

  TStnRunSummary*  rs;
  
  TString* tt;

  float lum_sum_onl    = 0;
  float lum_sum_onl_gr = 0;

  float lum_sum_ofl    = 0;
  float lum_sum_ofl_gr = 0;

  float lum_sum_tap    = 0;
  float lum_sum_tap_gr = 0;

  printf("-------------------------------------------------------------------");
  printf("------------------------------------------------------------\n");
  printf("                                          ------------ online -------   --- offline --------- \n");
  printf(" run number    trigger  table             C L L L C C C C C S C P S I L   C C C C C S C P S I L");
  printf(" good   onl_lumi    onl_sum   off_lumi    ofl_sum \n");      	    
  printf("                                          L 1 2 3 A O M M M M E E V S 0   A 0 M M M M E E V S 0\n");
  printf("                                          C       L T U P X X S S X L 0   L T U P X X S S X L 0\n");
  printf("-------------------------------------------------------------------");
  printf("------------------------------------------------------------\n");

  for (int run=RMin; run<=RMax; run++) {
    rs       = grl->GetRunSummary(run);
    //    printf(" run, rs = %7i %08x\n",run,rs);
    if (rs) {
      tt       = &rs->TriggerTableName();

      good_run = grl->GoodRun(run);
//-----------------------------------------------------------------------------
// do something, calculate luminosity, for example
//-----------------------------------------------------------------------------
      lum_sum_tap += rs->LumiTape();
      lum_sum_onl += rs->OnlineLumiRS();
      lum_sum_ofl += rs->OfflineLumiRS();
      if (good_run || PrintAll) {
	lum_sum_tap_gr += rs->LumiTape();
	lum_sum_onl_gr += rs->OnlineLumiRS();
	lum_sum_ofl_gr += rs->OfflineLumiRS();

	printf(" %7i   %-30s",run,tt.Data());
	printf(" %1i %1i %1i %1i %1i %1i %1i %1i %1i %1i %1i %1i %1i %1i %1i   %1i %1i %1i %1i %1i %1i %1i %1i %1i %1i %1i    %1i",
	       rs->ClcStatusBit(),
	       rs->L1tStatusBit(),
	       rs->L2tStatusBit(),
	       rs->L3tStatusBit(),
	       rs->CalStatusBit(),
	       rs->CotStatusBit(),
	       rs->CmuStatusBit(),
	       rs->CmpStatusBit(),
	       rs->CmxStatusBit(),
	       rs->SmxStatusBit(), 
	       rs->CesStatusBit(), 
	       rs->PesStatusBit(), 
	       rs->SvxStatusBit(),
	       rs->IslStatusBit(),
	       rs->L00StatusBit(),
	       
	       rs->CalOfflineBit(),
	       rs->CotOfflineBit(),
	       rs->CmuOfflineBit(),
	       rs->CmpOfflineBit(),
	       rs->CmxOfflineBit(),
	       rs->SmxOfflineBit(),
	       rs->CesOfflineBit(),
	       rs->PesOfflineBit(),
	       rs->SvxOfflineBit(),
	       rs->IslOfflineBit(),
	       rs->L00OfflineBit(),
	       good_run            );
	printf(" %10.3f %10.3f",rs->OnlineLumiRS (),lum_sum_onl_gr);
	printf(" %10.3f %10.3f",rs->OfflineLumiRS(),lum_sum_ofl_gr);
	printf("\n");
      }
    }
  NEXT_RUN:;
  }
  printf(" TAPE    lumi in all  runs = %10.5f pb^-1\n",lum_sum_tap/1000.);
  printf(" TAPE    lumi in good runs = %10.5f pb^-1\n",lum_sum_tap_gr/1000.);
  printf(" ONLINE  lumi in all  runs = %10.5f pb^-1\n",lum_sum_onl/1000.);
  printf(" ONLINE  lumi in good runs = %10.5f pb^-1\n",lum_sum_onl_gr/1000.);
  printf(" OFFLINE lumi in all  runs = %10.5f pb^-1\n",lum_sum_ofl/1000.);
  printf(" OFFLINE lumi in good runs = %10.5f pb^-1\n",lum_sum_ofl_gr/1000.);
}
