//-----------------------------------------------------------------------------
// file: pbar2m/scripts/display_phi.C
//
// test macro : fill and display histogram of N(combo hits) in the event
// no includes needed
// hide variables - dos it have any effect ?
//-----------------------------------------------------------------------------
namespace {
  TH1F*               h_ch       (nullptr);	       // to be accessed interactively
  TH1F*               h_tc[100];	       // to be accessed interactively
  TCanvas*            c          (nullptr);
}

mu2e::MuHitDisplay* m_disp (nullptr);
//-----------------------------------------------------------------------------
// Mode = 0: begin job (run)
// Mode = 1: event
// Mode = 2: end job (run)
//-----------------------------------------------------------------------------
void display_phi(int Mode, TModule* Module) {
  printf("display_001 called: Mode = %i, Module = %p\n",Mode,Module);

  TStnVisManager* vm = TStnVisManager::Instance();

  m_disp = (mu2e::MuHitDisplay*) Module;

  if (Mode == 0) {
//-----------------------------------------------------------------------------
// book histograms , perform other initializations
//-----------------------------------------------------------------------------
    c = new TCanvas;
    c->cd();
    // h_ch    = new TH1F ("hist","hist",32,0,6.4);
    h_ch    = new TH1F ("hist","hist",21,0,6.3);
    h_ch->Draw();

    for (int i=0; i<100; i++) {
      // h_tc[i] = new TH1F(Form("h_tc_%02i",i),Form("Time cluster %02i",i),32,0,6.4);
      h_tc[i] = new TH1F(Form("h_tc_%02i",i),Form("Time cluster %02i",i),21,0,6.3);
      h_tc[i]->SetLineColor(kRed+2);
      h_tc[i]->SetFillColor(kRed+2);
      h_tc[i]->SetFillStyle(3004);
    }
  }
  else if (Mode == 1) {
//-----------------------------------------------------------------------------
// fill histograms
//-----------------------------------------------------------------------------
    TEvdTimeClusterVisNode* vn = (TEvdTimeClusterVisNode*) vm->FindNode("TimeClusterVisNode");
    const mu2e::ComboHitCollection*    chColl = vn->ChColl();
    const mu2e::TimeClusterCollection* tcColl = vn->TcColl();

    int nh = chColl->size();
    printf("display_phi: nh:  %i\n",nh);

    c->cd();

    h_ch->Reset();
    for (int i=0; i<nh; i++) {
      const mu2e::ComboHit* ch = &chColl->at(i);
      float phi = ch->phi();
      if (phi < 0) phi = phi+2*M_PI;
      h_ch->Fill(phi);
    }

    h_ch->Draw();

    int ntc = tcColl->size();
    printf("display_phi: ntc:  %i\n",ntc);

    for (int itc=0; itc<ntc; itc++) {
      const mu2e::TimeCluster* tc = &tcColl->at(itc);

      h_tc[itc]->Reset();
      int n2 = tc->nhits();
      printf("display_phi: TC n2:  %i\n",n2);

      for (int ih=0; ih<n2; ih++) {
        const mu2e::ComboHit* ch = &chColl->at(tc->hits().at(ih));
        float phi = ch->phi();
        if (phi < 0) phi = phi+2*M_PI;
        h_tc[itc]->Fill(phi);
      }

      h_tc[itc]->Draw("sames");
    }

    c->Modified();
    c->Update();
  }
  else if (Mode == 2) {
//-----------------------------------------------------------------------------
// end run
//-----------------------------------------------------------------------------
  }
}
