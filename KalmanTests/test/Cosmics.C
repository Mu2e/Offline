void Cosmics(TTree* cr, const char* page="prod") {
  TString tpage(page);

  if(tpage=="prod"){
    TH2F* pxyep = new TH2F("pxyep","Electron production position;x(mm); y(mm)",100,-800,800,100,-800,800);
    TH2F* pxyem = new TH2F("pxyem","Electron production position;x(mm); y(mm)",100,-800,800,100,-800,800);
  
    pxyep->SetMarkerStyle(8);
    pxyep->SetMarkerColor(kBlue);
    pxyep->SetStats(0);
    pxyem->SetMarkerStyle(8);
    pxyem->SetMarkerColor(kRed);
    pxyem->SetStats(0);
    
    cr->Project("pxyem","ey:ex","umcpdgid==11");
    cr->Project("pxyep","ey:ex","umcpdgid==-11");

    TH2F* pxzep = new TH2F("pxzep","Electron production position;z(mm); x(mm)",100,1700,3600,100,-800,800);
    TH2F* pxzem = new TH2F("pxzem","Electron production position;z(mm); x(mm)",100,1700,3600,100,-800,800);
  
    pxzep->SetMarkerStyle(8);
    pxzep->SetMarkerColor(kBlue);
    pxzep->SetStats(0);
    pxzem->SetMarkerStyle(8);
    pxzem->SetMarkerColor(kRed);
    pxzem->SetStats(0);
  

    cr->Project("pxzem","ex:ez","umcpdgid==11");
    cr->Project("pxzep","ex:ez","umcpdgid==-11");

    TH2F* pxzmp = new TH2F("pxzmp","Projected #mu production position y=0;z(mm); x(mm)",100,0,4500,100,-2000,2000);
    TH2F* pxzmm = new TH2F("pxzmm","Projected #mu production position y=0;z(mm); x(mm)",100,0,4500,100,-2000,2000);
  
    pxzmp->SetMarkerStyle(8);
    pxzmp->SetMarkerColor(kCyan);
    pxzmp->SetStats(0);
    pxzmm->SetMarkerStyle(8);
    pxzmm->SetMarkerColor(kOrange);
    pxzmm->SetStats(0);
  

    cr->Project("pxzmm","ppx:ppz","umcpdgid==13");
    cr->Project("pxzmp","ppx:ppz","umcpdgid==-13");

    TH1F* mcmomem = new TH1F("mcmomem","True momentum;MeV",50,50,230);
    TH1F* mcmomep = new TH1F("mcmomep","True momentum;MeV",50,50,230);
    TH1F* mcmommm = new TH1F("mcmommm","True momentum;MeV",50,50,230);
    TH1F* mcmommp = new TH1F("mcmommp","True momentum;MeV",50,50,230);
    mcmomem->SetStats(0);
    mcmomem->SetLineColor(kRed);
    mcmomep->SetStats(0);
    mcmomep->SetLineColor(kBlue);
    mcmommm->SetStats(0);
    mcmommm->SetLineColor(kOrange);
    mcmommp->SetStats(0);
    mcmommp->SetLineColor(kCyan);

    cr->Project("mcmomem","umcmom","umcpdgid==11");
    cr->Project("mcmomep","umcmom","umcpdgid==-11");
    cr->Project("mcmommm","umcmom","umcpdgid==13");
    cr->Project("mcmommp","umcmom","umcpdgid==-13");
/*
    TH1F* rtandem = new TH1F("rtandem","Reco tan#lambda",50,-1.2,-0.3);
    TH1F* rtandep = new TH1F("rtandep","Reco tan#lambda",50,-1.2,-0.3);
    TH1F* rtandmm = new TH1F("rtandmm","Reco tan#lambda",50,-1.2,-0.3);
    TH1F* rtandmp = new TH1F("rtandmp","Reco tan#lambda",50,-1.2,-0.3);
    rtandem->SetStats(0);
    rtandem->SetLineColor(kRed);
    rtandep->SetStats(0);
    rtandep->SetLineColor(kBlue);
    rtandmm->SetStats(0);
    rtandmm->SetLineColor(kOrange);
    rtandmp->SetStats(0);
    rtandmp->SetLineColor(kCyan);

    cr->Project("rtandem","utd","umcpdgid==11");
    cr->Project("rtandep","utd","umcpdgid==-11");
    cr->Project("rtandmm","utd","umcpdgid==13");
    cr->Project("rtandmp","utd","umcpdgid==-13");

    TH1F* rd0em = new TH1F("rd0em","Reco d0;mm",50,-400,400);
    TH1F* rd0ep = new TH1F("rd0ep","Reco d0;mm",50,-400,400);
    TH1F* rd0mm = new TH1F("rd0mm","Reco d0;mm",50,-400,400);
    TH1F* rd0mp = new TH1F("rd0mp","Reco d0;mm",50,-400,400);

    rd0em->SetStats(0);
    rd0em->SetLineColor(kRed);
    rd0ep->SetStats(0);
    rd0ep->SetLineColor(kBlue);
    rd0mm->SetStats(0);
    rd0mm->SetLineColor(kOrange);
    rd0mp->SetStats(0);
    rd0mp->SetLineColor(kCyan);

    cr->Project("rd0em","ud0","umcpdgid==11");
    cr->Project("rd0ep","ud0","umcpdgid==-11");
    cr->Project("rd0mm","ud0","umcpdgid==13");
    cr->Project("rd0mp","ud0","umcpdgid==-13");

    TH1F* dummy = new TH1F("dummy","",10,-1,1);
    dummy->SetStats(0);
*/
    TLegend* leg = new TLegend(0.7,0.5,0.9,0.9);
    leg->AddEntry(mcmomem,"e^{-}","L");
    leg->AddEntry(mcmomep,"e^{+}","L");
    leg->AddEntry(mcmommm,"#mu^{-}","L");
    leg->AddEntry(mcmommp,"#mu^{+}","L");

    TCanvas* pcan = new TCanvas("pcan","Production",1200,800);

    pcan->Divide(2,2);
    pcan->cd(1);
    mcmommp->Draw();
    mcmomem->Draw("same");
    mcmomep->Draw("same");
    mcmommm->Draw("same");
    leg->Draw();

    pcan->cd(2);
    pxyem->Draw();
    pxyep->Draw("same");

    pcan->cd(3);
    pxzem->Draw();
    pxzep->Draw("same");

/*    pcan->cd(4);
    rtandem->Draw();
    rtandep->Draw("same");
    rtandmm->Draw("same");
    rtandmp->Draw("same");

    pcan->cd(5);
    rd0em->Draw();
    rd0ep->Draw("same");
    rd0mm->Draw("same");
    rd0mp->Draw("same");

    pcan->cd(6);
    dummy->Draw();
*/
    pcan->cd(4);
    pxzmm->Draw();
    pxzmp->Draw("same");


  } else if(tpage =="select") {
    TH2F* cmom = new TH2F("cmom","Downstream vs Upstream momentum;P_{d} (MeV);P_{u} (MeV)",100,60,200,100,60,200);
    TH2F* ctand = new TH2F("ctand","Downstream vs Upstream tan#lambda;tan(#lambda)_{d};tan(#lambda)_{u}",100,-1.25,-0.25,100,0.25,1.25);
    TH2F* cp0 = new TH2F("cp0","Downstream vs Upstream #phi_{0};radians;radians",100,-3.15,3.15,100,-3.15,3.15);
    TH2F* cd0 = new TH2F("cd0","Downstream vs Upstream d_{0};mm;mm",100,-400,400,100,-400,400);
    TH1F* dt0 = new TH1F("dt0","Downstream - Upstream t_{0};nsec",100,50,120);
    cmom->SetStats(0);
    ctand->SetStats(0);
    cd0->SetStats(0);
    cp0->SetStats(0);
    dt0->SetStats(0);
    
    cr->Project("cmom","dmom:umom");
    cr->Project("ctand","dtd:utd");
    cr->Project("cp0","dp0:up0");
    cr->Project("cd0","dd0:ud0");
    cr->Project("dt0","dt0-ut0");

    TCanvas* scan = new TCanvas("scan","Selection",1200,800);
    scan->Divide(3,2);
    scan->cd(1);
    cmom->Draw();
    scan->cd(2);
    ctand->Draw();
    scan->cd(3);
    cp0->Draw();
    scan->cd(4);
    cd0->Draw();
    scan->cd(5);
    dt0->Draw();


  }



}
