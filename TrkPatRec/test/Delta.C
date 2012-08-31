void Delta(TTree* ddiag, const char* page="rho") {
  TString spage(page);
  TCut con("mcmom>100&&mctd<1.0&&(nconv>=10||nconv/nsh>0.5)");
  TCut bkg("nconv<=2||nconv/nsh<0.05");
  TCut rhocut("prho>410&&prho<680");

  if(spage == "rho"){
    TH1F* prhocon = new TH1F("prhocon","Peak #rho;#rho (mm)",100,380,720);
    TH1F* prhobkg = new TH1F("prhobkg","Peak #rho;#rho (mm)",100,380,720);
    TH1F* mrhocon = new TH1F("mrhocon","Median #rho;#rho (mm)",100,380,720);
    TH1F* mrhobkg = new TH1F("mrhobkg","Median #rho;#rho (mm)",100,380,720);
    prhocon->SetLineColor(kRed);
    prhobkg->SetLineColor(kBlue);
    mrhocon->SetLineColor(kRed);
    mrhobkg->SetLineColor(kBlue);

    prhocon->SetStats(0);
    prhobkg->SetStats(0);
    mrhocon->SetStats(0);
    mrhobkg->SetStats(0);

    ddiag->Project("prhocon","prho",con);
    ddiag->Project("prhobkg","prho",bkg);
    ddiag->Project("mrhocon","rmed",con);
    ddiag->Project("mrhobkg","rmed",bkg);

    prhocon->Scale(10);
    mrhocon->Scale(10);

    TLegend* rleg = new TLegend(0.2,0.7,0.8,0.9);
    rleg->AddEntry(prhobkg,"Background peaks","L");
    rleg->AddEntry(prhocon,"Conversion peaks (X10)","L");

    TCanvas* rhocan = new TCanvas("rhocan","rhocan",800,600);
    rhocan->Divide(2,1);
    rhocan->cd(1);
    prhobkg->Draw();
    prhocon->Draw("same");
    rleg->Draw();
    rhocan->cd(2);
    mrhobkg->Draw();
    mrhocon->Draw("same");

  } else if(spage == "spread"){

    TH1F* srhocon = new TH1F("srhocon","Sigma of hit #rho distribution;#sigma #rho (mm)",100,0,200);
    TH1F* srhobkg = new TH1F("srhobkg","Sigma of hit #rho distribution;#sigma #rho (mm)",100,0,200);
    TH1F* sphicon = new TH1F("sphicon","Sigma of hit #rho*#phi distribution;#sigma #rho*#phi (mm)",100,0,200);
    TH1F* sphibkg = new TH1F("sphibkg","Sigma of hit #rho*#phi distribution;#sigma #rho*#phi (mm)",100,0,200);

    srhocon->SetLineColor(kRed);
    srhobkg->SetLineColor(kBlue);
    sphicon->SetLineColor(kRed);
    sphibkg->SetLineColor(kBlue);

    srhocon->SetStats(0);
    srhobkg->SetStats(0);
    sphicon->SetStats(0);
    sphibkg->SetStats(0);

    ddiag->Project("srhocon","sqrt(srho)",con);
    ddiag->Project("srhobkg","sqrt(srho)",bkg);
    ddiag->Project("sphicon","prho*sqrt(sphi)",con);
    ddiag->Project("sphibkg","prho*sqrt(sphi)",bkg);

    srhocon->Scale(10);
    sphicon->Scale(10);

    TLegend* sleg = new TLegend(0.2,0.7,0.8,0.9);
    sleg->AddEntry(srhobkg,"Background peaks","L");
    sleg->AddEntry(srhocon,"Conversion peaks (X10)","L");

    TCanvas* scan = new TCanvas("scan","scan",800,600);
    scan->Divide(2,1);
    scan->cd(1);
    srhobkg->Draw();
    srhocon->Draw("same");
    sleg->Draw();
    scan->cd(2);
    sphibkg->Draw();
    sphicon->Draw("same");
  
  } else if(spage == "z"){

    TH1F* zmincon = new TH1F("zmincon","peak zmin;zmin (mm)",100,-1600,1600);
    TH1F* zminbkg = new TH1F("zminbkg","peak zmin;zmin (mm)",100,-1600,1600);
    TH1F* zmaxcon = new TH1F("zmaxcon","peak zmax;zmax (mm)",100,-1600,1600);
    TH1F* zmaxbkg = new TH1F("zmaxbkg","peak zmax;zmax (mm)",100,-1600,1600);
    TH1F* zgapcon = new TH1F("zgapcon","peak zgap;zgap (mm)",100,0,1600);
    TH1F* zgapbkg = new TH1F("zgapbkg","peak zgap;zgap (mm)",100,0,1600);

    zmincon->SetLineColor(kRed);
    zminbkg->SetLineColor(kBlue);
    zmaxcon->SetLineColor(kRed);
    zmaxbkg->SetLineColor(kBlue);
    zgapcon->SetLineColor(kRed);
    zgapbkg->SetLineColor(kBlue);

    zmincon->SetStats(0);
    zminbkg->SetStats(0);
    zmaxcon->SetStats(0);
    zmaxbkg->SetStats(0);
    zgapcon->SetStats(0);
    zgapbkg->SetStats(0);

    ddiag->Project("zmincon","zmin",con);
    ddiag->Project("zminbkg","zmin",bkg);
    ddiag->Project("zmaxcon","zmax",con);
    ddiag->Project("zmaxbkg","zmax",bkg);
    ddiag->Project("zgapcon","zgap",con);
    ddiag->Project("zgapbkg","zgap",bkg);

    zmincon->Scale(10);
    zmaxcon->Scale(10);
    zgapcon->Scale(10);

    TLegend* zleg = new TLegend(0.2,0.7,0.8,0.9);
    zleg->AddEntry(zminbkg,"Background peaks","L");
    zleg->AddEntry(zmincon,"Conversion peaks (X10)","L");

    TCanvas* zcan = new TCanvas("zcan","zcan",1200,800);
    zcan->Divide(2,2);
    zcan->cd(1);
    zminbkg->Draw();
    zmincon->Draw("same");
    zleg->Draw();
    zcan->cd(2);
    zmaxbkg->Draw();
    zmaxcon->Draw("same");
    zcan->cd(3);
    zgapbkg->Draw();
    zgapcon->Draw("same");

  } else if(spage == "stations"){

    TH1F* smincon = new TH1F("smincon","peak smin;smin (mm)",36,-0.5,35.5);
    TH1F* sminbkg = new TH1F("sminbkg","peak smin;smin (mm)",36,-0.5,35.5);
    TH1F* smaxcon = new TH1F("smaxcon","peak smax;smax (mm)",36,-0.5,35.5);
    TH1F* smaxbkg = new TH1F("smaxbkg","peak smax;smax (mm)",36,-0.5,35.5);
    TH1F* nsmisscon = new TH1F("nsmisscon","peak nsmiss;nsmiss (mm)",36,-0.5,35.5);
    TH1F* nsmissbkg = new TH1F("nsmissbkg","peak nsmiss;nsmiss (mm)",36,-0.5,35.5);
    TH1F* nscon = new TH1F("nscon","peak ns;ns (mm)",36,-0.5,35.5);
    TH1F* nsbkg = new TH1F("nsbkg","peak ns;ns (mm)",36,-0.5,35.5);

    smincon->SetLineColor(kRed);
    sminbkg->SetLineColor(kBlue);
    smaxcon->SetLineColor(kRed);
    smaxbkg->SetLineColor(kBlue);
    nsmisscon->SetLineColor(kRed);
    nsmissbkg->SetLineColor(kBlue);
    nscon->SetLineColor(kRed);
    nsbkg->SetLineColor(kBlue);

    smincon->SetStats(0);
    sminbkg->SetStats(0);
    smaxcon->SetStats(0);
    smaxbkg->SetStats(0);
    nsmisscon->SetStats(0);
    nsmissbkg->SetStats(0);
    nscon->SetStats(0);
    nsbkg->SetStats(0);

    ddiag->Project("smincon","smin",con);
    ddiag->Project("sminbkg","smin",bkg);
    ddiag->Project("smaxcon","smax",con);
    ddiag->Project("smaxbkg","smax",bkg);
    ddiag->Project("nsmisscon","nsmiss",con);
    ddiag->Project("nsmissbkg","nsmiss",bkg);
    ddiag->Project("nscon","ns",con);
    ddiag->Project("nsbkg","ns",bkg);

    smincon->Scale(10);
    smaxcon->Scale(10);
    nsmisscon->Scale(10);
    nscon->Scale(10);

    TLegend* stleg = new TLegend(0.2,0.7,0.8,0.9);
    stleg->AddEntry(sminbkg,"Background peaks","L");
    stleg->AddEntry(smincon,"Conversion peaks (X10)","L");

    TCanvas* stcan = new TCanvas("stcan","stcan",1200,800);
    stcan->Divide(2,2);
    stcan->cd(1);
    sminbkg->Draw();
    smincon->Draw("same");
    stleg->Draw();
    stcan->cd(2);
    smaxbkg->Draw();
    smaxcon->Draw("same");
    stcan->cd(3);
    nsmissbkg->Draw();
    nsmisscon->Draw("same");
    stcan->cd(4);
    nsbkg->Draw();
    nscon->Draw("same");
  }

}
