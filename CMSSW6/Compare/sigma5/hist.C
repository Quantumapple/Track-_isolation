void hist()
{
    TFile *fInput1 = new TFile("plot_signal.root","OPEN");
    TFile *fInput2 = new TFile("plot_bkg.root","OPEN");

    TH1F *signal = (TH1F*)fInput1->Get("h1");
    TH1F *minbias = (TH1F*)fInput2->Get("h1");
    
    signal->SetName("signal");
    minbias->SetName("minbias");

    signal->Scale(1./signal->Integral(1,100));
    minbias->Scale(1./minbias->Integral(1,100));
   
    cout << "==============================================" << endl << endl;
    cout << "Considering all tracks" << endl << endl;

    cout << "When under 0.08" << endl;
    cout << "Signal: " << signal->Integral(1,8) << ", background: " << 1.-minbias->Integral(1,8) << endl;
    cout << endl;
    
    cout << "When under 0.1" << endl;
    cout << "Signal: " << signal->Integral(1,10) << ", background: " << 1.-minbias->Integral(1,10) << endl;
    cout << endl;
    
    cout << "When under 0.15" << endl;
    cout << "Signal: " << signal->Integral(1,15) << ", background: " << 1.-minbias->Integral(1,15) << endl;
    cout << endl;
     
    cout << "When under 0.2" << endl;
    cout << "Signal: " << signal->Integral(1,20) << ", background: " << 1.-minbias->Integral(1,20) << endl;
    cout << endl;
    cout << "==============================================" << endl;
    cout << endl;

    signal->SetFillStyle(3004);
    signal->SetFillColor(kBlue);
    signal->SetLineColor(kBlue);
    signal->SetLineWidth(2);
    signal->SetXTitle("Isolation");
    signal->GetXaxis()->CenterTitle(true);
    signal->GetXaxis()->SetTitleSize(0.05);
    signal->SetYTitle("Normalized number of L1 e/#gamma");
    signal->GetYaxis()->CenterTitle(true);
    signal->GetYaxis()->SetTitleSize(0.05);
    
    minbias->SetFillStyle(3005);
    minbias->SetFillColor(kRed);
    minbias->SetLineColor(kRed);
    minbias->SetLineWidth(2);

    TCanvas *cc = new TCanvas("cc","",1366,768);
    cc->SetLeftMargin(0.12);
    cc->SetBottomMargin(0.12);
    cc->SetLogy();

    signal->Draw("HIST");
    minbias->Draw("HIST sames");

    TLegend *lgd1 = new TLegend(0.485,0.68,0.705,0.85);
    lgd1->SetTextSize(0.045);
    lgd1->AddEntry(signal,"signal","lp");
    lgd1->AddEntry(minbias,"minbias","lp");
    lgd1->Draw();

    TH2F *signal2 = (TH2F*)fInput1->Get("h2");
    TH2F *minbias2 = (TH2F*)fInput2->Get("h2");

    signal2->SetXTitle("Number of tracks");
    signal2->GetXaxis()->CenterTitle(true);
    signal2->GetXaxis()->SetTitleSize(0.05);
    signal2->SetYTitle("Isolation");
    signal2->GetYaxis()->CenterTitle(true);
    signal2->GetYaxis()->SetTitleSize(0.05);

    TCanvas *c1 = new TCanvas("c1","",800,800);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    signal2->Draw("colz");

    TLegend *l1 = new TLegend(0.5,0.5,0.7,0.7);
    l1->SetTextSize(0.05);
    l1->AddEntry(signal2,"signal");
    l1->Draw();

    minbias2->SetXTitle("Number of tracks");
    minbias2->GetXaxis()->CenterTitle(true);
    minbias2->GetXaxis()->SetTitleSize(0.05);
    minbias2->SetYTitle("Isolation");
    minbias2->GetYaxis()->CenterTitle(true);
    minbias2->GetYaxis()->SetTitleSize(0.05);

    TCanvas *c2 = new TCanvas("c2","",800,800);
    c2->SetLeftMargin(0.12);
    c2->SetBottomMargin(0.12);
    minbias2->Draw("colz");
    
    TLegend *l2 = new TLegend(0.5,0.5,0.75,0.7);
    l2->SetTextSize(0.05);
    l2->AddEntry(minbias2,"minbias");
    l2->Draw();
    
    TH1F *signal3 = (TH1F*)fInput1->Get("a1");
    TH1F *minbias3 = (TH1F*)fInput2->Get("a1");
    
    signal3->SetName("signal3");
    minbias3->SetName("minbias3");

    signal3->Scale(1./signal3->Integral(1,100));
    minbias3->Scale(1./minbias3->Integral(1,100));
    
    cout << "Considering when 2 or more tracks" << endl << endl;
    
    cout << "When under 0.08" << endl;
    cout << "Signal: " << signal3->Integral(1,8) << ", background: " << 1.-minbias3->Integral(1,8) << endl;
    cout << endl;
    
    cout << "When under 0.1" << endl;
    cout << "Signal: " << signal3->Integral(1,10) << ", background: " << 1.-minbias3->Integral(1,10) << endl;
    cout << endl;
    
    cout << "When under 0.15" << endl;
    cout << "Signal: " << signal3->Integral(1,15) << ", background: " << 1.-minbias3->Integral(1,15) << endl;
    cout << endl;
     
    cout << "When under 0.2" << endl;
    cout << "Signal: " << signal3->Integral(1,20) << ", background: " << 1.-minbias3->Integral(1,20) << endl;
    cout << endl;

    signal3->SetFillStyle(3004);
    signal3->SetFillColor(kBlue);
    signal3->SetLineColor(kBlue);
    signal3->SetLineWidth(2);
    signal3->SetXTitle("Isolation");
    signal3->GetXaxis()->CenterTitle(true);
    signal3->GetXaxis()->SetTitleSize(0.05);
    
    minbias3->SetFillStyle(3005);
    minbias3->SetFillColor(kRed);
    minbias3->SetLineColor(kRed);
    minbias3->SetLineWidth(2);

    TCanvas *ccc = new TCanvas("ccc","",1366,768);
    ccc->SetBottomMargin(0.12);
    ccc->SetLogy();

    signal3->Draw("HIST");
    minbias3->Draw("HIST sames");
}
