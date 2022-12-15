{
gROOT ->Reset();
gStyle ->SetOptStat(1001110);
double x1, x2;
TCanvas*c1=new TCanvas("sim","",640,512);
// TH1F*h1=new TH1F("scinti1","",128,0,768);
// TH1F*h2=new TH1F("scinti2","",128,0,768);
// TH1F*h3=new TH1F("sum","",128,0,768);
TH2S*hist2d=new TH2S("hist","sim-150deg",128,0,768,128,0,768);
ifstream ifs("./scinti_150deg.dat");
while(ifs>>x1>>x2)
{
    // if (661-50 < x1+x2 && x1+x2 < 661+50 /*&& 80 < x1 && 80 < x2*/)
    // {
        // h1->Fill(x1);
        // h2->Fill(x2);
        // h3->Fill(x1+x2);
        hist2d->Fill(x1, x2);
    // }
}
// h1->SetLineColor(2);
// h2->SetLineColor(4);
// h3->SetLineColor(6);
// h1->Draw();
// h2->Draw("same");
// h3->Draw("same");
hist2d->Draw("colz");
// gPad->SetLogx(1);
// gPad->SetLogy(1);
// h1->SetStats(0);
// h2->SetStats(0);
// h3->SetStats(0);
hist2d->SetStats(0);
hist2d->GetXaxis()->SetTitle("Scintillator 1 energy [keV]");
hist2d->GetXaxis()->CenterTitle();
hist2d->GetYaxis()->SetTitle("Scintillator 2 energy [keV]");
hist2d->GetYaxis()->CenterTitle();
// TLegend *legend = new TLegend( 0.8, 0.68, 0.99, 0.78) ;
// legend->AddEntry( h1, "scinti1" , "l") ;
// legend->AddEntry( h2, "scinti2" , "l") ;
// legend->AddEntry( h3, "sum" , "l" ) ;
// legend->SetFillColor(0);
// legend->Draw();
}
