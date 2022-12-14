{
gROOT ->Reset();
gStyle ->SetOptStat(1001110);
double x1, x2;
TCanvas*c1=new TCanvas("sim");
TH1F*h1=new TH1F("scinti1","",128,0,768);
TH1F*h2=new TH1F("scinti2","",128,0,768);
// TH1F*h3=new TH1F("sum","",128,0,768);
ifstream ifs1("./scinti_all.dat");
ifstream ifs2("./scinti_first.dat");
while(ifs1>>x1)
{
        h1->Fill(x1);
}
while(ifs2>>x2)
{
        h2->Fill(x2);
}
h1->SetLineColor(2);
h2->SetLineColor(4);
// h3->SetLineColor(6);
h1->Draw();
h2->Draw("same");
// h3->Draw("same");
gPad->SetLogy(1);
h1->SetStats(0);
h2->SetStats(0);
// h3->SetStats(0);
TLegend *legend = new TLegend( 0.8, 0.68, 0.99, 0.78) ;
legend->AddEntry( h1, "all" , "l") ;
legend->AddEntry( h2, "first" , "l") ;
// legend->AddEntry( h3, "sum" , "l" ) ;
legend->SetFillColor(0);
legend->Draw();
}
