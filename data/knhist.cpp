{
gROOT ->Reset();
gStyle ->SetOptStat(1001110);
double x1, x2;
TCanvas*c1=new TCanvas("sim");
TH1F*h1=new TH1F("scinti1","",768,0,768);
TH1F*h2=new TH1F("scinti2","",768,0,768);
ifstream ifs("./scinti.dat");
while(ifs>>x1>>x2)
{
h1->Fill(x1);
h2->Fill(x2);
}
h1->SetLineColor(2);
h2->SetLineColor(4);
h1->Draw();
h2->Draw("same");
gPad->SetLogy(1);
TLegend *legend = new TLegend( 0.8, 0.68, 0.99, 0.78) ;
legend->AddEntry( h1, "scinti1" , "l") ;
legend->AddEntry( h2, "scinti2" , "l") ;
// legend->SetFillColor(0);
legend->Draw();
}
