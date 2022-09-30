{
gROOT ->Reset();
gStyle ->SetOptStat(1001110);
double x1, x2, y;
TCanvas*c1=new TCanvas("sim");
TH1F*h1=new TH1F("first","",768,0,768);
TH1F*h2=new TH1F("all","",768,0,768);
ifstream ifs1("./scinti_first.dat");
while(ifs1>>x1)
h1->Fill(x1);
h1->SetLineColor(2);
h1->Draw();
ifstream ifs2("./scinti_all.dat");
while(ifs2>>x2)
h2->Fill(x2);
h2->SetLineColor(4);
h2->Draw("same");
gPad->SetLogy(1);
TLegend *legend = new TLegend( 0.8, 0.68, 0.99, 0.78) ; //（）の中は位置の指定（左下の x , y 、右上の x , y ）
legend->AddEntry( h1, "first" , "l") ; // AddEntry( pointer , "interpretation" , "option" )
legend->AddEntry( h2, "all" , "l") ; // option は　"f"=box, "l"="L"=line, "p"=marker
// legend->SetFillColor(0);
legend->Draw() ;
}
