{
gROOT ->Reset();
gStyle ->SetOptStat(1001110);
double x1, x2, sum, pmt1slope = (1)/1, pmt1intersec = -(0)/1, pmt2slope = (1)/1, pmt2intersec = -(0)/1;
// double x1, x2, sum, pmt1slope = 1, pmt1intersec = 0, pmt2slope = 1, pmt2intersec = 0;
TCanvas*c1=new TCanvas("data","",640,512);
// TH1F*h1=new TH1F("scinti1","",128,0,768);
// TH1F*h2=new TH1F("scinti2","",128,0,768);
// TH1F*h3=new TH1F("sum","",128,0,768);
// TH1F*h1=new TH1F("scinti1","",128,0,2048);
// TH1F*h2=new TH1F("scinti2","",128,0,2048);
// TH1F*h3=new TH1F("sum","",128,0,2048);
// TH1F*h1=new TH1F("scinti1","",128,0,300);
// TH1F*h2=new TH1F("scinti2","",128,0,300);
// TH1F*h3=new TH1F("sum","",128,0,300);
TH2S*hist2d=new TH2S("hist","sim",768,0,768,768,0,768);
double deg1 = 2*TMath::Pi()/6, la = 11.17, lb = 20, lc = sqrt((la*la)+(lb*lb)-(2*la*lb*cos(deg1)));
double deg2 = TMath::Pi()-acos(((lb*lb)+(lc*lc)-(la*la))/(2*lb*lc)), cutrange = 50;
double ene1 = 661.6-(661.6/(1+(661.6/510.9)*(1-cos(deg1)))), ene2 = (661.6/(1+(661.6/510.9)*(1-cos(deg1))));
double ene3 = (661.6/(1+(661.6/510.9)*(1-cos(deg2)))), ene4 = 661.6-(661.6/(1+(661.6/510.9)*(1-cos(deg2))));
ifstream ifs("./scinti_120deg.dat");
while(ifs>>x1>>x2)
{
    if (661.6-cutrange < (x1*pmt1slope+pmt1intersec)+(x2*pmt2slope+pmt2intersec) && (x1*pmt1slope+pmt1intersec)+(x2*pmt2slope+pmt2intersec) < 661.6+cutrange /*&& 100 < x1 && 100 < x2*/)
    // if (std::abs(184.3 - (x2*pmt2slope+pmt2intersec)) < 100)
    {
        // h1->Fill(x1*pmt1slope+pmt1intersec);
        // h2->Fill(x2*pmt2slope+pmt2intersec);
        // h3->Fill((x1*pmt1slope+pmt1intersec)+(x2*pmt2slope+pmt2intersec));
        hist2d->Fill(x1*pmt1slope+pmt1intersec, x2*pmt2slope+pmt2intersec);
    }
}
// h1->SetStats(0);
// h1->SetLineColor(2);
// h2->SetStats(0);
// h2->SetLineColor(4);
// h3->SetStats(0);
// h3->SetLineColor(6);
// h1->Draw();
// h2->Draw("same");
// h3->Draw("same");
hist2d->Draw("colz");
// gPad->SetLogx(1);
// gPad->SetLogy(1);
// h1->Fit("gaus", "", "", 1600, 2048);
// h2->Fit("gaus", "", "", 1400, 2048);
hist2d->SetStats(0);
hist2d->GetXaxis()->SetTitle("Scintillator 1 energy [keV]");
hist2d->GetXaxis()->CenterTitle();
hist2d->GetYaxis()->SetTitle("Scintillator 2 energy [keV]");
hist2d->GetYaxis()->CenterTitle();
double one2two = hist2d->Integral(ene1-cutrange, ene1+cutrange, ene2-cutrange, ene2+cutrange);
double two2one = hist2d->Integral(ene3-cutrange, ene3+cutrange, ene4-cutrange, ene4+cutrange);
std::cout << "deg1: " << deg1 << endl << "deg2: " << deg2 << endl << "ene1: " << ene1 << endl << "ene2: " << ene2 << endl << "ene3: " << ene3 << endl << "ene4: " << ene4 << endl;
std::cout << "1->2: " << one2two << std::endl;
std::cout << "2->1: " << two2one << std::endl;
std::cout << "ratio: " << two2one/one2two << std::endl;
// TLegend *legend = new TLegend( 0.8, 0.68, 0.99, 0.78) ;
// legend->AddEntry( h1, "scinti1" , "l") ;
// legend->AddEntry( h2, "scinti2" , "l") ;
// legend->AddEntry( h3, "sum", "l");
// legend->SetFillColor(0);
// legend->Draw();
}
