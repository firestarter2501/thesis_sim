{
gROOT ->Reset();
gStyle ->SetOptStat(1001110);
double x,y;
TCanvas*c1=new TCanvas("sim");
TH1S*hist=new TH1S("sim","sim",756,0,756);
ifstream ifs("./scinti_0.dat");
while(ifs>>x)
hist->Fill(x);
// c1->GetXaxis()->SetRangeUser(1,110000);;
hist->Draw();
gPad->SetLogy(1);
}
