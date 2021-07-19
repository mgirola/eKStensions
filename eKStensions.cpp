/* Analisi Statistica dei Dati A.A 2020/2021
 * Authors: Massimo Girola, Roberto Moretti, Emanuele Mazzola 
 * Date: 15 July 2021
 * compile with the following line:
   g++ -o eKStensions eKStensions.cpp `root-config --cflags --glibs`
*/

//STL libraries
#include <iostream>
#include <chrono>   //execution time
#include <cmath>    //pow
#include <stdexcept> //runtime error
#include <vector>
#include <utility> //std::make_pair
#include <string>
#include <iomanip>
#include <ctime> //std::time
#include <algorithm>

//ROOT libraries
#include <TRandom3.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TColor.h>
#include <TString.h>
#include <TApplication.h>
#include <TAxis.h>
#include <TH1.h>
#include <TF1.h>
#include <TF2.h>
#include <TH2.h>
#include <TSystem.h>
#include <Math/GoFTest.h>

#define PRINTCONTAINER(x) std::cout<<"{ "; for(const auto &el : x ) std::cout<<el<<", "; std::cout<<" }"<<std::endl

//STL usage
using std::cout, std::endl, std::to_string;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds, std::chrono::seconds, std::chrono::minutes;
using std::vector;
using std::string;
using std::setw, std::setprecision;
using std::srand, std::rand;
using ROOT::Math::GoFTest;

vector<double> GenData1D( TF1* pdf, int N ){
	vector<double> data(N);
	for( double &el: data )
		el = pdf->GetRandom();
	return data;
}

TF1* GenerateEmpiricalCDF_1D( vector<double> data, double x0, double x1, TString name = "Empirical CDF" ){
	auto IndicatorFunc = [](double x, double X){ return (X<=x) ? 1 : 0; }; //Indicator function, is 1 when X<=x
	auto EmpFunc = [IndicatorFunc, data](double* x, double* p){ //definition of empirical distribution function (cumulative)
		double partial_sum = 0;
		for( const auto & val : data )
			partial_sum += IndicatorFunc( x[0], val );
		return (1./data.size())*partial_sum;
	};
	
	return new TF1( name, EmpFunc,  x0, x1, 0 );
}

void PlotCDFvsEmpCDF_1D(vector<double> data, TF1* myPdf, TString name = "CDF vs Data"){
	auto myC = new TCanvas(name, name, 0, 0, 500, 500);
	myC->cd();
	myPdf->DrawIntegral();
	double xMin, xMax;
	myPdf->GetRange(xMin,xMax);
	auto EmpCDF = GenerateEmpiricalCDF_1D(data, xMin, xMax, name + " Empirical CDF");
	EmpCDF->SetNpx(1000);
	EmpCDF->Draw("same");
	myC->BuildLegend();
}

//simple example where fitted function chi2 accepts hypotesis while KS rejects it
//uses GoFTest in the ROOT toolkit to evaluate 1D KS
void KSvsChi21D(){
	
	double N = 1, mu = 0, sigma = 1;
	double bkg_fraction = 0.0005;

	TF1* myPdf = new TF1("myPdf", "gausn(0) + [3]*(1 + cos(2*x))", -10, 10);
	myPdf->SetParameters(N,mu,sigma, bkg_fraction);
	myPdf->SetNpx(10000);

	auto data1D = GenData1D( myPdf, 1e4 );
	
	auto h_KSvsChi2 = new TH1D("h_KSvsChi2","h_KSvsChi2",100,-10,10);
	for( const auto &val : data1D )
		h_KSvsChi2->Fill(val);
	h_KSvsChi2->Sumw2();
	h_KSvsChi2->Scale(1./h_KSvsChi2->GetEntries(),"width");

	int ndf = h_KSvsChi2->GetNbinsX();
	double chi2 = h_KSvsChi2->Chisquare(myPdf);
	double p = TMath::Prob(chi2, ndf);
	double chi2BakerCousins = h_KSvsChi2->Chisquare(myPdf, "L");
	double p2 = TMath::Prob(chi2BakerCousins, ndf);
	cout<<"chi2 results with original function used to generate data: no fit was performed"<<endl;
	cout<<"\"normal chi2\":\n\t p = "<<p<<", chi2 = "<<chi2<<", ndf = "<<ndf<<endl;
	cout<<" Baker Cousins chi2:\n\t p = "<<p2<<", chi2 = "<<chi2BakerCousins<<", ndf = "<<ndf<<endl<<endl;
	
	TF1* myGausF = new TF1("myGausn","gausn",-10,10);
	myGausF->SetParameters(N, mu, sigma);
	myGausF->SetNpx(10000);
	chi2 = h_KSvsChi2->Chisquare(myGausF);
	p = TMath::Prob(chi2, ndf);
	chi2BakerCousins = h_KSvsChi2->Chisquare(myGausF, "L");
	p2 = TMath::Prob(chi2BakerCousins, ndf);
	cout<<"chi2 results with pure gaussian function used to generate data: no fit was performed"<<endl;
	cout<<"\"normal chi2\":\n\t p = "<<p<<", chi2 = "<<chi2<<", ndf = "<<ndf<<endl;
	cout<<" Baker Cousins chi2:\n\t p = "<<p2<<", chi2 = "<<chi2BakerCousins<<", ndf = "<<ndf<<endl<<endl;
	
	
	auto test = GoFTest(data1D.size(), &data1D[0], GoFTest::EDistribution::kGaussian);
	double t, D;
	test.KolmogorovSmirnovTest(t, D);
	cout<<"KS for gaussian function:"<<endl;
	cout<<"t = "<<t<<", D = "<<D<<endl;

	auto test2 = GoFTest(data1D.size(), &data1D[0], *myPdf, GoFTest::EUserDistribution::kPDF, -10, 10);
	double t2, D2;
	test2.KolmogorovSmirnovTest(t2,D2);
	cout<<"KS for originial function"<<endl;
	cout<<"t2 = "<<t2<<", D2 = "<<D2<<endl;

	auto myC = new TCanvas( "KSvsChi2", "KS vs Chi2", 500, 500 );
	myC->cd();
	h_KSvsChi2->Draw();
	myGausF->SetLineColor(kRed);
	myGausF->Draw("l same");
	myPdf->SetLineColor(kBlue);
	myPdf->Draw("same"); 
	myC->BuildLegend();

	PlotCDFvsEmpCDF_1D(data1D, myPdf);	
}

void eKStensions(){}

int main(){
	gStyle->SetOptFit(111111);
	auto myApp = new TApplication("myApp", NULL, NULL);
	
	KSvsChi21D();
	
	eKStensions();
	myApp->Run();
	return 0;
}
