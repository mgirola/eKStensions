/* Analisi Statistica dei Dati A.A 2020/2021
 * Author: Massimo Girola 
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
#include <utility>

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
using std::pair;
using ROOT::Math::GoFTest;

//simple function to extract data from PDF
//uses TRandom3 mersenne generator
vector<double> GenData1D( TF1* pdf, int N ){
	auto ran_gen = new TRandom3();
	ran_gen->SetSeed();
	vector<double> data(N);
	for( double &el: data )
		el = pdf->GetRandom(ran_gen);
	return data;
}

//returns pointer to the empirical CDF (continuous steps-like function which increases by 1/n at each real value in vector data)
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

//just a funtion to plot the CDF and the empirical CDF on the same canvas
void PlotCDFvsEmpCDF_1D(vector<double> data, TF1* myPdf, TString name = "CDF vs Data"){
	auto myC = new TCanvas(name, name, 0, 0, 500, 500);
	myC->cd();
	myPdf->SetLineColor(kBlue);
	myC->SetGridx(); myC->SetGridy();
	myPdf->DrawIntegral();
	double xMin, xMax;
	myPdf->GetRange(xMin,xMax);
	auto EmpCDF = GenerateEmpiricalCDF_1D(data, xMin, xMax, "Empirical CDF");
	EmpCDF->SetNpx(1000);
	EmpCDF->SetLineColor(kRed);
	EmpCDF->Draw("same");
	myC->BuildLegend();
	myC->Update();
}

//1D implementation of MC simulation to generate
//the distribution of the minimum distances
pair<TF1*, TH1*> MCStatistcalKS_1D( TF1* myPdfGen, TF1* myPdfTest, int nDataGen, int N ){
	cout<<"Running MC simulation to build KS statistics..."<<endl;

	double xMin, xMax;
	myPdfGen->GetRange(xMin, xMax);
	myPdfTest->SetRange(xMin, xMax); //set the same range
	vector<double> minimum_Dists(N);
	TString name = "Minimum KS Distances";
	auto hMinDist = new TH1D(name, name, 100, -0.0001, 0.04);
	int i = 0;
	for( auto & D : minimum_Dists ){	
		//---------generate data according to myPdfGen------------
		auto data1D = GenData1D( myPdfGen, nDataGen );	
		
		//---------KS test for original myPdfTest------------
		GoFTest test = GoFTest(data1D.size(), &data1D[0], *myPdfTest, GoFTest::EUserDistribution::kPDF, xMin, xMax);
		double t;
		test.KolmogorovSmirnovTest(t,D);
		hMinDist->Fill(D);
		cout<<100.*(double)(i++)/(double)N<<"\tD = "<<D<<endl;
	}
	auto myC = new TCanvas(name, name, 500, 0, 1000, 500);
	myC->cd();
	hMinDist->Scale(1./N);
	hMinDist->Draw();

	//auto myC2 = new TCanvas(name+" emp CDF", name+" emp CDF", 0, 500, 500, 1000);
	//myC2->cd();
	auto minDistCDF = GenerateEmpiricalCDF_1D(minimum_Dists, -0.001, 0.04, "KS Distances empirical CDF");
	minDistCDF->SetNpx(10000);
	minDistCDF->Draw("same");
	
	return {minDistCDF, hMinDist};
}

//the following function implements a case where chi2 behaves differently from KS.
//It is a simple example where fitted function chi2 accepts hypotesis while KS rejects it
//It uses GoFTest in the ROOT toolkit to evaluate 1D KS
void KSvsChi2_1D(){
	cout<<"Running test KS versus chi2..."<<endl;

	//---------generate data------------
	int nDataGen = 1e4;
	double N = 1, mu = 0, sigma = 1;
	double bkg_fraction = 0.0005, const_bkg_component = 1;
	TF1* myPdf = new TF1("myPdf", "gausn(0) + [3]*((1 + [4]) + cos(2*x))", -10, 10);
	myPdf->SetParameters(N,mu,sigma, bkg_fraction, const_bkg_component);
	myPdf->SetNpx(10000);
	auto data1D = GenData1D( myPdf, nDataGen );	
	
	//---------KS test for gaus------------
	cout<<endl<<"------------------Kolmogorov-Smirnov Results--------------------"<<endl;
	auto test = GoFTest(data1D.size(), &data1D[0], GoFTest::EDistribution::kGaussian);
	double t, D;
	test.KolmogorovSmirnovTest(t, D);
	cout<<"KS for gaussian function:"<<endl;
	cout<<"\tt = "<<t<<", D = "<<D<<endl;

	//---------KS test for original function------------
	auto test2 = GoFTest(data1D.size(), &data1D[0], *myPdf, GoFTest::EUserDistribution::kPDF, -10, 10);
	double t2, D2;
	test2.KolmogorovSmirnovTest(t2,D2);
	cout<<"KS for originial function"<<endl;
	cout<<"\tt2 = "<<t2<<", D2 = "<<D2<<endl;
	
	//---------plot empirical CDF and CDF------------
	PlotCDFvsEmpCDF_1D(data1D, myPdf);	
	
	//---------group data in bins to perform chi2 calculation------------
	cout<<endl<<"------------------chi2 Results--------------------"<<endl;
	auto h_KSvsChi2 = new TH1D("h_KSvsChi2","h_KSvsChi2",100,-10,10);
	for( const auto &val : data1D )
		h_KSvsChi2->Fill(val);
	h_KSvsChi2->Sumw2();
	h_KSvsChi2->Scale(1./h_KSvsChi2->GetEntries(),"width");

	//---------chi2 test for original function------------
	int ndf = h_KSvsChi2->GetNbinsX(); //notice that no fit is performed so that ndf = number of bins
	double chi2 = h_KSvsChi2->Chisquare(myPdf);
	double p = TMath::Prob(chi2, ndf);
	double chi2BakerCousins = h_KSvsChi2->Chisquare(myPdf, "L");
	double p2 = TMath::Prob(chi2BakerCousins, ndf);
	cout<<"chi2 results with original function used to generate data: no fit was performed"<<endl;
	cout<<"\t\"normal\" chi2:\n\t\t p = "<<p<<", chi2 = "<<chi2<<", ndf = "<<ndf<<endl;
	cout<<"\tBaker Cousins chi2:\n\t\t p = "<<p2<<", chi2 = "<<chi2BakerCousins<<", ndf = "<<ndf<<endl;
	
	//---------chi2 test for gaus------------
	TF1* myGausF = new TF1("myGausn","gausn(0) + [3]*(1+[4])",-10,10);
	myGausF->SetParameters(N, mu, sigma, bkg_fraction, const_bkg_component);
	myGausF->SetNpx(10000);
	chi2 = h_KSvsChi2->Chisquare(myGausF);
	p = TMath::Prob(chi2, ndf);
	chi2BakerCousins = h_KSvsChi2->Chisquare(myGausF, "L");
	p2 = TMath::Prob(chi2BakerCousins, ndf);
	cout<<"chi2 results with pure gaussian function used to generate data: no fit was performed"<<endl;
	cout<<"\t\"normal\" chi2:\n\t\t p = "<<p<<", chi2 = "<<chi2<<", ndf = "<<ndf<<endl;
	cout<<"\tBaker Cousins chi2:\n\t\t p = "<<p2<<", chi2 = "<<chi2BakerCousins<<", ndf = "<<ndf<<endl<<endl;
	
	//---------plot generated and binned data vs the two functions------------
	auto myC = new TCanvas( "KSvsChi2", "KS vs Chi2", 500, 500 );
	myC->cd();
	myC->SetLogy();
	myC->SetGridx(); myC->SetGridy();
	h_KSvsChi2->Draw();
	myGausF->SetLineColor(kRed);
	myGausF->Draw("l same");
	myPdf->SetLineColor(kBlue);
	myPdf->Draw("same"); 
	myC->BuildLegend();

	cout<<endl<<"the test has been completed!"<<endl;

	//Perform MC simulation to build test statistics to compare with the theoretical t-value returned by GoFTest
	cout<<endl<<"building test statistics of minimum KS distances from MC simulation..."<<endl;
	int MC_generation;
	auto MCresults = MCStatistcalKS_1D( myPdf, myPdf, nDataGen, 1000 );
	double t2_MC = MCresults.first->Integral(D2, 10);
	double t2_MC_hist, err = MCresults.second->IntegralAndError(MCresults.second->GetBin(D2), MCresults.second->GetNbinsX(), err);
	double diff_t2 = t2-t2_MC;
	double diff_t2_hist = t2-t2_MC_hist;
	cout<<"MC results: "<<endl;
	cout<<"\t t2_MC = "<<t2_MC<<endl;
	cout<<"\t Difference between t2 from MC and from GoFTest theroretical values = "<<diff_t2<<endl;
	cout<<"\t t2_MC hist = "<<t2_MC_hist<<endl;
	cout<<"\t Difference between t2 from MC and from GoFTest theroretical values = "<<diff_t2_hist<<endl;
	cout<<"MC statistics completed"<<endl<<endl;

}



void eKStensions(){}

int main(){
	gStyle->SetOptFit(111111);
	auto myApp = new TApplication("myApp", NULL, NULL);
	
	KSvsChi2_1D();
	
	eKStensions();
	myApp->Run();
	return 0;
}
