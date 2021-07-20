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
#include <numeric> //std::accumulate

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
using std::chrono::duration;
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
	TF1* myPdf1 = (TF1*)myPdf->Clone("integral "+name);
	myPdf1->SetLineColor(kBlue);
	myC->SetGridx(); myC->SetGridy();
	myPdf1->DrawIntegral();
	double xMin, xMax;
	myPdf1->GetRange(xMin,xMax);
	auto EmpCDF = GenerateEmpiricalCDF_1D(data, xMin, xMax, name+" Empirical CDF");
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
	myPdfTest->SetRange(xMin, xMax); //set the same range in the pdfs
	vector<double> minimum_Dists(N);
	TString name = "Minimum KS Distances";
	cout<<"Repeating "<<N<<" KS tests, generating "<< nDataGen <<" data points each time....."<<endl;
	if(N*nDataGen>1000*1e04) cout<<"Warning: this might take a long time............"<<endl;
	auto t1 = high_resolution_clock::now();
	//generate data and repeat KS test N times
	for( auto & D : minimum_Dists ){	
		//---------generate data according to myPdfGen------------
		auto data1D = GenData1D( myPdfGen, nDataGen );	
		
		//---------KS test for original myPdfTest------------
		GoFTest test = GoFTest(data1D.size(), &data1D[0], *myPdfTest, GoFTest::EUserDistribution::kPDF, xMin, xMax);
		double t;
		test.KolmogorovSmirnovTest(t,D); //evaluate minimum distance (KS GoF test)
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<seconds>( t2 - t1 ).count();
	cout<<"statistics for KS minimum distances has been generated in "<<duration<<" seconds (";
	cout<<1000.*(double)duration/N<<" ms for one iteration)"<<endl<<endl;
	double Min = *std::min_element(minimum_Dists.begin(), minimum_Dists.end());	
	double Max = *std::max_element(minimum_Dists.begin(), minimum_Dists.end());	
	auto hMinDist = new TH1D(name, name, 100, Min, Max); //prepare histogram to store the minimum distances
	for( const auto & D : minimum_Dists )
		hMinDist->Fill(D); //fill histogram
	auto myC = new TCanvas(name, name, 500, 0, 1000, 500);
	myC->cd();
	hMinDist->Scale(1./hMinDist->Integral());

	//evaluate the empirical CDF of the minimum distances found in the MC simulation
	//this will be used along with the foundamental theorem of calculus 
	//where it will be interpreted as the primitive of the KS minimum distances distribution
	//and it will be used to extract the t-value by evaluating the difference between [F(x_max)-F(D)]=(integral of the distrib from D to x_max)
	//where D is the minimum distance returend by KS test and x_max should be the upper limit range such that (F(x_max)=1).
	auto minDistCDF = GenerateEmpiricalCDF_1D(minimum_Dists, Min, Max, "KS Distances empirical CDF");
	minDistCDF->SetNpx(10000);
	minDistCDF->Draw();	
	hMinDist->Draw("same");

	return {minDistCDF, hMinDist};
}

//the following function implements a case where chi2 behaves differently from KS.
//It is a simple example where fitted function chi2 accepts hypotesis while KS rejects it
//It uses GoFTest in the ROOT toolkit to evaluate 1D KS
void Example_KSvsChi2_1D(){
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
	int MC_generation = 1000;
	auto MCresults = MCStatistcalKS_1D( myPdf, myPdf, nDataGen, MC_generation );
	double X0, X1;
	MCresults.first->GetRange(X0, X1);
	double t2_MC = MCresults.first->Eval(X1)-MCresults.first->Eval(D2);
	double err_t2MChist, t2_MC_hist = MCresults.second->IntegralAndError(MCresults.second->FindBin(D2), MCresults.second->FindBin(X1), err_t2MChist);
	double diff_t2 = t2-t2_MC;
	double diff_t2_hist = t2-t2_MC_hist;
	cout<<"MC results: "<<endl;
	cout<<"\t t2_MC = "<<t2_MC<<endl;
	cout<<"\t Difference between t2 from MC and from GoFTest theroretical values = "<<diff_t2<<endl;
	cout<<"\t t2_MC hist = "<<t2_MC_hist<<" +- "<<err_t2MChist<<endl;
	cout<<"\t Difference between t2 from MC and from GoFTest theroretical values = "<<diff_t2_hist<<endl;
	cout<<"MC statistics completed"<<endl<<endl;
}

//the following function implement an extension of the KS test
//for when some of the parameters of the distribution are 
//estimated from the data themselves.
//According to various references this test makes sense only if
//the estimated parameters are parameters of scale or location.
//When the data follows a normal distribution but
//the mean and the sigma are estimated from the data
//we are in the special case known in literature
//as the Lilliefors Test, that we can use to check for consistency.
pair<int, vector<double>> KS_test_extended_to_fitted_1D( vector<double> data, TF1* myFittedPdfTest, int nFittedPars, double nSamples_forMC, double x0, double x1 ){
	cout<<"performing KS test extendended to fitted 1D data..."<<endl;
		
	int ndf = data.size() - nFittedPars; //extension of degrees of freedom to the KS test in analogy with chi2

	//1. build distribution of KS minimum distances using MC method with data generated from the estimated (fitted) pdf
	//FIXME auto MCresults = MCStatistcalKS_1D( myFittedPdfTest, (TF1*)myFittedPdfTest->Clone(), data.size(), nSamples_forMC );
		
	
	//2. evaluate the KS distance on the original data	
	auto test = GoFTest(data.size(), &data[0], *myFittedPdfTest, GoFTest::EUserDistribution::kPDF, x0, x1);
	double t2, D;
	test.KolmogorovSmirnovTest(t2,D);
	cout<<"D = "<<D<<endl;
	
	//3. evaluate the t-value from the MC built distribution in two different ways for the fitted PDF
	double X0, X1;
	MCresults.first->GetRange(X0, X1);
	double t_MC = MCresults.first->Eval(X1)-MCresults.first->Eval(D); //uses empirical CDF
	double err_tMChist, t_MC_hist = MCresults.second->IntegralAndError(MCresults.second->FindBin(D), MCresults.second->FindBin(X1), err_tMChist); //groups data in histogram
	double diff_t = t2-t_MC;
	double diff_t_hist = t2-t_MC_hist;
	cout<<"GoFTest t results: "<<endl;
	cout<<"\t t = "<<t2<<endl;
	cout<<"MC results: "<<endl;
	cout<<"\t t2_MC = "<<t_MC<<endl;
	cout<<"\t t2_MC hist = "<<t_MC_hist<<" +- "<<err_tMChist<<endl;
	cout<<"Difference between the two:"<<endl;
	cout<<"\t Difference between t2 from MC (est CDF) and from GoFTest theroretical values = "<<diff_t<<endl;
	cout<<"\t Difference between t2 from MC (hist   ) and from GoFTest theroretical values = "<<diff_t_hist<<endl;

	//4. build a table specific for this pdf and for this ndf
	//   t-values: .20, .15, .10, .05, .01
	vector<double> table_row = {.20, .15, .10, .05, .01};
	cout<<"if the value of D exceeds the critical values in the table one rejects the hypotesis that the data are from the estimated pdf"<<endl;
	cout<<"in other words, any value of D larger than the ones indicated in the table is significant at the indicated level of significance"<<endl;
	cout<<"in other words, for a fixed value of D, bigger t-values are better in the sense that they indicate better agreement with the data"<<endl;
	cout<<"ndf = "<<ndf<<" "; PRINTCONTAINER(table_row);
	for( auto & val : table_row ){
		//have to find the D corresponding to val and substiute val with it
		auto myF = MCresults.first;
		auto sup = myF->Eval(X1);
		auto y = 1 - val;
		val = myF->GetX(y); //let do the inversion to ROOT and hope it works
	}
	cout<<"ndf = "<<ndf<<" "; PRINTCONTAINER(table_row);
	cout<<endl<<"completed KS test extendended to fitted 1D data..."<<endl;
	return {ndf, table_row};
}

void Lilliefors_Consistency_Test(){
	
	int N = 10000; //how many times repeat the test in MC simulation

	//---------generate gausn data------------
	int SampleSize = 10;
	double norm = 1, mu = 2, sigma = 1;
	double X0 = mu - 5*sigma, X1 = mu + 5 *sigma;
	TF1* myPdf = new TF1("myPdf", "gausn", X0, X1);
	myPdf->SetParameters(norm,mu,sigma);
	myPdf->SetNpx(10000);
	auto data1D = GenData1D( myPdf, SampleSize );
	PlotCDFvsEmpCDF_1D(data1D, myPdf, "lilliefors before fit");
	
	//-------estimate parameters from data-------
	double sum = std::accumulate(std::begin(data1D), std::end(data1D), 0.0);
	double m =  sum / data1D.size();
	double accum = 0.0;
	std::for_each (std::begin(data1D), std::end(data1D), [&](const double d) {
		accum += (d - m) * (d - m);
	});
	double stdev = sqrt(accum / (data1D.size()-1));
	myPdf->SetParameters(norm, m, stdev);
	cout<<"estimated mean = "<<m<<endl<<"estimated stddev = "<<stdev<<endl;

	//----perform MC KS test----
	auto results = KS_test_extended_to_fitted_1D(data1D, myPdf, 2, N, X0, X1);

	//PlotCDFvsEmpCDF_1D(data1D, myPdf, "lilliefors after fit");
	

}



void eKStensions(){}

int main(){
	gStyle->SetOptFit(111111);
	auto myApp = new TApplication("myApp", NULL, NULL);
	
	//Example_KSvsChi2_1D();
	Lilliefors_Consistency_Test();
	
	eKStensions();
	myApp->Run();
	return 0;
}
