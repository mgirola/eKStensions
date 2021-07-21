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
#include <tuple>

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
// For Fitting
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/UnBinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/FitResult.h"
#include "Fit/DataOptions.h"
#include "Fit/FitConfig.h"


// For defining the functions
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"

// For plotting
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"


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
using std::tuple;
using ROOT::Math::GoFTest;

//simple function to extract data from PDF
//uses TRandom3 mersenne generator
vector<double> GenData1D( TF1* pdf, int N ){
	auto ran_gen = new TRandom3();
	ran_gen->SetSeed();
	if(pdf->GetNpx()<5000) pdf->SetNpx(5000);
	vector<double> data(N);
	for( double &el: data )
		el = pdf->GetRandom(ran_gen);
	delete ran_gen;
	return data;
}

//simple funct to extract from 2d pdf
vector<pair<double,double>> GenData2D( TF2* pdf, int N ){
	auto ran_gen = new TRandom3();
	ran_gen->SetSeed();
	pdf->SetNpx(1000);
	pdf->SetNpy(1000);
	vector<pair<double,double>> data(N);
	double x,y;
	for( pair<double,double> &el : data ){
		pdf->GetRandom2(x,y,ran_gen);
		el = {x,y};
	}
	delete ran_gen;
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

double QKS( double lambda, int trunc ){
	double sum = 0;
	for( int i = 1; i<trunc; i++ )
		sum += pow(-1,i-1)*exp(-2.*pow((double)i*lambda,2));
	return 2*sum;
}


TF1 QKS_func( double lambda_coeff, int trunc, double xMin, double xMax ){
	return TF1("QKS",[trunc,lambda_coeff](double * x, double * p){
				return QKS(lambda_coeff*x[0],trunc);
			},xMin,xMax,0);
}

double correlationCoefficient(vector<pair<double,double>> data){ //FIXME: this does not work and give also values greater than 1 and smaller than -1
    int sum_X = 0, sum_Y = 0, sum_XY = 0;
    int squareSum_X = 0, squareSum_Y = 0;
    for (const auto & p : data){
	double X = p.first, Y = p.second;
        sum_X += X;
        sum_Y += Y;
        sum_XY += X * Y;
        squareSum_X += X * X;
        squareSum_Y += Y * Y;
    }
    // use formula for calculating correlation coefficient.
    int n = data.size();
    double corr = (double)(n * sum_XY - sum_X * sum_Y) / sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
    return corr;
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
//meant to be used for pdf indipendent from the data
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
void unbinnedSimpleFit(vector<double> x, TF1* myFunction, double* origin_pars, double xMin, double xMax){
	gROOT->Reset();
	ROOT::Fit::DataOptions opt;
	ROOT::Fit::DataRange range(xMin,xMax); 
	ROOT::Fit::UnBinData data(opt, range, x.size());
  	for( auto it:x )
    		data.Add(it);
	ROOT::Math::WrappedMultiTF1 fitFunction( *myFunction, myFunction->GetNdim() );
	ROOT::Fit::Fitter fitter;
	fitter.SetFunction( fitFunction, false );
	double* initialParams = origin_pars;
	fitter.Config().SetParamsSettings(myFunction->GetNpar(),initialParams);
	fitter.Config().SetUpdateAfterFit();

	fitter.LikelihoodFit(data); //unbinned likelihood fit

	//ROOT::Fit::FitResult r=fitter.Result();
	//r.Print(std::cout);

	return;
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
pair<int, vector<double>> KS_test_extended_to_fitted_1D( vector<double> data, TF1* myPdfGen, TF1* myFittedPdfTest, double* origin_pars, int nFittedPars, double nSamples_forMC, double x0, double x1 ){
	cout<<"performing KS test extendended to fitted 1D data..."<<endl;
		
	int ndf = data.size() - nFittedPars; //extension of degrees of freedom to the KS test in analogy with chi2

	//-----------------------------------------------------------------------------------------------------------------
	//1. build distribution of KS minimum distances using MC method with data following the (true) myPdfGen PDF
	//-----------------------------------------------------------------------------------------------------------------
	cout<<"Running MC simulation to build corrected KS statistics..."<<endl;
	vector<double> minimum_Dists(nSamples_forMC);
	TString name = "Minimum KS Distances";
	cout<<"Repeating "<<nSamples_forMC<<" KS tests, generating "<< data.size() <<" data points each time....."<<endl;
	if(nSamples_forMC*data.size()>10000*30) cout<<"Warning: this might take a long time............"<<endl;
	auto time1 = high_resolution_clock::now();
	//generate data and repeat KS test N times
	auto myPdf = (TF1*)myPdfGen->Clone("original pdf");
	for( auto & D : minimum_Dists ){	
		//---------generate data according to myPdf------------
		auto data1D = GenData1D( myPdfGen, data.size() );	
		
		//---------unbinned fit of the fit funct to the data---------
		unbinnedSimpleFit(data1D, myPdf, origin_pars, x0, x1);

		//---------KS test for original myPdfGen ------------
		GoFTest test = GoFTest(data1D.size(), &data1D[0], *myPdf, GoFTest::EUserDistribution::kPDF, x0, x1);
		double t;
		test.KolmogorovSmirnovTest(t,D); //evaluate minimum distance (KS GoF test)
	}
	delete myPdf;
	//print some infos
	auto time2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<seconds>( time2 - time1 ).count();
	cout<<"statistics for KS minimum distances has been generated in "<<duration<<" seconds (";
	cout<<1000.*(double)duration/nSamples_forMC*data.size()<<" ms for one iteration)"<<endl<<endl;
	//find minimum and maximum in the vector of minimum KS distances
	double Min = *std::min_element(minimum_Dists.begin(), minimum_Dists.end());	
	double Max = *std::max_element(minimum_Dists.begin(), minimum_Dists.end());	
	auto hMinDist = new TH1D(name, name, 100, Min, Max); //prepare histogram to store the minimum distances
	for( const auto & D : minimum_Dists )
		hMinDist->Fill(D); //fill histogram
	//prepare canvas and normalize histogram
	auto myC = new TCanvas(name, name, 500, 0, 1000, 500);
	myC->cd();
	hMinDist->Scale(1./hMinDist->Integral());
	//generate the empirical CDF from the vector of the minimum distances, this will be used to estimate the t-values by reversing it
	auto minDistCDF = GenerateEmpiricalCDF_1D(minimum_Dists, Min, Max, "KS Distances empirical CDF");
	minDistCDF->SetNpx(10000);
	minDistCDF->Draw();	
	hMinDist->Draw("same");
	//store results in a pair
	pair<TF1*, TH1*> MCresults = {minDistCDF, hMinDist};
		
	
	//-----------------------------------------------------------------------------------------------------------------
	//2. evaluate the KS distance on the original data fitted with the (possibly wrong) pdf	
	//-----------------------------------------------------------------------------------------------------------------
	auto test = GoFTest(data.size(), &data[0], *myFittedPdfTest, GoFTest::EUserDistribution::kPDF, x0, x1);
	double t2, D;
	test.KolmogorovSmirnovTest(t2,D);
	cout<<"D = "<<D<<endl;
	
	//-----------------------------------------------------------------------------------------------------------------
	//3. evaluate the t-value from the MC built distribution in two different ways for the original pdf
	//-----------------------------------------------------------------------------------------------------------------
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

	//-----------------------------------------------------------------------------------------------------------------
	//4. build a table specific for this pdf and for this ndf
	//-----------------------------------------------------------------------------------------------------------------
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
		auto y = sup - val;
		val = myF->GetX(y); //let do the inversion to ROOT and hope it works
	}
	cout<<"ndf = "<<ndf<<" "; PRINTCONTAINER(table_row);
	cout<<endl<<"completed KS test extendended to fitted 1D data..."<<endl;
	return {ndf, table_row};
}

//function meant to give results (in particular t-values) meant to 
//compare with the ones reported by original lilliefors article 
//to check the consistency of the current implemetation
void Lilliefors_Consistency_Test( int N = 1e04, int SampleSize = 10 ){

	//---------generate gausn data------------
	double mu = 5, sigma = 2;
	double X0 = mu - 5*sigma, X1 = mu + 5 *sigma;
	TF1* myPdf = new TF1("myPdf", "ROOT::Math::normal_pdf(x,[1],[0])", X0, X1);
	myPdf->SetParameters(mu,sigma);
	double * original_pars = myPdf->GetParameters();
	myPdf->SetNpx(10000);
	auto data1D = GenData1D( myPdf, SampleSize );
	//PlotCDFvsEmpCDF_1D(data1D, (TF1*)myPdf->Clone("lilliefors before fit func"), "lilliefors before fit");
	
	//-------estimate parameters from data-------
	TF1* myPdfFitted = (TF1*)(myPdf->Clone("myPdfFitted"));
	unbinnedSimpleFit(data1D, myPdfFitted, original_pars, X0, X1);

	//----perform MC KS test----
	auto results = KS_test_extended_to_fitted_1D(data1D, myPdf, myPdfFitted, original_pars, 2, N, X0, X1);
}

void Extended_Lilliefors_Test( int N = 1e04, int SampleSize = 10 ){
	
	//---------generate gausn data------------
	double mu = 5, sigma = 2;
	double X0 = mu - 5*sigma, X1 = mu + 5 *sigma;
	double bkg = 0.01;
	TF1* myPdf = new TF1("myPdf", "ROOT::Math::normal_pdf(x,[1],[0])", X0, X1);
	TF1* myPdfGen = new TF1("myPdfFit", "ROOT::Math::normal_pdf(x,[1],[0])+[2]*(1.+(1./3.)*sin(2*x))", X0, X1);
	myPdf->SetParameters(mu,sigma);
	myPdfGen->SetParameters(mu, sigma, bkg);
	myPdfGen->SetLineColor(kBlue);
	myPdfGen->Draw();
	myPdf->Draw("same");
	double * original_pars = myPdf->GetParameters();
	myPdf->SetNpx(10000);
	
	//generate data according to custom pdf
	auto data1D = GenData1D( myPdfGen, SampleSize );
	//PlotCDFvsEmpCDF_1D(data1D, (TF1*)myPdf->Clone("lil before fit func"), "lil before fit");
	
	//-------estimate parameters from data-------
	//fit data with normal pdf (which should not agree with the generated data)
	TF1* myPdfFitted = (TF1*)(myPdf->Clone("myPdfFitted")); //fit data with the wrong pdf
	unbinnedSimpleFit(data1D, myPdfFitted, original_pars, X0, X1);
	
	//----perform MC KS test with normal pdf----
	auto results = KS_test_extended_to_fitted_1D(data1D, myPdf, myPdfFitted, original_pars, 2, N, X0, X1);
}

//extension of the KS test to 2 dimensions based on ROOT method
//in this function there is no fit of the data
//we just generate data according to a 2D distribution
//and we use KS test to establish agreement
void KS_Test2D_ROOT(int SampleSize = 50){
	//prepare data
	double muX = -4, muY = +3;
	double sigmaX = 2, sigmaY = 2.8;
	double rho = 0.8; //correlation btw -1 and 1
	double x0, y0, x1, y1;
	x0 = std::max( abs(std::max(muX+5*sigmaX,muY+5*sigmaY)), abs(std::min(muX-5*sigmaX,muY-5*sigmaY)) );
	y0 = -x0; x1 = x0; y1 = x0; x0 = -x0;
	TF2* myPdfGen = new TF2("myPdfGen", "ROOT::Math::bigaussian_pdf(x,y,[0],[1],[2],[3],[4])", x0, x1, y0, y1);
	myPdfGen->SetParameters(sigmaX, sigmaY, rho, muX, muY);
	auto data = GenData2D(myPdfGen, SampleSize);
	TH2D* h1 = new TH2D("h1","h1",1000,x0,x1,1000,y0,y1);
	for( const auto & p : data )
		h1->Fill(p.first, p.second);
	
	vector<double> epss = {0, .0001, .001, .01, .05, .08, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .66, .7, .75, .8, .9, 1 };
	vector<double> ts = epss;
	for( auto & eps : ts ){
		TF2 myPdf = TF2("myPdfModified", "ROOT::Math::bigaussian_pdf(x,y,[0],[1],[2],[3],[4])", x0, x1, y0, y1);
		myPdf.SetParameters(sigmaX-eps, sigmaY-eps, rho, muX+eps, muY-eps);
		TH2D h2 = TH2D("h2","h2",1000,x0,x1,1000,y0,y1);
		h2.FillRandom("myPdfModified",1e07);
		
		double t = h1->KolmogorovTest(&h2,"D");
		eps = t;
	}
	auto myG = new TGraph(epss.size(), &epss[0], &ts[0]);
	myG->SetMarkerStyle(20);
	myG->Draw();
}

//function to estimate the maximum distance from TF2 and 2D data in a customized way
//for each point in the dataset we divide the domain in four quadrants and we 
//evaluate the fraction of points in each quadrant and the integral of the function in each quadrant
//then we find the distances between each fraction and we take the maximum distance
double KolmogorovSmirnovDistance2DCustom(TF2 myPdf, vector<pair<double,double>> data){
	double xmin, xmax, ymin, ymax;
	myPdf.GetRange( xmin, ymin, xmax, ymax );
	double x, y;
	double pdfNorm = myPdf.Integral(xmin,xmax,ymin,ymax);
	double nData = (double)data.size();
	double f1=0,f2=0,f3=0,f4=0; //fraction func in the four quadrants
	double d1=0,d2=0,d3=0,d4=0; //fraction points in the four quadrants
	double D1=0,D2=0,D3=0,D4=0;
	vector<double> Dmaxs;
	for( const auto & p : data ){
		x = p.first; y = p.second;
		if( x<xmin || x>xmax || y<ymin || y>ymax ) 
			throw std::runtime_error("values out of range");
		
		//frazione dati nei quattro quadranti
		for( const auto & P : data ){
			double X = P.first, Y = P.second;
			if( X==x, Y == y )
				continue;
			else if( X>x && Y>y )
				d1 += 1.;
			else if( X<x && Y>y )
				d2 += 1.;
			else if( X<x && Y<y )
				d3 += 1.;
			else if( X>x && Y<y )
				d4 += 1.;
		}
		d1 /= nData;
		d2 /= nData;
		d3 /= nData;
		d4 /= nData;
		
		//frazione funzione nei quattro quadranti
		f1 = myPdf.Integral(x,xmax,y,ymax)/pdfNorm;
		f2 = myPdf.Integral(xmin,x,y,ymax)/pdfNorm;
		f3 = myPdf.Integral(xmin,x,ymin,y)/pdfNorm;
		f4 = myPdf.Integral(x,xmax,ymin,y)/pdfNorm;
		
		//distanze
		D1 = abs(d1-f1);
		D2 = abs(d2-f2);
		D3 = abs(d3-f3);
		D4 = abs(d4-f4);

		//find maximum distance of the four 
		double max12 = D1 > D2 ? D1 : D2;
		double max34 = D3 > D4 ? D3 : D4;
		double maxD = max12 > max34 ? max12 : max34;

		Dmaxs.push_back(maxD);
	}

	return *std::max_element(Dmaxs.begin(), Dmaxs.end());	
}

//generate MC CDF distribution of KS statistics
//returns CDF meant to be used to extract t-value for a certain D
pair<TF1*,TH1*> MC_CDFfor2DCustom( TF2 myPdf, int SampleSize, int N ){
	cout<<"Running MC simulation for 2D with "<<N<<" iterations and "<<SampleSize<<" data in each sample...."<<endl; 
	vector<double> Ds(N);
	int ii = 0;
	for( auto & D : Ds ){
		cout<<ii++<<endl;
		auto data = GenData2D(&myPdf, SampleSize); //gen data with following myPdf
		D = KolmogorovSmirnovDistance2DCustom(myPdf,data); //find D under for data generated from myPdf
	}
	double minD = *std::min_element(Ds.begin(), Ds.end());
	double maxD = *std::max_element(Ds.begin(), Ds.end());
	auto KolmCDF2D = GenerateEmpiricalCDF_1D(Ds,minD,maxD,"empiricalCDF2D");
	KolmCDF2D->SetNpx(10000);
	auto KolmCDF2DCompl = new TF1("2DQKS_MC",[KolmCDF2D](double*x, double *p){return 1.-KolmCDF2D->Eval(x[0]);},minD,maxD,0);

	TString name = "HistoMinDist2D";
	auto hMinDist = new TH1D(name, name, 100, minD, maxD); //prepare histogram to store the minimum distances
	for( const auto & D : Ds )
		hMinDist->Fill(D); //fill histogram
	hMinDist->Scale(1./hMinDist->Integral());
	//store results in a pair
	cout<<"MC simulation 2D completed"<<endl;
	pair<TF1*, TH1*> MCresults = {KolmCDF2DCompl, hMinDist};
	return MCresults;
}



tuple<double,double,double> KolmogorovSmirnov2DCustom( TF2 myPdf, vector<pair<double,double>> data, int N_MC = 100, TString drawOpt = "draw" ){
	
	int SampleSize = data.size();

	//MC based evaluation of t-value
	static bool AlreadyExecuted = false, AlreadyDrawn = false;
	static pair<TF1*, TH1*> MCresults;
	TF1* myFunc;
	if ( not AlreadyExecuted ){
		MCresults = MC_CDFfor2DCustom(*((TF2*)(myPdf.Clone("myPdfTest"))),SampleSize,N_MC);
		AlreadyExecuted = true;
	}else
		cout<<"Warning, already executed, using old statistics, ensure you have not changed myPdf since the first run"<<endl;
	
	//Theoretical QKS in 2D case:
	double rho = myPdf.GetParameter(2);//correlationCoefficient(data);
	cout<<"rho used: "<<rho<<endl<<"estimated rho = "<<correlationCoefficient(data)<<" rho from bigaussian pdf = "<<myPdf.GetParameter(2)<<endl;
	double lambda_coeff = sqrt(SampleSize)/(1.+sqrt(1.-(pow(rho,2)*(0.25-0.75/sqrt(SampleSize)))));
	int trunc = 1000;
	double DminRange, DmaxRange;
	MCresults.first->GetRange(DminRange,DmaxRange);
	auto QKSf = (TF1*)(QKS_func(lambda_coeff, trunc,DminRange,DmaxRange).Clone("theoreticalQKS"));
	if(SampleSize<20)
		cout<<"WARNING: this expression for QKS in the 2D case is valid only for NSamples>20"<<endl;
	//draw
	if( drawOpt.Contains("draw") and not AlreadyDrawn){
		TString name ="2D MC and QKS distrib";
		auto myC = new TCanvas(name, name, 500, 0, 1000, 500);
		myC->cd();
		//generate the empirical CDF from the vector of the minimum distances, this will be used to estimate the t-values by reversing it
		MCresults.first->Draw();	
		MCresults.second->Draw("same");
		QKSf->SetLineColor(kTeal);
		QKSf->Draw("same");
		myC->BuildLegend();
		myC->Update();
		AlreadyDrawn = true;
	}
	
	//now evaluate distance and t-value for our data to see if they match or not with myPdf
	//using both MC and theoretical QKS
	double D = KolmogorovSmirnovDistance2DCustom(myPdf, data);
	double t = QKSf->Eval(D), tMC=MCresults.first->Eval(D);
	return {t,tMC,D};
}
//extension of the KS test to 2 dimensions based on custom method
//in this function there is no fit of the data
//we just generate data according to a 2D distribution
//and we use KS test to establish agreement
void KS_Test2D_Custom(int SampleSize = 100){
	
	int N_MC = 100;
	
	//prepare data
	double muX = -4, muY = +3;
	double sigmaX = 2, sigmaY = 2.8;
	double rho = 0.8; //correlation btw -1 and 1
	double x0, y0, x1, y1;
	x0 = std::max( abs(std::max(muX+5*sigmaX,muY+5*sigmaY)), abs(std::min(muX-5*sigmaX,muY-5*sigmaY)) );
	y0 = -x0; x1 = x0; y1 = x0; x0 = -x0;
	TF2* myPdfGen = new TF2("myPdfGen", "ROOT::Math::bigaussian_pdf(x,y,[0],[1],[2],[3],[4])", x0, x1, y0, y1);
	myPdfGen->SetParameters(sigmaX, sigmaY, rho, muX, muY);
	auto data = GenData2D(myPdfGen, SampleSize);
	TH2D* h1 = new TH2D("h1","h1",1000,x0,x1,1000,y0,y1);
	for( const auto & p : data )
		h1->Fill(p.first, p.second);
	
	vector<double> epss = {0, .0001, .001, .01, .05, .08, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .66, .7, .75, .8, .9, 1 };
	vector<double> ts = epss;
	vector<double> tsMC = epss;
	vector<double> Ds = epss;
	{
		int i = 0;
		for( auto & eps : ts ){
			TF2 myPdf = TF2("myPdfModified", "ROOT::Math::bigaussian_pdf(x,y,[0],[1],[2],[3],[4])", x0, x1, y0, y1);
			myPdf.SetParameters(sigmaX-eps, sigmaY-eps, rho, muX+eps, muY-eps);
			
			auto result = KolmogorovSmirnov2DCustom(*((TF2*)myPdf.Clone("pdf_compare")),data, N_MC); //KS
			eps = std::get<0>(result);
			tsMC[i] = std::get<1>(result);
			Ds[i] = std::get<2>(result);
			i++;
		}
	}
	
	TString namemg = "2DcustomResults;eps;prob";
	TCanvas* myC = new TCanvas(namemg,namemg,700,700);
	myC->Divide(2,1);
	auto mySubC = myC->cd(1);
	auto mg = new TMultiGraph();
	mg->SetNameTitle(namemg,namemg);
	auto myG = new TGraph(epss.size(), &epss[0], &ts[0]);
	myG->SetTitle("prob from QKS;prob;eps");
	myG->SetMarkerStyle(20);
	myG->SetMarkerColor(kTeal);
	mg->Add(myG);
	auto myG1 = new TGraph(epss.size(), &epss[0], &tsMC[0]);
	myG1->SetTitle("prob from MC;prob;eps");
	myG1->SetMarkerStyle(4);
	myG1->SetMarkerColor(kRed);
	mg->Add(myG1);
	mg->Draw("ap");
	mySubC->BuildLegend();
	//myC->SetLogx();
	myC->cd(2);
	auto myG2 = new TGraph(epss.size(), &epss[0], &Ds[0]);
	myG2->SetTitle("distance;eps;D");
	myG2->SetMarkerStyle(20);
	myG2->SetMarkerColor(kBlack);
	myG2->Draw("ap");	
	myC->Update();

}

void eKStensions(){}

int main(){
	gStyle->SetOptFit(111111);
	auto myApp = new TApplication("myApp", NULL, NULL);

	//Example_KSvsChi2_1D();
	int N = 1e04, SampleSize = 1000;
	//Lilliefors_Consistency_Test(N, SampleSize);
	//Extended_Lilliefors_Test(N, SampleSize);

	//KS_Test2D_ROOT();
	KS_Test2D_Custom();

	eKStensions();
	myApp->Run();
	return 0;
}
