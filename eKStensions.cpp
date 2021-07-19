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
#include <numeric> //std::iota

//ROOT libraries
#include <TRandom3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TColor.h>
#include <TString.h>
#include <TApplication.h>
#include <TAxis.h>
#include <TH1.h>

//STL usage
using std::cout, std::endl, std::to_string;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds, std::chrono::seconds, std::chrono::minutes;
using std::vector;
using std::string;
using std::setw, std::setprecision;
using std::srand, std::rand;



void es4(){} //utility function to open the file as a ROOT MACRO

int main(){
	auto myApp = new TApplication("myApp", NULL, NULL);
	//es4();
	myApp->Run();
	return 0;
}
