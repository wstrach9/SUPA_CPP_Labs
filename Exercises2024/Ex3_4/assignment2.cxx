/*
William Strachan
02/12/2024
Objectives: 
	1) Fit mystery dataset to a set of finite distribution functions
	2) Determine the best-fit distribution as a guess of which distribution the data is sampled from
	3) Create an additional sample from this best-fit distribution using Metropolis Method
	4) Comment on accuracy of the Metropolis method, improve by playing with the Normal standard distribution?
	5) Estimate pi to 10dp by calling a radius and a number of pseudo-random numbers.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <numbers>
#include <algorithm>
#include <numeric>

#include "FiniteFunctions.h"

using std::numbers::pi;

int main(){
	//data reading
	std::string path = "Outputs/data/";
	std::string file = "MysteryData00204.txt";
	std::string pathfile = path + file;
	std::vector<double> mysterydata = readmysterydata(pathfile);
	//data in a vector of floats.
	
	//output filenames (will be appended .png or .data)
	std::string Normalout = "NormalModel";
	std::string CauchyLorentzout = "CauchyLorentzModel";
	std::string RevCrysBallout = "RevCrysBallModel";
	
	bool askforrange = true;
	bool askforsigmas = true;
	bool askforgammas = true;
	bool askforNChi2bins = false;
	bool plotNormal = true;
	bool calcNormal = true;
	bool plotCauchy = true;
	bool calcCauchy = true;
	bool optimiseCauchyLorentz = true;
	bool optimiseNormal = true;
	bool calcinvxsqr = true;
	bool plotinvxsqr = true;
	bool calcRevCrysBall = true;
	bool askforsigmasR = true;
	bool askforns = true;
	bool askforalphas = true;
	bool askforRevCrysBall_optimisecycles = true;
	bool plotRevCrysBall = true;
	bool optimiseRevCrysBall = true;
	int RevCrysBall_optimisecycles;

	//global plotting variables
	std::pair<double,double> GlobalRange = RangeMinMax(askforrange);
	int Nbins = 30; //number of bins for mysterydata histogram
	double mysterymean = std::accumulate(mysterydata.begin(), mysterydata.end(), 0u)/mysterydata.size();
	
	//Plot invxsquared
	FiniteFunction invxsquared;
	invxsquared.setRangeMin(GlobalRange.first);
	invxsquared.setRangeMax(GlobalRange.second);
	invxsquared.plotData(mysterydata, Nbins);
	/*
	//Plot Normal
	if (calcNormal){
		Normal Final_Normal(GlobalRange.first, GlobalRange.second, Normalout);
		if (optimiseNormal){
			// we will take the mean as the sample mean for brevity.
			std::vector<double> sigmas;
			sigmas = parametersweep("sigma", askforsigmas, Final_Normal);
			std::vector<Normal> Normals;
			for (double sigma : sigmas){
				Normals.push_back(Normal(GlobalRange.first, GlobalRange.second, mysterymean, sigma, Normalout));
			}
			std::vector<double> chi2sweep;
			for (Normal item : Normals){
				chi2sweep.push_back(item.CalcChi2(&mysterydata));
			}
			std::pair<int,double> optimumsigma = findmin(chi2sweep);
			std::cout << "Best-fit standard deviation in sweep = " << sigmas[optimumsigma.first] << std::endl;
			std::cout << "With corresponding R-Squared value = " << optimumsigma.second << std::endl;
			std::cout << "Neighbouring iterations: ";
			if (optimumsigma.first > 0){
				std::cout << "sigma = " << sigmas[optimumsigma.first - 1] << ", R-Squared = " << chi2sweep[optimumsigma.first - 1] << ";";
			}
			else{
				std::cout << "no smaller iterations of sigma; "
			}
			if (optimumsigma.first < sigmas.size() - 1){
				std::cout << "sigma = " << sigmas[optimumsigma.first + 1] << ", R-Squared = " << chi2sweep[optimumsigma.first + 1] << ".";
			}
			else{
				std::cout << "no larger iterations of sigma.";
			}
			Final_Normal.setmu(mysterymean);
			Final_Normal.setsigma(optimumsigma.second);
		}
		else{
			std::cout << "Setting defaults for Normal Distribution mu = 0; sigma = 1.";
			Final_Normal.setmu();
			Final_Normal.setsigma();
		}
		if (plotNormal){
			Final_Normal.plotData(mysterydata, Nbins, true);
			Final_Normal.plotFunction();
			Final_Normal.printInfo();
			~Final_Normal();
		}
	}
	
	//Plot Cauchy-Lorentz
	if (calcCauchy){
		CauchyLorentz Final_CauchyLorentz(GlobalRange.first, GlobalRange.second, CauchyLorentzout);
		if (optimiseCauchyLorentz){
			//Assuming the sample mean is a good-enough estimate of x0. For brevity, we will not optimise for x0.
			std::vector<double> gammas;
			gammas = parametersweep("gamma", askforgammas);
			std::vector<CauchyLorentz> CauchyLorentzs;
			for (double gamma : gammas){
				CauchyLorentzs.push_back(CauchyLorentz(GlobalRange.first, GlobalRange.second, mysterymean, gamma, CauchyLorentzout));
			}
			std::vector<double> chi2sweep;
			for (CauchyLorentz item : CauchyLorentzs){
				chi2sweep.push_back(item.CalcChi2(&mysterydata));
			}
			std::pair<int,double> optimumgamma = findmin(chi2sweep);
			std::cout << "Best-fit gamma in sweep = " << gammas[optimumgamma.first] << std::endl;
			std::cout << "With corresponding R-Squared value = " << optimumgamma.second << std::endl;
			std::cout << "Neighbouring iterations: ";
			if (optimumgamma.first > 0){
				std::cout << "gamma = " << gammas[optimumgamma.first - 1] << ", R-Squared = " << chi2sweep[optimumgamma.first - 1] << ";";
			}
			else{
				std::cout << "no smaller iterations of gamma; "
			}
			if (optimumgamma.first < gammas.size() - 1){
				std::cout << "gamma = " << gammas[optimumgamma.first + 1] << ", R-Squared = " << chi2sweep[optimumgamma.first + 1] << ".";
			}
			else{
				std::cout << "no larger iterations of gamma.";
			}
			Final_CauchyLorentz.setmu(mysterymean);
			Final_CauchyLorentz.setgamma(optimumgamma.second);
		}
		else{
			std::cout << "Setting defaults for Cauchy-Lorentz Distribution mu = 0; gamma = 1.";
			Final_CauchyLorentz.setmu();
			Final_CauchyLorentz.setgamma();
		}
		if (plotCauchyLorentz){
			Final_CauchyLorentz.plotData(mysterydata, Nbins, true);
			Final_CauchyLorentz.plotFunction();
			Final_CauchyLorentz.printInfo();
			~Final_CauchyLorentz();
		}
	}
	*/
	
	
	//plot Reverse-Crystal-Ball
	if (calcRevCrysBall){
		RevCrysBall Final_RevCrysBall(GlobalRange.first, GlobalRange.second, RevCrysBallout);
		if (optimiseRevCrysBall){
			if (askforRevCrysBall_optimisecycles){
				std::cout << "Parameters mu, sigma, n, alpha will be optimised one by one over a cyclical pattern: [ABCDADCBDBAC]..." << std::endl;
				std::string prompt = "Enter a number of cycles to optimise over. (integer > 0)";
				RevCrysBall_optimisecycles = intcin(prompt);
				if (RevCrysBall_optimisecycles < 1){
					std::cout << "Must run at least 1 cycle. Running 1 cycle...";
					RevCrysBall_optimisecycles = 1;
				}
			}
			else{
				RevCrysBall_optimisecycles = 3;
			}
			bool askforguesses = true;
			std::vector<std::string> paramnames = {"mu", "sigma", "n", "alpha"};
			std::vector<std::pair<double,double>> guessrange = getguesses(askforguesses, 4, paramnames, mysterymean);
			std::vector<std::vector<std::pair<double,double>>> optimisationhistory;
			optimisationhistory = TourDeSequence(paramnames, guessrange, 30, mysterydata, RevCrysBall_optimisecycles,true);
			Final_RevCrysBall.setmu(optimisationhistory[-1][0].first);
			Final_RevCrysBall.setsigma(optimisationhistory[-1][1].first);
			Final_RevCrysBall.setn(optimisationhistory[-1][2].first);
			Final_RevCrysBall.setalpha(optimisationhistory[-1][3].first);
		}
		
		/*
		else{
			std::cout << "Setting defaults for Cauchy-Lorentz Distribution mu = 0; gamma = 1.";
			Final_RevCrysBall.setmu(mysterymean);
			Final_RevCrysBall.setsigma();
			Final_RevCrysBall.setchi2(Final_RevCrysBall.CalcChi2(mysterydata));
		}
		*/
		if (plotRevCrysBall){
			Final_RevCrysBall.plotData(mysterydata, Nbins, true);
			Final_RevCrysBall.plotFunction();
			Final_RevCrysBall.printInfo();
			Final_RevCrysBall.~RevCrysBall();
		}
	}

	
	return 0;
}




