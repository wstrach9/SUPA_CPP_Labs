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
#include <string>
#include <vector>
#include "FiniteFunctions.h"
#include <filesystem> //To check extensions in a nice way

#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 


int main(){
	std::string path = "Outputs/data/";
	std::string mysteryfile = "MysteryData12005.txt";
	
	std::string invx_out = "invxsquared";
	std::string Normal_out = "Normal";
	std::string CL_out = "CauchyLorentz";
	std::string RCB_out = "RevCrysBall";
	
	int Nbins = 30;
	bool askforrange = true;
	bool askforparams = false;
	bool loud = true;
	
	std::vector<double> mysterydata = readmysterydata(path + mysteryfile);
	
	/*For each class, member variables to set:
		m_rmin
		m_rmax
		m_outfile
		m_params	(for Normal, C-L, RCB)*/
	
	//create instances (_inst for instance)
	std::cout << "Plotting 1/(1+x^2) distribution." << std::endl;
	FiniteFunction invx_inst;
	invx_inst.setRangeMinMax(mysterydata, askforrange, loud);
	invx_inst.setOutfile(invx_out);
	invx_inst.plotData(mysterydata, Nbins, true);
	invx_inst.~FiniteFunction();

	std::cout << "Plotting Normal distribution." << std::endl;
	Normal Normal_inst;
	Normal_inst.setRangeMinMax(mysterydata, askforrange, loud);
	Normal_inst.setOutfile(Normal_out);
	Normal_inst.setparams(askforparams, mysteryfile, loud, path);
	Normal_inst.plotData(mysterydata, Nbins, true);
	Normal_inst.plotFunction();
	Normal_inst.~Normal();
	
	/*std::cout << "Plotting Cauchy-Lorentz distribution." << std::endl;
	CauchyLorentz CL_inst;
	CL_inst.setRangeMinMax(mysterydata, askforrange, loud);
	CL_inst.setOutfile(CL_out);
	CL_inst.setparams(askforparams, mysteryfile, loud, path);
	CL_inst.plotData(mysterydata, Nbins, true);
	CL_inst.~CauchyLorentz();

	std::cout << "Plotting Reverse Crystal Ball distribution." << std::endl;
	RevCrysBall RCB_inst;
	RCB_inst.setRangeMinMax(mysterydata, askforrange, loud);
	RCB_inst.setOutfile(RCB_out);
	RCB_inst.setparams(askforparams, mysteryfile, loud, path);
	RCB_inst.plotData(mysterydata, Nbins, true);
	RCB_inst.~RevCrysBall();*/
	
	return 0;
}