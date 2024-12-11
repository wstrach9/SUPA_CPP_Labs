/*
William Strachan
started: 02/12/2024
last updated: 11/12/2024
Objectives: 
	1) Fit mystery dataset to a set of finite distribution functions
		- implemented by:
		a) asking for user input for fit parameters
		b) plotting and saving fit parameters to an archive document, along with a Chi-Squared value
		c) on subsequent plotting, user can choose to automatically choose the best-fit that exists in archive file
	2) Determine the best-fit distribution as a guess of which distribution the data is sampled from
		- user will have to do this manually
			- can run a bash macro (I assume) to populate the archive fit and run analysis on the Chi-squared values to be sure.
			- or just visually I guess.
	3) Create an additional sample from this best-fit distribution using Metropolis Method
		- done
	4) Comment on accuracy of the Metropolis method, improve by playing with the Normal standard distribution?
		- Accuracy is dependent on the choice of metropolis Normal standard deviation (sigma) and sample size.
			- Larger sample size is always better
			- sigma too large -> acceptance rate low -> spikey histogram
			- sigma too small -> won't pan the whole function domain -> results skewed towards start point
	5) avoid std::organ<brain> myBrain.insane(true);
		- mission failed, we'll get em next time.
*/

#include <iostream>
#include <string>
#include <vector>
#include "FiniteFunctions.h"
#include <filesystem> //To check extensions in a nice way

#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 


int main(int argc, char *argv[]){
	std::vector<std::string> inputs(argv,argv+argc);
	
	std::string mysteryfile = "MysteryData00204.txt";
	
	//no need to chage these
	std::string path = "Outputs/data/";
	std::string invx_out = "invxsquared";
	std::string Normal_out = "Normal";
	std::string CL_out = "CauchyLorentz";
	std::string RCB_out = "RevCrysBall";
	
	int Nbins = 30;
	bool askforrange = true;
	bool askedforrange = false;
	std::pair<double,double> GlobalRange;
	bool askforparams = true;
	bool loud = false;
	double metropolis_stddev = 30;
	
	if (argc > 1){
		if (inputs[1] == "--help" || inputs[1] == "-help" || inputs[1] == "help"){
			std::cout << "Expecting 1 to 5 additional arguements:" << std::endl;
			std::cout << "1) (required) MysteryData filename. (eg. MysteryData12005.txt)" << std::endl;
			std::cout << "2) Number of histogram bins for data and sample plotting. Default = 30" << std::endl;
			std::cout << "3) Input new fit parameters (mu, sigma, etc.) to try? (y/n) Default = y" << std::endl;
			std::cout << "3) y -> asks for user input, n -> searches for the best previous fit for the data." << std::endl;
			std::cout << "4) Input custom plotting range? (y/n) Default = y" << std::endl;
			std::cout << "4) y -> asks for user input, n -> sets range to integer floor/ceiling of the data." << std::endl;
			std::cout << "5) Print extra messages to terminal? (y/n) Default = n" << std::endl;
			return 0;
		}
	}

	switch(argc){
		case 1:{
			std::cout << "Please specify a file within " << path <<  std::endl;
			return -1;
			break;
		}
		case 2:{
			mysteryfile = inputs[1];
			std::cout << "Processing " + mysteryfile + " with " << Nbins << " histogram bins." <<std::endl;
			break;
		}
		case 3:{
			mysteryfile = inputs[1];
			try {
			Nbins = std::stoi(inputs[2]);
			}
			catch (std::invalid_argument& e) {
			std::cout << "Invalid input for number of histogram bins. Expected integer." <<std::endl;
			return -1;
			}
			std::cout << "Processing " + mysteryfile + " with " << Nbins << " histogram bins." <<std::endl;
			break;
		}
		case 4:{
			mysteryfile = inputs[1];
			try {
			Nbins = std::stoi(inputs[2]);
			}
			catch (std::invalid_argument& e) {
			std::cout << "Invalid input for number of histogram bins. Expected integer." <<std::endl;
			return -1;
			}
			if (inputs[3] == "y" || inputs[3] == "Y"){
				askforparams = true;
			}
			else if(inputs[3] == "n" || inputs[3] == "N"){
				askforparams = false;
			}
			else{
				std::cout << "Invalid input. Expected y or n." << std::endl;
				return -1;
			}
			std::cout << "Processing " + mysteryfile + " with " << Nbins << " histogram bins." <<std::endl;
			break;
		}
		case 5:{
			mysteryfile = inputs[1];
			try {
			Nbins = std::stoi(inputs[2]);
			}
			catch (std::invalid_argument& e) {
			std::cout << "Invalid input for number of histogram bins. Expected integer." <<std::endl;
			return -1;
			}
			if (inputs[3] == "y" || inputs[3] == "Y"){
				askforparams = true;
			}
			else if(inputs[3] == "n" || inputs[3] == "N"){
				askforparams = false;
			}
			else{
				std::cout << "Invalid input. Expected y or n." << std::endl;
				return -1;
			}
			if (inputs[4] == "y" || inputs[4] == "Y"){
				askforrange = true;
			}
			else if(inputs[4] == "n" || inputs[4] == "N"){
				askforrange = false;
			}
			else{
				std::cout << "Invalid input. Expected y or n." << std::endl;
				return -1;
			}
			std::cout << "Processing " + mysteryfile + " with " << Nbins << " histogram bins." <<std::endl;
			break;
		}
		case 6:{
			mysteryfile = inputs[1];
			try {
			Nbins = std::stoi(inputs[2]);
			}
			catch (std::invalid_argument& e) {
			std::cout << "Invalid input for number of histogram bins. Expected integer." <<std::endl;
			return -1;
			}
			//askforparams
			if (inputs[3] == "y" || inputs[3] == "Y"){
				askforparams = true;
			}
			else if(inputs[3] == "n" || inputs[3] == "N"){
				askforparams = false;
			}
			else{
				std::cout << "Invalid input. Expected y or n." << std::endl;
				return -1;
			}
			//askforrange
			if (inputs[4] == "y" || inputs[4] == "Y"){
				askforrange = true;
			}
			else if(inputs[4] == "n" || inputs[4] == "N"){
				askforrange = false;
			}
			else{
				std::cout << "Invalid input. Expected y or n." << std::endl;
				return -1;
			}
			//loud
			if (inputs[5] == "y" || inputs[5] == "Y"){
				loud = true;
			}
			else if(inputs[5] == "n" || inputs[5] == "N"){
				loud = false;
			}
			else{
				std::cout << "Invalid input. Expected y or n." << std::endl;
				return -1;
			}
			std::cout << "Processing " + mysteryfile + " with " << Nbins << " histogram bins." <<std::endl;
			break;
		}
		default:{
			std::cout << "Too many arguments given. Aborting..." << std::endl;
			return -1;
			break;
		}
	}

	std::vector<double> mysterydata = readmysterydata(path + mysteryfile);
	
	{
	std::cout << std::endl << "Plotting 1/(1+x^2) distribution." << std::endl;
	FiniteFunction invx_inst;
	invx_inst.setRangeMinMax(mysterydata, askedforrange, GlobalRange, askforrange, loud);
	invx_inst.setOutfile(invx_out);
	invx_inst.plotData(mysterydata, Nbins, true);
	invx_inst.plotFunction();
	std::vector<double> sampledata = invx_inst.metropolis(mysterydata.size(), metropolis_stddev);
	invx_inst.plotData(sampledata, Nbins, false);
	}

	{
	std::cout << std::endl << "Plotting Normal distribution." << std::endl;
	Normal Normal_inst;
	Normal_inst.setRangeMinMax(mysterydata, askedforrange, GlobalRange, askforrange, loud);
	Normal_inst.setOutfile(Normal_out);
	bool NormParams = askforparams;
	Normal_inst.setparams(NormParams, mysteryfile, loud, "previous_fits/");
	Normal_inst.plotData(mysterydata, Nbins, true);
	Normal_inst.plotFunction();
	if (NormParams){ //only add to history if using a new set of distribution parameters
		Normal_inst.printtofithist(true, mysteryfile, loud, "previous_fits/");
	}
	std::vector<double> sampledata = Normal_inst.metropolis(mysterydata.size(), metropolis_stddev);
	Normal_inst.plotData(sampledata, Nbins, false);
	}

	{
	std::cout << std::endl << "Plotting Cauchy-Lorentz distribution." << std::endl;
	CauchyLorentz CL_inst;
	CL_inst.setRangeMinMax(mysterydata, askedforrange, GlobalRange, askforrange, loud);
	CL_inst.setOutfile(CL_out);
	bool CLParams = askforparams;
	CL_inst.setparams(CLParams, mysteryfile, loud, "previous_fits/");
	CL_inst.plotData(mysterydata, Nbins, true);
	CL_inst.plotFunction();
	if (CLParams){//only add to history if using a new set of distribution parameters
		CL_inst.printtofithist(true, mysteryfile, loud, "previous_fits/");
	}
	std::vector<double> sampledata = CL_inst.metropolis(mysterydata.size(), metropolis_stddev);
	CL_inst.plotData(sampledata, Nbins, false);
	}

	{
	std::cout << std::endl << "Plotting Reverse Crystal Ball distribution." << std::endl;
	RevCrysBall RCB_inst;
	RCB_inst.setRangeMinMax(mysterydata, askedforrange, GlobalRange, askforrange, loud);
	RCB_inst.setOutfile(RCB_out);
	bool RCBParams = askforparams;
	RCB_inst.setparams(RCBParams, mysteryfile, loud, "previous_fits/");
	RCB_inst.plotData(mysterydata, Nbins, true);
	RCB_inst.plotFunction();
	if (RCBParams){//only add to history if using a new set of distribution parameters
		RCB_inst.printtofithist(true, mysteryfile, loud, "previous_fits/");
	}
	std::vector<double> sampledata = RCB_inst.metropolis(mysterydata.size(), metropolis_stddev);
	RCB_inst.plotData(sampledata, Nbins, false);
	}
	return 0;
}