#include <iostream>
#include <string>
#include <vector>
#include "FiniteFunctions.h"
#include <filesystem> //To check extensions in a nice way
#include <numbers>
#include <algorithm>

#include "/workspaces/SUPA_CPP_Labs/GNUplot/gnuplot-iostream.h" //Needed to produce plots (not part of the course) 

using std::numbers::pi;

//Global functions:
std::vector<double> readmysterydata(std::string filepath){
	std::vector<double> points;
	std::ifstream input_file;
	input_file.open(filepath);  //opens file as input stream
	if (input_file.fail()){		//checks file opened
		std::cout<<"Couldn't open file: "<< filepath << std::endl;}
	std::string line;
	while (std::getline(input_file,line)){ //append all \n (or whitespace) separated lines to output vector
		points.push_back(std::stod(line));
	}
	input_file.close(); //file fully read, can close
	return points;
}

int intcin(std::string prompt){
	int input;
	std::cout << prompt << std::endl;
	std::cin >> input;
	while (!std::cin){
		std::cin.clear();
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cout << "Invalid. " + prompt << std::endl;
		std::cin >> input;
	}
	return input;
}

double doublecin(std::string prompt){
	double input;
	std::cout << prompt << std::endl;
	std::cin >> input;
	while (!std::cin){
		std::cin.clear();
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cout << "Invalid. " + prompt << std::endl;
		std::cin >> input;
	}
	return input;
}

bool yncin(std::string prompt){
	std::cout << prompt + " (y/n)" << std::endl;
	std::string input;
	std::cin >> input;
	while (input != "y" && input != "Y" && input != "n" && input != "N"){
		std::cout << "Invalid. " + prompt + " (y/n)" << std::endl;
		std::cin >> input;
	}
	if (input == "Y" || input == "y"){
		return true;
	}
	else{
		return false;
	}
}

// <Reverse Crystal Ball optimisation functions>
std::vector<std::pair<double,double>> getguesses(bool askforguesses, std::vector<std::string> paramnames, double mysterymean){
	//general
	std::string prompt;
	std::string warning;
	int cutoff;
	double tempguess;
	double temprange;
	std::vector<std::pair<double,double>> output;
	for (std::string paramname : paramnames){
		if (askforguesses){
			if (paramname == "mu" || paramname == "x0"){
				prompt = "Use sample mean (" + std::to_string(mysterymean) + ") as estimate of " + paramname + "?";
				if (yncin(prompt)){
					tempguess = mysterymean;
				}
				else{
					prompt = "Enter an estimate of "+ paramname +". (double)";
					tempguess = doublecin(prompt);
				}
				prompt = "How far (+/-) from this estimate to test fit? (double)";
				temprange = doublecin(prompt);
			}
			else{
				if(paramname == "alpha"){
					cutoff = 1;
				}
				else{cutoff = 0;}
				std::cout << paramname + " cannot be less than " << cutoff << ". If indicated range crosses " << cutoff << ", the estimate will be increased accordingly." << std::endl;
				prompt = "Enter an estimate of "+ paramname +". (double)";
				tempguess = doublecin(prompt);
				prompt = "How far (+/-) from this estimate to test fit? (double)";
				temprange = doublecin(prompt);
				if (!(tempguess-temprange > cutoff)){
					// add a little to get around zero division errors
					tempguess = cutoff + temprange + 1e-5;
				}
			}
		}
		else{//default guesses
			if(paramname == "mu" || paramname == "x0"){
				tempguess = mysterymean;
				temprange = 5;
			}
			else{
				tempguess = 2;
				temprange = 1-1e-5;
			}
		}
		output.push_back(std::make_pair(tempguess,temprange));
	}
	return output;
}

void changeguess(double& guess, double& range, double best, std::string paramname){
    double max = guess + range;
    double min = guess - range;
    double tomax = max - best;
    double tomin = best - min;
    //if tomax or tomin too small, reduce range by a little
    //if tomax or tomin small, redule by more.
    double newrange = std::max(tomax,tomin) - range;
    //if in centre of range already, new range is zero
    //if on the edge of range, range remains the same
    range = newrange;
    guess = best;
	checkzerox(guess, range, paramname);
}

void checkzerox(double& guess, double range, std::string paramname){
	if (paramname == "alpha"){
		if (!(guess-range > 1)){
			guess = 1 + range + 1e-5;
		}
	}
	else if (paramname == "sigma" || paramname == "gamma" || paramname == "n"){
		if (!(guess-range > 1)){
			guess = 1e-5 + range;
		}
	}
}

std::vector<std::vector<std::pair<double,double>>> TourDeSequence(std::vector<std::string> paramnames, int Nbins, std::pair<double,double> GlobalRange, std::string RevCrysBallout, std::vector<std::pair<double, double>>& guessrange, int halfResolution, std::vector<double> data, int n_cycles, bool loud = false){ 
	//get it? because we're cycling now, hahaha...
	//returns history of chi2 history
	std::vector<std::vector<std::pair<double,double>>> output;
	//calculate starting Chi2:
	RevCrysBall StartDistribution = RevCrysBall(GlobalRange.first, GlobalRange.second, guessrange[0].first, guessrange[1].first, guessrange[2].first, guessrange[3].first, RevCrysBallout);
	int Parambins = (halfResolution-1)/2; //halfResolution is the number of additional points created above/below guessrange[:].first
	double startchi2 = StartDistribution.CalcChi2(data, Parambins);
	//starting Chi2 calculated
	if (loud){
		std::cout << "Starting optimisation of Reverse Crystal Ball Distribution with:" <<std::endl;
		std::cout << "mu = " << guessrange[0].first << " +/- " << guessrange[0].second << ";  ";
		std::cout << "sigma = " << guessrange[1].first << " +/- " << guessrange[1].second << ";  ";
		std::cout << "n = " << guessrange[2].first << " +/- " << guessrange[2].second << ";  ";
		std::cout << "alpha = " << guessrange[3].first << " +/- " << guessrange[3].second << std::endl;
		std::cout << "and starting Chi-squared fit value of: " << startchi2 << ".";
	}
	double chi2improvement; // smaller is better, Chi2 dropping by more per cycle. We will shrink the error for precision accordingly.
	for (int i = 0; i < n_cycles; i++){
		std::vector<std::pair<double,double>> onecycle = cycleoptimise(Nbins, RevCrysBallout, GlobalRange, guessrange,  paramnames, halfResolution, data);
		if (i == 0){
			chi2improvement = onecycle[-1].second/startchi2;
		}
		else{
			chi2improvement = onecycle[-1].second/onecycle[-2].second;
		}
		for (int j = 0; j < 4; j++){
			guessrange[j].second = guessrange[j].second*chi2improvement;
		}
		if (loud){
			std::cout << "after " << i+1 << " cycles, we have optimised to:"<< std::endl; 
			std::cout << "mu = " << guessrange[0].first << " +/- " << guessrange[0].second << ";  ";
			std::cout << "sigma = " << guessrange[1].first << " +/- " << guessrange[1].second << ";  ";
			std::cout << "n = " << guessrange[2].first << " +/- " << guessrange[2].second << ";  ";
			std::cout << "alpha = " << guessrange[3].first << " +/- " << guessrange[3].second << std::endl;
			std::cout << "and starting Chi-squared fit value of: " << startchi2 << ".";
		}
		output.push_back(onecycle);
	}
	return output;
}

std::vector<std::pair<double,double>> cycleoptimise(int Nbins, std::string RevCrysBallout, std::pair<double,double> GlobalRange, std::vector<std::pair<double, double>>& guessrange, std::vector<std::string> paramnames, int halfResolution, std::vector<double> data){
	//outputs array (length 12) of optimised parameter values and chi2 values after each optimisation run in cycle

	//guessrange contains pairs of value guesses and error guesses (respectively) for each mu, sigma, n, alpha.
	//guessrannge[#] and param_spaces[#] corresponds to: # = 0 mu; # = 1 sigma; # = 2 n; # = 3 alpha
	//use these to create a basis of each to iterate off.
	std::vector<std::vector<double>> param_spaces;
	std::vector<std::pair<double,double>> onecycle;
	for (int i = 0; i < 4; i++){
		param_spaces.push_back(parameterspace(guessrange, halfResolution, paramnames, i));
		//param_spaces[0:3][0:(2*halfResolution+1)]
	}
	int paramindex;
	bool cyclecomplete = false;
	int place = 0;
	while (!cyclecomplete){
		paramindex = fourcycle(place, cyclecomplete);
		onecycle.push_back(singleoptimise(Nbins, GlobalRange, RevCrysBallout, param_spaces, data, paramindex));
		guessrange[paramindex].first = onecycle[-1].first;
	}
	return onecycle;
}

std::pair<double,double> singleoptimise(int Nbins, std::pair<double,double> GlobalRange, std::string RevCrysBallout, std::vector<std::vector<double>> param_spaces, std::vector<double> data, int paramindex){
	//return optimised pair[parameter, chi2]
	int Parambins = param_spaces[0].size();
	int midpoint = (Parambins - 1)/2;
	std::vector<RevCrysBall> RevCrysBalls;
	std::vector<double> chi2vals;
	for (int i = 0; i < param_spaces[0].size(); i++){
		RevCrysBalls.push_back(RevCrysBall());
		RevCrysBalls[-1].setRangeMin(GlobalRange.first);
		RevCrysBalls[-1].setRangeMax(GlobalRange.second);
		RevCrysBalls[-1].setOutfile(RevCrysBallout);
		switch(paramindex){
			case 0:
				RevCrysBalls[-1].setmu(param_spaces[0][i]);
				RevCrysBalls[-1].setsigma(param_spaces[1][midpoint]);
				RevCrysBalls[-1].setn(param_spaces[2][midpoint]);
				RevCrysBalls[-1].setalpha(param_spaces[3][midpoint]);
				break;
			case 1:
				RevCrysBalls[-1].setmu(param_spaces[0][midpoint]);
				RevCrysBalls[-1].setsigma(param_spaces[1][i]);
				RevCrysBalls[-1].setn(param_spaces[2][midpoint]);
				RevCrysBalls[-1].setalpha(param_spaces[3][midpoint]);
			case 2:
				RevCrysBalls[-1].setmu(param_spaces[0][midpoint]);
				RevCrysBalls[-1].setsigma(param_spaces[1][midpoint]);
				RevCrysBalls[-1].setn(param_spaces[2][i]);
				RevCrysBalls[-1].setalpha(param_spaces[3][midpoint]);
			case 3:
				RevCrysBalls[-1].setmu(param_spaces[0][midpoint]);
				RevCrysBalls[-1].setsigma(param_spaces[1][midpoint]);
				RevCrysBalls[-1].setn(param_spaces[2][midpoint]);
				RevCrysBalls[-1].setalpha(param_spaces[3][i]);
		}
		chi2vals.push_back(RevCrysBalls[-1].CalcChi2(data, Nbins));
	}
	std::pair<int, double> argmin = findmin(chi2vals);
	return std::make_pair(param_spaces[paramindex][argmin.first], argmin.second);
}

int fourcycle(int &place, bool &cyclecomplete){ 
	//outputs useable paramindex
	cyclecomplete = false;
	std::vector<int> sequence = {0,1,2,3,0,3,2,1,3,1,0,2};
	int turn;
	if ((place > -1 && place < 12)){
		turn = sequence[place];
	}
	if ((place > -1 && place < 11)){
		place++;
	}
	else{
		place = 0;
		cyclecomplete = true;
	}
	return turn;
}

std::vector<double> parameterspace(std::vector<std::pair<double, double>> guessrange, int halfResolution, std::vector<std::string> paramnames, int paramindex){
	std::vector<double> param_space;
	checkzerox(guessrange[paramindex].first, guessrange[paramindex].second, paramnames[paramindex]);
	double stepwidth = guessrange[paramindex].second/halfResolution;
	for (int i = -halfResolution; i < halfResolution + 1; i++){
		param_space.push_back(guessrange[paramindex].first + i*stepwidth);
	}
	return param_space;
}

std::pair<int, double> findmin(std::vector<double> chi2vals){	
	std::vector<double>::iterator minimum = std::min_element(chi2vals.begin(), chi2vals.end());
	int minimum_index = std::distance(chi2vals.begin(), minimum);
	return std::make_pair(minimum_index, *minimum);
}

std::vector<double> parametersweep(bool askforparams, std::string paramname, double& Lowbound, double& Uppbound, int& Niterations){
	//outputs a parameter
	if (askforparams){
		std::cout << "Enter Lower bound to optimise " + paramname + ". (double > 0)" << std::endl;
		std::cin >> Lowbound;
		while (!std::cin || !(Lowbound > 0)){ // bool flag if input i'n't int or double trouble
			std::cin.clear(); //clear flag
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //clear cin input
			std::cout << "Invalid. Enter Lower bound to optimise " + paramname + ". (double > 0)" << std::endl;
			std::cin >> Lowbound;
		}
		std::cout << "Enter Upper bound to optimise " + paramname + ". (double)" << std::endl;
		std::cin >> Uppbound;
		while (!std::cin){ // bool flag if input i'n't int or double trouble
			std::cin.clear(); //clear flag
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //clear cin input
			std::cout << "Invalid. Enter Upper bound to optimise " + paramname + ". (double > 0)" << std::endl;
			std::cin >> Lowbound;
		}
		std::cout << "Enter number of iterations to optimise " + paramname + " over. (integer)" << std::endl;
		std::cin >> Niterations;
		while (!std::cin){ // bool flag if input i'n't int or double trouble
			std::cin.clear(); //clear flag
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //clear cin input
			std::cout << "Invalid. Enter number of iterations to optimise " + paramname + " over. (integer)" << std::endl;
			std::cin >> Niterations;
		}
	}
	if (Niterations < 2){
		std::cout << "Must iterate at least twice! Iterating twice..." << std::endl;
		Niterations = 2;
	}
	if (paramname == "n" && !(Lowbound >1)){
		std::cout << "n must be greater than 1. Setting lower bound to 2..." << std::endl;
		Lowbound = 2;
	}
	if (!(Uppbound > Lowbound)){
		std::cout << "Upper bound must be greater than lower bound. Setting Upper bound to lower bound + 1..." << std::endl;
		Uppbound = Lowbound + 1;
	}
	std::vector<double> parameters;
	double step = (Uppbound - Lowbound)/(Niterations-1);
	for (int i = 0; i < Niterations; i++){
		parameters.push_back(Lowbound + i*step);
	}
	return parameters;
}
// </Reverse Crystal Ball optimisations functions>

std::pair<double,double> RangeMinMax(bool askforrange, double rmin, double rmax){
	if (askforrange){
		std::cout << "Enter global Lower bound for fitting and plotting. (double)" << std::endl;
		std::cin >> rmin;
		while (!std::cin){ // bool flag if input i'n't int or double trouble
			std::cin.clear(); //clear flag
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //clear cin input
			std::cout << "Not numeric. Enter global Lower bound for fitting and plotting. (double)" << std::endl;
			std::cin >> rmin;
		}
		//std::cout << "Lower bound is: "<<rmin<<std::endl;
		std::cout << "Enter global Upper bound for fitting and plotting. (double)" << std::endl;
		std::cin >> rmax;
		while (!std::cin){ // bool flag if input i'n't int or double trouble
			std::cin.clear(); //clear flag
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //clear cin input
			std::cout << "Not numeric. Enter global Upper bound for fitting and plotting. (double)" << std::endl;
			std::cin >> rmax;
		}
		//std::cout << "Upper bound is: "<< rmax << std::endl;
	}
	return std::make_pair(rmin, rmax);
}


//Normal class definitions
//constructor
Normal::Normal(){
	FiniteFunction();
	this->checkPath("NormalPDF");
}
Normal::Normal(double range_min, double range_max, double mu, double sigma, std::string outfile){
	FiniteFunction(range_min,range_max,outfile);
		this->setmu(mu);
		this->setsigma(sigma);
}
//getters
double Normal::mu(){
return m_mu;}
double Normal::sigma(){
return m_sigma;}
//setters
void Normal::setmu(double MU = 0){
m_mu = MU;}
void Normal::setsigma(double SIGMA = 1){
m_sigma = SIGMA;}
//distribution function
double Normal::Normalpdf(double x){
	double normalisation = 1/(m_sigma*std::sqrt(2*pi));
	double exponential = std::exp(-0.5*std::pow((x-m_mu)/m_sigma,2));
	return (normalisation * exponential);
}
double Normal::callFunction(double x){
	return this->Normalpdf(x);
}
//end Normal definitions

//C-L definitions
//constructor
CauchyLorentz::CauchyLorentz(){
	FiniteFunction();
	this->checkPath("CauchyLorentzPDF");
}
CauchyLorentz::CauchyLorentz(double range_min, double range_max, double x0, double gamma,std::string outfile){
	FiniteFunction(range_min,range_max,outfile);
	this->setx0(x0);
	this->setgamma(gamma);
}
//getters
double CauchyLorentz::x0(){
return m_x0;}
double CauchyLorentz::gamma(){
return m_gamma;}
//setters
void CauchyLorentz::setx0(double X0){
m_x0 = X0;}
void CauchyLorentz::setgamma(double GAMMA){
	if (!(GAMMA > 0)){
		std::cout<<"Invalid entry of Cauchy Lorentz Distribution gamma (must be greater than 0). Setting to 1 instead..."<<std::endl;
		m_gamma = 1;
	}
	else{
		m_gamma = GAMMA;
	}
}
//pdf
double CauchyLorentz::CauchyLorentzpdf(double x){
	return 1/(pi*m_gamma*(1 + std::pow((x-m_x0)/m_gamma,2)));
}
double CauchyLorentz::callFunction(double x){
	return this->CauchyLorentzpdf(x);
}
//end C-L definitions

//RevCrysBall definitions
//constructor
RevCrysBall::RevCrysBall(){
	FiniteFunction();
	this->checkPath("RevCrysBallPDF");
}
/*
RevCrysBall::RevCrysBall(double range_min, double range_max, std::string outfile){
	FiniteFunction(range_min,range_max,outfile);
}*/
RevCrysBall::RevCrysBall(double range_min, double range_max, double mu, double sigma, double n, double alpha, std::string outfile){
	FiniteFunction(range_min,range_max,outfile);
	this->setmu(mu);
	this->setsigma(sigma);
	this->setn(n);
	this->setalpha(alpha);
}
//getters
double RevCrysBall::mu(){
return m_mu;}
double RevCrysBall::sigma(){
return m_sigma;}
double RevCrysBall::alpha(){
return m_alpha;}
double RevCrysBall::n(){
return m_n;}
//setters
void RevCrysBall::setmu(double MU){
m_mu = MU;}
void RevCrysBall::setsigma(double SIGMA){
m_sigma = SIGMA;}
void RevCrysBall::setalpha(double ALPHA){
	if (!(ALPHA > 0)){
		std::cout<<"Invalid entry of Reverse Crystal Ball Distribution alpha (must be greater than 0). Setting to 1 instead..."<<std::endl;
		m_alpha = 1;
	}
	else{
		m_alpha = ALPHA;
	}
}
void RevCrysBall::setn(double N){
	if (!(N > 0)){
		std::cout<<"Invalid entry of Reverse Crystal Ball Distribution n (must be greater than 1). Setting to 2 instead..."<<std::endl;
		m_n = 2;
	}
	else{
		m_n = N;
	}
}
//intermediates A,B,C,D
double RevCrysBall::calcA(){
return std::pow(m_n/m_alpha,m_n)*std::exp(-0.5*m_alpha*m_alpha);}
double RevCrysBall::calcB(){
return (m_n/m_alpha)-m_alpha;}
double RevCrysBall::calcC(){
return (m_n/(m_alpha*(m_n - 1)))*std::exp(-0.5*m_alpha*m_alpha);}
double RevCrysBall::calcD(){
return std::sqrt(0.5*pi)*(1 + std::erf(m_alpha/std::sqrt(2)));}
double RevCrysBall::calcN(){
return 1/(m_sigma*(this->calcC() + this->calcD()));}
//distribution function
double RevCrysBall::RevCrysBallpdf(double x){
	double RevCrysBallout;
	double discriminant = ((x - m_mu)/m_sigma);
	if (discriminant > - m_alpha){
		RevCrysBallout = this->calcN()*std::exp(-0.5*discriminant*discriminant);
	}
	else{
		RevCrysBallout = this->calcN()*this->calcA()*std::pow((this->calcB() - discriminant),-m_n);
	}
	return RevCrysBallout;
}
double RevCrysBall::callFunction(double x){
	return this->RevCrysBallpdf(x);
}

//end RevCrysBall definitions

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Original file contents below 
	//exceptions:
	//integrate function written.
	//Chi Squared test to improve fits
	//chi2 member variable to compare goodness of fit
using std::filesystem::path;

//Empty constructor
FiniteFunction::FiniteFunction(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
}

//initialised constructor
FiniteFunction::FiniteFunction(double range_min, double range_max, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
}

//Plots are called in the destructor
//SUPACPP note: They syntax of the plotting code is not part of the course
FiniteFunction::~FiniteFunction(){
  Gnuplot gp; //Set up gnuplot object
  this->generatePlot(gp); //Generate the plot and save it to a png using "outfile" for naming 
}

/*
###################
//Setters
###################
*/ 
void FiniteFunction::setRangeMin(double RMin) {m_RMin = RMin;};
void FiniteFunction::setRangeMax(double RMax) {m_RMax = RMax;};
void FiniteFunction::setchi2(double CHI2) {m_Chi2 = CHI2;};
void FiniteFunction::setOutfile(std::string Outfile) {this->checkPath(Outfile);};

/*
###################
//Getters
###################
*/ 
double FiniteFunction::rangeMin() {return m_RMin;};
double FiniteFunction::rangeMax() {return m_RMax;};
double FiniteFunction::chi2() {return m_Chi2;};

/*
###################
//Function eval
###################
*/ 
double FiniteFunction::invxsquared(double x) {return 1/(1+x*x);};
double FiniteFunction::callFunction(double x) {return this->invxsquared(x);}; //(overridable)

double FiniteFunction::CalcChi2(std::vector<double>& points, int Nbins){
	//first x; second y in pairs.
	//compare observed weighted freq. to expected weighted freq.
	std::vector<std::pair<double,double>> m_data0 = makeHist(points, Nbins); //want something like m_data but can't redeclare!
	double cumulative = 0;
	for (int i = 0; i < Nbins; i++){
		cumulative += (m_data0[i].second - this->callFunction(m_data0[i].first))*(m_data0[i].second - this->callFunction(m_data0[i].first))/this->callFunction(m_data0[i].first);
	}
	return cumulative;
}

/*
###################
Integration by hand (output needed to normalise function when plotting)
###################
*/ 
double FiniteFunction::integrate(int Ndiv){ //private
  //ToDo write an integrator
  //Writing a Riemann sum algorithm.
  double stripwidth = (m_RMax - m_RMin)/Ndiv;
  double cumulative = 0;
  for (int i = 0; i < Ndiv; i++){
	cumulative += stripwidth*this->callFunction(m_RMin + stripwidth*(0.5 + i));
  }
  return cumulative;
}
double FiniteFunction::integral(int Ndiv) { //public
  if (Ndiv <= 0){
    std::cout << "Invalid number of divisions for integral, setting Ndiv to 1000" <<std::endl;
    Ndiv = 1000;
  }
  if (m_Integral == NULL || Ndiv != m_IntDiv){
    m_IntDiv = Ndiv;
    m_Integral = this->integrate(Ndiv);
    return m_Integral;
  }
  else return m_Integral; //Don't bother re-calculating integral if Ndiv is the same as the last call
}

/*
###################
//Helper functions 
###################
*/
// Generate paths from user defined stem
void FiniteFunction::checkPath(std::string outfile){
 	path fp = outfile;
 	m_FunctionName = fp.stem(); 
 	m_OutData = m_FunctionName+".data";
 	m_OutPng = m_FunctionName+".png";
}

//Print (overridable)
void FiniteFunction::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

/*
###################
//Plotting
###################
*/

//Hack because gnuplot-io can't read in custom functions, just scan over function and connect points with a line... 
void FiniteFunction::plotFunction(){
  m_function_scan = this->scanFunction(10000);
  m_plotfunction = true;
}

//Transform data points into a format gnuplot can use (histogram) and set flag to enable drawing of data to output plot
//set isdata to true (default) to plot data points in black, set to false to plot sample points in blue
void FiniteFunction::plotData(std::vector<double> &points, int Nbins, bool isdata){
  if (isdata){
    m_data = this->makeHist(points,Nbins);
    m_plotdatapoints = true;
  }
  else{
    m_samples = this->makeHist(points,Nbins);
    m_plotsamplepoints = true;
  }
}


/*
  #######################################################################################################
  ## SUPACPP Note:
  ## The three helper functions below are needed to get the correct format for plotting with gnuplot
  ## In theory you shouldn't have to touch them
  ## However it might be helpful to read through them and understand what they are doing
  #######################################################################################################
 */

//Scan over range of function using range/Nscan steps (just a hack so we can plot the function)
std::vector< std::pair<double,double> > FiniteFunction::scanFunction(int Nscan){
  std::vector< std::pair<double,double> > function_scan;
  double step = (m_RMax - m_RMin)/(double)Nscan;
  double x = m_RMin;
  //We use the integral to normalise the function points
  if (m_Integral == NULL) {
    std::cout << "Integral not set, doing it now" << std::endl;
    this->integral(Nscan);
    std::cout << "integral: " << m_Integral << ", calculated using " << Nscan << " divisions" << std::endl;
  }
  //For each scan point push back the x and y values 
  for (int i = 0; i < Nscan; i++){
    function_scan.push_back( std::make_pair(x,this->callFunction(x)/m_Integral));
    x += step;
  }
  return function_scan;
}

//Function to make histogram out of sampled x-values - use for input data and sampling
std::vector< std::pair<double,double> > FiniteFunction::makeHist(std::vector<double> &points, int Nbins){

  std::vector< std::pair<double,double> > histdata; //Plottable output shape: (midpoint,frequency)
  std::vector<int> bins(Nbins,0); //vector of Nbins ints with default value 0 
  int norm = 0;
  for (double point : points){
    //Get bin index (starting from 0) the point falls into using point value, range, and Nbins
    int bindex = static_cast<int>(floor((point-m_RMin)/((m_RMax-m_RMin)/(double)Nbins)));
    if (bindex<0 || bindex>=Nbins){
      continue;
    }
    bins[bindex]++; //weight of 1 for each data point
    norm++; //Total number of data points
  }
  double binwidth = (m_RMax-m_RMin)/(double)Nbins;
  for (int i=0; i<Nbins; i++){
    double midpoint = m_RMin + i*binwidth + binwidth/2; //Just put markers at the midpoint rather than drawing bars
    double normdata = bins[i]/((double)norm*binwidth); //Normalise with N = 1/(Ndata*binwidth)
    histdata.push_back(std::make_pair(midpoint,normdata));
  }
  return histdata;
}

//Function which handles generating the gnuplot output, called in destructor
//If an m_plot... flag is set, the we must have filled the related data vector
//SUPACPP note: They syntax of the plotting code is not part of the course
void FiniteFunction::generatePlot(Gnuplot &gp){

  if (m_plotfunction==true && m_plotdatapoints==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotdatapoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
  }
  else if (m_plotfunction==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title 'function'\n";
    gp.send1d(m_function_scan);
  }

  else if (m_plotdatapoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_data);
  }

  else if (m_plotsamplepoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_samples);
  }
}