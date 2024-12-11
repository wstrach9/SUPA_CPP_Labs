//FiniteFunctions backup

#include <iostream>
#include <string>
#include <vector>
#include "FiniteFunctions.h"
#include <filesystem> //To check extensions in a nice way
#include <numbers>
#include <random>

#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 


//global functions

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

bool yncin(std::string prompt){
	std::cout << prompt + " (y/n)" << std::endl;
	char yn;
	std::cin >> yn;
	int count = 0;
	while (yn != 'y' && yn != 'n' && yn != 'Y' && yn != 'N'){
		std::cout << "Invalid input. " << prompt + " (y/n)" << std::endl;
		std::cin >> yn;
		count++;
		if (count == 5){
			std::cout << "Too many invalid inputs. Defaulting to yes..." << std::endl;
			yn = 'y';
		}
	}
	if (yn == 'y' || yn == 'Y'){return true;}
	else{return false;};	
}

int intcin(std::string prompt){
	std::cout << prompt + " (integer)" << std::endl;
	int input;
	std::cin >> input;
	while (!std::cin){ // bool flag if input i'n't int
		std::cin.clear(); //clear flag
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //clear cin input
		std::cout << "Invalid input. " << prompt + " (integer)" << std::endl;
		std::cin >> input;
	}
	return input;
}

double doublecin(std::string prompt){
	std::cout << prompt + " (double)" << std::endl;
	double input;
	std::cin >> input;
	while (!std::cin){ // bool flag if input i'n't int
		std::cin.clear(); //clear flag
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //clear cin input
		std::cout << "Invalid input. " << prompt + " (double)" << std::endl;
		std::cin >> input;
	}
	return input;
}

std::pair<int, double> findargmin(std::vector<double> vector){	
	std::vector<double>::iterator minimum = std::min_element(vector.begin(), vector.end());
	int minimum_index = std::distance(vector.begin(), minimum);
	return std::make_pair(minimum_index, *minimum);
}

double findmin(std::vector<double> vector){	
	std::vector<double>::iterator minimum = std::min_element(vector.begin(), vector.end());
	return *minimum;
}

double findmax(std::vector<double> vector){
	std::vector<double>::iterator maximum = std::max_element(vector.begin(), vector.end());
	return *maximum;
}	

//class functions
	//added FiniteFunction functions
bool FiniteFunction::readfithistory(std::vector<double>& bestfit, std::string mysteryfile, bool loud, std::string path){
	// defaults: path - ""; loud - false
	// outputs distribution parameters for historical best fit.
	// open + read file to 2D vector
	std::string mysterynumbers; //create full history file name
	int check;
	for (size_t i = 0; i < mysteryfile.size(); i++){
		check = isdigit(mysteryfile[i]);
		if (check){
			mysterynumbers += mysteryfile[i];
		}
	}
	std::string historyfile = path + m_historyfile + "_fithistory_" + mysterynumbers + ".txt";
	//history file name generated
	if (!(std::filesystem::exists(historyfile))){
		std::cout << "No history of running this distribution before. Please input fit distribution parameters." << std::endl;
		return true;
	}
	else{
		std::ifstream input_file;
		input_file.open(historyfile);
		if (input_file.fail()){		//checks file opened
			std::cout<< "Couldn't open file: " << historyfile << std::endl;
		}
		std::vector<std::vector<double>> allfitvals;
		std::string line;
		while (std::getline(input_file,line)){
			std::stringstream ls(line);
			std::vector<double> linefitvals;
			while (ls.good()){
				std::string entry;
				std::getline(ls, entry, ',');
				linefitvals.push_back(std::stod(entry));
			}
			allfitvals.push_back(linefitvals);
		}
		input_file.close();
		std::vector<double> chi2vals;
		for (int i = 0; i < allfitvals.size(); i++){
			chi2vals.push_back(allfitvals[i][0]);
		}
		std::pair<int,double> historicalbest = findargmin(chi2vals);
		if (loud){
			std::cout << "Using historical best fit, with Chi Squared value of " << historicalbest.second;
			std::cout << " and parameters: " << std::endl;
		}
		int j = 1;
		for (std::string paramname : m_paramnames){
			if (loud){
				std::cout << paramname + ": " << allfitvals[historicalbest.first][j] << std::endl;
			}
			bestfit.push_back(allfitvals[historicalbest.first][j]);
			j++;
		}
	}
	return false;
}

double FiniteFunction::calcChi2(bool isdata){
	std::vector<std::pair<double,double>> hist_data;
	if (isdata){
		hist_data = m_data;
	}
	else{
		hist_data = m_samples;
	}
	double cumulative = 0;
	for (size_t i = 0; i < hist_data.size(); i++){
		cumulative += (hist_data[i].second - this->callFunction(hist_data[i].first))*(hist_data[i].second - this->callFunction(hist_data[i].first))/this->callFunction(hist_data[i].first);
	}
	return cumulative;
}





std::vector<double> FiniteFunction::metropolis(int Npoints, double stddev){
	std::vector<double> generatedsample;
	double x_i;
	double random_y;
	std::vector<double> to_test;
	double A;
	double T;
	std::random_device rd;
	std::mt19937 mtEngine{rd()};
	std::uniform_real_distribution<double> uniform_PDF{0, 1};
	x_i = uniform_PDF(mtEngine);
	x_i = x_i*(m_RMax-m_RMin) + m_RMin;
	generatedsample.push_back(x_i);
	std::normal_distribution<double> Normal_y_PDF{x_i, stddev};
	while (generatedsample.size() < Npoints){
		x_i = generatedsample.back();
		random_y = Normal_y_PDF(mtEngine);
		to_test = {this->callFunction(random_y)/this->callFunction(x_i) ,1};
		A = findmin(to_test);
		T = uniform_PDF(mtEngine);
		if (T < A && random_y <= m_RMax && random_y >= m_RMin){
			generatedsample.push_back(random_y);
		}
		else{
			generatedsample.push_back(x_i);
		}
	}
	return generatedsample;
}

















void FiniteFunction::printtofithist(bool isdata, std::string datafile, bool loud, std::string histpath){
	std::string mysterynumbers;
	int check;
	for (size_t i = 0; i < datafile.size(); i++){
		check = isdigit(datafile[i]);
		if (check){
			mysterynumbers += datafile[i];
		}
	}
	std::string historyfile = histpath + m_historyfile + "_fithistory_" + mysterynumbers + ".txt";
	std::filesystem::create_directory(histpath); 		//create directory (returns false if already exists)
	std::ofstream historyfilestream;
	if (!(std::filesystem::exists(historyfile))){
		std::ofstream tempfstream;
		tempfstream.open(historyfile);
		tempfstream.close();
	}
	historyfilestream.open(historyfile, std::ios::app);
	if (!historyfilestream.is_open()){	//check opened correctly
		std::cout << "Error: failed to open history file: " << m_historyfile + "_fithistory_" + mysterynumbers + ".txt" << std::endl;
	}
	else if (loud){
		std::cout << "History file: " << m_historyfile + "_fithistory_" + mysterynumbers + ".txt" << " open." <<std::endl;
	}
	//datawriting here
	double Chi2 = this->calcChi2(isdata);	//calculate chi squared
	std::string line = std::to_string(Chi2);	
	for (size_t i = 0; i < m_paramnames.size(); i++){
		line += (',' + std::to_string(m_params[i]));
	}
	historyfilestream << line << std::endl;
	historyfilestream.close();
	if (loud){
		std::cout << m_historyfile + "_fithistory_" + mysterynumbers + ".txt" << " closed." << std::endl;
	}
}

void FiniteFunction::setRangeMinMax(std::vector<double> mysterydata, bool& askedforrange, std::pair<double,double>& GlobalRange, bool askforrange, bool loud){
	if (askedforrange == true){
		m_RMin = GlobalRange.first;
		m_RMax = GlobalRange.second;
		if (loud){
			std::cout << "Plotting between " << m_RMin << " and " << m_RMax << std::endl;
		}
	}
	else{
	double minimum = findmin(mysterydata);
	double maximum = findmax(mysterydata);
	double rmin;
	double rmax;
	if (askforrange){
		std::string prompt = "Enter lower bound for plotting.";
		rmin = doublecin(prompt);
		prompt = "Enter upper bound for plotting.";
		rmax = doublecin(prompt);
		if (rmin > minimum || rmax < maximum){
			std::cout << "Data exists outside of given range. Increasing to include..." <<std::endl;
			if (rmin > minimum){rmin = floor(minimum);};
			if (rmax < maximum){rmax = ceil(maximum);};
		}			
	}
	else{
		if (loud){
			std::cout << "Smallest data point: " << minimum;
			std::cout << "; largest: " << maximum << std::endl;
		}
		rmin = floor(minimum);
		rmax = ceil(maximum);
	}
	if (loud){
		std::cout << "Plotting between " << rmin << " and " << rmax << std::endl;
	}
	m_RMin = rmin;
	m_RMax = rmax;
	askedforrange = true;
	GlobalRange = std::make_pair(rmin,rmax);
	}
}
	
void FiniteFunction::selectparams(bool& askforparams, std::string mysteryfile, bool loud, std::string path){
	std::vector<double> params;
	if (!(askforparams)){
		askforparams = this->readfithistory(params, mysteryfile, loud, path);
	}
	if (askforparams){
		std::string prompt;
		double param;
		for (std::string paramname : m_paramnames){
			prompt = "Enter a value for " + paramname + ".";
			param = doublecin(prompt);
			if ((paramname == "sigma" || paramname == "alpha" || paramname == "gamma") && !(param > 0) ){
				std::cout << paramname + " must be greater than 0. Setting to 1..." << std::endl;
				param = 1;
			}
			else if (paramname == "n" && !(param > 1)){
				std::cout << paramname + " must be greater than 1. Setting to 2..." << std::endl;
				param = 2;
			}
			params.push_back(param);
		}
	}
	m_params = params;
}
	//
	
	//Normal functions
void Normal::setparams(bool& askforparams, std::string mysteryfile, bool loud, std::string path){
	this->selectparams(askforparams, mysteryfile, loud, path);
	this->setmu(m_params[0]);
	this->setsigma(m_params[1]);
}

double Normal::Normalpdf(double x){
	double normalisation = 1/(m_sigma*std::sqrt(2*std::numbers::pi));
	double exponential = std::exp(-0.5*std::pow((x-m_mu)/m_sigma,2));
	double Normalout = normalisation * exponential;
	return Normalout;
}

double Normal::callFunction(double x){
	return this->Normalpdf(x);
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


	//Cauchy-Lorentz functions
void CauchyLorentz::setparams(bool& askforparams, std::string mysteryfile, bool loud, std::string path){
	this->selectparams(askforparams, mysteryfile, loud, path);
	this->setx0(m_params[0]);
	this->setgamma(m_params[1]);
}

double CauchyLorentz::CauchyLorentzpdf(double x){
	return 1/(std::numbers::pi*m_gamma*(1 + std::pow((x-m_x0)/m_gamma,2)));
}

double CauchyLorentz::callFunction(double x){
	return this->CauchyLorentzpdf(x);
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
m_gamma = GAMMA;}


	//Reverse Crystal Ball functions
void RevCrysBall::setparams(bool& askforparams, std::string mysteryfile, bool loud, std::string path){
	this->selectparams(askforparams, mysteryfile, loud, path);
	this->setmu(m_params[0]);
	this->setsigma(m_params[1]);
	this->setn(m_params[2]);
	this->setalpha(m_params[3]);
}

double RevCrysBall::calcA(){
return std::pow(m_n/m_alpha,m_n)*std::exp(-0.5*m_alpha*m_alpha);}
double RevCrysBall::calcB(){
return (m_n/m_alpha)-m_alpha;}
double RevCrysBall::calcC(){
return (m_n/(m_alpha*(m_n - 1)))*std::exp(-0.5*m_alpha*m_alpha);}
double RevCrysBall::calcD(){
return std::sqrt(0.5*std::numbers::pi)*(1 + std::erf(m_alpha/std::sqrt(2)));}
double RevCrysBall::calcN(){
return 1/(m_sigma*(this->calcC() + this->calcD()));}

double RevCrysBall::RevCrysBallpdf(double x){
	double RevCrysBalloutput;
	double discriminant = ((x - m_mu)/m_sigma);
	if (discriminant > - m_alpha){
		RevCrysBalloutput = this->calcN()*std::exp(-0.5*discriminant*discriminant);
	}
	else{
		RevCrysBalloutput = this->calcN()*this->calcA()*std::pow((this->calcB() - discriminant),-m_n);
	}
	return RevCrysBalloutput;
}

double RevCrysBall::callFunction(double x){
	return this->RevCrysBallpdf(x);
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
void RevCrysBall::setn(double N){
m_n = N;}
void RevCrysBall::setalpha(double ALPHA){
m_alpha = ALPHA;}


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
  m_historyfile = outfile;
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
void FiniteFunction::setOutfile(std::string Outfile){
	m_historyfile = Outfile;
	this->checkPath(Outfile);
};

/*
###################
//Getters
###################
*/ 
double FiniteFunction::rangeMin() {return m_RMin;};
double FiniteFunction::rangeMax() {return m_RMax;};

/*
###################
//Function eval
###################
*/ 
double FiniteFunction::invxsquared(double x) {return 1/(1+x*x);};
double FiniteFunction::callFunction(double x) {return this->invxsquared(x);}; //(overridable)

/*
###################
Integration by hand (output needed to normalise function when plotting)
###################
*/ 
double FiniteFunction::integrate(int Ndiv){ //private
  //ToDo write an integrator
  double stripwidth = (m_RMax - m_RMin)/Ndiv;
  double cumulative =0;
  for (int i = 0; i < Ndiv; i++){
	cumulative += stripwidth*this->callFunction(m_RMin + stripwidth*(0.5+i));
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
  std::cout << "fit parameters for " << m_historyfile << " probability distribution:" << std::endl;
  for (size_t i = 0; i < m_paramnames.size(); i++){
	  std::cout << m_paramnames[i] << ": " << m_params[i] << std::endl;
  }
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
    if (bindex<0 || bindex>Nbins){
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
