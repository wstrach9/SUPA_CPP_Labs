#include <string>
#include <vector>
#include "/workspaces/SUPA_CPP_Labs/GNUplot/gnuplot-iostream.h"

#pragma once //Replacement for IFNDEF

std::vector<double> readmysterydata(std::string);
int intcin(std::string prompt);
bool yncin(std::string prompt);
double doublecin(std::string prompt);
std::vector<std::pair<double, double>> getguesses(bool askforguesses, int n_params, std::vector<std::string> paramnames, double mysterymean);
std::pair<int, double> findmin(std::vector<double>);
std::vector<double> parametersweep(std::string, bool, double&, double&, int&);
std::pair<double,double> RangeMinMax(bool, double rmin = -5, double rmax = 5);
//Functions to optimise fit for Reverse Crystal Ball
void checkzerox(double& guess, double range, std::string paramname);
void changeguess(double& guess, double& range, double best, std::string paramname);
std::vector<std::vector<std::pair<double,double>>> TourDeSequence(std::vector<std::string>, std::vector<std::pair<double, double>>&, int, std::vector<double>, int, bool);
std::vector<std::pair<double,double>> cycleoptimise(int, std::string, std::pair<double,double>, std::vector<std::pair<double, double>>&, std::vector<std::string> , int, std::vector<double>);
std::pair<double,double> singleoptimise(int, std::pair<double,double> GlobalRange, std::string, std::vector<std::vector<double>>, std::vector<double>, int);
int fourcycle(int&, bool&);
std::vector<double> parameterspace(std::vector<std::pair<double, double>>, int, std::vector<std::string>, int);

class FiniteFunction{

public:
  FiniteFunction(); //Empty constructor
  FiniteFunction(double range_min, double range_max, std::string outfile); //Variable constructor
  ~FiniteFunction(); //Destructor
  double rangeMin(); //Low end of the range the function is defined within
  double rangeMax(); //High end of the range the function is defined within
  double chi2(); //quality of data fit
  double integral(int Ndiv = 1000); 
  std::vector< std::pair<double,double> > scanFunction(int Nscan = 1000); //Scan over function to plot it (slight hack needed to plot function in gnuplot)
  void setRangeMin(double RMin);
  void setRangeMax(double RMax);
  void setOutfile(std::string outfile);
  void setchi2(double CHI2);
  void plotFunction(); //Plot the function using scanFunction
  double CalcChi2(std::vector<double> &points, int);
  //Plot the supplied data points (either provided data or points sampled from function) as a histogram using NBins
  void plotData(std::vector<double> &points, int NBins, bool isdata=true); //NB! use isdata flag to pick between data and sampled distributions
  virtual void printInfo(); //Dump parameter info about the current function (Overridable)
  virtual double callFunction(double x); //Call the function with value x (Overridable)
  
  
  //Protected members can be accessed by child classes but not users
protected:
  double m_RMin;
  double m_RMax;
  double m_Integral;
  int m_IntDiv = 0; //Number of division for performing integral
  std::string m_FunctionName;
  std::string m_OutData; //Output filename for data
  std::string m_OutPng; //Output filename for plot
  std::vector< std::pair<double,double> > m_data; //input data points to plot
  std::vector< std::pair<double,double> > m_samples; //Holder for randomly sampled data 
  std::vector< std::pair<double,double> > m_function_scan; //holder for data from scanFunction (slight hack needed to plot function in gnuplot)
  bool m_plotfunction = false; //Flag to determine whether to plot function
  bool m_plotdatapoints = false; //Flag to determine whether to plot input data
  bool m_plotsamplepoints = false; //Flag to determine whether to plot sampled data 
  double integrate(int Ndiv);
  std::vector< std::pair<double, double> > makeHist(std::vector<double> &points, int Nbins); //Helper function to turn data points into histogram with Nbins
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp);
  double m_Chi2;
  
private:
  double invxsquared(double x); //The default functional form
};

/* parent class FiniteFunction. Contains:
- constructor FiniteFunction w/ arguments (and individual get/setters):
	- rangeMin = lower bound of domain
	- rangeMax = upper bound of domain
	- outfile  = filename for outputs of printing
- public functions:
	- scanFunction(Nscan = 1000) scan function to plot
	- plotFunction() (using scanFunction)
	- plotData(&points, NBins, bool isdata) (isdata flag to choose data or sampled distributions))
	- virtual printInfo() (return parameter information for current function)
	- vitual callFunction(x) (calculate function for value x)
- protected variables:
	- m_RMin, m_RMax protected range variables
	- m_Integral
	- m_IntDiv (division for calculating integral)
	- m_FunctionName
	- m_OutData (filename for output data)
	- m_OutPng (filename for output plot)
	- m_data (input data)
	- m_samples (holder for random sample data)
	- m_function_scan (holder for scanFunction)
	- bool m_plotfunction (?)
	- bool m_plotdatapoints (?)
	- bool m_plotsamplepoints (?)
	- integrate(Ndiv) (indegrates function)
	- makeHist(&points, Nbins) - convert data to histogram with Nbins bins
	- checkPath(outstring) - check data and png paths are correct
	- generatePlot(Gnuplot &gp) - call Gnuplot to create plot
private function:
	- invxsquared(x) - basic function 1/(1+x**2) */


class Normal: public FiniteFunction{
public:
	//constructors
	Normal(); //empty
	Normal(double range_min, double range_max, std::string outfile) : FiniteFunction(range_min,range_max,outfile){}; 
	//overload parent (?)
	Normal(double range_min, double range_max, double mu, double sigma, std::string outfile); // new arguments
	//getters for new arguments
	double mu();
	double sigma();
	//setters
	void setmu(double MU);
	void setsigma(double SIGMA);
private:
	double m_mu;
	double m_sigma;
	double Normalpdf(double x);
	double callFunction(double x);
};

class CauchyLorentz: public FiniteFunction{
public:
	//constructors
	CauchyLorentz();//empty
	CauchyLorentz(double range_min, double range_max, std::string outfile) : FiniteFunction(range_min,range_max,outfile){}; 
	//overload parent (?)
	CauchyLorentz(double range_min, double range_max, double x0, double gamma,std::string outfile); // new arguments
	//getters for new arguments
	double x0();
	double gamma();
	//setters
	void setx0(double);
	void setgamma(double);
private:
	double m_x0;
	double m_gamma;
	double CauchyLorentzpdf(double x);
	double callFunction(double x);
};

class RevCrysBall: public FiniteFunction{
public:
	//constructors
	RevCrysBall(); //empty
	RevCrysBall(double range_min, double range_max, std::string outfile) : FiniteFunction(range_min,range_max,outfile){}; 
	//overload parent (?)
	RevCrysBall(double range_min, double range_max, double mu, double sigma, double n, double alpha, std::string outfile); //new arguments
	//getters for new arguments
	double mu();
	double sigma();
	double n();
	double alpha();
	//setters
	void setmu(double MU);
	void setsigma(double SIGMA);
	void setn(double N);
	void setalpha(double ALPHA);
private:
	double m_mu;
	double m_sigma;
	double m_n;
	double m_alpha;
	double calcA();
	double calcB();
	double calcC();
	double calcD();
	double calcN();
	double RevCrysBallpdf(double x);
	double callFunction(double x);
};