//FiniteFunctions.h backup

#include <string>
#include <vector>
#include "gnuplot-iostream.h"

#pragma once //Replacement for IFNDEF

//global functions

std::vector<double> readmysterydata(std::string filepath);
bool yncin(std::string prompt);
int intcin(std::string prompt);
double doublecin(std::string prompt);
std::pair<int, double> findargmin(std::vector<double> vector);
double findmin(std::vector<double> vector);
double findmax(std::vector<double> vector);



class FiniteFunction{

public:
  //added
  void setRangeMinMax(std::vector<double>, bool askforrange = false, bool loud = false);
  void printtofithist(std::vector<double>, int, std::string, std::string path = "");
  double calcChi2(std::vector<double>, int);
  //original
  FiniteFunction(); //Empty constructor
  FiniteFunction(double range_min, double range_max, std::string outfile); //Variable constructor
  ~FiniteFunction(); //Destructor
  double rangeMin(); //Low end of the range the function is defined within
  double rangeMax(); //High end of the range the function is defined within
  double integral(int Ndiv = 1000); 
  std::vector< std::pair<double,double> > scanFunction(int Nscan = 1000); //Scan over function to plot it (slight hack needed to plot function in gnuplot)
  void setRangeMin(double RMin);
  void setRangeMax(double RMax);
  void setOutfile(std::string outfile);
  void plotFunction(); //Plot the function using scanFunction
  
  //Plot the supplied data points (either provided data or points sampled from function) as a histogram using NBins
  void plotData(std::vector<double> &points, int NBins, bool isdata=true); //NB! use isdata flag to pick between data and sampled distributions
  virtual void printInfo(); //Dump parameter info about the current function (Overridable)
  virtual double callFunction(double x); //Call the function with value x (Overridable)

  //Protected members can be accessed by child classes but not users
protected:
  //added:
  bool readfithistory(std::vector<double>&, std::string, bool loud = false, std::string path = "");
  void selectparams(bool&, std::string, bool loud = false, std::string path = "");
  std::vector<double> m_params;
  std::vector<std::string> m_paramnames;
  //original:
  double m_RMin;
  double m_RMax;
  double m_Integral;
  int m_IntDiv = 0; //Number of division for performing integral
  std::string m_FunctionName;
  std::string m_OutData; //Output filename for data
  std::string m_OutPng; //Output filename for plot
  std::string m_historyfile;
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
  
  
private:
  double invxsquared(double x); //The default functional form
};

class Normal : public FiniteFunction {
public:
	Normal() : FiniteFunction(){m_paramnames = {"mu","sigma"};};
	Normal(double range_min, double range_max, std::string outfile) : FiniteFunction(range_min,range_max,outfile){};
	void setparams(bool&, std::string, bool loud = false, std::string path = ""); //sets mu and sigma according to best previous fit (bool f)or user input (t)
	double mu();
	double sigma();
	void setmu(double MU);
	void setsigma(double SIGMA);
	virtual double callFunction(double x);
private:
	double m_mu;
	double m_sigma;
	double Normalpdf(double x);
};

class CauchyLorentz : public FiniteFunction {
public:
	CauchyLorentz() : FiniteFunction(){m_paramnames = {"x0","gamma"};};
	CauchyLorentz(double range_min, double range_max, std::string outfile) : FiniteFunction(range_min,range_max,outfile){};
	void setparams(bool&, std::string, bool loud = false, std::string path = "");
	double x0();
	double gamma();
	void setx0(double X0);
	void setgamma(double GAMMA);
	virtual double callFunction(double x);
private:
	double m_x0;
	double m_gamma;
	double CauchyLorentzpdf(double x);
};

class RevCrysBall : public FiniteFunction { 
public:
	RevCrysBall() : FiniteFunction(){m_paramnames = {"mu","sigma", "n", "alpha"};};
	RevCrysBall(double range_min, double range_max, std::string outfile) : FiniteFunction(range_min,range_max,outfile){};
	void setparams(bool&, std::string, bool loud = false, std::string path = "");
	double mu();
	double sigma();
	double n();
	double alpha();
	void setmu(double MU);
	void setsigma(double SIGMA);
	void setn(double N);
	void setalpha(double ALPHA);
	virtual double callFunction(double x);
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
};