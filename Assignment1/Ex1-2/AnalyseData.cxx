
// compile with:
// g++ -std=c++20 -w /workspaces/SUPA_CPP_Labs/Assignment1/Ex1-2/AnalyseData.cxx /workspaces/SUPA_CPP_Labs/Assignment1/Ex1-2/myFunctions.cxx -o Analyse

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include "/workspaces/SUPA_CPP_Labs/Assignment1/Ex1-2/myFunctions.h"

int main(){
	/* Problem 1
	Print n lines of input2d_float.txt as plain text to terminal
	*/
	std::string path = "/workspaces/SUPA_CPP_Labs/Exercises2024/Ex1_2/";
	std::string file = "input2D_float.txt";
	std::string pathfile = path + file;			//set the file to read
	
	std::vector<std::string> datastr;
	datastr = read(pathfile);
	
	printn(datastr, file, path);	//prints first n (user input) lines to terminal
	
	/*Problem 2
	Calculate magnitudes of each datapoint, ask user if they want to print to terminal
	*/
	std::vector<std::vector<float>> data2D;
	data2D = float2D(datastr, file);				//float the data
	
	std::vector<float> magnitudes;
	magnitudes = allmag2D(data2D);			// calculates magnitudes with entry indices correlating to data2D
	
	printmags(magnitudes, path); //prints first n or all (user input) magnitudes to terminal
	
	/*Problem 3
	Calculate least squares fit line coefficients, ask to print to terminal
	Calculate chi-squared test, ask to print to terminal
	*/
	
	// p and q coefficients
	float p = pcoeff(data2D);
	float q = qcoeff(data2D);
	
	//Chi-2 test
	file = "error2D_float.txt";
	pathfile = path + file;
	
	std::vector<std::string> errorstr;
	errorstr = read(pathfile);
	
	std::vector<std::vector<float>> error2D;
	error2D = float2D(errorstr, file); //need to input number of headerlines again.
	
	float Chi2;
	Chi2 = chi2(data2D, error2D);
	
	printfit(p,q,Chi2, path);
	
	/* Problem 4
	Calculate x^int(y) without using std::pow() or a for/while loop
	*/
	
	std::vector<float> XtotheY;
	XtotheY = allxtothey(data2D);
	printxtothey(XtotheY, path);
    return 0;
}