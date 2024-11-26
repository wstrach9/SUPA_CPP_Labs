#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

std::vector<std::string> read(std::string pathfile){ // argument should be on string of a path + file name in folder
	std::ifstream input_file;
	input_file.open(pathfile);  //opens file as input stream
	if (input_file.fail()){		//checks file opened
		std::cout<<"Couldn't open file: "<< pathfile << std::endl;}
	std::string line;
	std::vector<std::string> datastr; //declare output vector and vector entries

	while (std::getline(input_file,line)){ //append all \n (or whitespace) separated lines to output vector
		datastr.push_back(line);
	}
	input_file.close(); //file fully read, can close
	return datastr;		//return output - lines of [pathfile]
}

void printn(std::vector<std::string> datastr, std::string file){
	int lendatastr = datastr.size();
	std::cout << "How many (integer) lines of " << file << " would you like to print?" << std::endl;
	int n;
	std::cin >> n;
	while (!std::cin){
		std::cin.clear();
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cout << "Not an integer. How many (integer) lines of " << file << " would you like to print?" << std::endl;
		std::cin >> n;
	}
	if (n>lendatastr){
	std::cout << "Sorry, there are only " << lendatastr << " lines in " << file << ", printing 5 lines instead.";
	n = 5;}
	std::cout << "The first " << n << " lines of " << file << " are:" << std::endl;
	for (int i = 0; i < n; i++){
		std::cout << datastr[i] << std::endl;
	}
}

std::vector<std::vector<float>> float2D(std::vector<std::string> datastr, std::string file){ //floats datastr from read(pathfile)
	std::vector<std::vector<float>> data2D;
	int i = 0;
	int headerlines;
	std::cout << "Input (integer) number of header lines in " << file << std::endl;
	std::cin >> headerlines; //user input: number of headerlines from plain text printout after read(pathfile)
	while (!std::cin){
		std::cin.clear();
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cout << "Not an integer. Input (integer) number of header lines in " << file << std::endl;
		std::cin >> headerlines;
	}
	for (std::string line : datastr){ // for each (','-separated) entry of datastr
		if (i < headerlines){ 
			i++;} // ignore headerlines
		else{
			std::stringstream linestream(line); 		// create a string stream for each line of datastr
			std::vector<float> datapoint;
			while (linestream.good()){					// repeat while there is string remaining
				std::string entry;				
				std::getline(linestream,entry,','); 	// entry = string up to a ','
				datapoint.push_back(std::stof(entry));}	// append entry to a vector. Once each for {x,y} entries.
			data2D.push_back(datapoint);				// append {x,y} float vector to 2D vector
		}				
	}
	return data2D;
}

float mag2D(float a,float b){
	return std::sqrt(std::pow(a,2)+std::pow(a,2));
}

std::vector<float> allmag2D(std::vector<std::vector<float>> data2D){
	std::vector<float> Allmag2D;
	for (std::vector<float> line : data2D){
		Allmag2D.push_back(mag2D(line[0],line[1]));}
	return Allmag2D;
}

void printmags(std::vector<float> Allmag2D){
	std::cout << "Print magnitudes of all x-y datapoints? (y/n/some)" <<std::endl;
	std::string yn;
	std::cin >> yn;
	while (yn != "y" && yn != "n" && yn != "some" && yn != "Y" && yn != "N" && yn != "Some"){
		std::cout << "Invalid response. Print magnitudes of all x-y datapoints? (y/n/some)" << std::endl;
		std::cin >> yn;
	}
	if (yn != "n" && yn != "N"){
		int n;
		if (yn == "y" || yn == "Y"){
			n = Allmag2D.size();
			std::cout << "The magnitudes of all the x-y datapoints are:" << std::endl;}
		else{
			std::cout << "How many magnitudes to print? (integer)" <<std::endl;
			std::cin >> n;
			while (!std::cin){
				std::cin.clear();
				std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				std::cout << "Not an integer. How many magnitudes to print? (integer)" <<std::endl;
				std::cin >> n;
			}
			std::cout << "The magnitudes of the first " << n << " x-y datapoints are:" << std::endl;
		}
		for (int i = 0; i<n; i++){
			std::cout << Allmag2D[i] << std::endl;
		}
	}
}

float sumxy(std::vector<std::vector<float>> xyvec){
	float Sumxy = 0;
	for (std::vector<float> xy : xyvec){
		Sumxy += xy[0]*xy[1];
	}
	return Sumxy;
}

float sumcolumn(std::vector<std::vector<float>> xyvec, int n){
	float Sumcolumn = 0;
	for (std::vector<float> xy : xyvec){
		Sumcolumn += xy[n];
	}
	return Sumcolumn;
}

float sumcolumn2(std::vector<std::vector<float>> xyvec, int n){
	float Sumcolumn2 = 0;
	for (std::vector<float> xy : xyvec){
		Sumcolumn2 += std::pow(xy[n],2);
	}
	return Sumcolumn2;
}

float pcoeff(std::vector<std::vector<float>> data2D){
	float N = data2D.size();
	float numerator = N*sumxy(data2D) - sumcolumn(data2D,0)*sumcolumn(data2D,1);
	float denominator = N*sumcolumn2(data2D,0)-std::pow(sumcolumn(data2D,0),2);
	return numerator/denominator;	
}

float qcoeff(std::vector<std::vector<float>> data2D){
	float N = data2D.size();
	float numerator = sumcolumn2(data2D,0)*sumcolumn(data2D,1) - sumxy(data2D)*sumcolumn(data2D,0);
	float denominator = N*sumcolumn2(data2D,0) - std::pow(sumcolumn(data2D,0),2);
	return numerator/denominator;
}

float expval(std::vector<std::vector<float>> data2D, float x){
	float expy = pcoeff(data2D)*x + qcoeff(data2D);
	return expy;
}

float chi2(std::vector<std::vector<float>> data2D, std::vector<std::vector<float>> error2D){
	int lendat = data2D.size();
	int lenerr = error2D.size();
	if (lendat != lenerr){ //escape sequence  if more/less data than errors
		std::cout << "Data and Error lists not the same the length. Could not compute chi squared test, defaulting to 0.0" << std::endl;
		return 0; //float 0 hopefully...
	}
	float Chi2 = 0;
	for (int i=0; i<lendat; i++){
		/*float component = std::pow(data2D[i][1] - expval(data2D,data2D[i][0]),2)/(std::pow(error2D[i][0],2)+std::pow(error2D[i][1],2)); 
		//assuming expected errors in x and y added in quadrature
		Chi2 += component; */
		Chi2 += std::pow(data2D[i][1] - expval(data2D,data2D[i][0]),2)/(std::pow(error2D[i][0],2)+std::pow(error2D[i][1],2)) ;
	}
	return Chi2;
}

void printfit(float p,float q,float Chi2){
	std::cout << "Print Least Squares fit coefficients? (y/n)" << std::endl;
	std::string yn;
	std::cin >> yn;
	while (yn != "y" && yn != "n" && yn != "Y" && yn != "N"){
		std::cout << "Invalid response. Print Least Squares fit coefficients? (y/n)" << std::endl;
		std::cin >> yn;
	}
	if (yn == "y" || yn == "Y"){
		std::cout << "The x-y data is fitted to the line y = px + q , with coefficients:" <<std::endl;
		std::cout << "p = " << p << std::endl;
		std::cout << "q = " << q << std::endl;
	}
	std::cout << "Print Chi Squared test value? (y/n)" << std::endl;
	std::cin >> yn;
	while (yn != "y" && yn != "n"&& yn != "Y" && yn != "N"){
		std::cout << "Invalid response. Print Chi Squared test value? (y/n)" << std::endl;
		std::cin >> yn;
	}
	if (yn == "y" || yn == "Y"){
		std::cout << "The x-y data is correlated with a Chi Squared value of:" << std::endl;
		std::cout << Chi2 << std::endl;
	}
}

float xtothey(float x, float y){ //for 0<y<5.4999..
	int inty = std::round(y);
	float Xtothey;
	switch(inty){
	case 0:{
		Xtothey = 1;
		break;}
	case 1:{
		Xtothey = x;
		break;}
	case 2:{
		Xtothey = x*x;
		break;}
	case 3:{
		Xtothey = x*x*x;
		break;}
	case 4:{
		Xtothey = x*x*x*x;
		break;}
	case 5:{
		Xtothey = x*x*x*x*x;
		break;}
	}
	return Xtothey;
}

void recurseiveexponent(std::vector<float>& XtotheY, std::vector<std::vector<float>> data2D, int& index){ //recursively raises x to the y and appends to XtotheY vector
	if (XtotheY.size() < data2D.size()){						//as long as we haven't done x^y for all datapoints
		float Xtothey = xtothey(data2D[index][0],data2D[index][1]); //individual x^y
		XtotheY.push_back(Xtothey);									//adds individual x^y to a vector of x^y
		index++;													//next datapoint
		recurseiveexponent(XtotheY,data2D,index);					//repeat
	}
} // a glorified while statement...

std::vector<float> allxtothey(std::vector<std::vector<float>> data2D){ //calls recursive function with required initial inputs
	int index = 0;
	std::vector<float> XtotheY;
	recurseiveexponent(XtotheY, data2D, index);
	return XtotheY;
} //

void printxtothey(std::vector<float> XtotheY){
	std::cout << "Print x^y (y rounded to integer) for all x-y datapoints? (y/n/some)" <<std::endl;
	std::string yn;
	std::cin >> yn;
	while (yn != "y" && yn != "n" && yn != "some" && yn != "Y" && yn != "N" && yn != "Some"){
		std::cout << "Invalid response. Print x^y (y rounded to integer) for all x-y datapoints? (y/n/some)" << std::endl;
		std::cin >> yn;
	}
    int n;
	if (yn != "n" && yn != "N"){
		if (yn == "y" || yn == "Y"){
			n = XtotheY.size();
			std::cout << "The values x^y (y rounded to integer) for all x-y datapoints are:" << std::endl;
		}
		else{
			std::cout << "How many values x^y to print? (integer)" <<std::endl;
			std::cin >> n;
			while (!std::cin){
				std::cin.clear();
				std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				std::cout << "Not an integer. How many values x^y to print? (integer)" <<std::endl;
				std::cin >> n;
			}
			std::cout << "The values x^y (y rounded to integer) for the first " << n << " x-y datapoints are:" << std::endl;
		}
		for (int i = 0; i<n; i++){
		std::cout << XtotheY[i] << std::endl;
		}
	}
}