#pragma once
std::vector<std::string> read(std::string);
void printn(std::vector<std::string>, std::string);
std::vector<std::vector<float>> float2D(std::vector<std::string>, std::string);
float mag2D(float ,float);
std::vector<float> allmag2D(std::vector<std::vector<float>>);
void printmags(std::vector<float>);
float sumxy(std::vector<std::vector<float>>);
float sumcolumn(std::vector<std::vector<float>>, int);
float sumcolumn2(std::vector<std::vector<float>>, int);
float pcoeff(std::vector<std::vector<float>>);
float qcoeff(std::vector<std::vector<float>>);
float expval(std::vector<std::vector<float>>, float);
float chi2(std::vector<std::vector<float>>,std::vector<std::vector<float>>);
void printfit(float,float,float);
float xtothey(float, float);
std::vector<float> allxtothey(std::vector<std::vector<float>>);
void recurseiveexponent(std::vector<float>&, std::vector<std::vector<float>>, int&);
void printxtothey(std::vector<float>);