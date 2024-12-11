/*
compile as:
g++ -std=c++20 -w pi_estimation.cxx -o pi

William Strachan
started: 11/12/2024
last updated: 11/12/2024
Objectives: 
	1) Estimate pi
        - Monte Carlo integration:
        1) Understand geometry
        2) Be very bad at darts
        3) Throw a million darts
        4) Inspect the damage
        5) Replace wall (optional)
*/

#include <random>
#include <vector>
#include <iostream>


std::vector<std::pair<double,double>> generateNxy(double radius, int n_random){
    //
    std::random_device rd;
	std::mt19937 mtEngine{rd()};
	std::uniform_real_distribution<double> uniform_PDF{0,radius};
	   
    std::vector<std::pair<double,double>> xy_points;
    std::pair<double,double> xy_single;
    for (int i = 0; i < n_random; i++){
        xy_single = std::make_pair(uniform_PDF(mtEngine),uniform_PDF(mtEngine));
        xy_points.push_back(xy_single);
    }
    return xy_points;
}

bool checkincircle(std::pair<double,double> xy_single, double radius){
    if (xy_single.first*xy_single.first + xy_single.second*xy_single.second < radius*radius){
        return true;
    }
    else{
        return false;
    }
}

int numberincircle(std::vector<std::pair<double,double>> xy_points, double radius){
    int cumulative = 0;
    for (std::pair<double,double> xy_single : xy_points){
        if (checkincircle(xy_single, radius)){
            cumulative += 1;
        }
    }
    return cumulative;
}

double pi_estimate(int incircle, int total_points){
    double pi_est = 4 * (double)incircle/total_points;
    return pi_est;
}

int main(int argc, char *argv[]){
    //It shall be assumed that the user can figure out the command line arguments with a bit of trial and error!

    std::vector<std::string> inputs(argv,argv+argc);
    double radius = 1;
    int n_random = 10000;

    switch(argc){
        case 1:{
            std::cout << "No arguments given. Using 10,000 random points with radius 1." << std::endl;
            break;
        }
        case 2:{
            std::cout << "No radius given. Using " << inputs[1] << " random points with radius 1." << std::endl;
            n_random = std::stoi(inputs[1]);
            break;
        }
        case 3:{
            std::cout << "Using " << inputs[1] << " random points with radius " << inputs[2] << "." << std::endl;
            n_random = std::stoi(inputs[1]);
            radius = std::stod(inputs[2]);
            break;
        }
    }

    std::vector<std::pair<double,double>> n_random_xy;
    n_random_xy = generateNxy(radius,n_random);
    int n_incircle = numberincircle(n_random_xy,radius);
    double pi_est = pi_estimate(n_incircle, n_random);
    std::cout.setf(std::ios::fixed);
    std::cout.precision(10);
    std::cout << pi_est << std::endl << std::endl;

    return 0;
}