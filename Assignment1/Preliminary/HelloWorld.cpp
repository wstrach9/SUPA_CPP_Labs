// Preliminary Exercises

#include <iostream>
#include <vector>
#include <math.h>

float magnitude(); //prototype (function defined after main) 

int main(){

std::cout << "Hello World!" << std::endl;


float x = 2.3;
float y = 4.5;
std::vector<float> x_y2D;   //empty vector
x_y2D.push_back(x);         //adds x to end of vector
x_y2D.push_back(y);         //adds y to end
float magnitude_xy;
magnitude_xy = sqrt(pow(x_y2D[0],2)+pow(x_y2D[1],2));
std::cout << "The magnitude of vector 2D vector is " << magnitude_xy << std::endl;

float magnitude_xyfn = magnitude(); //magnitude calculated by magnitude()
std::cout << "The magnitude of vector 2D vector is " << magnitude_xyfn << std::endl;
return 0;
}

float magnitude(){
//function to calculate magnitude of 2D vector (x,y)

    std::string user_input;
    std::cout << "Input value of vector x component:" << std::endl;
    std::cin >> user_input;
    float x = std::stof(user_input); //string->float
    //std::cout<<x<<std::endl;
    std::cout << "Input value of y component:" << std::endl;
    std::cin >> user_input;
    float y = std::stof(user_input);
    //std::cout<<y<<std::endl;
    return sqrt(pow(x,2)+pow(y,2));
    }