#include "thomas.cpp"
#include <iostream>
#include <vector>
using namespace std;

// For brief compiling and run:
// g++ -std=c++17 -O2 -Wall -o thomas_test thomas_test.cpp && ./thomas_test

int main()
{
    vector<double> a = {-1,-1,-1}; // sub-diagonal
    vector<double> b = {2,2,2,2}; // main diagonal
    vector<double> c = {-1,-1,-1}; // super-diagonal
    vector<double> f = {1,0,0,1}; // right-hand side
    vector<double> x; // solution vector
    thomas(a,b,c,f,x);

    std::cout << "Solution: ";
    for (const auto& xi : x) {
        std::cout << xi << " ";
    }
    std::cout << std::endl;

    return 0;
}