//Compile with g++ -o boost_installed boost_installed.cpp
//g++ -g -I/home/ibsimu/boost_1_73_0 
//  -L/home/ibsimu/boost_1_73_0/stage/lib -o temp_option_groups temp_option_groups.cpp
// -H to check for included libraries

// g++ -lboost_program_options -o boost_installed boost_installed.cpp

#include <iostream>
#include <boost/array.hpp>

using namespace std;
int main(){
  boost::array<int, 4> arr = {{1,2,3,4}};
  cout << "If boost is installed, it prints 1: " << arr[0];
  return 0;
}