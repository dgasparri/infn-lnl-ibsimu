#include <iostream>
#include <stdlib.h>

int main() {
    std::string a = "-2.9119e-16";
    double b = strtod(a.c_str(), nullptr);
    std::cout<<"Originale: "<<a<<" valore: "<<b<<std::endl;
    return 0;

}