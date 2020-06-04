//g++ -o boost_options boost_options.cpp -lboost_program_options

/*
command line:

    ./boost_options --include-path abcd --include-value 3 --include-path bcde --include-value 5


output:
    0 ip: abcd v: 31 ip: bcde v: 5

*/



// https://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <string>
#include <iostream>

namespace bpo = boost::program_options;

using namespace std;


int main(int argc, char *argv[]) {

    int opt;
    bpo::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("optimization", bpo::value<int>(&opt)->default_value(10), 
    "optimization level")
        ("include-path,I", bpo::value< vector<string> >(), 
    "include path")
        ("include-value,V", bpo::value< vector<int> >(), 
    "include value")

        ("input-file", bpo::value< vector<string> >(), "input file")
    ;
    bpo::positional_options_description p;
    p.add("input-file", -1);

    bpo::variables_map vm;
    bpo::store(bpo::command_line_parser(argc, argv).
          options(desc).positional(p).run(), vm);
    bpo::notify(vm);
    const vector<string> &ip = vm["include-path"].as<vector<string>>();
    const vector<int> &iv = vm["include-value"].as<vector<int>>();
    const std::vector<string>::size_type include_size = ip.size(); 
    for(std::vector<string>::size_type i=0; i < include_size; i++) {
        cout<<i<<" ip: "<<ip[i]<<" v: "<<iv[i];
    }

}
