// for rand and srand
#include <stdlib.h>

#include <vector>

#include <utility>
using std::pair;

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;

/*
The parent chooser can be implemented in C or C++. The C++ way of using a class
is presented. An object of this class will be created from Python, and be used
like a regular Python class.

Using this approach, and the -shadow option of SWIG, a Python shadow class will be 
created in demoNonRandomMating.py. The C++ code is put to module
_demoNonRandomMating.so (or .pyd under windows), and will be imported into
demoNonRandomMating.py.

*/
class parentsChooser_cpp
{

public:
	// a constructor takes all locations of male and female.
	parentsChooser_cpp(const std::vector<double> & mx, const std::vector<double> & my,
		const std::vector<double> & fx, const std::vector<double> & fy,
		unsigned int seed) 
		: male_x(mx), male_y(my), female_x(fx), female_y(fy)
	{
		if (male_x.size() != male_y.size() || female_x.size() != female_y.size()) {
			cout << "Wrong male or female locations";
			// it is recommended to throw an exception here.
		}
		srand(seed);
	}
	
	// choose a male, and find a female in his vincinity.
	// Of course, some of probability density function can be used.
	// In this example, I am using the system rand() function,
	// just for demonstration. Any serious simulation should use a qualified
	// random number generator library. You may want to have a look at
	// 1. GNU Scientific Library (what simuPOP uses)
	// 2. boost C++ library
	// 3. Numeric Recipes (C or C++ version)
	pair<unsigned long, unsigned long> chooseParents()
	{
		// index of male.
		unsigned long male = rand() % male_x.size();
		
		// calculate distance between this male, and all female, and choose
		// the closest one.
		double dist;
		double mx = male_x[male];
		double my = male_y[male];
		double mindist = -1;
		unsigned long female = 0;
		for (size_t i = 0; i < female_x.size(); ++i) {
			dist = (female_x[i]-mx)*(female_x[i]-mx) + (female_y[i]-my)*(female_y[i]-my);
			if (mindist == -1 || mindist > dist) {
			    mindist = dist;
				female = i;
			}
		}
		return std::make_pair(male, female);
	}

private:
	std::vector<double> male_x;
	std::vector<double> male_y;
	std::vector<double> female_x;
	std::vector<double> female_y;
};

