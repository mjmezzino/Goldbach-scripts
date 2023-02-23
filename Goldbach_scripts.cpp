// A C++ program which can be used to compute
// data supporting the three conjectures in this paper.
// Each option generates a user-specified output file
// or options to send the output to the concole.
//
// Input Options and output formats
//
// Option 1: Prime number generation on [a,b]
//	 Each output line consists of the pattern 
//	  {p,s}:
//	 	p is the prime
//		s = 1 -> odd symmetry
//		s = 2 -> even symmetry
//   Output summary:
//		# primes found
//		# primes missed
//		prob of a hit
//
// Option 2: Goldbach solutions on [a,b]
//	 Each output line consists of the patterns 
//	  {e,{p,e-p},{s_p,s_e-p}}:
//	 	e is an even number greater than 8
//	 	p is the prime
//		s_p = 1 -> odd symmetry for prime p
//		s_p = 2 -> even symmetry for prime p
//		s_e-p = 1 -> odd symmetry for prime e - p
//		s_e-p = 2 -> even symmetry for prime e - p
//
// Option 3: n random samples on [a,b]
//	 Each output line consists of the patterns 
//	  {e,{p,e-p},{s_p,s_e-p}}:
//	 	e is a sampled even number in [a,b]
//	 	p is the prime
//		s_p = 1 -> odd symmetry for prime p
//		s_p = 2 -> even symmetry for prime p
//		s_e-p = 1 -> odd symmetry for prime e - p
//		s_e-p = 2 -> even symmetry for prime e - p
//
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iterator>
#include <vector>
#include <ctime>
#include <Windows.h>
#include <Lmcons.h>

using namespace std;

ofstream fout;
ifstream fin;

string filename = "~test.txt";
double filtered_count = 0, total_count = 0;

using namespace std;

// Check for prime number
bool is_prime(long long n) {
	// Assumes that n is a positive natural number
	// We know 1 is not a prime number
	if (n == 1) {
		return false;
	}

	long long i = 2;
	// This will loop from 2 to int(sqrt(x))
	while (i * i <= n) {
		// Check if i divides x without leaving a remainder
		if (n % i == 0) {
			// This means that n has a factor in between 2 and sqrt(n)
			// So it is not a prime number
			return false;
		}
		i += 1;
	}
	// If we did not find any factor in the above loop,
	// then n is a prime number
	return true;
}

// Compute GCD
long long findGCD(long long a, long b) {
	if (b == 0)
		return a;
	return findGCD(b, a % b);
}

// Separate numerator and denominator
long long* reduce(long long num, long long den)
{
	long long fract[2]{};

	long long gcd = findGCD(num, den);
	num /= gcd;
	den /= gcd;
	fract[0] = num;
	fract[1] = den;
	//	std::cout << "Lowest Fraction : " << num << "/" << den;
	return fract;
}

// Compute Extended GCD
long long xGCD(long long a, long long b, long long& x, long long& y) {
	if (b == 0) {
		x = 1;
		y = 0;
		return a;
	}

	long long x1, y1, gcd = xGCD(b, a % b, x1, y1);
	x = y1;
	y = x1 - (a / b) * y1;
	return gcd;
}

//Compute elliptic curve group from y^2 = x^3 + b x - b 
// where b = 2p - 3 starting with the tangent at (1,1)
vector<long> group(long long p)
{
	long long* fract;
	vector <long> coords = {};
	vector <vector<long>>  points = {};
	vector <long> result = {};

	bool trace = false;
	bool turn = false;
	int loop = 0;

	long long x, y, f, g, x2, y2, x3, y3, num, den;

	int index = 0;
	int x1 = 1;
	int y1 = 1;

	if (trace) std::cout << "x1= " << x1 << endl;
	if (trace) std::cout << "y1= " << y1 << endl;

	coords.clear();
	coords.push_back(x1);
	coords.push_back(y1);
	points.push_back(coords);

// Compute 2P
	xGCD(1, p, f, y);
	if (f < 0) f = p + f;
	if (trace) std::cout << "f= " << f << endl;
	g = (long long)(p * f) % p;
	if (trace) std::cout << "g= " << g << endl;
	x2 = (long long)(pow(g, 2) - 2 * x1) % p;
	if (x2 < 0) x2 = p + x2;
	if (trace) std::cout << "x2= " << x2 << endl;
	y2 = (long long)(g * (x1 - x2) - y1) % p;
	if (y2 < 0) y2 = p + y2;
	if (trace) std::cout << "y2= " << y2 << endl;
	coords.clear();
	coords.push_back(x2);
	coords.push_back(y2);
	points.push_back(coords);

	if (x2 == 1 && y2 + 1 == p)
	{
		result.push_back(1);
		result.push_back(2);
		result.push_back(1);
		result.push_back(1);
		return result;
	}

// Compute kP, 2 < k 
	long long i;
	for (i = 1; i <= 2 * p; i++)
	{
		fract = reduce((y2 - y1), (x2 - x1));
		num = fract[0];
		den = fract[1];
		if (den < 0) den = p + den;
		xGCD(den, p, f, y);
		if (trace) std::cout << "#1= " << den << endl;
		if (f < 0) f = p + f;
		if (trace) std::cout << "f= " << f << endl;
		g = (long long)(num * f) % p;
		if (trace) std::cout << "g= " << g << endl;
		x3 = (long long)(pow(g, 2) - x2 - x1) % p;
		if (x3 < 0) x3 = p + x3;
		if (trace) std::cout << "x" << i + 2 << " = " << x3 << endl;
		y3 = (long long)(g * (x1 - x3) - y1) % p;
		if (y3 < 0) y3 = p + y3;
		if (trace) std::cout << "y" << i + 2 << " = " << y3 << endl;

		//  Test for even symmetric group
		if (!turn
			&& (points[i].at(0) == points[i - 1].at(0))
			&& (points[i].at(1)) + points[i - 1].at(1) == p)
		{
			if (trace)
			{
				std::cout << i << " - " << points[i].at(0) << " - " << points[i - 1].at(0) << endl;
				std::cout << i << " - " << (points[i].at(1) + points[i - 1].at(1)) << endl;
			}
			i++;
			loop = i;
			break;
		}

		//  Test for odd symmetric group
		if (turn
			&& (points[i].at(0) == points[i - 2].at(0))
			&& (points[i].at(1)) + points[i - 2].at(1) == p)
		{
			if (trace)
			{
				std::cout << i << " - " << points[i].at(0) << " - " << points[i - 2].at(0) << endl;
				std::cout << i << " - " << (points[i].at(1) + points[i - 2].at(1)) << endl;
			}
			loop = i;
			break;
		}

		x2 = x3;
		y2 = y3;
		coords.clear();
		coords.push_back(x3);
		coords.push_back(y3);
		points.push_back(coords);

		if (x3 == 1 && y3 + 1 == p)
		{
			loop = i;
			break;
		}

		// kP off the curve y^2 = x^3 + b x - b 
		if ((f * den % p) != 1)
		{
			result.push_back(0);
			result.push_back(0);
			result.push_back(0);
			result.push_back(i + 2);
			points.clear();
			return result;
		}

		// Potential odd symmetric group 
		if (y3 == 0) turn = true;
	}

	if (!loop)
	{
		result.push_back(0);
		result.push_back(0);
	}
	if (!turn) // Even group order
	{
		result.push_back(2);
		result.push_back(2 * (i - 1));
	}
	else //Odd group order
	{
		result.push_back(1);
		result.push_back(2 * i - 1);
	}

	result.push_back(turn);
	result.push_back(loop);

	points.clear();
	return result;
}

// Find Goldbach prime pair
long long goldbach(long long e)
{
	bool trace = false;
	double l2;
	long long* fract;
	long long rtn = 0, us, p;

	for (long long b = 3; b <= e; b += 4)
	{
		total_count++;
		p = (3 + b) / 2;

		// Test from D-set analysis added
		if (e - 2 * p > 0 && findGCD(e, e - 2 * p) != 2)
		{
			filtered_count++;
			continue;
		}

		us = 3 * p / 2 + 2 * sqrt(p);
		vector <long> g = group(p);
		int symmetric1 = g.at(0);
		l2 = pow(g.at(1), 2);
		if (trace)
			std::cout << "Group p: " << p << " : " << symmetric1 << " - " << us << " < " << l2 << endl;
		// Conjectur 3 test 
		if (symmetric1 && us < l2) 
		{
			if (is_prime(p) == false)
			{
				std::cout << "\n\n***** Error: p = " << p << " is not prime.\n" << endl;
				rtn = 0;
				break;
			}

			std::cout << e << ":            \r";
			std::cout << e << ": " << p << "\r";

			us = 3 * (e - p) / 2 + 2 * sqrt(e - p);
			g = group(e - p);
			int symmetric2 = g.at(0);
			l2 = pow(g.at(1), 2);
			if (trace)
				std::cout << "Group e-p: "<< e-p << " : " << symmetric2 << " - " << us << " < " << l2 << endl;
			if (symmetric2 && us < l2)
			{
				if (is_prime(e - p))
				{
				  if (fout.is_open())
					fout << "{" << e << " , {" << p << "," << e - p << "} , {" << symmetric1 << "," << symmetric2 << "}}" << endl;
				  else if (toupper(filename[0]) == 'C')
					cout << "{" << e << " , {" << p << "," << e - p << "} , {" << symmetric1 << "," << symmetric2 << "}}" << endl;
				  rtn = e;
				  break;
				}
				else
				{
					std::cout << "\n\n***** Error: e - p = " << e - p << " is not prime.\n" << endl;
					rtn = 0;
					break;
				}
			}
		}
	}
	return rtn;
}

double diffclock(long clock1, long clock2)
{
	double diffticks = clock1 - clock2;
	double diffms = (diffticks) / (CLOCKS_PER_SEC / 1000);
	return diffms;
}

void format(double t, char endline)
{
	long millis = t;

	long seconds = millis / 1000;
	millis %= 1000;

	long minutes = seconds / 60;
	seconds %= 60;

	long hours = minutes / 60;
	minutes %= 60;

	long days = hours / 24;
	hours %= 24;

	if (days > 1)
		std::cout << days << " Days : ";
	else if (days == 1)
		std::cout << days << " Day : ";
	if (hours > 1)
		std::cout << hours << " Hours : ";
	else if (hours == 1)
		std::cout << hours << " Hour : ";
	if (minutes > 1)
		std::cout << minutes << " Minutes : ";
	else if (minutes == 1)
		std::cout << minutes << " Minute : ";
	if (seconds > 1)
		std::cout << seconds << " Seconds : ";
	else if (seconds == 1)
		std::cout << seconds << " Second : ";
	if (millis > 1)
		std::cout << millis << " MilliSeconds";
	else if (millis == 1)
		std::cout << millis << " MilliSecond";
	if (endline) std::cout << endline;
}

string username()
{
	char buffer[UNLEN + 1]{};
	DWORD len = UNLEN + 1;
	if (GetUserNameA(buffer, &len)) return buffer;
	else return "";
}

int main()
{
	string  report = "1";
	double hit = 0, miss = 0, l2;
	long long us = 0, count = 0, rtn , start = 0, end = 0, lower = 0, upper = 0, p;

	std::cout << "\nSelect option (1 = Primes, 2 = Goldbach solution, 3 = Random sample): ";
		std::cin >> report;
	if (report == "") report = "1";

	std::cout << "\nEnter a fully qualified file name or " << endl;
	std::cout << "'c' for console or 's' for summary only" << endl;
	std::cout << "  [e.g., c:\\\\Users\\\\username\\\\Desktop\\\\filename.txt]" << endl;
	std::cout << "  [  or, ~filename creates filename on the Desktop]" << endl << endl << ":";

	std::cin >> filename;

	if (toupper(filename[0]) == 'C')
		cout << "Console output." << endl;

	size_t found = filename.find("~");
	if (found == 0)
		filename = "C:\\Users\\" + username() + "\\Desktop\\" + filename.substr(found + 1);

	if (report == "3") // Random range
	{
		std::cout << "\nEnter lower and upper limits: ";
		std::cin >> lower;
		std::cin >> upper;
		if (lower > upper)
		{
			std::cout << "\nLower limit > Upper limit" << endl;
			std::cin.get();
			return 0;
		}
		std::cout << "\nEnter number of samples: ";
		std::cin >> count;
	}
	else
	{
		std::cout << "\nEnter start and end: ";
		std::cin >> start;
		std::cin >> end;
		if (start > end)
		{
			std::cout << "\nStart > End" << endl;
			std::cin.get();
			return 0;
		}
	}
	std::cout << endl;

	long start_time = clock();

	if (found == 0)
		fout.open(filename);

	if (report == "1") // Primes
	{
		if (start % 2 == 0) start++;
		start = max(3, start);

		for (p = start; p <= end; p += 2)
		{
			if (p > 5 && p % 5 == 0) continue;

			std::cout << p << "\r";

			vector <long> g = group(p);
			int symmetric = g.at(0);
			if (symmetric)
			{
				us = 3 * p / 2 + 2 * sqrt(p);
				l2 = pow(g.at(1), 2);
				if (symmetric && us < l2)
				{
					if (is_prime(p))
					{
						hit++;
						if (fout.is_open())
							fout << "{" << p << " , " << g.at(0) << "}" << endl; 
						else if (toupper(filename[0]) == 'C')
							cout << "{" << p << " , " << g.at(0) << "}" << endl;
					}
					else
					{
						std::cout << "\n\n***** Error: " << p << " is not prime.\n" << endl;
						break;
					}
				}
				if (is_prime(p) && l2 <= us)
				{
					miss++;
					if (toupper(filename[0]) != 'S')
						std::cout << "Missed prime = " << p << endl;
				}

			}
		}
		std::cout << "            " << endl;
		std::cout << "Found: " << hit << endl;
		std::cout << "Missed: " << miss << endl;
		std::cout << "Hit Rate: " << 1.0 - miss/(hit + miss) << endl << endl;
	}
	else if (report == "2") // Goldbach
	{
		if (start % 2 == 1) start++;

		for (long long e = start; e <= end; e += 2)
		{
			std::cout << e << "\r";

			if (e != goldbach(e))
				std::cout << "***** Error ***** - No solution: " << e << endl;
		}
		std::cout << "Filter efficiency = " << (int) (100*filtered_count / total_count) << "%" << endl << endl;

	}
	else // Random sample in range
	{
		srand(time(0));
		while (count-- > 0)
		{
			long long e = ((double)rand() / (RAND_MAX)) * (upper - lower) + lower;
			if (e % 2 == 1) e++;

			std::cout << e << "\r";

			if (e != goldbach(e))
			{
				std::cout << "***** Error ***** - No solution: " << e << endl;
				break;
			}
		}
		std::cout << "Filter efficiency = " << (int)(100 * filtered_count / total_count) << "%" << endl << endl;
	}

	if (fout.is_open())
		fout.close();

	long end_time = clock();
	format(diffclock(end_time, start_time), '\n');

	std::cin.get();

	return 0;
}
