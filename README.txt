Title

  Goldbach-scripts.cpp

Description

A C++ program which can be used to compute
data supporting the three conjectures in this paper.
Each option generates a user-specified output file 
or options to send the output to the console.

Input Options and output formats

 Option 1: Prime number generation on [a,b] from Conjecture 1.
	 Each output line consists of the pattern 
	  {p , s}:
	 	p is the prime
		s = 1 -> odd symmetry
		s = 2 -> even symmetry
   Output summary:
		# primes found
		# primes missed
		prob of a hit

 Option 2: Goldbach solutions on [a,b] from Conjecture 3.
	 Each output line consists of the patterns 
	  {e , {p , e-p}, {s_p , s_e-p}}:
	 	e is an even number greater than 8
	 	p is the prime
		s_p = 1 -> odd symmetry for prime p
		s_p = 2 -> even symmetry for prime p
		s_e-p = 1 -> odd symmetry for prime e - p
		s_e-p = 2 -> even symmetry for prime e - p 

 Option 3: n random samples on [a,b] from Conjecture 3.
	 Each output line consists of the patterns 
	  {e , {p , e-p}, {s_p , s_e-p}}:
	 	e is a sampled even number in [a, b]
	 	p is the prime
		s_p = 1 -> odd symmetry for prime p
		s_p = 2 -> even symmetry for prime p
		s_e-p = 1 -> odd symmetry for prime e - p
		s_e-p = 2 -> even symmetry for prime e - p

Dependencies

  Self contained. Compile with any C++ compiler.

Author

  Michael Mezzino - michaelmezzino@comcast.net

Version History

  01 Initial release