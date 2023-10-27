[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/IIHfOF4Q)
# project1-quicksort  

## Compiling the code  
A Makefile is given in the directory, simply use ``make`` to compile the code. If you're compiling the code on a M-series Mac, add the ``MXMAC=1`` option:  
```bash
make MXMAC=1  
```
## Running the code  
Enter the course environment by typing ``cs214`` in command line. Then you can run the code with:  
```bash
numactl -i all ./quicksort [num_elements] [num_rounds] [random_seed]  
```
If not specified, the default values are ``num_elements=100000000``, ``num_rounds=3`` and ``random_seed=1``.  
``numactl`` can help distribute memory pages evenly on multiple sockets. More information can be found in [Linux man page](https://linux.die.net/man/8/numactl).  

## Getting started  
Please implement your parallel quicksort in the file ``quicksort.h``. All other files will be replaced in our testing.  
For simplicity, you can implement distribution generators in ``quicksort.cpp``, but we will use our own generator when testing for fair comparison.  

## About generators  
uniform generator: https://cplusplus.com/reference/random/uniform_int_distribution/  
exponential distribution: https://cplusplus.com/reference/random/exponential_distribution/  
Zipfian distribution: https://en.wikipedia.org/wiki/Zipf%27s_law  

## Optimization hints  
+ Coarsening.  
+ Use ``std::sort`` for your base case. It's already well-optimized.  
+ Pre-allocate memory. Memory allocation can not be done in parallel.  
+ Use a hash function to replace built-in ``rand()`` since it uses system call and cannot be parallelized.  
+ Use a random element as the pivot. Further optimized it with multiple pivots.  
+ Optimize your **reduce/scan** code first because it's an important building block in this algorithm. You can even have a separate unit test for it.  

