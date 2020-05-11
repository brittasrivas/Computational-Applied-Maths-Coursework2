For Q2 I have designed a class system. There is an 'Abstract2DFunction'
class, from which I have derived the class: 'SourceFunction2' which implements the
source function, boundary function and u_exact function specified in the question.

The important class is 'AbstractPDE', from which I have derived the class:
'PoissonPDE' which implements the PDE specified in the question.
In order to solve the PDE you must pass in an 'AbstractFunction' then go through
the steps of constructing the matrix system of equations AU=F and then solving
to find the U approximation and its corresponding error norm.

'Driver.cpp' file is our main cpp file.

Compile and run all the files by using the makefile. Navigate to this directory
using the terminal and then type 'make' then './Driver'.
Use 'make clean' to restore the directory back to its original state.

Results are saved into the csv files which were then plotted using Google Sheets.
