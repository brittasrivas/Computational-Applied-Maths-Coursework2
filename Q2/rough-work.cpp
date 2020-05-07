#include <iostream>
#include <vector>
#include <cmath>

int main(int argc, char* argv[])
{
  int m = 3;

  std::vector< std::vector<double> > A;

  A.resize(m, std::vector<double>(m,1));

  // So can resize and keep original elements :D - use this to make the gauss augmented matrix
  A.resize(m, std::vector<double>(m+1));

  int size = A.size();

  for(int i=0; i<size; i++)
  {
    for(int j=0; j<size+1; j++)
    {
      std::cout << A[i][j] << "\t";
    }
    std::cout << "\n";
  }


//MOD CHECK
  for(int i=0; i<20; i++)
  {
    if ((i%5) != 0) //got rid of each multiple of 5 including 0
    {
      std::cout << i << "\n";
    }
  }


  return 0;
}
