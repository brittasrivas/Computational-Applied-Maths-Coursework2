#include <iostream>
#include <stdlib.h> //used for rand()
#include <vector>


void fillMatrix(std::vector< std::vector<double> > &M)
{
  srand((unsigned)time(NULL));

  // M[0][0] = (double)(rand()) / (double)(RAND_MAX + 1.0);
  // M[0][1] = (double)(rand()) / (double)(RAND_MAX + 1.0);
  // M[1][0] = (double)(rand()) / (double)(RAND_MAX + 1.0);
  // M[1][1] = (double)(rand()) / (double)(RAND_MAX + 1.0);

  for(int i=0; i<2; i++)
  {
    for(int j=0; j<2; j++)
    {
      M[i][j] = (double)(rand()) / (double)(RAND_MAX + 1.0);
    }
  }

}


int main(int argc, char* argv[])
{

  //double M[2][2];

  std::vector< std::vector<double> > M(2, std::vector<double>(2));

  fillMatrix(M);

  for (int i=0; i<2; i++)
  {
    for (int j=0; j<2; j++)
    {
      std::cout << M[i][j] << "\n";
    }
  }

  std::cout << "Size of matrix: " << M.size() << "\n";

  int number = 1e6;
  std::cout << number << "\n";


  return 0;
}
