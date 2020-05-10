#include <iostream>
#include <stdlib.h> //used for rand()


void fillMatrix(double M[2][2])
{
  M[0][0] = rand() / (RAND_MAX + 1.0);
  M[0][1] = rand() / (RAND_MAX + 1.0);
  M[1][0] = rand() / (RAND_MAX + 1.0);
  M[1][1] = rand() / (RAND_MAX + 1.0);

}


int main(int argc, char* argv[])
{

  double M[2][2];

  fillMatrix(M);

  for (int i=0; i<2; i++)
  {
    for (int j=0; j<2; j++)
    {
      std::cout << M[i][j] << "\n";
    }
  }

  return 0;
}
