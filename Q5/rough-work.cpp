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
    for(int j=0; j<3; j++)
    {
      //M[i][j] = (double)(rand()) / (double)(RAND_MAX + 1.0);
      M[i][j] = i+1;
    }
  }

}


double* matvec_mult(const std::vector< std::vector<double> > M, const double* v)
/* Calculates Matrix-Vector multiplication: M*v
Matrix M is size mxn. */
{
  int m = M.size();
  int n = M[0].size();
  double *ans;
  ans = new double[m];

  for(int i=0; i<m; i++)
  {
    double row_sum = 0.0;
    for (int j=0; j<n; j++)
    {
      row_sum += M[i][j] * v[j];
    }
    ans[i] = row_sum;
  }

  return ans;
}



int main(int argc, char* argv[])
{
  int number = 1e6;
  std::cout << number << "\n";

  //double M[2][2];

  std::vector< std::vector<double> > M(2, std::vector<double>(3));

  fillMatrix(M);

  for (int i=0; i<2; i++)
  {
    for (int j=0; j<3; j++)
    {
      std::cout << M[i][j] << "\t";
    }
    std::cout << "\n";
  }

  std::cout << "No of rows of matrix: " << M.size() << "\n";
  std::cout << "No of cols of matrix: " << M[0].size() << "\n";



  double *vector;
  vector = new double[3];
  vector[0] = vector[1] = vector[2] = 1.0;
  //std::cout << "Size of vector: " << vector.size() << "\n";
  //CAN NOT GET SIZE OF A POINTER ARRAY

  double *ans;
  ans = new double[2];

  ans = matvec_mult(M, vector);

  std::cout << "Mat vec mult ans: \n";

  for (int i=0; i<2; i++)
  {
    std::cout << ans[i] << "\n";
  }

  delete[] vector;
  delete[] ans;
  return 0;
}
