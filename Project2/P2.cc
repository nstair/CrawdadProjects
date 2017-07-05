#include <iostream>
#include <fstream>
#include <cassert>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
#include <cmath>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;


using namespace std;

int main()
//1.) reading files
{
  // open filename
  std::ifstream coord("h2o_geom.txt");
  assert(coord.good());

  int natom;
  // read the number of atoms from filename
  coord >> natom;
  cout << "the number of atoms is: " << natom << endl;

  double geom[natom][natom];
  double zval[natom];


  /* This wasn't working for me, Porque?
  for(int i=0; i < natom; i++){
    coord >> zval[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];
    cout << geom[i][0] << endl;
    cout << geom[i][2] << endl;
    cout << geom[i][3] << endl;
  }*/

  for(unsigned int i=0; i < natom; i++)
    coord >> zval[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];

  coord.close();

  std::ifstream hessn("h2o_hessian.txt");
  assert(hessn.good());

  int hs0;
  int hs1;

  hessn >> hs0;
  hessn >> hs1;

  cout << "first item in hessn: " << hs0 << endl;
  cout << "second item in hessn: " << hs1 << endl;
  cout << "The number of atoms is:\n" << natom << endl;

  /*
  double H[natom*3][natom*3];
  for(int i=0; i < 3*natom; i++){
    for(int j=0; j < natom; j++){
      hessn >> H[i][j*3] >> H[i][j*3+1] >> H[i][j*2+1];
      cout << hessn << endl;
    }
  }
  */
  hessn.close();
  /*
  //2.) finding mass weighted hessian matrix
  double amas[natom];
  amas[0] = 16.01;
  amas[1] = amas[2] = 1.01;

  double Hm[natom][natom];
  Matrix Hmm(9,9);

  for (int i=0; i < natom*3; i++){
    for (int j=0; j < natom*3; j++){
      Hm[i][j] = H[i][j]/sqrt(amas[i]*amas[j]);
      Hmm(i,j) = Hm[i][j];
    }
  }

  cout << "  " << endl;
  cout << "The mass weighted hessian is:\n" << Hmm << endl;
  cout << "  " << endl;
  for (int i=0; i < natom*3; i++){
    for (int j=0; j < natom*3; j++){
      cout << "Hm value: " << H[i][j] << "   >> i " << i << " j " << j << endl;
    }
  }

  */
  return 0;
}
