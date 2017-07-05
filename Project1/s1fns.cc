#include "molecule.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <cassert>
#include <cmath>


//define print_geom function
void Molecule::print_geom()
{
  for(int i=0; i < natom; i++)
    printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
}

//define translate funcion (won't need for step 1)
void Molecule::translate(double x, double y, double z)
{
  for(int i=0; i < natom; i++) {
     geom[i][0] += x;
     geom[i][1] += y;
     geom[i][2] += z;
  }
}

double Molecule::bond(int atom1, int atom2)
{
  //allocate space here if not in main
  //rij = new double* [natom];
  //for(int i=0; i < natom; i++)
    //rij[i] = new double[natom];

  //for(int i=0; i < atom1; i++) {
    //for(int j=0; i < atom2; j++) {
      return sqrt(
        (geom[atom1][0]-geom[atom2][0])*(geom[atom1][0]-geom[atom2][0])
        + (geom[atom1][1]-geom[atom2][1])*(geom[atom1][1]-geom[atom2][1])
        + (geom[atom1][2]-geom[atom2][2])*(geom[atom1][2]-geom[atom2][2])
                        );
    //}
  //}

}

double Molecule::uv()
{
  ex = new double* [natom];
  ey = new double* [natom];
  ez = new double* [natom];
  for(int i = 0; i < natom; i++){
    ex[i] = new double[natom];
    ey[i] = new double[natom];
    ez[i] = new double[natom];
  }
  for(int i=0; i < natom; i++){
    for(int j=0; j < natom; j++){
      if( j == i ) {
        ex[i][j] = ex[j][i] = 0.0;
        ey[i][j] = ey[j][i] = 0.0;
        ez[i][j] = ez[j][i] = 0.0;
      } else {
        ex[i][j] = ex[j][i] = -(geom[i][0] - geom[j][0])/bond(i,j);
        ey[i][j] = ey[j][i] = -(geom[i][1] - geom[j][1])/bond(i,j);
        ez[i][j] = ez[j][i] = -(geom[i][2] - geom[j][2])/bond(i,j);
      }
    }
  }
  return 0;
}

double Molecule::angle(int atom1, int atom2, int atom3)
{
  double val;
  val = acos(
    ex[atom1][atom2]*ex[atom2][atom3]
    + ey[atom1][atom2]*ey[atom2][atom3]
    + ez[atom1][atom2]*ez[atom2][atom3]
                );
    // still unsure as to why one of the vectors is flipping :(
  if (val > 0.0005) {
    return (acos(-1.0) - val);
  } else {
    return 0.0;
  }
}

double Molecule::oop(int a, int b, int c, int d)
{
  double ebcd_x = ey[c][b] * ez[c][d] - ez[c][b] * ey[c][d];
  double ebcd_y = ez[c][b] * ex[c][d] - ex[c][b] * ez[c][d];
  double ebcd_z = ex[c][b] * ey[c][d] - ey[c][b] * ex[c][d];

  double exx = ebcd_x * ex[c][a];
  double eyy = ebcd_y * ey[c][a];
  double ezz = ebcd_z * ez[c][a];

  double theta = (exx + eyy + ezz)/sin(angle(b,c,d));

  if(theta < -1.0) theta = asin(-1.0);
  else if(theta > 1.0) theta = asin(1.0);
  else theta = asin(theta);

return theta;
}

double Molecule::tor(int a, int b, int c, int d)
{
  double eabc_x = ey[b][a] * ez[b][c] - ez[b][a] * ey[b][c];
  double eabc_y = ez[b][a] * ex[b][c] - ex[b][a] * ez[b][c];
  double eabc_z = ex[b][a] * ey[b][c] - ey[b][a] * ex[b][c];

  double ebcd_x = ey[c][b] * ez[c][d] - ez[c][b] * ey[c][d];
  double ebcd_y = ez[c][b] * ex[c][d] - ex[c][b] * ez[c][d];
  double ebcd_z = ex[c][b] * ey[c][d] - ey[c][b] * ex[c][d];

  double e_xx = eabc_x * ebcd_x;
  double e_yy = eabc_y * ebcd_y;
  double e_zz = eabc_z * ebcd_z;

  double theta = (e_xx + e_yy + e_zz)/sin(angle(b,c,d))/sin(angle(a,b,c));

  if(theta < -1.0) theta = asin(-1.0);
  else if(theta > 1.0) theta = asin(1.0);
  else theta = asin(theta);

  return theta;
}

// funciton for finding com
double Molecule::com()
{
  mass_sum = 0.0;
  x_com = 0.0;
  y_com = 0.0;
  z_com = 0.0;

  for(int i=0; i < natom; i++)
    mass_sum += zvals[i];

  for(int i=0; i < natom; i++){
    x_com += (geom[i][0]*zvals[i])/mass_sum;
    y_com += (geom[i][1]*zvals[i])/mass_sum;
    z_com += (geom[i][2]*zvals[i])/mass_sum;
  }
  /* ? what retun stmt do I use if all I want the
  funciton to do is create valued variables for
  later use? */
  return 0;
}



//define class constructor function
Molecule::Molecule(std::string filename, int q)
{
  charge = q;

  // open filename
  std::ifstream input(filename);
  assert(input.good());

  // read the number of atoms from filename
  input >> natom;

  // allocate space
  zvals = new int[natom];
  // for geom, it's making [natom] dimensional dp colum vecors OF 3 dim row vectors!
  geom = new double* [natom];
  for(int i=0; i < natom; i++)
    geom[i] = new double[3];

  for(unsigned int i=0; i < natom; i++)
    input >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];

  input.close();
}

//define class deconstructor funcion (for removal of memory once used)
Molecule::~Molecule()
{
  delete[] zvals;
  for(int i=0; i < natom; i++)
    delete[] geom[i]; //delete[] ex[i]; delete[] ey[i]; delete[] ez[i];
  delete[] geom; //delete ex; delete ey; delete ez;
  //I think I will also need to delete the bond distance pointers?
}
