
#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

using namespace std;

int main()

{
  Molecule mol("geom.dat",0);

  //1.) reading the .dat file and reprinting the values using a class fn
  cout << "Number of atoms:" << mol.natom << endl;
  cout << "Input Cartesian coordinates:\n";
  mol.print_geom();

  // 2.) finding bond distances
  /*? is there a preference beteen the below
  and allocating memory vie douubel rij[mol.natom][mol.natom],
  in terms of defining a 'pointer to double' typer object?   */

  double **rij;
  rij = new double* [mol.natom];
  for(int i=0; i < mol.natom; i++)
    rij[i] = new double[mol.natom];

  cout << "the bond distance is:" << endl;
  for(int i=0; i < mol.natom; i++){
    for(int j=0; j < mol.natom; j++){
      rij[i][j] = mol.bond(i,j);
      cout << " i= " << i << " j= " << j << " bond length: " << mol.bond(i,j) << endl;
    }
  }
  // when should I delete this, or is that even nesicary in the using std lib?

  cout << "  " << endl;

  // 3.) finding bond angles
  double ***anijk = new double** [mol.natom];
  for(int i = 0; i < mol.natom; i++){
    anijk[i] = new double*[mol.natom];
    for(int j = 0; j < mol.natom; j++){
      anijk[i][j] = new double[mol.natom];
    }
  }

  mol.uv();
  //cout << "test val for uv matrix:" << mol.ex[4][5] << endl;
  for(int i=0; i < mol.natom; i++){
    for(int j=0; j < i; j++){
      for(int k=0; k < j; k++){
        anijk[i][j][k] = mol.angle(i,j,k);
        if (i!=j && j!=k && i!=k){
          if(mol.angle(i,j,k)*(180/acos(-1.0)) != 90.0 && mol.angle(i,j,k)*(180/acos(-1.0)) != 0.0)
            cout << "i= " << i;
            cout << "  j= " << j;
            cout << "  k= " << k << endl;
            cout << " the bond angle is:" << mol.angle(i,j,k)*(180/acos(-1.0)) << endl;
        };
        //cout << "i=:"<< i << endl;
        //cout << "j=:"<< j << endl;
        //cout << "k=:"<< k << endl;
      }
    }
  }

  cout << "  " << endl;

  // 4. & 5.) finding out of plane angles and torsion angles
  for(int i=0; i < mol.natom; i++){
    for(int j=0; j < mol.natom; j++){
      for(int k=0; k < mol.natom; k++){
        for(int l=0; l < j; l++){
          if(i!=j && i!=k && i!=l && j!=k && j!=l && k!=l){
            if(rij[i][k] < 4.0 && rij[k][j] < 4.0 && rij[k][l] < 4.0){
              cout << "i= " << i;
              cout << "  j= " << j;
              cout << "  k= " << k;
              cout << "  l= " << l  << endl;
              cout << " the out of plane angle is:" << (180/acos(-1.0))*mol.oop(i,j,k,l) << endl;
              /* there may be more 'interesing' torsion
                angles but I threw them in
                this loop just for fun! */
              cout << " the dihedral angle is:" << (180/acos(-1.0))*mol.tor(i,j,k,l) << endl;
              cout << "  " << endl;
            }
          }
        }
      }
    }
  }

  cout << "  " << endl;

  // 6.) center of mass translation
  mol.com();
  // testing odd values
  cout << "mass_sum:" << mol.mass_sum << endl;
  cout << "the center of mass in X is:" << mol.x_com << endl;
  cout << "the center of mass in Y is:" << mol.y_com << endl;
  cout << "the center of mass in Z is:" << mol.z_com << endl;

  // translating c.o.m.
  mol.translate(-mol.x_com, mol.y_com, -mol.z_com);
  mol.com();
  cout << "the com after translation in X is:" << mol.x_com << endl;
  cout << "the com after translation in Y is:" << mol.y_com << endl;
  cout << "the com after translation in Z is:" << mol.z_com << endl;

  cout << "  " << endl;

  // 7.) principle moments of inertia
  // Jeff & Sam said not to allocate memory, so I'm no longer going to do that :)

  double I[3][3] = {{0.0}};

  for(int i=0; i < mol.natom; i++){
    I[0][0] += 2.0*mol.zvals[i]*(mol.geom[i][1]+mol.geom[i][2])*(mol.geom[i][1]+mol.geom[i][2]);
    I[1][1] += 2.0*mol.zvals[i]*(mol.geom[i][0]+mol.geom[i][2])*(mol.geom[i][0]+mol.geom[i][2]);
    I[2][2] += 2.0*mol.zvals[i]*(mol.geom[i][0]+mol.geom[i][1])*(mol.geom[i][0]+mol.geom[i][1]);

    I[1][0] = I[0][1] += 2.0*mol.zvals[i]*mol.geom[i][1]*mol.geom[i][0];
    I[2][0] = I[0][2] += 2.0*mol.zvals[i]*mol.geom[i][2]*mol.geom[i][0];
    I[2][1] = I[1][2] += 2.0*mol.zvals[i]*mol.geom[i][2]*mol.geom[i][1];
  };

  cout << "Ixx: " << I[0][0] << " Ixy: " << I[0][1] << " Ixz: " << I[0][2] << endl;
  cout << "Iyx: " << I[1][0] << " Iyy: " << I[1][1] << " Iyz: " << I[1][2] << endl;
  cout << "Izx: " << I[2][0] << " Izy: " << I[2][1] << " Izz: " << I[2][2] << endl;

  cout << "  " << endl;

  // canned program for finding eiganvalues of this system?

  double Ia = 51.5092;
  double Ib = 211.263;
  double Ic = 349.616;

  double m_brad = 5.2918e-11;
  double kg_amu = 1.6605e-27;

  cout << "principel moments of inertia (amu*br^2):\n";
  cout << "  " << Ia << "  " << Ib << "  " << Ic << endl;
  cout << "  " << endl;
  cout << "principel moments of inertia (kg*m^2):\n";
  cout << "  " << Ia * m_brad * m_brad * kg_amu;
  cout << "  " << Ib * m_brad * m_brad * kg_amu;
  cout << "  " << Ic * m_brad * m_brad * kg_amu << endl;

  cout << "  " << endl;

  Ia *= m_brad * m_brad * kg_amu;
  Ib *= m_brad * m_brad * kg_amu;
  Ic *= m_brad * m_brad * kg_amu;

  double h = 6.626070e-34;
  double c = 3.000e+8;
  double A; double B; double C;

  A = h/(100.0*8.0*acos(-1.0)*acos(-1.0)*c*Ia);
  B = h/(100.0*8.0*acos(-1.0)*acos(-1.0)*c*Ib);
  C = h/(100.0*8.0*acos(-1.0)*acos(-1.0)*c*Ic);

  cout << "The rontational constants are (cm^-1):\n";
  cout << "  " << A << "  " << B << "  " << C << endl;
  cout << "  " << endl;

  A *= 3.0e+4;
  B *= 3.0e+4;
  C *= 3.0e+4;

  cout << "The rontational constants are (Mhz):\n";
  cout << "  " << A << "  " << B << "  " << C << endl;
  cout << "  " << endl;


  // 8.) finding rotational constants (big conversion value problem...)


  /*(I belive this will be done by the deconstructor (but does it need to be called?)
  //~Molecule mol()?
  delete[] zval;
  delete[] x;  delete[] y;  delete[] z; */


  for(int i=0; i < mol.natom; i++)
    delete[] rij[i];
    //delete unit vec matricies
  delete[] rij;

  return 0;
}
