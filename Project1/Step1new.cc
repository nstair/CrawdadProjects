#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include "molecule.h"

using namespace std;

int main(int argc, char *argv[])
{
  Molecule mol("geom.dat",0);

  cout << "the number of atoms is:" << mol.natom << endl;
  cout << "the atomic number and goemetry of mol is:\n";
  mol.print_geom();

  return 0;
}
