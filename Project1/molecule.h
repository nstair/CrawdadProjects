#include <string>

using namespace std;

//define class iteslf
class Molecule
{
  public:
    int natom;
    int charge;
    int *zvals;
    double **geom;
    double **ex;
    double **ey;
    double **ez;
    double mass_sum;
    double x_com;
    double y_com;
    double z_com;
    string point_group;

    void print_geom();
    void rotate(double phi);
    void translate(double x, double y, double z);
    double bond(int atom1, int atom2);
    double uv();
    double com();
    double oop(int a, int b, int c, int d);
    double tor(int a, int b, int c, int d);
    double angle(int atom1, int atom2, int atom3);
    double torsion(int atom1, int atom2, int atom3, int atom4);

    Molecule(const std::string filename, int q);
    ~Molecule();
};
