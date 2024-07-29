#ifndef MY_FUNCTION
#define MY_FUNCTION

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <map>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cmath>

#define my_epsilon 0.0000001
using namespace std;

// TODO
void find_max_span(vector<double *> &result);

void recenter(vector<double *> &a,int Lknot);

void rotate_array(vector<double *> &x,int n);
void rotate_array_vector(vector<double *> &x,int n);

void write_map(map<int,int> &hist_count,fstream &write);

int write_core_x_tail(vector<int> &knotSize,vector<double *> &point,int x_tail);
int write_core_tail(vector<int> &knotSize,vector<double *> &point,int n_tail);
int write_core(vector<int> &knotSize,vector<double *> &point);//把相同纽结大小的构象写入同一个文件

int get_segment(vector<int> &segment,int offset,int point);
bool pair_compare(pair<int,double>a,pair<int,double>b);

int judge_triangle(double *a,double *b,double *c,double *plain,double *line_1,double *line_2);

double cal_interSection(double *line1_1,double *line1_2,double *line2_1,double *line2_2,int *up_down);

vector<double> cross(const vector<double> a,const vector<double> b);

void cross_product(double a[3],double b[3], double res[3]);

void cal_normals(double I[3], double J[3], double K[3], double planeijk[4]);

void read_data(vector<double *> & x,char *s);

int read_data_lammps(vector<double *> &x,fstream &read);

int read_data_pdb(vector<double *> &x,fstream &read);

int read_data_cpp(vector<double *> &x,fstream &read);

void write_data_cpp(vector<double *> &x,string s);

int read_data_cpp_oxdna(vector<double *> &x,fstream &read,int NB);

void clean_pointer(vector<double *> &x);
void print_knot(vector<double *> &points3D);

void error_out(vector<double *> &x,string s="error.txt");

#endif