//  从vector<double *> 输入，先用一遍kmt算法简化，再计算Alexander多项式

#ifndef KNOT_TYPE
#define KNOT_TYPE
#include<iostream>
#include<vector>
#include<string>
#include<map>
#include"myfunction.h"
using namespace std;

string get_knottype(vector<double *> &points);
void get_alxeander_map(fstream &read);

void print_alexander_map();

vector<double *> KMT_open_chain(vector<double *> points);

string get_knottype_ring_faster(vector<double *> &points);
string get_knottype_ring(vector<double *> &points);
string get_knottype_open_faster(vector<double *> &points);

#endif