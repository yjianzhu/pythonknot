//  从vector<double *> 输入，先用一遍kmt算法简化，再计算Alexander多项式

#ifndef KNOT_TYPE
#define KNOT_TYPE
#include<iostream>
#include<vector>
#include<string>
#include<map>
#include"myfunction.h"
using namespace std;

vector<int> get_gauss_notation(vector<double *> &points);
string get_knottype(vector<double *> &points);

vector<double *> KMT_open_chain(vector<double *> points);
vector<double *> KMT(vector<double*> points);

string get_knottype_ring_faster(vector<double *> &points);
string get_knottype_ring(vector<double *> &points);
string get_knottype_open_faster(vector<double *> &points);

#endif