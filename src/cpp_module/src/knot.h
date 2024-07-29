#ifndef KNOT
#define KNOT

#include "myfunction.h"
#include "knottype.h"

class knot
{
private:
    /* data */
public:
    //vector<double*> Simplify_points;
    vector<double*> points3D;
    vector<int> ends;
    //vector<vector<double>> interSection_Matrix;
    //vector<vector<int>> interSection_Matrix_up_down;

    knot(vector<double*> &a,vector<int> b/* args */);
    knot(vector<double*> &a,int b/* args */);
    knot(knot *source)
    {
        points3D=source->points3D;
        ends=source->ends;
    }
    knot(vector<double*> &a)
    {
        points3D=a;
    }
    
    void print_knot();
    ~knot();


    void KMT();

    vector<int> knot_size(string &return_knot_type,string &knottype_desired);
    vector<int> knot_size_debug(string &return_knot_type,string &knottype_desired);
    vector<int> knot_size_ring(string &return_knot_type,string &knottype_desired);
    vector<int> knot_size_sonknot(string &mother_knot_type,const string &knottype_desired);
    vector<int> knot_size_sonknot2(string &mother_knot_type,const string &knottype_desired);
    //void get_interSection_Matrix();

    void print_IM();

    void find_max_span(vector<double *> & result);

    //void get_gauss_notation();
};

#endif