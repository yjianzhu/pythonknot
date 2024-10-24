#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include "homfly.h"
#include <getopt.h>
#include <chrono>
#include <tuple>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

/* comment
2021/11/24
*/
using namespace std;
using namespace Eigen;
#define epsilon 0.0000001
int cross_num=0;

int flag_KMT =0; // 0 for no KMT, 1 for KMT

bool pair_compare(pair<int,double>a,pair<int,double>b)
{
    return a.second<b.second;
}
int judge_triangle(double *a,double *b,double *c,double *plain,double *line_1,double *line_2)
{
    // a b c 三角形三个顶点。 line_1 and line_2 是向量
    // 输出1 表示有重叠
    double dis1=line_1[0]*plain[0]+line_1[1]*plain[1]+line_1[2]*plain[2]+plain[3];
    double dis2=line_2[0]*plain[0]+line_2[1]*plain[1]+line_2[2]*plain[2]+plain[3];
    if( dis1*dis2>epsilon)
        return 0;
    //上面三行代码判断是否为平面同侧。
    if(fabs(dis1)<epsilon and fabs(dis2)<epsilon  )
    {
        //两个点都在面上，暂时先不给用
        return 1;
    }
    Vector3d T,d,E1,E2,M,K;
    T<<line_1[0]-a[0],line_1[1]-a[1],line_1[2]-a[2];
    d<<line_2[0]-line_1[0],line_2[1]-line_1[1],line_2[2]-line_1[2];
   
    E1<<b[0]-a[0],b[1]-a[1],b[2]-a[2];
    E2<<c[0]-a[0],c[1]-a[1],c[2]-a[2];
    M=d.cross(E2);
    double det=M.dot(E1);
    K=T.cross(E1);
    double t=K.dot(E2)/det;
    double u=M.dot(T)/det;
    double v=K.dot(d)/det;
    
    if(u<=0 or v<=0 or u+v>=1 or t>=1 or t<=0)
        return 0;
    //cout<<det<<"\tuvt\t"<<u<<'\t'<<v<<'\t'<<t<<endl;
    return 1;    
}

double cal_interSection(double *line1_1,double *line1_2,double *line2_1,double *line2_2,int *up_down)
{
    Vector2d k,b;
    b<<line2_1[0]-line1_1[0],line2_1[1]-line1_1[1];
    Matrix2d x;
    x<<line1_2[0]-line1_1[0],line2_1[0]-line2_2[0],line1_2[1]-line1_1[1],line2_1[1]-line2_2[1];
    if(x.determinant()==0)  return 0;//平行
    k=x.inverse()*b;
    if(k(0,0)<0+epsilon or k(0,0)>1-epsilon or k(1,0)<0+epsilon or k(1,0)>1-epsilon)    return 0;//用一个小量解决可能出现的共点情况

    if((k(0,0)*(line1_2[2]-line1_1[2])+line1_1[2])>(k(1,0)*(line2_2[2]-line2_1[2])+line2_1[2]))
        {//第一个片段在上方
            *up_down=1;
            if(x(0,0)*(-x(1,1))+x(1,0)*x(0,1)>0)
                return k(0,0);
            else
                return -1*k(0,0);
        }
    else
        {//第二个片段在下方
            *up_down=-1;
            if(x(0,0)*(-x(1,1))+x(1,0)*x(0,1)>0)
                return -1*(k(0,0));
            else
                return k(0,0);
        }
}
vector<double> cross(const vector<double> a,const vector<double> b)
{
    if(a.size()!=3 or b.size()!=3) printf("error, vector product should be two vectors with 3 elements.\n");
    vector<double> result{0,0,0};

    result[0]=a[1]*b[2]-a[2]*b[1];
    result[1]=-a[0]*b[2]+a[2]*b[0];
    result[2]=a[0]*b[1]-a[1]*b[0];

    return result;
}
void cross_product_hull(double a[3],double b[3], double res[3]) {
   res[0] = a[1]*b[2] - a[2]*b[1];
   res[1] = a[2]*b[0] - a[0]*b[2];
   res[2] = a[0]*b[1] - a[1]*b[0];
}
void cal_normals(double I[3], double J[3], double K[3], double planeijk[4]) {
   int d;
   double  vij[3], vik[3], vtmp[3], vabs;

   for(d=0;d<3;d++) vij[d] = J[d] - I[d];
   for(d=0;d<3;d++) vik[d] = K[d] - I[d];
   cross_product_hull(vij, vik, vtmp);
   vabs = sqrt( vtmp[0]*vtmp[0] + vtmp[1]*vtmp[1] + vtmp[2]*vtmp[2] );
   if(fabs(vabs)<epsilon) return;
   for(d=0;d<3;d++) planeijk[d] = vtmp[d]/vabs;
   planeijk[3]  = 0.;
   for(d=0;d<3;d++)  planeijk[3] += I[d]*planeijk[d];
   planeijk[3]  = -planeijk[3];

}
// class definition
class knot
{
private:
    /* data */
public:
    vector<double*> points3D;
    vector<int> ends;
    vector<vector<double>> interSection_Matrix;
    vector<vector<int>> interSection_Matrix_up_down;
    knot(vector<double*> a,vector<int> b/* args */);
    knot(vector<double*> a)
    {
        points3D=a;
    }
    void print_knot(){
        cout<<points3D.size()<<endl;
        for(int i=0;i<points3D.size();i++)
        {
            printf("%10.5lf\t%10.5lf\t%10.5lf\t\n",points3D[i][0],points3D[i][1],points3D[i][2]);
        }
        if(ends.size()==0)  printf("knot without ends\n");
        else
        {
            for(int i=0;i<ends.size();i++)
                printf("%d\t",ends[i]);
            printf("\n");
        }
    }
    ~knot();


    void KMT()
    {
        while (true)
        {
            int number=points3D.size();
            int flag;
            for(int i=1;i<=points3D.size();i++)
            {
                //cout<<i<<'\t'<<points3D.size()<<endl;
                double plain[4]={0};
                flag=0;
                cal_normals(points3D[i-1],points3D[i%points3D.size()],points3D[(i+1)%points3D.size()],plain);

                if(fabs(plain[0])<epsilon and fabs(plain[1])<epsilon and fabs(plain[2])<epsilon)
                    {//如果三个点在一条线上，可以省去
                    points3D.erase(points3D.begin()+i%points3D.size());
                    i--;
                    //cout<<i<<endl;
                    //cout<<points3D.size()<<endl;
                    continue;
                    }

                for(int j=i+1;j<=points3D.size()+i-2;j++)
                    {
                        if(judge_triangle(points3D[i-1],points3D[i%points3D.size()],points3D[(i+1)%points3D.size()],plain,points3D[j%points3D.size()],points3D[(j+1)%points3D.size()]))
                            {
                                flag=1;
                                break;
                            }
                    }
                if(flag==0)
                {
                    points3D.erase(points3D.begin()+i%points3D.size());
                    i--;
                    //cout<<i<<endl;
                    //cout<<points3D.size()<<endl;
                }

            }
            if(number==points3D.size())
                break;
        }   
    }


    void get_interSection_Matrix()
    {
        int n=points3D.size();
        int up_down;
        vector<vector<double>> i_Matrix(n,vector<double> (n,0));
        vector<vector<int>> i_Matrix_up_down(n,vector<int> (n,0));
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {   
                if(i==j)    continue;
                up_down=0;
                i_Matrix[i][j]=cal_interSection(points3D[i],points3D[(i+1)%n],points3D[j],points3D[(j+1)%n],&up_down);
                i_Matrix_up_down[i][j]=up_down;
            }
        }
        interSection_Matrix=i_Matrix;
        interSection_Matrix_up_down=i_Matrix_up_down;
    }

    void print_IM()
    {
        printf("\n Intersection Matrix:\n");
        int n=points3D.size();
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                printf("%5.10f\t",interSection_Matrix[i][j]);
            }
            printf("\n");
        }

        printf("\n Intersection Matrix up or down:\n");
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                printf("%5d\t",interSection_Matrix_up_down[i][j]);
            }
            printf("\n");
        }
    }

    std::tuple<int, std::string> get_gauss_notation()
    {
        string notation;
        notation+=to_string(ends.size()+1)+' ';
        map<pair<int,int>,pair<int,int>> cross;
        int n=points3D.size();
        int count=0;

        for(int i=0;i<n;i++)
        {
            for(int j=i+1;j<n;j++)
            {
                if(fabs(interSection_Matrix[i][j])>epsilon)
                    {
                        int isOver0=interSection_Matrix[i][j]>0.?1:-1;//大于0的记为1，小于0-1
                        cross.insert(map<pair<int,int>,pair<int,int>>::value_type (pair<int,int>(i,j), pair<int,int>(count,isOver0)));
                        count++;
                    }
            }
        }
        // for(auto it=cross.begin();it!=cross.end();it++)
        // {
        //     printf("%d\t%d\n",it->second.first,it->second.second);
        // }
        notation+=to_string(cross.size()*2)+' ';
        cross_num+=cross.size();
        for(int i=0;i<n;i++)
        {   
            vector<pair<int,double>> oneLine;
            for(int j=0;j<n;j++)
                {
                    if(fabs(interSection_Matrix[i][j])>epsilon)
                        oneLine.push_back(pair<int,double>(j,fabs(interSection_Matrix[i][j])));
                }
            sort(oneLine.begin(),oneLine.end(),pair_compare);
            for(auto it=oneLine.begin();it!=oneLine.end();it++)
            {   
                pair<int,int> temp;
                if(i<it->first) {temp.first=i;temp.second=it->first;}
                else {temp.first=it->first;temp.second=i;}
                notation+=to_string(cross[temp].first)+' ';
                notation+=to_string(interSection_Matrix_up_down[i][it->first])+' ';
            }
        }
        for(auto i:cross)
        {
            notation+=to_string(i.second.first)+' '+to_string(i.second.second)+' ';
        }
        //cout<<notation<<endl;

        char * writable = new char[notation.size() + 1];
        copy(notation.begin(), notation.end(), writable);
        writable[notation.size()] = '\0';
        auto result = homfly_str(writable);

        delete[] writable;

        // result 是char *, 需要转换为 std::string
        return std::make_tuple(cross.size(),std::string(result));
    }
};


knot::knot(vector<double*> a,vector<int> b/* args */)
{
    points3D=a;
    ends=b;
}

knot::~knot()
{
    //printf("class has been killed!\n");
}




void move3D(vector<double *> a,vector<int> b)
{//不能用
    int N_components=b.size()+1;
    b.push_back(a.size());
    vector<double *> single_knot;
    single_knot.assign(a.begin(),a.begin()+b[0]);
    // knot temp(single_knot,b);
    // temp.print_knot();


    //接口留给links，目前只用单一闭合链。
    for(int component=1;component<N_components;component++)
    {
        vector<double *> single_knot;
        printf("%d\n",b[component]);
        single_knot.assign(a.begin()+b[component-1],a.begin()+b[component]);
        knot temp(single_knot);
        temp.print_knot();
    }


}

int read_data_cpp(vector<double *> &x,fstream &read)
{
    double *temp_point;
    int i_molecule=0,NB;
    string aline;
    read>>NB;
    if(read.fail()) return 0;
    getline(read,aline);
    getline(read,aline);

    for(int i=0;i<NB;i++)
    {   
        temp_point=new double[3];
        read>>i_molecule>>temp_point[0]>>temp_point[1]>>temp_point[2];
        x.push_back(temp_point);
    }
    return 1;
}

void write_data_cpp(vector<double *> &x,string s)
{
    string s1="simplify"+s;

    fstream write;
    write.open(s1,ios::app);
    write<<x.size()<<endl;
    write<<'\n';
    for(int i=0;i<x.size();i++)
    {
        write<<left<<setw(12)<<1<<setw(12)<<x[i][0]<<'\t'<<setw(12)<<x[i][1]<<'\t'<<setw(12)<<x[i][2]<<'\t'<<endl;
    }

    write.close();
}

// 将NumPy数组转换为std::vector<double*>
std::vector<double *> numpy_to_vector(py::array_t<double> numpy_array) {
    py::buffer_info buf_info = numpy_array.request();
    auto ptr = static_cast<double *>(buf_info.ptr);
    size_t num_elements = buf_info.shape[0];
    size_t dimension = buf_info.shape[1];

    std::vector<double *> point;
    for (size_t i = 0; i < num_elements; ++i) {
        point.push_back(ptr + i * dimension);
    }
    return point;
}

py::tuple homfly_str_from_numpy(py::array_t<double> numpy_array, bool flag_KMT) {
    // 将输入的NumPy数组转换为vector<double *>
    std::vector<double *> point = numpy_to_vector(numpy_array);

    // 初始化 knot 对象
    knot a(point);

    // 根据 flag_KMT 调用不同的方法
    if (flag_KMT) {
        a.KMT();
    }

    // 计算交叉矩阵并获得Gauss符号
    a.get_interSection_Matrix();
    auto result = a.get_gauss_notation();
    return py::make_tuple(std::get<0>(result), std::get<1>(result));

    // 调用homfly_str，输出homfly多项式
    // char *out;
    // char input[] = " 1 6 0 1 1 -1 2 1 0 -1 1 1 2 -1 0 1 1 1 2 1 "; // 假设这是需要传递的数据
    //out = homfly_str(input);
}

PYBIND11_MODULE(homfly, m) {
    m.def("homfly_str", &homfly_str_from_numpy, "A function that computes the homfly polynomial",
          py::arg("numpy_array"), py::arg("flag_KMT") = false);
}


