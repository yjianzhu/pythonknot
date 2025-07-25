#include "knottype.h"
#include <ginac/ginac.h>
#include <mutex>
#include "knot_alex_table.h"
#include "hull2.h"

using namespace std;
extern std::mutex ginac_mutex;

// extern GiNaC::symbol t;
// extern GiNaC::lst syms;
// extern std::map<GiNaC::ex,string,GiNaC::ex_is_less>  alexander_polynomial;  //查阅官方文档给出的给ex排序的方式

// only for open chain
vector<int> get_gauss_notation(vector<double *> &points)
{
    int n=points.size();
    if(n<3)
    {   
        return vector<int>(0);
    }

    int up_down=0;
    vector<vector<double>> i_Matrix(n-1,vector<double> (n-1,0));
    vector<vector<int>> i_Matrix_up_down(n-1,vector<int> (n-1,0));
    for(int i=0;i<n-1;i++)
    {
        for(int j=0;j<n-1;j++)
        {   
            if(i==j)    continue;
            up_down=0;
            i_Matrix[i][j]=cal_interSection(points[i],points[(i+1)],points[j],points[(j+1)],&up_down);
            i_Matrix_up_down[i][j]=up_down;
        }
    }

    map<pair<int,int>,pair<int,int>> crossing;
    vector<int> segment; //存片段信息
    // crossing 前两个存哪两个线段交叉，后面的是交叉点计数，以及正负。
    int count_crossing=0,count_segement=0;//从0开始

    for(int i=0;i<n-1;i++)
    {
        vector<pair<int,double>> oneLine;//记录这个线段与哪些片段交叉

        for(int j=0;j<n-1;j++)
        {
            if(fabs(i_Matrix[i][j])>my_epsilon)
            {
                oneLine.push_back(pair<int,double>(j,fabs(i_Matrix[i][j])));
            }
        }

        sort(oneLine.begin(),oneLine.end(),pair_compare);//排序，知道这个线段上的交叉点信息。

        for(auto it=oneLine.begin();it!=oneLine.end();it++)
        {
            pair<int,int> temp;
            if(i<it->first) {temp.first=i;temp.second=it->first;}
            else {temp.first=it->first;temp.second=i;}

            if(crossing.count(temp)==0)
            {
                int isOver0=i_Matrix[i][it->first]>0.?1:-1;
                crossing.insert(pair<pair<int,int>,pair<int,int>>(temp,pair<int,int>(count_crossing,isOver0)));
                count_crossing++;
            }
            //确保新的交叉点已经加入了map；

            if(i_Matrix_up_down[i][it->first]==-1)
            {
                count_segement++;
                segment.push_back(i);
            }
        }
    }

    // 计算extended gauss notation
    std::vector<int> notation;
    std::vector<int> viewed(count_crossing,0);  // 记录当前交叉点是否访问过

    for(int index=0;index<n-1;index++)
    {
        vector<pair<int,double>> oneLine;
        for(int j=0;j<n-1;j++)
        {
            if(fabs(i_Matrix[index][j])>my_epsilon)
                oneLine.push_back(pair<int,double>(j,fabs(i_Matrix[index][j])));
        }
        sort(oneLine.begin(),oneLine.end(),pair_compare);
        for(auto it=oneLine.begin();it!=oneLine.end();it++)
        {   
            pair<int,int> temp;
            if(index<it->first) {temp.first=index;temp.second=it->first;}
            else {temp.first=it->first;temp.second=index;}

            if(viewed[crossing[temp].first]==0)
            {
                viewed[crossing[temp].first]=1;
                notation.push_back((crossing[temp].first+1)* i_Matrix_up_down[index][it->first]);
            }
            else
            {
                // 第二次访问，存入手性信息
                notation.push_back((crossing[temp].first+1)* crossing[temp].second);
            }
        }
    }
    return notation;
}

vector<double *> KMT_open_chain(vector<double *> points)
{
    int flag=0;
    
    while (true)
    {
        int number=points.size();//作为判断是否结束循环的标志，当没有vertex可以删除的时候，就退出。

        for(int i=1;i<points.size()-1;i++)
        {
            //cout<<i<<'\t'<<points.size()<<endl;
            double plain[4]={0};
            flag=0;
            cal_normals(points[i-1],points[i],points[(i+1)],plain);

            if(fabs(plain[0])<my_epsilon and fabs(plain[1])<my_epsilon and fabs(plain[2])<my_epsilon)
                {//如果三个点在一条线上，可以省去
                points.erase(points.begin()+i%points.size());
                i--;
                //cout<<i<<endl;
                //cout<<points.size()<<endl;
                continue;
                }

            for(int j=0;j<points.size()-1;j++)
                {
                    if(j==i-1 or j==i) continue;
                    if(judge_triangle(points[i-1],points[i],points[(i+1)],plain,points[j],points[(j+1)]))
                        {
                            flag=1;
                            break;
                        }
                }
            if(flag==0)
            {
                points.erase(points.begin()+i%points.size());
                i--;
                //cout<<i<<endl;
                //cout<<points.size()<<endl;
            }

        }
        if(number==points.size())
            break;
    }
    return points;
}

vector<double *> KMT(vector<double*> points)
{
    int flag=0;
    
    while (true)
    {
        int number=points.size();//作为判断是否结束循环的标志，当没有vertex可以删除的时候，就退出。

        for(int i=1;i<=points.size();i++)
        {
            //cout<<i<<'\t'<<points.size()<<endl;
            double plain[4]={0};
            flag=0;
            cal_normals(points[i-1],points[i%points.size()],points[(i+1)%points.size()],plain);

            if(fabs(plain[0])<my_epsilon and fabs(plain[1])<my_epsilon and fabs(plain[2])<my_epsilon)
                {//如果三个点在一条线上，可以省去
                points.erase(points.begin()+i%points.size());
                i--;
                //cout<<i<<endl;
                //cout<<points.size()<<endl;
                continue;
                }

            for(int j=i+1;j<=points.size()+i-2;j++)
                {
                    if(judge_triangle(points[i-1],points[i%points.size()],points[(i+1)%points.size()],plain,points[j%points.size()],points[(j+1)%points.size()]))
                        {
                            flag=1;
                            break;
                        }
                }
            if(flag==0)
            {
                points.erase(points.begin()+i%points.size());
                i--;
                //cout<<i<<endl;
                //cout<<points.size()<<endl;
            }

        }
        if(number==points.size())
            break;
    }
    return points;

}

//void get_interSection_Matrix(vector<double *> &points)
//*************************存交叉点矩阵
//*************************     i_Matrix 存的是交叉的位置，范围从0-1.即第一个线段，第一个点到第二个点的线段上的位置
//*************************     从

string get_knottype_by_matrix_open(vector<double *> &points)
{
    //先算矩阵元素，再用ginac库计算行列式
    
    int n=points.size();
    if(n<3)
    {   
        return "1";
    }

    int up_down=0;
    vector<vector<double>> i_Matrix(n-1,vector<double> (n-1,0));
    vector<vector<int>> i_Matrix_up_down(n-1,vector<int> (n-1,0));
    for(int i=0;i<n-1;i++)
    {
        for(int j=0;j<n-1;j++)
        {   
            if(i==j)    continue;
            up_down=0;
            i_Matrix[i][j]=cal_interSection(points[i],points[(i+1)],points[j],points[(j+1)],&up_down);
            i_Matrix_up_down[i][j]=up_down;
        }
        //这里可以优化一下，j>i开始，但是函数要重写，保证得到两个交叉点的位置。
    }

    map<pair<int,int>,pair<int,int>> crossing;
    vector<int> segment; //存片段信息
    // crossing 前两个存哪两个线段交叉，后面的是交叉点计数，以及正负。
    int count_crossing=0,count_segement=0;//从0开始

    for(int i=0;i<n-1;i++)
    {
        vector<pair<int,double>> oneLine;//记录这个线段与哪些片段交叉

        for(int j=0;j<n-1;j++)
        {
            if(fabs(i_Matrix[i][j])>my_epsilon)
            {
                oneLine.push_back(pair<int,double>(j,fabs(i_Matrix[i][j])));
            }
        }

        sort(oneLine.begin(),oneLine.end(),pair_compare);//排序，知道这个线段上的交叉点信息。

        for(auto it=oneLine.begin();it!=oneLine.end();it++)
        {
            pair<int,int> temp;
            if(i<it->first) {temp.first=i;temp.second=it->first;}
            else {temp.first=it->first;temp.second=i;}

            if(crossing.count(temp)==0)
            {
                int isOver0=i_Matrix[i][it->first]>0.?1:-1;
                crossing.insert(pair<pair<int,int>,pair<int,int>>(temp,pair<int,int>(count_crossing,isOver0)));
                count_crossing++;
            }
            //确保新的交叉点已经加入了map；

            if(i_Matrix_up_down[i][it->first]==-1)
            {
                count_segement++;
                segment.push_back(i);
            }
        }
    }
    
    //print for test
    // cout<<"segement:";
    // for(int i=0;i<segment.size();i++)
    //     cout<<segment[i]<<'\t';
    // cout<<endl;

    // cout<<"crossing number:"<<crossing.size()<<'\t'<<count_crossing<<endl;

    //end test
    std::unique_lock<std::mutex> lock(ginac_mutex);
    GiNaC::matrix m(count_crossing,count_crossing);
    //GiNaC::ex a=m.determinant();
    
    for(int i=0;i<n-1;i++)
    {
        vector<pair<int,double>> oneLine;//记录这个线段与哪些片段交叉

        for(int j=0;j<n-1;j++)
        {
            if(fabs(i_Matrix[i][j])>my_epsilon)
            {
                oneLine.push_back(pair<int,double>(j,fabs(i_Matrix[i][j])));
            }
        }

        sort(oneLine.begin(),oneLine.end(),pair_compare);

        int offset=0; //同一个线段上有两个片段。
        for(auto it=oneLine.begin();it!=oneLine.end();it++)
        {
            pair<int,int> temp;
            if(i<it->first) {temp.first=i;temp.second=it->first;}
            else {temp.first=it->first;temp.second=i;}

            if(crossing[temp].second==1)//positive crossing
            {
                if(i_Matrix_up_down[i][it->first]==-1)
                    {
                        //cout<<crossing[temp].first<<endl;
                        int i_segment=get_segment(segment,offset,i);
                        m.operator()(crossing[temp].first,(i_segment)%count_crossing)-=1;
                        m.operator()(crossing[temp].first,(i_segment+1)%count_crossing)+=t;
                        offset++;
                        
                    }
                    else
                    {
                        //cout<<crossing[temp].first<<endl;
                        int i_segment=get_segment(segment,offset,i);
                        //cout<<i_segment<<endl;
                        m.operator()(crossing[temp].first,(i_segment)%count_crossing)+=1-t;
                    }

            }
            else    //negative crossing
            {
                if(i_Matrix_up_down[i][it->first]==-1)
                    {
                        //cout<<crossing[temp].first<<endl;
                        int i_segment=get_segment(segment,offset,i);
                        m.operator()(crossing[temp].first,(i_segment)%count_crossing)+=t;
                        m.operator()(crossing[temp].first,(i_segment+1)%count_crossing)-=1;
                        offset++;
                    }
                else
                    {
                        //cout<<crossing[temp].first<<endl;
                        int i_segment=get_segment(segment,offset,i);
                        //cout<<i_segment<<endl;
                        m.operator()(crossing[temp].first,(i_segment)%count_crossing)+=1-t;
                    }
            }
        }
    }

    if(count_crossing==0)   
    {   
        return "1";
    }
    //*********官方文档没有写清楚这里怎么获得子矩阵，sub_matrix, reduced_matrix 返回的是ex.手写子矩阵
    GiNaC::matrix m_sub(count_crossing-1,count_crossing-1);
    for(int i=0;i<count_crossing-1;i++)
    {
        for(int j=0;j<count_crossing-1;j++)
            m_sub.operator()(i,j)=m.operator()(i,j);
    }

    
    //用ldegree 和degree两个函数，返回系数，保证最低系数为0；
    GiNaC::ex poly_now=m_sub.determinant();
    
    if(poly_now.ldegree(t)!=0)
    {
        poly_now=(poly_now*pow(t,-poly_now.ldegree(t))).expand();
    }

    string a;

    if(alexander_polynomial.count(poly_now)==1)   
    {
        a=alexander_polynomial[poly_now];
    }
    else if ( alexander_polynomial.count(-poly_now)==1)
    {
        a=alexander_polynomial[-poly_now];
    }
    else
        cout<<"knot type not found, with polynomial as:"<<poly_now<<endl;
    
    return a;
}

//闭合链的alexander polynomial
string get_knottype_by_matrix(vector<double *> &points)
{
    //先算矩阵元素，再用ginac库计算行列式
    int n=points.size();
    int up_down;
    vector<vector<double>> i_Matrix(n,vector<double> (n,0));
    vector<vector<int>> i_Matrix_up_down(n,vector<int> (n,0));
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {   
            if(i==j)    continue;
            up_down=0;
            i_Matrix[i][j]=cal_interSection(points[i],points[(i+1)%n],points[j],points[(j+1)%n],&up_down);
            i_Matrix_up_down[i][j]=up_down;
        }
        //这里可以优化一下，j>i开始，但是函数要重写，保证得到两个交叉点的位置。
    }

    map<pair<int,int>,pair<int,int>> crossing;
    vector<int> segment; //存片段信息
    // crossing 前两个存那两个线段交叉，后面的是交叉点计数，以及正负。
    int count_crossing=0,count_segement=0;//从0开始

    for(int i=0;i<n;i++)
    {
        vector<pair<int,double>> oneLine;//记录这个线段与哪些片段交叉

        for(int j=0;j<n;j++)
        {
            if(fabs(i_Matrix[i][j])>my_epsilon)
            {
                oneLine.push_back(pair<int,double>(j,fabs(i_Matrix[i][j])));
            }
        }

        sort(oneLine.begin(),oneLine.end(),pair_compare);//排序，知道这个线段上的交叉点信息。

        for(auto it=oneLine.begin();it!=oneLine.end();it++)
        {
            pair<int,int> temp;
            if(i<it->first) {temp.first=i;temp.second=it->first;}
            else {temp.first=it->first;temp.second=i;}

            if(crossing.count(temp)==0)
            {
                int isOver0=i_Matrix[i][it->first]>0.?1:-1;
                crossing.insert(pair<pair<int,int>,pair<int,int>>(temp,pair<int,int>(count_crossing,isOver0)));
                count_crossing++;
            }
            //确保新的交叉点已经加入了map；

            if(i_Matrix_up_down[i][it->first]==-1)
            {
                count_segement++;
                segment.push_back(i);
            }
        }
    }
    //print for test
    // cout<<"segement:";
    // for(int i=0;i<segment.size();i++)
    //     cout<<segment[i]<<'\t';
    // cout<<endl;

    // cout<<"crossing number:"<<crossing.size()<<'\t'<<count_crossing<<endl;

    //end test

    count_segement=0;
    //GiNaC::ex poly=-1+t-t*t;
    GiNaC::matrix m(count_crossing,count_crossing);
    //GiNaC::ex a=m.determinant();

    for(int i=0;i<n;i++)
    {
        vector<pair<int,double>> oneLine;//记录这个线段与哪些片段交叉

        for(int j=0;j<n;j++)
        {
            if(fabs(i_Matrix[i][j])>my_epsilon)
            {
                oneLine.push_back(pair<int,double>(j,fabs(i_Matrix[i][j])));
            }
        }

        sort(oneLine.begin(),oneLine.end(),pair_compare);

        int offset=0; //同一个线段上有两个片段。
        for(auto it=oneLine.begin();it!=oneLine.end();it++)
        {
            pair<int,int> temp;
            if(i<it->first) {temp.first=i;temp.second=it->first;}
            else {temp.first=it->first;temp.second=i;}

            if(crossing[temp].second==1)//positive crossing
            {
                if(i_Matrix_up_down[i][it->first]==-1)
                    {
                        //cout<<crossing[temp].first<<endl;
                        int i_segment=get_segment(segment,offset,i);
                        m.operator()(crossing[temp].first,i_segment)+=-1;
                        m.operator()(crossing[temp].first,(i_segment+1)%count_crossing)+=t;
                        offset++;
                    }
                    else
                    {
                        //cout<<crossing[temp].first<<endl;
                        int i_segment=get_segment(segment,offset,i);
                        //cout<<i_segment<<endl;
                        m.operator()(crossing[temp].first,i_segment)+=1-t;
                    }

            }
            else    //negative crossing
            {
                if(i_Matrix_up_down[i][it->first]==-1)
                    {
                        //cout<<crossing[temp].first<<endl;
                        int i_segment=get_segment(segment,offset,i);
                        m.operator()(crossing[temp].first,i_segment)+=t;
                        m.operator()(crossing[temp].first,(i_segment+1)%count_crossing)+=-1;
                        offset++;
                    }
                    else
                    {
                        //cout<<crossing[temp].first<<endl;
                        int i_segment=get_segment(segment,offset,i);
                        //cout<<i_segment<<endl;
                        m.operator()(crossing[temp].first,i_segment)+=1-t;
                    }
            }
        }
    }

    //*********官方文档没有写清楚这里怎么获得子矩阵，sub_matrix, reduced_matrix 返回的是ex.手写子矩阵
    GiNaC::matrix m_sub(count_crossing-1,count_crossing-1);
    for(int i=0;i<count_crossing-1;i++)
    {
        for(int j=0;j<count_crossing-1;j++)
        m_sub.operator()(i,j)=m.operator()(i,j);
    }

    //用ldegree 和degree两个函数，返回系数，保证最低系数为0；
    GiNaC::ex poly_now=m_sub.determinant();
    cout<<poly_now<<endl;
    
    if(poly_now.ldegree(t)!=0)
    {
        poly_now=(poly_now*pow(t,-poly_now.ldegree(t))).expand();
    }
    cout<<poly_now<<endl;
    string a;

    if(alexander_polynomial.count(poly_now)==1)   
    {
        a=alexander_polynomial[poly_now];
    }
    else if ( alexander_polynomial.count(-poly_now)==1)
    {
        a=alexander_polynomial[-poly_now];
    }
    return a;
}

string get_knottype(vector<double *> &points)
{
    vector<double *> simplify_points=KMT(points);//无论有无简化过，都走一遍KMT
    //print_knot(simplify_points);
    string Alexander=get_knottype_by_matrix(simplify_points);
    return Alexander;
}

string get_knottype_ring(vector<double *> &points)
{
    vector<double *> simplify_points=points;
    //************************************************************  给points添加尾部
    double *start_point;
    start_point=new double[3];
    for(int i=0;i<3;i++)
        start_point[i]=simplify_points[0][i];
    simplify_points.push_back(start_point);
    simplify_points=KMT(simplify_points);   // 适用于一开始没有kmt，放在这里加速

    string  result=get_knottype_by_matrix_open(simplify_points);
    delete[] start_point;

    return  result;
}

string get_knottype_ring_faster(vector<double *> &points)
{
    vector<double *> simplify_points=KMT(points);
    //************************************************************  给points添加尾部
    double *start_point;
    start_point=new double[3];
    for(int i=0;i<3;i++)
        start_point[i]=simplify_points[0][i];
    simplify_points.push_back(start_point);


    string  result=get_knottype_by_matrix_open(simplify_points);
    delete[] start_point;

    return  result;
}

string get_knottype_open_faster(vector<double *> &points)
{
    vector<double *> simplify_points=KMT_open_chain(points);
    //************************************************************  给points添加尾部
    double *start_point,*end_point;
    start_point=new double[3];
    end_point=new double[3];

    try
    {
        if(simplify_points.size()>4) 
        {
            int flag_add_end=hull_ends(simplify_points,simplify_points.size(),start_point,end_point);

            //cout<<flag_add_end<<endl;

            if(flag_add_end == 1) 
            {       
                simplify_points.insert(simplify_points.begin(),start_point);
                simplify_points.push_back(end_point);

            } 
            else 
            { // it is shorter to close two ends by a line
                for(int i=0;i<3;i++)
                    start_point[i]=simplify_points[0][i];
                simplify_points.push_back(start_point);
            }
        }
    }
    catch(int e)
    {

    };


    string result=get_knottype_by_matrix_open(simplify_points);

    delete[] start_point,end_point;
    return  result;
}