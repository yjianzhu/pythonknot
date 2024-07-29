#include "knot.h"
#include "knottype.h"
#include <iostream>
#include <mutex>

using namespace std;
void knot::find_max_span(vector<double *> &result)
{
    int n=this->points3D.size();

    double max_x=-1e6,max_y=-1e6,max_z=-1e6,min_x=1e6,min_y=1e6,min_z=1e6;
    for(int i=0;i<n;i++)
    {
        max_x=max(this->points3D[i][0],max_x);
        max_y=max(this->points3D[i][1],max_y);
        max_z=max(this->points3D[i][2],max_z);

        min_x=min(this->points3D[i][0],min_x);
        min_y=min(this->points3D[i][1],min_y);
        min_z=min(this->points3D[i][2],min_z);
    }
    double dis_x=max_x-min_x,
            dis_y=max_y-min_y,
            dis_z=max_z-min_z;
    
    
    if(dis_x<dis_z and dis_x<dis_y)
    {
        for(int i=0;i<n;i++)
        {
            result[i][2]=points3D[i][0];
            result[i][1]=points3D[i][1];
            result[i][0]=points3D[i][2];
        }
        //cout<<"x-z"<<endl;
    }
    else if(dis_y<dis_x and dis_y<dis_z)
    {
        for(int i=0;i<n;i++)
        {
            result[i][2]=points3D[i][1];
            result[i][0]=points3D[i][0];
            result[i][1]=points3D[i][2];
        }
        //cout<<"y-z"<<endl;
    }
    else
    {
        for(int i=0;i<n;i++)
        {
            result[i][1]=points3D[i][1];
            result[i][0]=points3D[i][0];
            result[i][2]=points3D[i][2];
        }
    }
}

knot::knot(vector<double*> &a,vector<int> b/* args */)
{
    points3D=a;   
    ends=b;
}

knot::knot(vector<double*> &a,int b/* args */)  //new 需要释放
{
    double * temp;
    for(int i=0;i<a.size();i++)
    {   
        temp=new double[3];
        temp[0]=a[i][0];temp[1]=a[i][1];temp[2]=a[i][2];
        points3D.push_back(temp);
    }
}

knot::~knot()
{
    //printf("class has been killed!\n");
}

void knot::print_knot(){
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

//*********************************************************implementation of KMT algorithm***********************************************************
void knot::KMT()
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

                if(fabs(plain[0])<my_epsilon and fabs(plain[1])<my_epsilon and fabs(plain[2])<my_epsilon)
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

