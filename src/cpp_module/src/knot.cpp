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

vector<int> knot::knot_size(string &return_knot_type,string &knottype_desired)
{
    vector<int> result(3,0);//result 0,1 存起始点终点，result2存大小。
    string knottype;

    //vector<double *> new_points(points3D.size(),new double[3]{0,0,0});      //big mistake 怎么老是犯这种。

    vector<double *> new_points;
    double *temp_point;
    for(int i=0;i<points3D.size();i++)
    {
        temp_point=new double[3]{0};
        new_points.push_back(temp_point);
    }

    find_max_span(new_points);

    knottype=get_knottype_open_faster(new_points);
    return_knot_type=knottype;
    if(knottype!=knottype_desired)
    {
        cout<<"error, current knottype:"<<knottype<<endl;
        //error_out(new_points);
        clean_pointer(new_points);
        return result;
    }
    //cout<<knottype<<endl;
    //knottype=get_knottype(copy_knot.points3D);
    //********************get knot core************** for ring

    //************************折半查找法
    int left=0,right=new_points.size()-1;
    int knot_left=0,knot_right=new_points.size()+1;
    //int knot_left_2=0,knot_right_2=new_points.size()+1;
    int medium=(left+right)/2;

    while(left<right)
    {
        medium=(left+right)/2;
        vector<double *> temp;
        temp.assign(new_points.begin(),new_points.begin()+medium+1);
        if(get_knottype_open_faster(temp)==knottype)
            right=medium;
        else
            left=medium+1;
    }
    knot_right=left+1;

    left=0;right=new_points.size()-1;
    while(left<right-1)
    {
        //cout<<"left and right"<<left<<'\t'<<right<<endl;
        medium=(left+right)/2;
        vector<double *> temp;
        temp.assign(new_points.begin()+medium,new_points.end());
        if(get_knottype_open_faster(temp)==knottype)
            left=medium;
        else
            right=medium;
    }
    knot_left=left+1;
   


    // //second
    // left=0;right=new_points.size()-1;
    // while(left<right)
    // {
    //     medium=(left+right)/2;
    //     vector<double *> temp;
    //     temp.assign(new_points.begin()+medium,new_points.end());
    //     if(get_knottype_open(temp)==knottype)
    //         left=medium+1;
    //     else
    //         right=medium-1;
    // }
    // knot_left_2=left+1;

    // left=knot_left_2-1;right=new_points.size()-1;
    // while(left<right)
    // {
    //     medium=(left+right)/2;
    //     vector<double *> temp;
    //     temp.assign(new_points.begin()+knot_left_2-1,new_points.begin()+medium+1);
    //     if(get_knottype_open(temp)==knottype)
    //         right=medium;
    //     else
    //         left=medium+1;
    // }
    // knot_right_2=left+1;

    // knot_left=max(knot_left,knot_left_2);
    // knot_right=min(knot_right,knot_right_2);

    // 错误原因是开链会产生knot，将原本不是knotcore的部分也能看成是knotcore。slipknot
    if(knot_left>=knot_right)
    {
        // 以knot_left 为左边界，找knot_right_temp;
        left=knot_left-1;right=new_points.size()-1;
        while(left<right-1)
        {
            //cout<<"left and right"<<left<<'\t'<<right<<endl;
            medium=(left+right)/2;
            vector<double *> temp;
            temp.assign(new_points.begin()+knot_left-1,new_points.begin()+medium+1);
            if(get_knottype_open_faster(temp)==return_knot_type)
                right=medium;
            else
                left=medium+1;
        }
        int knot_right_temp=left+1;

        // 以knot_right 为右边界，找knot_left_temp;
        left=0;right=knot_right-1;
        while(left<right-1)
        {
            //cout<<"left and right"<<left<<'\t'<<right<<endl;
            medium=(left+right)/2;
            vector<double *> temp;
            temp.assign(new_points.begin()+medium,new_points.begin()+knot_right);
            if(get_knottype_open_faster(temp)==return_knot_type)
                left=medium;
            else
                right=medium;
        }
        int knot_left_temp=left+1;

        if(knot_right_temp-knot_left>knot_right-knot_left_temp)
        {
            knot_left=knot_left_temp;
        }
        else
        {
            knot_right=knot_right_temp;
        }
    }


    int flag_left=0;
    while(knot_left>1 and knot_right<new_points.size() and knot_left<knot_right)
    {
        vector<double *> temp;
        temp.assign(new_points.begin()+knot_left-1,new_points.begin()+knot_right);
        if(get_knottype_open_faster(temp)==knottype)
            break;
        else 
        {
            if(flag_left==1)
            {
                knot_right++;
                flag_left=0;
            }    
            else
            {
                knot_left--;
                flag_left=1;
            }
        }
    }

    result[0]=knot_left-1;
    result[1]=knot_right-1;
    result[2]=knot_right-knot_left+1;

    clean_pointer(new_points);
    return result;
}

vector<int> knot::knot_size_ring(string &return_knot_type,string &knottype_desired)
{
    vector<int> result(3,0);
    int n=points3D.size();

    vector<double *> new_points;
    double *temp_point;
    for(int i=0;i<n;i++)
    {
        temp_point=new double[3]{0};
        new_points.push_back(temp_point);
    }

    find_max_span(new_points);

    int count=0;
    vector<vector<int>> knotSize5(4,vector<int>(3,0));
    for(int i=0;i<4;i++)
    {
        return_knot_type=get_knottype_ring_faster(new_points);
        if(return_knot_type!=knottype_desired)
        {
            count+=n/4;
            rotate_array_vector(new_points,n/5);
            continue;
        }

        
        int left=0,right=new_points.size()-1;
        int knot_left=0,knot_right=new_points.size()+1;
        int medium=(left+right)/2;

        while(left<right)
        {
            medium=(left+right)/2;
            vector<double *> temp;
            temp.assign(new_points.begin(),new_points.begin()+medium+1);
            if(get_knottype_open_faster(temp)==return_knot_type)
                right=medium;
            else
                left=medium+1;
        }
        knot_right=left+1;

        left=0;right=new_points.size()-1;
        while(left<right-1)
        {
            //cout<<"left and right"<<left<<'\t'<<right<<endl;
            medium=(left+right)/2;
            vector<double *> temp;
            temp.assign(new_points.begin()+medium,new_points.end());
            if(get_knottype_open_faster(temp)==return_knot_type)
                left=medium;
            else
                right=medium;
        }
        knot_left=left+1;

        
        // 错误原因是开链会产生knot，将原本不是knotcore的部分也能看成是knotcore。
        if(knot_left>=knot_right)
        {
            // 以knot_left 为左边界，找knot_right_temp;
            left=knot_left-1;right=new_points.size()-1;
            while(left<right-1)
            {
                //cout<<"left and right"<<left<<'\t'<<right<<endl;
                medium=(left+right)/2;
                vector<double *> temp;
                temp.assign(new_points.begin()+knot_left-1,new_points.begin()+medium+1);
                if(get_knottype_open_faster(temp)==return_knot_type)
                    right=medium;
                else
                    left=medium+1;
            }
            int knot_right_temp=left+1;

            // 以knot_right 为右边界，找knot_left_temp;
            left=0;right=knot_right-1;
            while(left<right-1)
            {
                //cout<<"left and right"<<left<<'\t'<<right<<endl;
                medium=(left+right)/2;
                vector<double *> temp;
                temp.assign(new_points.begin()+medium,new_points.begin()+knot_right);
                if(get_knottype_open_faster(temp)==return_knot_type)
                    left=medium;
                else
                    right=medium;
            }
            int knot_left_temp=left+1;

            if(knot_right_temp-knot_left>knot_right-knot_left_temp)
            {
                knot_left=knot_left_temp;
            }
            else
            {
                knot_right=knot_right_temp;
            }
        }

        // 最后的检查
        int flag_left=0;
        while(knot_left>1 and knot_right<new_points.size() and knot_left<=knot_right)
        {
            vector<double *> temp;
            temp.assign(new_points.begin()+knot_left-1,new_points.begin()+knot_right);
            if(get_knottype_open_faster(temp)==return_knot_type)
                break;
            else 
            {
                if(flag_left==1)
                {
                    knot_right++;
                    flag_left=0;
                }    
                else
                {
                    knot_left--;
                    flag_left=1;
                }
            }
        }

        // 可能会出现带尾巴的,检查右边是否存在尾巴,删去几个点的knottype更复杂，就操作

        left=knot_left-1;right=knot_right-1;
        while (left<right-1)
        {
            medium=(left+right)/2;
            vector<double *> temp;
            temp.assign(new_points.begin()+medium,new_points.begin()+knot_right);
            string temp_knottype=get_knottype_open_faster(temp);
            if(temp_knottype==return_knot_type)
                left=medium;
            else
            {
                if(temp_knottype[0]>return_knot_type[0])
                    left=medium;
                else
                    right=medium;
            }
        }
        knot_left=left+1;

        // 可能会出现带尾巴的,检查左边是否存在尾巴,如果删去几个点的knottype更复杂，就操作
        left=knot_left-1;right=knot_right-1;
        while (left<right-1)
        {
            medium=(left+right)/2;
            vector<double *> temp;
            temp.assign(new_points.begin()+knot_left-1,new_points.begin()+medium+1);
            string temp_knottype=get_knottype_open_faster(temp);
            if(temp_knottype==return_knot_type)
                right=medium;
            else
            {
                if(temp_knottype[0]>return_knot_type[0])
                    right=medium;
                else
                    left=medium;
            }
        }
        knot_right=right+1;

        knotSize5[i][0]=(knot_left-1+count);
        knotSize5[i][1]=(knot_right-1+count);
        knotSize5[i][2]=(knot_right-knot_left+1);
        rotate_array(new_points,n/4);
        count+=n/4;
    }


    //从数组里找最小的
    int knot_left=0,knot_right=0;
    int min_size=n;
    for(int i=0;i<4;i++)
    {   
        //debug
        //cout<<knotSize5[i][0]<<'\t'<<knotSize5[i][1]<<'\t'<<knotSize5[i][2]<<endl;

        if(knotSize5[i][2]>0 and knotSize5[i][2]<min_size)
        {
            knot_left=knotSize5[i][0];
            knot_right=knotSize5[i][1];
            min_size=knotSize5[i][2];
        }
    }

    result[0]=knot_left%n;
    result[1]=knot_right%n;
    result[2]=knot_right-knot_left+1;

    clean_pointer(new_points);
    return result;

}