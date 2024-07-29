#include "myfunction.h"

/* comment
2021/12/23
*/
using namespace std;
using namespace Eigen;
extern int flag_ring_open;

void find_max_span(vector<double *> &result)
{
    int n=result.size();
    double max_x=-1e6,max_y=-1e6,max_z=-1e6,min_x=1e6,min_y=1e6,min_z=1e6;

    for(int i=0;i<n;i++)
    {
        max_x = max_x>result[i][0]?max_x:result[i][0];
        max_y = max_y>result[i][1]?max_y:result[i][1];
        max_z = max_z>result[i][2]?max_z:result[i][2];
        min_x = min_x<result[i][0]?min_x:result[i][0];
        min_y = min_y<result[i][1]?min_y:result[i][1];
        min_z = min_z<result[i][2]?min_z:result[i][2];
    }
    double dis_x = max_x-min_x;
    double dis_y = max_y-min_y;
    double dis_z = max_z-min_z;

    if(dis_x<dis_z and dis_x<dis_y)
    {
        for(int i=0;i<n;i++)
        {
            // 交换x和z
            double temp=result[i][0];
            result[i][0]=result[i][2];
            result[i][2]=temp;
        }
    }
    else if(dis_y<dis_x and dis_y<dis_z)
    {
        for(int i=0;i<n;i++)
        {
            // 交换y和z
            double temp=result[i][1];
            result[i][1]=result[i][2];
            result[i][2]=temp;
        }
    }
}

void rotate_array_vector(vector<double *> &x,int n)
{
    double *temp;
    int NB=x.size();
    if(n<0)
    {
        n=-n;
        for(int i=0;i<n;i++)
        {
            temp=new double[3];
            temp[0]=x[NB-1][0];temp[1]=x[NB-1][1];temp[2]=x[NB-1][2];
            delete[] x[NB-1];
            x.pop_back();
            x.insert(x.begin(),temp);
        }
    }

    else
    {
        for(int i=0;i<n;i++)
        {
            temp=new double[3];
            temp[0]=x[0][0];temp[1]=x[0][1];temp[2]=x[0][2];
            delete[] x[0];
            x.erase(x.begin());
            x.push_back(temp);
        }
    }
}


//n<0向右转 ,n>0 向左转
void rotate_array(vector<double *> &x,int n)
{   
    if(n<0)
    {
        n=-n;
        double temp[n][3];
        int NB=x.size();
        for(int i=0;i<n;i++)
        {  
            for(int d=0;d<3;d++)
                temp[i][d]=x[NB-n+i][d];
        }
        for(int i=NB-n-1;i>=0;i--)
        {
            for(int d=0;d<3;d++)
                x[i+n][d]=x[i][d];
        }
        for(int i=0;i<n;i++)
        {
            for(int d=0;d<3;d++)
                x[i][d]=temp[i][d];
        }
    }
    else
    {
        double temp[n][3];
        int NB=x.size();
        for(int i=0;i<n;i++)
        {  
            for(int d=0;d<3;d++)
                temp[i][d]=x[i][d];
        }
        for(int i=0;i<NB-n;i++)
        {
            for(int d=0;d<3;d++)
                x[i][d]=x[i+n][d];
        }
        for(int i=0;i<n;i++)
        {
            for(int d=0;d<3;d++)
                x[i+NB-n][d]=temp[i][d];
        }
    }
}

void write_map(map<int,int> &hist_count,fstream &write)//输出map到文件
{
    for(auto iter=hist_count.begin();iter!=hist_count.end();iter++)
    {
        write<<iter->first<<'\t'<<iter->second<<endl;
    }
}

/**
 * @brief 写入纽结core文件，在前后都添加tail个珠子，最好使用于开链
 * //TODO:
 * @param knotSize 
 * @param point 
 * @param tail ：添加的珠子数目
 * @return int 
 */
int write_core_tail(vector<int> &knotSize,vector<double *> &point,int n_tail)
{
    if(flag_ring_open==0 and knotSize[0]-n_tail<0 or knotSize[1]+n_tail>=point.size())
    {
        cout<<"knotcore 带尾巴，数组越界"<<endl;
        return 1;
    }
    fstream write_for_rmsd;
    int n=point.size();
    char  a[100];
    sprintf(a,"result/%d_tail_%d.txt",n_tail,knotSize[2]);
    write_for_rmsd.open(a,ios::app);

    write_for_rmsd<<knotSize[2]+2*n_tail<<endl<<endl;   // debug
    if(knotSize[0]>knotSize[1])
    {
        for(int i=knotSize[0]-n_tail;i<=knotSize[1]+n+n_tail;i++)
        {
            write_for_rmsd<<1<<'\t'<<point[i%n][0]<<'\t'<<point[i%n][1]<<'\t'<<point[i%n][2]<<endl;
        }
    }
    else
    {
        for(int i=knotSize[0]-n_tail;i<=knotSize[1]+n_tail;i++)
        {
            if(flag_ring_open==0 and (i<0 or i>=point.size()))
            {
                cout<<"knotcore 带尾巴，数组越界"<<endl;
                return 1;
            }
            write_for_rmsd<<1<<'\t'<<point[(i+n)%n][0]<<'\t'<<point[(i+n)%n][1]<<'\t'<<point[(i+n)%n][2]<<endl;
        }
    }
    write_for_rmsd.close();
    return 1;
}

// save x shape of knot core
int write_core_x_tail(vector<int> &knotSize,vector<double *> &point,int x_tail)
{
    if(flag_ring_open==0 and knotSize[0]-x_tail<0 or knotSize[1]+x_tail>=point.size())
    {
        cout<<"knotcore 带尾巴，数组越界"<<endl;
        return 1;
    }
    fstream write_for_rmsd;
    char  a[100];
    sprintf(a,"result/x_tail.txt");
    write_for_rmsd.open(a,ios::app);

    write_for_rmsd<<4*x_tail+2<<endl<<endl;
    if(knotSize[0]>knotSize[1])
    {   //必定是环形链
        int n=point.size();
        for(int i=knotSize[0]-x_tail;i<=knotSize[0]+x_tail;i++)
        {
            write_for_rmsd<<1<<'\t'<<point[(i<0?i+n:i)%n][0]<<'\t'<<point[(i<0?i+n:i)%n][1]<<'\t'<<point[(i<0?i+n:i)%n][2]<<endl;
        }
        for( int i=knotSize[1]-x_tail;i<=knotSize[1]+x_tail;i++)
        {
            write_for_rmsd<<1<<'\t'<<point[(i<0?i+n:i)%n][0]<<'\t'<<point[(i<0?i+n:i)%n][1]<<'\t'<<point[(i<0?i+n:i)%n][2]<<endl;
        }
    }
    else
    {
        int n=point.size();
        for(int i=knotSize[0]-x_tail;i<=knotSize[0]+x_tail;i++)
        {
            if(flag_ring_open==0 and (i<0 or i>=point.size()))
            {
                cout<<"knotcore 带尾巴，数组越界"<<endl;
                return 1;
            }
            write_for_rmsd<<1<<'\t'<<point[(i<0?i+n:i)%n][0]<<'\t'<<point[(i<0?i+n:i)%n][1]<<'\t'<<point[(i<0?i+n:i)%n][2]<<endl;
        }
        for(int i=knotSize[1]-x_tail;i<=knotSize[1]+x_tail;i++)
        {
            if(flag_ring_open==0 and (i<0 or i>=point.size()))
            {
                cout<<"knotcore 带尾巴，数组越界"<<endl;
                return 1;
            }
            write_for_rmsd<<1<<'\t'<<point[(i<0?i+n:i)%n][0]<<'\t'<<point[(i<0?i+n:i)%n][1]<<'\t'<<point[(i<0?i+n:i)%n][2]<<endl;
        }
    }
    write_for_rmsd.close();
    return 1;

}

int write_core(vector<int> &knotSize,vector<double *> &point)//输出的是knotSize是从0到N-1;
{
    fstream write_for_rmsd;
    char  a[100];
    sprintf(a,"result/core_%d.txt",knotSize[2]);
    write_for_rmsd.open(a,ios::app);

    write_for_rmsd<<knotSize[2]<<endl<<endl;
    if(knotSize[0]>knotSize[1])
    {
        int n=point.size();
        for(int i=knotSize[0];i<=knotSize[1]+n;i++)
        {
            write_for_rmsd<<1<<'\t'<<point[i%n][0]<<'\t'<<point[i%n][1]<<'\t'<<point[i%n][2]<<endl;
        }
    }
    else
    {
        for(int i=knotSize[0];i<=knotSize[1];i++)
        {
            write_for_rmsd<<1<<'\t'<<point[i][0]<<'\t'<<point[i][1]<<'\t'<<point[i][2]<<endl;
        }
    }
    write_for_rmsd.close();
    return 1;
}

int get_segment(vector<int> &segment,int offset,int line_segment)
{
    for(int i=0;i<segment.size();i++)
    {
        if(line_segment<=segment[i])
        {
            return offset+i;
        }
    }
    return 0;
}

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
    if( dis1*dis2>my_epsilon)
        return 0;
    //上面三行代码判断是否为平面同侧。
    if(fabs(dis1)<my_epsilon and fabs(dis2)<my_epsilon  )
    {
        //两个点都在面上，暂时先不给用
        return 1;
    }
    // if( fabs(dis1)<my_epsilon )
    // {
    //     Matrix2d ab;
    //     ab<<b[0]-a[0],c[0]-a[0],b[1]-a[1],c[1]-a[1];
    //     Vector2d b0,x0;
    //     b0<<line_1[0]-a[0],line_1[1]-a[1];
    //     x0=ab.inverse()*b0;
    //     if(x0(0,0)<0 or x0(1,0)<0 or x0(0,0)+x0(1,0)>1)
    //         return 0;
    //     else return 1;
    // }
    // if( fabs(dis2)<my_epsilon )
    // {
    //     Matrix2d ab;
    //     ab<<b[0]-a[0],c[0]-a[0],b[1]-a[1],c[1]-a[1];
    //     Vector2d b0,x0;
    //     b0<<line_2[0]-a[0],line_2[1]-a[1];
    //     x0=ab.inverse()*b0;
    //     if(x0(0,0)<-my_epsilon or x0(1,0)<-my_epsilon or x0(0,0)+x0(1,0)>1+my_epsilon)
    //         return 0;
    //     else return 1;
    // }
    
    

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
    if(k(0,0)<0+my_epsilon or k(0,0)>1-my_epsilon or k(1,0)<0+my_epsilon or k(1,0)>1-my_epsilon)    return 0;//用一个小量解决可能出现的共点情况

    if((k(0,0)*(line1_2[2]-line1_1[2])+line1_1[2])>(k(1,0)*(line2_2[2]-line2_1[2])+line2_1[2]))
        {//第一个片段在上方
            *up_down=1;
            if(x(0,0)*(-x(1,1))+x(1,0)*x(0,1)>0)
                return k(0,0);
            else
                return -1*k(0,0);
        }
    else
        {//第二个片段在上方
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

void cross_product(double a[3],double b[3], double res[3]) {
   res[0] = a[1]*b[2] - a[2]*b[1];
   res[1] = a[2]*b[0] - a[0]*b[2];
   res[2] = a[0]*b[1] - a[1]*b[0];
}
void cal_normals(double I[3], double J[3], double K[3], double planeijk[4]) {
   int d;
   double  vij[3], vik[3], vtmp[3], vabs;

   for(d=0;d<3;d++) vij[d] = J[d] - I[d];
   for(d=0;d<3;d++) vik[d] = K[d] - I[d];
   cross_product(vij, vik, vtmp);
   vabs = sqrt( vtmp[0]*vtmp[0] + vtmp[1]*vtmp[1] + vtmp[2]*vtmp[2] );
   if(fabs(vabs)<my_epsilon) return;
   for(d=0;d<3;d++) planeijk[d] = vtmp[d]/vabs;
   planeijk[3]  = 0.;
   for(d=0;d<3;d++)  planeijk[3] += I[d]*planeijk[d];
   planeijk[3]  = -planeijk[3];

}
// class definition



// void move3D(vector<double *> a,vector<int> b)
// {//不能用
//     int N_components=b.size()+1;
//     b.push_back(a.size());
//     vector<double *> single_knot;
//     single_knot.assign(a.begin(),a.begin()+b[0]);
//     // knot temp(single_knot,b);
//     // temp.print_knot();


//     //接口留给links，目前只用单一闭合链。
//     for(int component=1;component<N_components;component++)
//     {
//         vector<double *> single_knot;
//         printf("%d\n",b[component]);
//         single_knot.assign(a.begin()+b[component-1],a.begin()+b[component]);
//         knot temp(single_knot);
//         temp.print_knot();
//     }
// }


void read_data(vector<double *> & x,char *s)
{
    printf("当前有%d个元素\n",int(x.size()));
    
    double *temp_point;
    char read_line[200];
    FILE *fp_read;

    fp_read = fopen(s,"r");
    if(fp_read==NULL) {fprintf(stderr,"cannot find your file\n"); exit(-1);}

    int i_molecule=0,NB;
    fscanf(fp_read,"%d\n",&NB);

    if(fgetc(fp_read)!='\n')
        fscanf( fp_read,"%[^\n]", read_line);
    //fgets(read_line,20,fp_read);
    for(int i=0;i<NB;i++)  
    {
        temp_point=new double[3];
        fscanf(fp_read,"%d %lf %lf %lf\n",&i_molecule,&temp_point[0],&temp_point[1],&temp_point[2]);
        x.push_back(temp_point);
    }
    fclose(fp_read);
}
int read_data_lammps(vector<double *> &x,fstream &read)
{
    double *temp_point;
    int i_molecule=0,NB,atom_type;
    string aline;

    for(int i=0;i<3;i++)
        if(!getline(read,aline)) 
            return 0;
    
    if(!(read>>NB))
        return 0;
    
    for(int i=0;i<6;i++)
        if(!getline(read,aline)) 
            return 0;
    
    for(int i=0;i<NB;i++)
    {   
        temp_point=new double[3];
        read>>i_molecule>>atom_type>>temp_point[0]>>temp_point[1]>>temp_point[2];
        x.push_back(temp_point);
    }
    getline(read,aline);
    return 1;
}

int read_data_pdb(vector<double *> &x,fstream &read)
{
    double *temp_point;
    int i_molecule=0,NB;
    string aline;
    // 从pdb文件中读取数据
    // 一行一行读，知道遇到MODEL 开头的行，开始记录，遇到ENDMDL 结束，只有以ATOM开头的行，第7，8，9三个元素记为一个点
    // 读取的时候，把第7，8，9三个元素读入，然后把第7，8，9三个元素的值赋给temp_point[0],temp_point[1],temp_point[2]

    while(getline(read,aline))
    {
        if(aline.substr(0,5)=="MODEL")
        {
            while(getline(read,aline))
            {
                if(aline.substr(0,6)=="ENDMDL")
                    return 1;
                if(aline.substr(0,4)=="ATOM")
                {
                    temp_point=new double[3];
                    temp_point[0]=stod(aline.substr(30,8));
                    temp_point[1]=stod(aline.substr(38,8));
                    temp_point[2]=stod(aline.substr(46,8));
                    x.push_back(temp_point);
                }
            }
        }
    }
    return 0;
}

int read_data_cpp(vector<double *> &x,fstream &read)
{
    double *temp_point;
    int i_molecule=0,NB;
    string aline;

    if(!(read>>NB)) 
    {   
        //read.seekg(0,ios::end);
        return 0;
    }
    if(!getline(read,aline)) 
        return 0;
    if(!getline(read,aline)) 
        return 0;

    for(int i=0;i<NB;i++)
    {   
        temp_point=new double[3];
        read>>i_molecule>>temp_point[0]>>temp_point[1]>>temp_point[2];
        x.push_back(temp_point);
    }
    //read.close();
    return 1;
}

int read_data_cpp_oxdna(vector<double *> &x,fstream &read,int NB)
{
    double *temp_point;
    int i_molecule=0;
    string aline;

    if(!getline(read,aline)) return 0;
    if(!getline(read,aline)) return 0;
    if(!getline(read,aline)) return 0;

    for(int i=0;i<NB;i++)
    {   
        temp_point=new double[3];
        read>>temp_point[0]>>temp_point[1]>>temp_point[2];
        getline(read,aline);
        x.push_back(temp_point);
    }
    //*******OxDNA double strand***********
    for(int i=0;i<NB;i++)
    {
        getline(read,aline);
    }
    return 1;
}

//清楚new的对象
void clean_pointer(vector<double *> &x)
{
    for(int i=0;i<x.size();i++)
    {
        delete[] x[i];
    }
    x.clear();
}

void write_data_cpp(vector<double *> &x,string s)
{
    string s1=s;

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

void error_out(vector<double *> &x,string s)
{
    fstream write;
    write.open(s,ios::app);
    write<<x.size()<<endl;
    write<<'\n';
    for(int i=0;i<x.size();i++)
    {
        write<<left<<setw(12)<<1<<setw(12)<<x[i][0]<<'\t'<<setw(12)<<x[i][1]<<'\t'<<setw(12)<<x[i][2]<<'\t'<<endl;
    }

    write.close();
}

void print_knot(vector<double *> &points3D)
{
    cout<<points3D.size()<<endl<<endl;
    for(int i=0;i<points3D.size();i++)
    {
        printf("1\t%10.5lf\t%10.5lf\t%10.5lf\t\n",points3D[i][0],points3D[i][1],points3D[i][2]);
    }
}



void recenter(double a[][3],int Lknot)
{
	double center[3]={0,0,0};
	for(int i=0;i<Lknot;i++)
	{	
		for(int d=0;d<3;d++)
			center[d]+=a[i][d];
	}
	for(int i=0;i<Lknot;i++)
	{	
		for(int d=0;d<3;d++)
			a[i][d]-=center[d]/Lknot;
	}

}

void recenter(vector<double *> &a,int Lknot)
{
	double center[3]={0,0,0};
	for(int i=0;i<Lknot;i++)
	{	
		for(int d=0;d<3;d++)
			center[d]+=a[i][d];
	}
	for(int i=0;i<Lknot;i++)
	{	
		for(int d=0;d<3;d++)
			a[i][d]-=center[d]/Lknot;
	}
}