#include"hull2.h"
#ifndef QUICKHULL_IMPLEMENTATION
#define QUICKHULL_IMPLEMENTATION
#include "quickhull.h"

#endif
int hull_ends(vector<double *> &x,  int L, double point1[3], double point2[3] )
{
    qh_vertex_t vertices[L];
    for(int i=0;i<L;i++)
    {
        vertices[i].x=x[i][0];
        vertices[i].y=x[i][1];
        vertices[i].z=x[i][2];
    }

    qh_mesh_t mesh=qh_quickhull3d(vertices,L);

    //*************** 凸包建立完成

    double vext1[3]{0}, vext2[3]{0}, xcen[3]{0};
    for(int i=0;i<L;i++)
    {   
        for(int d=0;d<3;d++)
        {
            xcen[d]+=x[i][d];
        }
    }

    for(int d=0;d<3;d++) 
        xcen[d] /= L;

    int N_planes= mesh.nnormals;

    // 原包里面没有包含计算平面的第四个值的部分。用数组存。算点到平面距离用点向量乘以法向量减去这个值。
    double face_dist[N_planes]{0},point1_dist[N_planes]{0},point2_dist[N_planes]{0};

    int min1_face,min2_face;
    double min1_dist=1e5,min2_dist=1e5;
    for(int i=0;i<N_planes;i++)
    {
        
        face_dist[i]+=mesh.normals[i].x * mesh.vertices[3*i].x;
        face_dist[i]+=mesh.normals[i].y * mesh.vertices[3*i].y;
        face_dist[i]+=mesh.normals[i].z * mesh.vertices[3*i].z;
        
        point1_dist[i]=qh__dist_point_plane(&vertices[0],&mesh.normals[i],face_dist[i]);
        if(point1_dist[i]<min1_dist)
        {
            min1_face=i;
            min1_dist=point1_dist[i];
        }
        point2_dist[i]=qh__dist_point_plane(&vertices[L-1],&mesh.normals[i],face_dist[i]);
        if(point2_dist[i]<min2_dist)
        {
            min2_face=i;
            min2_dist=point2_dist[i];
        }
    }
    // todo list: 这里可以让二者判断两个点在面上就停止计算，可以减少计算量



    double dis_end_end=0;
    for(int d=0;d<3;d++)
        dis_end_end+=pow((x[L-1][d]-x[0][d]),2);
    dis_end_end=sqrt(dis_end_end);

    // if(dis_end_end<(min1_dist+min2_dist))
    //     return 0;//在凸包内两点，距离最近，能否直接用线连起来。

    //**********在面上的决策与在外面的决策都不太好，可能需要重新思考。
    if(min1_dist<my_epsilon_2)
    {
        for(int d=0;d<3;d++)
        {
            vext1[d]=x[0][d]-xcen[d];
        }
    }
    else
    {
        double dot=0;
        for(int d=0;d<3;d++)
        {
            dot+=x[0][d]*mesh.normals[min1_face].v[d];
        }
        dot-=face_dist[min1_face];
        for(int d=0;d<3;d++)
        {
            vext1[d]=-dot*mesh.normals[min1_face].v[d];
        }
    }

    if(min2_dist<my_epsilon_2)
    {
        for(int d=0;d<3;d++)
        {
            vext2[d]=x[L-1][d]-xcen[d];
        }
    }
    else
    {
        double dot=0;
        for(int d=0;d<3;d++)
        {
            dot+=x[L-1][d]*mesh.normals[min2_face].v[d];
        }
        dot-=face_dist[min2_face];
        for(int d=0;d<3;d++)
        {
            vext2[d]=-dot*mesh.normals[min2_face].v[d];
        }
    }


    //rescale的值需要决策
    double rescale_max1=20,rescale_max2=20;
    double max_distace=0;

    for(int i=0;i<x.size();i++)
    {
        for(int d=0;d<3;d++)
        {
            if(fabs(x[i][d]-xcen[d])>max_distace)
                max_distace=fabs(x[i][d]-xcen[d]);
        }
    }

    double min_xyz=1e5;
    for(int d=0;d<3;d++)
    {
        if(fabs(vext1[d])<min_xyz)
            min_xyz=fabs(vext1[d]);
    }
    rescale_max1=max_distace/min_xyz;

    min_xyz=1e5;
    for(int d=0;d<3;d++)
    {
        if(fabs(vext2[d])<min_xyz)
            min_xyz=fabs(vext2[d]);
    }
    rescale_max2=max_distace/min_xyz;


    for(int d=0;d<3;d++)  point1[d]  = x[0][d]   + rescale_max1 * vext1[d];
    for(int d=0;d<3;d++)  point2[d]  = x[L-1][d] + rescale_max2 * vext2[d];
    //std::cout<<min1_dist<<std::endl;

    qh_free_mesh(mesh);
    return 1;
}