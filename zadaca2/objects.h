#pragma once
#include "geometry.h"
#include "material.h"
#include "ray.h"
#include <iostream>
#include <cmath>
#include <limits>
using namespace std;

struct Object
{
    Material material;
    
    virtual bool ray_intersect(const Ray &ray, float &t, Vec3f &normal) const = 0;    
};

struct Sphere : Object
{
    Vec3f c; // centar
    float r; // radius
    
    Sphere(const Vec3f &c, const float &r, const Material &mat) : c(c), r(r)
    {
        Object::material = mat;
    }
    
    bool ray_intersect(const Ray &ray, float &t, Vec3f &normal) const
    {
        float d2 = ray.direction * ray.direction;   
        Vec3f e_minus_c = ray.origin - c;             
        
        float disc = pow(ray.direction * e_minus_c, 2) - d2 * (e_minus_c * e_minus_c - r * r);  
        
        bool ray_inside_sphere = e_minus_c * e_minus_c <= r * r;    
        
        if (disc < 0)        
        {
            return false;            
        }
        else
        {
            if (ray_inside_sphere)     
            {
                t = (-ray.direction * e_minus_c + sqrt(disc))/d2;   
            }
            else
            {
                t = (-ray.direction * e_minus_c - sqrt(disc))/d2;                  
            }
            
            Vec3f normal_origin = ray.origin + ray.direction * t;      
            normal = (normal_origin - c).normalize();                  
            
            return t > 0;   
        }
    }
};



vector<float> racunanje (const vector<Vec3f> &t11,const Ray &ray){
    vector<float> a;
    a.assign(4,0.f);
    Vec3f aminusorg = t11[0]-ray.origin;
    
    float det = determinant(ray.direction ,-t11[1],-t11[2]);
    float detkoef=determinant(aminusorg ,-t11[1],-t11[2]);
    float detbeta=determinant(ray.direction ,aminusorg,-t11[2]);
    float detgama=determinant(ray.direction ,-t11[1],aminusorg); 

    if (det!=0){
        a[1]=detbeta/det;
        a[2]=detgama/det;
        a[3]=detkoef/det;
    }
    a[0]=1-a[1]-a[2];
    return a;
} 

struct Cuboid : public Object{
 
    Vec3f c1; // donji kut (𝑥𝑚𝑖𝑛, 𝑦𝑚𝑖𝑛, 𝑧𝑚𝑖𝑛)
    Vec3f c2; // gornji kut (𝑥𝑚𝑎𝑥, 𝑦𝑚𝑎𝑥, 𝑧𝑚𝑎𝑥)
      

    Cuboid(const Vec3f &c1, const Vec3f &c2, const Material &mat) : c1(c1), c2(c2)
    {
        Object::material = mat;
    }

    bool ray_intersect(const Ray &ray, float &t, Vec3f &normal) const
    {   
        float x1 = c1.x;
        float y1 = c1.y;
        float z1 = c1.z;
        float x2 = c2.x;
        float y2 = c2.y;
        float z2 = c2.z;

        float mint=std::numeric_limits<float>::max();
        bool pogodak=false;
        
       
        Vec3f prvi ={x1 ,y2,z1};
        Vec3f drugi  ={x1,y1,z1};
        Vec3f treci = {x2,y1,z1};
        vector<Vec3f> t11 = {prvi,drugi-prvi,treci-prvi};  
        vector<float> rez = racunanje(t11,ray);  //racuna alfa beta gama i t
        //usporeduje parametre alfa ,beta,gama,t
        if(0<=rez[0] &&rez[0] <=1 && 0<=rez[1] &&rez[1] <=1 && 0<=rez[2] &&rez[2] <=1 && rez[3]>0){
            if (rez[3]<mint){
                mint=rez[3];
                normal= cross (drugi-prvi,treci-prvi);    
            }   
            pogodak=true;
        }
        prvi={x1,y2,z1};
        drugi={x2,y2,z1};
        treci={x2,y1,z1};
        t11={prvi,drugi-prvi,treci-prvi};
        rez = racunanje(t11,ray);
        if(0<=rez[0] &&rez[0] <=1 && 0<=rez[1] &&rez[1] <=1 && 0<=rez[2] &&rez[2] <=1 && rez[3]>0  ){
            if (rez[3]<mint){
                mint=rez[3];
                normal = cross (drugi-prvi,treci-prvi);    
            }   
            pogodak=true;
        }
        
        

        //vector<Vec3f> t21 = ({x1,y1,z2},{x1,y1,z1},{x1,y2,z2});
         prvi={x1,y2,z2};
        drugi={x1,y1,z2};
        treci={x2,y1,z2};
        t11={prvi,drugi-prvi,treci-prvi};
        rez = racunanje(t11,ray);
        if(0<=rez[0] &&rez[0] <=1 && 0<=rez[1] &&rez[1] <=1 && 0<=rez[2] &&rez[2] <=1 && rez[3]>0){
            if (rez[3]<mint){
                mint=rez[3];
                normal = cross (drugi-prvi,treci-prvi);    
            }   
            pogodak=true;
        }
        //vector<Vec3f> t22 = ({x2,y2,z2},{x1,y1,z1},{x1,y2,z2});
         prvi={x1,y2,z2};
        drugi={x2,y2,z2};
        treci={x2,y1,z2};
        t11={prvi,drugi-prvi,treci-prvi};
        rez = racunanje(t11,ray);
        if(0<=rez[0] &&rez[0] <=1 && 0<=rez[1] &&rez[1] <=1 && 0<=rez[2] &&rez[2] <=1 && rez[3]>0){
            if (rez[3]<mint){
                mint=rez[3];
                normal = cross (drugi-prvi,treci-prvi);    
            }   
            pogodak=true;
        }
        
        
      
        //vector<Vec3f> t31 = ({x1,y2,z1},{x1,y2,z2},{x2,y2,z2});
         prvi={x1,y2,z2};
        drugi={x1,y2,z1};
        treci={x2,y2,z1};
        t11={prvi,drugi-prvi,treci-prvi};
        rez = racunanje(t11,ray);
        if(0<=rez[0] &&rez[0] <=1 && 0<=rez[1] &&rez[1] <=1 && 0<=rez[2] &&rez[2] <=1 && rez[3]>0){
            if (rez[3]<mint){
                mint=rez[3];
                normal = cross (drugi-prvi,treci-prvi);    
            }   
            pogodak=true;
        }
        //vector<Vec3f> t32 = ({x2,y2,z1},{x1,y2,z2},{x2,y2,z2});
         prvi={x1,y2,z2};
        drugi={x2,y2,z2};
        treci={x2,y2,z1};
        t11={prvi,drugi-prvi,treci-prvi};
        rez = racunanje(t11,ray);
        if(0<=rez[0] &&rez[0] <=1 && 0<=rez[1] &&rez[1] <=1 && 0<=rez[2] &&rez[2] <=1 && rez[3]>0){
            if (rez[3]<mint){
                mint=rez[3];
                normal = cross (drugi-prvi,treci-prvi);    
            }   
            pogodak=true;
        }
        //vector<Vec3f> t41 = ({x1,y2,z1},{x1,y2,z2},{x1,y1,z1});
         prvi={x1,y1,z2};
        drugi={x1,y1,z1};
        treci={x2,y1,z1};
        t11={prvi,drugi-prvi,treci-prvi};
        rez = racunanje(t11,ray);
        if(0<=rez[0] &&rez[0] <=1 && 0<=rez[1] &&rez[1] <=1 && 0<=rez[2] &&rez[2] <=1 && rez[3]>0){
            if (rez[3]<mint){
                mint=rez[3];
                normal = cross (drugi-prvi,treci-prvi);    
            }   
            pogodak=true;
        }
        //vector<Vec3f> t42 = ({x1,y2,z1},{x1,y2,z2},{x1,y1,z2});
         prvi={x1,y1,z2};
        drugi={x2,y1,z2};
        treci={x2,y1,z1};
        t11={prvi,drugi-prvi,treci-prvi};
        rez = racunanje(t11,ray);
        if(0<=rez[0] &&rez[0] <=1 && 0<=rez[1] &&rez[1] <=1 && 0<=rez[2] &&rez[2] <=1 && rez[3]>0){
            if (rez[3]<mint){
                mint=rez[3];
                normal = cross (drugi-prvi,treci-prvi);    
            }   
            pogodak=true;
        }
        
        //vector<Vec3f> t51 = ({x2,y2,z1},{x2,y1,z2},{x2,y2,z2});
         prvi={x2,y2,z1};
        drugi={x2,y1,z1};
        treci={x2,y1,z2};
        t11={prvi,drugi-prvi,treci-prvi};
        rez = racunanje(t11,ray);
        if(0<=rez[0] &&rez[0] <=1 && 0<=rez[1] &&rez[1] <=1 && 0<=rez[2] &&rez[2] <=1 && rez[3]>0){
            if (rez[3]<mint){
                mint=rez[3];
                normal = cross (drugi-prvi,treci-prvi);    
            }   
            pogodak=true;
        }
        //vector<Vec3f> t52 = ({x2,y2,z1},{x2,y1,z1},{x2,y2,z2});
         prvi={x2,y2,z1};
        drugi={x2,y2,z2};
        treci={x2,y1,z2};
        t11={prvi,drugi-prvi,treci-prvi};
        rez = racunanje(t11,ray);
        if(0<=rez[0] &&rez[0] <=1 && 0<=rez[1] &&rez[1] <=1 && 0<=rez[2] &&rez[2] <=1 && rez[3]>0){
            if (rez[3]<mint){
                mint=rez[3];
                normal = cross  (drugi-prvi,treci-prvi);    
            }   
            pogodak=true;
        }
        
        //vector<Vec3f> t61 = ({x1,y2,z1},{x2,y2,z1},{x1,y1,z1});
         prvi={x1,y2,z1};
        drugi={x1,y1,z1};
        treci={x1,y1,z2};
        t11={prvi,drugi-prvi,treci-prvi};
        rez = racunanje(t11,ray);
        if(0<=rez[0] &&rez[0] <=1 && 0<=rez[1] &&rez[1] <=1 && 0<=rez[2] &&rez[2] <=1 && rez[3]>0){
            if (rez[3]<mint){
                mint=rez[3];
                normal = cross (drugi-prvi,treci-prvi);    
            }   
            pogodak=true;
        }
        //vector<Vec3f> t62 = ({x1,y2,z1},{x2,y2,z1},{x2,y1,z1});
         prvi={x1,y2,z1};
        drugi={x1,y2,z2};
        treci={x1,y1,z2};
        t11={prvi,drugi-prvi,treci-prvi};
        rez = racunanje(t11,ray);
        if(0<=rez[0] &&rez[0] <=1 && 0<=rez[1] &&rez[1] <=1 && 0<=rez[2] &&rez[2] <=1 && rez[3]>0){
            if (rez[3]<mint){
                mint=rez[3];
                normal = cross  (drugi-prvi,treci-prvi);    
            }   
            pogodak=true;
        }
        t=mint;
        
       return (pogodak);
    }
};
 
 