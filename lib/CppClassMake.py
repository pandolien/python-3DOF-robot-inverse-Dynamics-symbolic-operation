from .NewtonEuler import *
class CCM:
    def __init__(self):
        self.h_start = '''
#ifndef __NE_H__
#define __NE_H__
#include <stddef.h>
#include <string.h>
#define pi 3.141592653589793
#define sp(x) ((x) * (x))
class NE{
    public:
        NE();
        ~NE();
    public:
        void Initialization();
        void Init(float,float);
        void Close();
    public:
        float g[3];
        float *M,*C,*G,*t;
        float *q,*dq,*ddq;
        float col,row;
    public:
        void Setq(float*);
        void Setdq(float*);
        void Setddq(float*);
        void SetGravity(float,float,float);
    public:
        void SetM_();
        void SetC_();
        void SetG_();
    private:
        float M_(float,float);
        float C_(float,float);
        float G_(float);
    public:
        voif Update();
        
        '''
        self.h_end = '''
};
#endif
'''
        self.cpp_start = '''
#include "NE.h"
#include "Mat/lMat.h"
NE::NE(){
    Initialization();
}
NE::~NE(){
    Close();
}
void NE::Initialization(){
    memset(g,0,sizeof(float)*3)
    M = NULL;
    C = NULL;
    G = NULL;
    t = NULL;
}
void NE::Init(float i,float j){
    M = new float[i*j];
    C = new float[i*j];
    G = new float[n];
    t = new float[n];
    memset(M,0,sizeof(float)*i,j);
    memset(C,0,sizeof(float)*i,j);
    memset(G,0,sizeof(float)*i);
    memset(t,0,sizeof(float)*i);
    col = j;
    row = i;
    
}
void NE::Close(){
    if(M)delete[](M);
    if(C)delete[](C);
    if(G)delete[](G);
    if(t)delete[](t);
    M = NULL;
    C = NULL;
    G = NULL;
    t = NULL;
}
void NE::Setq(float *q){
    this->q = q;
}
void NE::Setdq(float *dq){
    this->dq = dq;
}
void NE::Setddq(float *ddq){
    this->ddq = ddq;
}
float NE::M_(float i,float j){
    return M[i*row+j];
}
float NE::C_(float i,float j){
    return C[i*row+j];
}
float NE::G_(float i){
    return G[i];
}
void NE::Update(){
    lMat M,C,G;
}
'''
        
        
