//
// Created by huacheng on 6/23/21.
//

#ifndef SALEC_INPUTPARSER_H
#define SALEC_INPUTPARSER_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Variable.h"
#include "Parallelism.h"
#include "Error.h"
#include "State.h"

#define MaxStrLen   100
#define MaxKeyId    200

typedef struct {
    char Key[MaxKeyId][MaxStrLen];
    char Value[MaxKeyId][MaxStrLen];
    int  Len;
} InputFile;

InputFile * OpenInputFile(const char fname[]);
InputFile * ParseInputFile(FILE * fp);
void * CloseInputFile(InputFile * ifp);
void SortInputFile(InputFile * ifp, int start, int end);
void SimpleSort(InputFile * ifp, int start, int end);
int SearchInput(InputFile * ifp,const char key[]);
void SearchTest(InputFile * ifp, const char key[]);
int GetValueI(InputFile * ifp,const char key[], char dvalue[]);
double GetValueD(InputFile * ifp,const char key[], char dvalue[]);
void GetTest(InputFile * ifp);
int GetValueS(InputFile * ifp,const char key[], char value[], char dvalue[]);
int GetStrk(char strlist[],int pos, char value[]);
double GetValueDk(InputFile * ifp,const char key[], int k,char dvalue[]);
int GetValueSk(InputFile * ifp,const char key[], char value[],int k,char dvalue[]);
int GetValueIk(InputFile * ifp,const char key[], int k,char dvalue[]);

void SetMeshInfo(InputFile * ifp, struct MeshInfo * _minfo, struct ProcessInfo *_pinfo, struct Mesh * _mesh);
void SRFInit(InputFile * ifp,struct StateReference * _s, int iMat);
void SetBC(InputFile * ifp, struct Mesh * _mesh);
void SetDamping(struct MeshInfo * _minfo, struct Mesh * _mesh);
void InitMeshANEOS(InputFile * ifp, struct Mesh * _mesh, struct ANEOSTable * ATB);
void InitMesh(InputFile * ifp, struct Mesh * _mesh, struct EosTable * _etb);
double Spacing(int k,double ext, int eL, int eR, int N);

double LoadTillEOS(struct TillotsonTable * _t, FILE * fp);
void EosInit(struct MeshInfo * _minfo, struct EosTable * _etb, struct StateReference * _srf);
void LoadStateRef(InputFile * ifp, struct MeshInfo * _minfo, struct StateReference * _srf);

void InitChk(InputFile * ifp, struct Mesh * _mesh, struct ProcessInfo * _pinfo, struct MeshInfo * _minfo);
void UpdateChk(struct Mesh * _mesh, struct MeshInfo * _minfo, struct ProcessInfo * _pinfo, int _step);
void CompareChk(struct Mesh * _mesh, struct MeshInfo * _minfo, struct ProcessInfo * _pinfo, int _step);
#endif //SALEC_INPUTPARSER_H
