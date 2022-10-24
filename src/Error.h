//
// Created by huacheng on 2020/12/22.
//

#ifndef SALEC_ERROR_H
#define SALEC_ERROR_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include "Variable.h"
#include "State.h"
#include "Output.h"
#include "Ghost.h"
#include "Setup.h"

void NullPointer(FILE *fp, const char *_message);
void CheckInt(FILE *fp, const char *_message, int x);
void RunPoint(FILE *fp, const char *_message);
void RunPoint_s(FILE *fp, const char *_message);
void CheckDouble(FILE *fp, const char *_message, double x);

void TestMeshInfo_air(struct MeshInfo *_info);
double StopWatch(int _control, int _k);
void CheckDoubleArray(FILE *fp, const char *_message, const double *x, int length);
void CheckDoubleArrayE(FILE *fp, const char *_message, const double *x, int length);
void CheckElement(const struct Element *_element);
double RandomDouble();
void CheckANEOSTable(struct ANEOSTable * _t);
void CheckEformat();
void CheckStateReference(struct StateReference * _s, struct ANEOSTable * _t);
void TestMesh(struct Mesh * E, const char * VtsName);
void TestVelStavlity(struct Vertex * iv);
void DataClean(struct ANEOSTable * ATB, struct Mesh * _mesh,int);

void Shuffle(int a[], int len);
#define __RUNPOINT__ 0
#endif //SALEC_ERROR_H
