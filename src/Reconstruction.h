//
// Created by huacheng on 6/17/21.
//

#ifndef SALEC_RECONSTRUCTION_H
#define SALEC_RECONSTRUCTION_H
#include "Variable.h"
#include "Error.h"
#define RecTOL      1e-6
#define MaxCutCycle 20
#define CutMaterial 0

double CutH(double x);
double CutVolume(double d,double _p[]);
double CutLength(double _vf,double _p[]);
double CutSf(double _p[]);
void InsertBoundaryBID(struct Boundary *ib);
void ConstructPlane(struct Element * ie);
void SetVnodeCheck(struct Element * ie);
double CutVOF(struct Boundary * ib, int tag);
void MapE2Vc(struct Vertex * iv);
void CutGradVOF(struct Element * ie, int nm);
void CropVelocity(struct Element * ie, double MaxVel);
void Sort3(double _p[]);

double CutLengthC(double _vf,double _p[]);
double CutLengthS(double _vf,double _p[]);
double CutLengthL(double _vf,double _p[]);
#endif //SALEC_RECONSTRUCTION_H
