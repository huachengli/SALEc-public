//
// Created by huacheng on 2020/12/25.
//

#ifndef SALEC_ITERATION_H
#define SALEC_ITERATION_H

#include "Variable.h"
#include "State.h"
#include "Reconstruction.h"
#include "Error.h"
#include <math.h>

void CalculateGradM(struct Element *_element, int nm);
void CalculateABV(struct Element *_element, double _artvis);
void CalculateABV2(struct Element *_element, double _artvis, double _artvis2);
void VelocityUpdate(struct Vertex * _vertex);
double CalculateL1b(struct Vertex * iv, double dt0, double MaxVel);
void CalculateL1c(struct Vertex * _vertex, double dt);

void CalculateBFV(struct Boundary *_boundary);
void CalculateCourant(struct Element * doe);
void CalculateVOFM(struct Boundary *_boundary, double CISSAMK1, int nm);
void EnergyUpdate(struct Element *_element, const struct MeshInfo * _minfo , double dt);

void CalculateL3a(struct Boundary *_boundary, int nm);
void CalculateL3b(struct Element *_e, const struct MeshInfo * _minfo);
void CalculateL3bOut(struct Element * _e, int nm);
void CalculateL3bIn(struct Element * _e, int nm);
void CalculateL3bAdjust(struct Element * _e, int nm);
void CalculateL3bSimple(struct Element * _e, int nm);
void ANEOSCalL3c(struct Element *_e, const struct ANEOSTable _a[], const struct StateReference _s[], int nm);
void CalculateL3c(struct Element *_e, const struct EosTable _et[], const struct StateReference _s[], const struct MeshInfo * _minfo);
void CalculateL3c_gas(struct Element *_e, const struct AirEOSTable _a[], int nm);
void FailureUpdate(struct Element *_element, double dt, const struct StateReference _s[], struct MeshInfo *_minfo);
void RheologyUpdate(struct Element *_element, double dt, const struct StateReference _s[], int nm);
void SyncForward(struct Element *_element);
void SyncBackward(struct Element * _e);

extern const double FlowWeight[BPE];

void MapE2Vb(struct Vertex * iv);
void MapV2Eb(struct Element * ie, const struct MeshInfo *_info);

void Refactor(const double * Uu, const double * Du, const double * Au, double * fu, int _len);
void CalculateVOFTVD(struct Boundary *_boundary, double CISSAMK1, int nm);

void CalculateL1cABF(struct Vertex * _vertex, double dt, const double ABFAn[]);
void CalculateL2ABF(struct Element *_element, int nm, double dt, const double ABFAn[]);
void CalculateABFAn(const double dt[], double ABFAn[]);

double vofMSTACS(double nVOFD, double Gf, double CourantD);
double vofCICSAM(double nVOFD, double Gf, double CourantD);
double vofSTACS(double nVOFD, double Gf);
double vofMHRIC(double nVOFD, double Gf);
#endif //SALEC_ITERATION_H
