//
// Created by huacheng on 9/14/21.
//

#ifndef SALEC_TRACER_H
#define SALEC_TRACER_H

#include "Variable.h"
#include "Parallelism.h"
#include "InputParser.h"
#define IdTracer (MaxField+1)

struct _Tracer;
typedef struct _Tracer
{
    int NeE;
    unsigned int Id;
    int matId;
    double Position[DIM];
    double PeakPre;
    double PeakTem;
    double DumpTmp;
    struct _Tracer * next;
} Tracer;

#define TRDLEN (DIM+3)
#define TRILEN (2)

struct _TracerInfo;
typedef struct _TracerInfo
{
    Tracer * head;
    Tracer * tail;
    unsigned int TracerNum;
    double ** Coordinate;
    unsigned int RTI[DIM3];
    unsigned int STI[DIM3];
    MPI_Request SendRecvReq[DIM3*2];
    int nReq;

    // buffer allocated before send/recv
    double *SendTracer[DIM3];
    int *SendTracerId[DIM3];
    double *RecvTracer[DIM3];
    int *RecvTracerId[DIM3];

    char Prefix[MaxStrLen];
    char Type[MaxStrLen];

    char EjectaModel[MaxStrLen];
    double Hthreshold;
    double Pcenter[DIM];
    int (*EjectaFunc)(Tracer*,struct _TracerInfo *, struct Element * _mesh);
    FILE * efp;
} TracerInfo;

void TailInsertList(TracerInfo * _tinfo,Tracer * _x,unsigned int _Len);

void TracerDel(TracerInfo * _tinfo, Tracer * _x);
void TailInsert(TracerInfo * _tinfo,Tracer * _x);
void InitTracer(TracerInfo * _tinfo, InputFile * ifp,struct Mesh * _mesh);
void OutputTracer(TracerInfo * _tinfo, struct Mesh * _mesh,int step,int rank);
void OutputTracerVtp(TracerInfo * _tinfo,struct Mesh * _mesh,int step,int rank);
void OutputEjecta(TracerInfo * _tinfo,struct Mesh * _mesh,double time,int rank);
void FlushEjecta(TracerInfo * _tinfo, int step, int rank);
void TracerMove(TracerInfo * _tinfo, struct Mesh * _mesh,double dt);

void TracerInfoExchangeBegin(TracerInfo * _tinfo, struct ProcessInfo * _pinfo, struct Mesh * _mesh);
void TracerInfoExchangeEnd(TracerInfo * _tinfo, struct ProcessInfo * _pinfo);
void TracerExchangeBegin(TracerInfo * _tinfo, struct ProcessInfo * _pinfo, struct Mesh * _mesh);
void TracerExchangeEnd(TracerInfo * _tinfo, struct ProcessInfo * _pinfo, struct Mesh * _mesh);
void TracerClean(TracerInfo * _tinfo);

void InitTracerInfo(TracerInfo * _tinfo, struct MeshInfo * _minfo);

int PlanetEjectaDetect(Tracer * _t, TracerInfo * _tinfo,struct Element *);
int FlatEjectaDetect(Tracer * _t, TracerInfo * _tinfo,struct Element *);
#endif //SALEC_TRACER_H
