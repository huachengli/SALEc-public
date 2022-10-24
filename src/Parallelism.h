//
// Created by huacheng on 2021/5/22.
//

#ifndef SALEC_PARALLELISM_H
#define SALEC_PARALLELISM_H
#include <mpi.h>
#include "Variable.h"
#include "Ghost.h"
#include <unistd.h>

#define IdPsiL3    0
#define IdPhiL3    1
#define IdMass     2
#define IdPres     3
#define IdPresABV  4
#define IdCsound   5
#define IdVOF      6
#define IdMomentum 7
#define IdStress   8
#define IdCourant  9
#define IdGradVOF 10
#define MaxField  11
#define COMM_ON    0

struct RouteTableE
{
    int priority;
    int destTag;
    int destProc;
    struct Element * local;
};

struct RequireMap
{
    int npcs;
    int prcs[10];
    int dest[10];
    int dist[10];
};

struct ProcessInfo
{
    int npgx;
    int npgy;
    int npgz;
    int npgs;
    int noffset;
    int ** BaseId;
    int ** nwp;
    int Nxx;
    int Nyy;
    int Nzz;

    int NeProc[DIM3];
    int nrmap;
    struct RequireMap * RMAP;
    int rank;
    char hostname[200];

    int nstb;
    int nrtb;
    struct RouteTableE * STB;
    struct RouteTableE * RTB;
    double ** PackSendCache;
    double ** PackRecvCache;

    MPI_Request ** SendReq;
    MPI_Request ** RecvReq;
    int nsr;
    int nrr;

    MPI_Request * SendPsi;
    MPI_Request * SendPhi;
    MPI_Request * SendMass;
    MPI_Request * SendPres;
    MPI_Request * SendPABV;
    MPI_Request * SendCsound;
    MPI_Request * SendVOF;
    MPI_Request * SendMoment;
    MPI_Request * SendStress;
    MPI_Request * SendCourant;

    MPI_Request * RecvPsi;
    MPI_Request * RecvPhi;
    MPI_Request * RecvMass;
    MPI_Request * RecvPres;
    MPI_Request * RecvPABV;
    MPI_Request * RecvCsound;
    MPI_Request * RecvVOF;
    MPI_Request * RecvMoment;
    MPI_Request * RecvStress;
    MPI_Request * RecvCourant;

    MPI_Request * OnSend[MaxField];
    MPI_Request * OnRecv[MaxField];
    MPI_Request * SendList[MaxField];
    MPI_Request * RecvList[MaxField];

    int NumSend;
    int NumRecv;
};
typedef double * (*DataFunc)(struct Element *);
extern DataFunc IdToData[MaxField];
extern int CacheLength[MaxField];
extern int SwitchR0Table[MaxField];

void InitProcessInfo(struct Mesh * _mesh, struct ProcessInfo * _pinfo);
void ParaSearchDomainE(struct Mesh *_mesh, struct ProcessInfo *_pinfo);
void ParaSearchDomainV(struct Mesh *_mesh, struct ProcessInfo *_pinfo);
void ParaBegin(struct ProcessInfo *_pinfo, int _control);
void ParaEndAll(struct ProcessInfo *_pinfo);
void ParaEnd(struct ProcessInfo *_pinfo, int _control);
void PIFClean(struct ProcessInfo *_pinfo);
int ParaRecvSendInit(struct ProcessInfo *_pinfo);
void ParaSortDomainE(struct Mesh *_mesh, struct ProcessInfo *_pinfo);

int ParaPackInit(struct ProcessInfo *_pinfo);
int ParaPackRecvInit(struct ProcessInfo *_pinfo);
int ParaPackSendInit(struct ProcessInfo *_pinfo);
int ParaPackId(struct ProcessInfo *_pinfo, int _control);
int ParaUnpackId(struct ProcessInfo *_pinfo, int _control);
int ParaPackBegin(struct ProcessInfo *_pinfo, int _control);
int ParaPackEnd(struct ProcessInfo *_pinfo, int _control);
int ParaPackEndAll(struct ProcessInfo *_pinfo);
int ParaPackRecvSendInit(struct ProcessInfo *_pinfo);
int PostPackEnd(struct ProcessInfo *_pinfo, int _control);
void ParaDebug();
void InfoClean(struct MeshInfo * _minfo, struct ProcessInfo *_pinfo);

extern MPI_Comm SALEC_CART_COMM;
void SalecCartCommCreate(int dim[], int period[]);

#endif //SALEC_PARALLELISM_H
