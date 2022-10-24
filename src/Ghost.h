//
// Created by huacheng on 2021/4/16.
//

#ifndef SALEC_GHOST_H
#define SALEC_GHOST_H

#include "Variable.h"
#include "Reconstruction.h"
#define MNGE 18
#define MNGV 18
#define IdFCEL1 0
#define IdFCEL2 1
#define IdFCEL3 2
#define IdFCEL4 3
#define IdFCEL5 4

struct GhostElement
{
    double WallNorm[DIM];
    struct Element * dest;
    struct Element * refe[MNGE];
    double Weight[MNGE];
    int Nrefe;
    int Type;
};

struct GhostVertex
{
    double WallNorm[DIM];
    struct Vertex * dest;
    struct Vertex * refe[MNGV];
    double Weight[MNGV];
    int Nrefe;
    int Type;
};

struct GhostElementMarker
{
    struct Element * el;
    int LongId;
    int ShortId[DIM];
    int bState;
    int dState;
    double ImagePosition[DIM];
    double WallNorm[DIM];
};

struct GhostVertexMarker
{
    struct Vertex * ve;
    int LongId;
    int ShortId;
    int bState;
    double ImagePosition[DIM];
    double WallNorm[DIM];
};

struct GhostMeshMarker
{
    int ngv;
    int nge;
//    int npe;
//    int npv;
//    int npb;
    struct GhostElement * ghtE;
    struct GhostVertex  * ghtV;
//    struct Element  ** phyE;
//    struct Vertex   ** phyV;
//    struct Boundary ** phyB;

    int ne;
    int nv;
    int * bitMapE;
    int * bitMapV;
    int (*CheckIn)(const double _pos[]);
    int (*ImageP)(const double _pos[], double _Ipos[]);
    int (*ImageNorm)(const double _pos[], double _Inorm[]);
};

int FindNeId3(int NeId[],int LongId, int depth,int nxp, int nyp, int nzp);
int SetGhostWeightE(struct GhostElement * ige, const double _Ipos[]);
int SetGhostWeightV(struct GhostVertex * igv, const double _Ipos[]);
void BuildGhostMeshMarker(struct Mesh *_mesh, struct GhostMeshMarker * _gmm);
void CutDomain(struct Mesh * _mesh, struct GhostMeshMarker * _gmm);
void FlushGhostVertex(struct GhostVertex * _igv);
void FlushGhostElement(struct GhostElement * _ige);
void FlushGhostMesh(struct GhostMeshMarker * _gmm, int _op);

void FlushGhostVertexa(struct GhostVertex * _igv);
void FlushGhostVertexb(struct GhostVertex * _igv);
void FlushGhostElementa(struct GhostElement * _ige);
void FlushGhostElementb(struct GhostElement * _ige);


/*
 * normal boundary condition implement.
 * SearchDomain in Setup will be removed.
 */
void _searchdomain(int * CondNote, int Idx, int Idy, int Idz, int Nxp, int Nyp, int Nzp, int offset);
void _searchdomain2(int * CondNote, int Idx, int Idy, int Idz, int Nxp, int Nyp, int Nzp, int offset);
void _SearchDomainV(struct Mesh *_mesh, int offset, struct MeshInfo * _meshinfo);
void SearchDomainV(struct Mesh *_mesh, struct MeshInfo * _meshinfo);
void _SearchDomainE(struct Mesh *_mesh, int offset, struct MeshInfo * _meshinfo);
void SearchDomainE2(struct Mesh *_mesh, struct MeshInfo * _meshinfo);
void FlushCondE(struct Mesh * _mesh);
void FlushCondV(struct Mesh * _mesh);
void FlushCondB(struct Mesh * _mesh);

void FlushVnode(struct Mesh * _mesh);
void FlushCondEId(struct Mesh * _mesh, int _control);
#endif //SALEC_GHOST_H
