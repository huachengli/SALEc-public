//
// Created by huacheng on 2021/4/16.
//

#include "Ghost.h"
#include "Iteration.h"


int FindNeId3(int NeId[],int LongId, int depth,int nxp, int nyp, int nzp)
{
    int Idx,Idy,Idz;
    ToShortId(&Idx,&Idy,&Idz,LongId,nxp,nyp,nzp);
    int kNeId = 0;
    for(int k=-depth;k<=depth;k++)
    {
        for(int j=-depth;j<=depth;j++)
        {
            for(int i=-depth;i<=depth;i++)
            {
                if(abs(k)!=depth && abs(j)!=depth && abs(i)!=depth) continue;
                NeId[kNeId++] = ToLongId(Idx+i,Idy+j,Idz+k,nxp,nyp,nzp);
            }
        }
    }
    return kNeId;
}

int SetGhostWeightE(struct GhostElement * ige, const double _Ipos[])
{
    double SumWeight = 0.0;
    for(int k=0;k<ige->Nrefe;k++)
    {
        ige->Weight[k] = 1.0/(1e-6 + Distance(ige->refe[k]->Center, _Ipos) * Distance(ige->refe[k]->Center, _Ipos));
        SumWeight += ige->Weight[k];
    }
    for(int k=0;k<ige->Nrefe;k++)
        ige->Weight[k] /= SumWeight;
    return ige->Nrefe;
}

int SetGhostWeightV(struct GhostVertex * igv, const double _Ipos[])
{
    double SumWeight = 0.0;
    for(int k=0;k<igv->Nrefe;k++)
    {
        igv->Weight[k] = 1.0/(1.0e-6 + Distance(igv->refe[k]->Position,_Ipos)*Distance(igv->refe[k]->Position,_Ipos));
        SumWeight += igv->Weight[k];
    }
    for(int k=0;k<igv->Nrefe;k++)
        igv->Weight[k] /= SumWeight;
}

void BuildGhostMeshMarker(struct Mesh *_mesh, struct GhostMeshMarker * _gmm)
{
    int * bitMapV;
    int * bitMapE;
    bitMapV = (int *) malloc(sizeof(int) * _mesh->nv);
    bitMapE = (int *) malloc(sizeof(int) * _mesh->ne);
    for(int k=0;k<_mesh->nv;k++) bitMapV[k] = 0;
    for(int k=0;k<_mesh->ne;k++) bitMapE[k] = 0;
    _gmm->nge = _gmm->ngv = 0;

    for(int k=0;k<_mesh->nv;k++)
    {
        if(_gmm->CheckIn(_mesh->Vertexes[k].Position))
            bitMapV[k] = 4; // in ghost domain
        else
            bitMapV[k] = 3; // in physical domain
    }

    for(int k=0;k<_mesh->ne;k++)
    {
        int bState = 0;
        for(int j=0;j<VPE;j++)
        {
            int offset = (int)(_mesh->Elements[k].NeV[j] - &_mesh->Vertexes[0]);
            if((4 == bitMapV[offset])||(2 == bitMapV[offset])) SetBit(&bState,j);
        }
        if(0==bState)
        {
            bitMapE[k] = 3; // in physical domain
        }
        else if(255==bState)
        {
            bitMapE[k] = 4; // in ghost domain
        }
        else
        {
            if(_gmm->CheckIn(_mesh->Elements[k].Center))
            {
                bitMapE[k] = 2; // in ghost domain and set as a ghost element
                _gmm->nge ++;
            }
            else
                bitMapE[k] = 1; // in ghost domain but not set as a ghost element

            for(int j=0;j<VPE;j++)
            {
                int offset = (int) (_mesh->Elements[k].NeV[j] - &_mesh->Vertexes[0]);
                if(offset >= _mesh->nv) exit(0);
                if(4 == bitMapV[offset])
                    bitMapV[offset] = 2; // in ghost domain and set as a ghost vertex
                else if(3 == bitMapV[offset])
                    bitMapV[offset] = 1; // near ghost domain but not set as a ghost vertex
            }
        }
    }


    // allocate the ghost elements and ghost vertex
    struct GhostElement * ghtE = (struct GhostElement *) malloc(sizeof(struct GhostElement)*_gmm->nge);
    int kghtE = -1;
    for(int k=0;k<_mesh->ne;k++)
    {
        if(2!=bitMapE[k])
            continue;
        else
            kghtE ++;
        ghtE[kghtE].dest = &_mesh->Elements[k];
        double ImagePos[DIM] = {0.0, 0.0, 0.0};
        _gmm->ImageP(_mesh->Elements[k].Center, ImagePos);
        int NeId[26];
        int rNeId = FindNeId3(NeId,k,1,_mesh->nxp-1,_mesh->nyp-1,_mesh->nzp-1);
        if(rNeId != 26)
        {
            fprintf(stdout,"error in FindNeId3\n");
            exit(0);
        }


        ghtE[kghtE].Nrefe = 0;
        for(int j=0;j<26;j++)
        {
            if(-1==NeId[j])
                continue;
            else if((bitMapE[NeId[j]]==1)||(bitMapE[NeId[j]]==3))
                ghtE[kghtE].refe[ghtE[kghtE].Nrefe++] = &_mesh->Elements[NeId[j]];
        }

        if(0==ghtE[kghtE].Nrefe)
        {
            int bitMapLocal[26] = {0};
            double CenterLocal[26][3] = {0.0};
            double DistanceLocal[26] = {0.0};
            double tmpX1[3] = {-30.0,0.0,0.0};
            int FlagLocal[26] = {-1};
            for(int i=0;i<26;i++)
            {
                if(-1!=NeId[i])
                {
                    bitMapLocal[i] = bitMapE[NeId[i]];
                    Copy(_mesh->Elements[NeId[i]].Center, CenterLocal[i]);
                    tmpX1[Z] = CenterLocal[i][Z];
                    DistanceLocal[i] = Distance(tmpX1,CenterLocal[i]);
                    FlagLocal[i] = _gmm->CheckIn(CenterLocal[i]);
                }
            }
            fprintf(stdout,"Error in FindNeId!\n");
        }
        else
            SetGhostWeightE(&ghtE[kghtE],ImagePos);
    }

    // allocate the ghost vertex
    for(int k=0;k<_mesh->nv;k++)
    {
        if(2 == bitMapV[k])
            _gmm->ngv ++;
    }
    struct GhostVertex * ghtV = (struct GhostVertex *) malloc(sizeof(struct GhostVertex)*_gmm->ngv);
    int kghtV = -1;
    for(int k=0;k<_mesh->nv;k++)
    {
        if(2 !=  bitMapV[k])
            continue;
        else
            kghtV ++;

        ghtV[kghtV].dest = &_mesh->Vertexes[k];
        double ImagePos[DIM] = {0.0, 0.0, 0.0};
        _gmm->ImageP(_mesh->Vertexes[k].Position,ImagePos);
        int NeId[26];
        int rNeId = FindNeId3(NeId,k,1,_mesh->nxp,_mesh->nyp,_mesh->nzp);
        if(rNeId != 26)
        {
            fprintf(stdout,"error in FindNeId3\n");
            exit(0);
        }
        ghtV[kghtV].Nrefe = 0;
        for(int j=0;j<rNeId;j++)
        {
            if(-1==NeId[j])
                continue;
            else if((3==bitMapV[NeId[j]])||(1==bitMapV[NeId[j]]))
            {
                ghtV[kghtV].refe[ghtV[kghtV].Nrefe++] = &_mesh->Vertexes[NeId[j]];
            }
        }

        if(0 == ghtV[kghtV].Nrefe)
        {
            fprintf(stdout,"cannot found around vertexes in FindNeId");
            exit(0);
        }
        else
            SetGhostWeightV(&ghtV[kghtV],ImagePos);

        _gmm->ImageNorm(_mesh->Vertexes[k].Position,ghtV[kghtV].WallNorm);
        ghtV[kghtV].Type = 0; // reversed boundary
    }

    _gmm->ne = _mesh->ne;
    _gmm->nv = _mesh->nv;
    _gmm->bitMapE = bitMapE;
    _gmm->bitMapV = bitMapV;
    _gmm->ghtE = ghtE;
    _gmm->ghtV = ghtV;
}

void CutDomain(struct Mesh * _mesh, struct GhostMeshMarker * _gmm)
{
    int npv = 0;
    int npe = 0;
    int npb = 0;
    int * NoteV = (int *) malloc(sizeof(int)*_mesh->ndv);
    int * NoteE = (int *) malloc(sizeof(int)*_mesh->nde);
    int * NoteB = (int *) malloc(sizeof(int)*_mesh->ndb);
    for(int k=0;k<_mesh->ndv;k++)
    {
        int offset = (int) (_mesh->DomainV[k] - &_mesh->Vertexes[0]);
        if((3 == _gmm->bitMapV[offset]) || (1 == _gmm->bitMapV[offset]))
        {
            npv++;
            NoteV[k] = 1;
        }
        else
            NoteV[k] = 0;
    }

    for(int k=0;k<_mesh->nde;k++)
    {
        int offset = (int) (_mesh->DomainE[k] - &_mesh->Elements[0]);
        if((3==_gmm->bitMapE[offset])||(1==_gmm->bitMapE[offset]))
        {
            npe++;
            NoteE[k] = 1;
        }
        else
            NoteE[k] = 0;
    }

    for(int k=0;k<_mesh->ndb;k++)
    {
        int offsetNe0 = (int) (_mesh->DomainB[k]->NeE[0] - &_mesh->Elements[0]);
        int offsetNe1 = (int) (_mesh->DomainB[k]->NeE[1] - &_mesh->Elements[0]);

        if((4!=_gmm->bitMapE[offsetNe0])&&(4!=_gmm->bitMapE[offsetNe1]))
        {
            npb ++;
            NoteB[k] = 1;
        }
        else
            NoteB[k] = 0;
    }

    struct Vertex ** phyV  = (struct Vertex **)   malloc(sizeof(struct Vertex *)   * npv);
    struct Element **phyE  = (struct Element **)  malloc(sizeof(struct Element *)  * npe);
    struct Boundary **phyB = (struct Boundary **) malloc(sizeof(struct Boundary *) * npb);
    int kphyV = -1;
    for(int k=0;k<_mesh->ndv;k++)
    {
        if(0==NoteV[k])
            continue;
        else
        {
            phyV[++kphyV] = _mesh->DomainV[k];
        }
    }

    int kphyE = -1;
    for(int k=0;k<_mesh->nde;k++)
    {
        if(0==NoteE[k])
            continue;
        else
        {
            phyE[++kphyE] = _mesh->DomainE[k];
        }
    }

    int kphyB = -1;
    for(int k=0;k<_mesh->ndb;k++)
    {
        if(0==NoteB[k])
            continue;
        else
        {
            phyB[++kphyB] = _mesh->DomainB[k];
            int offsetUp0 = (int) (_mesh->DomainB[k]->UpE[0] - &_mesh->Elements[0]);
            int offsetUp1 = (int) (_mesh->DomainB[k]->UpE[1] - &_mesh->Elements[0]);
            if(4 == _gmm->bitMapE[offsetUp0])
                _mesh->DomainB[k]->UpE[0] = _mesh->DomainB[k]->NeE[0];
            if(4 == _gmm->bitMapE[offsetUp1])
                _mesh->DomainB[k]->UpE[1] = _mesh->DomainB[k]->NeE[1];
        }
    }

    free(NoteV);
    free(NoteE);
    free(NoteB);
    free(_mesh->DomainV);
    free(_mesh->DomainE);
    free(_mesh->DomainB);
    _mesh->DomainV = phyV;
    _mesh->DomainE = phyE;
    _mesh->DomainB = phyB;
    _mesh->ndv = npv;
    _mesh->nde = npe;
    _mesh->ndb = npb;
}

void FlushGhostVertex(struct GhostVertex * _igv)
{
    _igv->dest->MassL1 = 0.0;
    ZeroEx(_igv->dest->VOF,NMT);
    Zero(_igv->dest->Velocity);
    Zero(_igv->dest->DistanceL1);

    for(int k=0;k<_igv->Nrefe;k++)
    {
        _igv->dest->MassL1 += _igv->Weight[k] * _igv->refe[k]->MassL1;
        ScalerAdditionEx(_igv->dest->VOF,_igv->refe[k]->VOF,_igv->Weight[k],NMT);
        ScalerAddition(_igv->dest->Velocity,_igv->refe[k]->Velocity,_igv->Weight[k]);
        ScalerAddition(_igv->dest->DistanceL1,_igv->refe[k]->DistanceL1,_igv->Weight[k]);
    }

    ScalerAddition(_igv->dest->DistanceL1,_igv->WallNorm,-2.0*Dot(_igv->WallNorm,_igv->dest->DistanceL1));
    ScalerAddition(_igv->dest->Velocity  ,_igv->WallNorm,-2.0*Dot(_igv->WallNorm,_igv->dest->Velocity  ));
}

void FlushGhostElement(struct GhostElement * _ige)
{
    for(int j=0;j<NMT;j++)
    {
        ZeroEx(_ige->dest->PhiL3[j],NPHI);
        Zero(_ige->dest->GradVOF[j]);
        for(int k=0;k<_ige->Nrefe;k++)
        {
            ScalerAdditionEx(_ige->dest->PhiL3[j],_ige->refe[k]->PhiL3[j],_ige->Weight[k],NPHI);
            ScalerAddition(_ige->dest->GradVOF[j],_ige->refe[k]->GradVOF[j],_ige->Weight[k]);
        }
    }
    ZeroEx(_ige->dest->PsiL3,NPSI);
    for(int k=0;k<_ige->Nrefe;k++)
        ScalerAdditionEx(_ige->dest->PsiL3,_ige->refe[k]->PsiL3,_ige->Weight[k],NPSI);

    // add pressure and ...
    _ige->dest->Pressure = 0.0;
    _ige->dest->ArtificialPressure = 0.0;
    _ige->dest->Density = 0.0;
    _ige->dest->Mass    = 0.0;
    for(int k=0;k<_ige->Nrefe;k++)
    {
        _ige->dest->Pressure += _ige->Weight[k] * _ige->refe[k]->Pressure;
        _ige->dest->ArtificialPressure += _ige->Weight[k] * _ige->refe[k]->ArtificialPressure;
        _ige->dest->Density  += _ige->Weight[k] * _ige->refe[k]->Density;
        _ige->dest->Mass     += _ige->Weight[k] * _ige->refe[k]->Mass;
    }
    SyncBackward(_ige->dest);
}

void FlushGhostMesh(struct GhostMeshMarker * _gmm, int _op)
{
    if(0 == _op)
    {
        // flush the vertex
        for(int k=0;k<_gmm->ngv;k++)
        {
            FlushGhostVertex(_gmm->ghtV+k);
        }
    }
    else if(1 == _op)
    {
        for(int k=0;k<_gmm->nge;k++)
            FlushGhostElement(_gmm->ghtE + k);
    }
}

void FlushGhostVertexa(struct GhostVertex * _igv)
{
    _igv->dest->MassL1 = 0.0;
    _igv->dest->CsoundL1 = 0.0;
    ZeroEx(_igv->dest->VOF,NMT);
    Zero(_igv->dest->Velocity);

    for(int k=0;k<_igv->Nrefe;k++)
    {
        _igv->dest->MassL1   += _igv->Weight[k] * _igv->refe[k]->MassL1;
        _igv->dest->CsoundL1 += _igv->Weight[k] * _igv->refe[k]->CsoundL1;
        ScalerAdditionEx(_igv->dest->VOF,_igv->refe[k]->VOF,_igv->Weight[k],NMT);
        ScalerAddition(_igv->dest->Velocity,_igv->refe[k]->Velocity,_igv->Weight[k]);
    }

    ScalerAddition(_igv->dest->Velocity  ,_igv->WallNorm, -2.0*Dot(_igv->WallNorm,_igv->dest->Velocity  ));
}

void FlushGhostVertexb(struct GhostVertex * _igv)
{
    Zero(_igv->dest->Velocity);
    Zero(_igv->dest->DistanceL1);

    for(int k=0;k<_igv->Nrefe;k++)
    {
        ScalerAddition(_igv->dest->Velocity,_igv->refe[k]->Velocity,_igv->Weight[k]);
        ScalerAddition(_igv->dest->DistanceL1,_igv->refe[k]->DistanceL1,_igv->Weight[k]);
    }

    ScalerAddition(_igv->dest->DistanceL1,_igv->WallNorm,-2.0*Dot(_igv->WallNorm,_igv->dest->DistanceL1));
    ScalerAddition(_igv->dest->Velocity  ,_igv->WallNorm,-2.0*Dot(_igv->WallNorm,_igv->dest->Velocity  ));
}

void FlushGhostElementa(struct GhostElement * _ige)
{
    for(int j=0;j<NMT;j++)
    {
        ZeroEx(_ige->dest->PhiL3[j],NPHI);
        ZeroEx(_ige->dest->sPhiL3[j],NPHI);
        for(int k=0;k<_ige->Nrefe;k++)
        {
            ScalerAdditionEx(_ige->dest->PhiL3[j],_ige->refe[k]->PhiL3[j],_ige->Weight[k],NPHI);
            ScalerAdditionEx(_ige->dest->sPhiL3[j],_ige->refe[k]->sPhiL3[j],_ige->Weight[k],NPHI);
        }
    }
}

void FlushGhostElementb(struct GhostElement * _ige)
{
    for(int j=0;j<NMT;j++)
    {
        ZeroEx(_ige->dest->PhiL3[j],NPHI);
        Zero(_ige->dest->GradVOF[j]);
        for(int k=0;k<_ige->Nrefe;k++)
        {
            ScalerAdditionEx(_ige->dest->PhiL3[j],_ige->refe[k]->PhiL3[j],_ige->Weight[k],NPHI);
            ScalerAddition(_ige->dest->GradVOF[j],_ige->refe[k]->GradVOF[j],_ige->Weight[k]);
        }
    }
    ZeroEx(_ige->dest->PsiL3,NPSI);
    for(int k=0;k<_ige->Nrefe;k++)
        ScalerAdditionEx(_ige->dest->PsiL3,_ige->refe[k]->PsiL3,_ige->Weight[k],NPSI);

    // add pressure and ...
    _ige->dest->Pressure = 0.0;
    _ige->dest->Density = 0.0;
    _ige->dest->Mass    = 0.0;
    _ige->dest->Csound  = 0.0;
    _ige->dest->Temperature = 0.0;
    for(int k=0;k<_ige->Nrefe;k++)
    {
        _ige->dest->Pressure += _ige->Weight[k] * _ige->refe[k]->Pressure;
        _ige->dest->Density  += _ige->Weight[k] * _ige->refe[k]->Density;
        _ige->dest->Mass     += _ige->Weight[k] * _ige->refe[k]->Mass;
        _ige->dest->Csound   += _ige->Weight[k] * _ige->refe[k]->Csound;
        _ige->dest->Temperature += _ige->Weight[k] * _ige->refe[k]->Temperature;
    }
    SyncBackward(_ige->dest);
}

void _searchdomain(int * CondNote, int Idx, int Idy, int Idz, int Nxp, int Nyp, int Nzp, int offset)
{
    int State1 = 0, State2 = 0;
    if(offset < 0) offset = abs(offset);

    if(offset>Idx) SetBit(&State1,X);
    if(offset>Idy) SetBit(&State1,Y);
    if(offset>Idz) SetBit(&State1,Z);
    if(Idx>Nxp-1-offset) SetBit(&State1,DIM+X);
    if(Idy>Nyp-1-offset) SetBit(&State1,DIM+Y);
    if(Idz>Nzp-1-offset) SetBit(&State1,DIM+Z);

    if(Idx<(offset-1))   SetBit(&State2,X);
    if(Idx>(Nxp-offset)) SetBit(&State2,DIM+X);
    if(Idy<(offset-1))   SetBit(&State2,Y);
    if(Idy>(Nyp-offset)) SetBit(&State2,DIM+Y);
    if(Idz<(offset-1))   SetBit(&State2,Z);
    if(Idz>(Nzp-offset)) SetBit(&State2,DIM+Z);

    if(0==State1)
        *CondNote = -1;
    else if((0!=State1)&&(0==State2))
    {
        int tIdx = 0, tIdy = 0, tIdz=0;
        if((offset-1)==Idx)   tIdx = 1;
        if((Nxp-offset)==Idx) tIdx = -1;

        if((offset-1)==Idy)   tIdy = 1;
        if((Nyp-offset)==Idy) tIdy = -1;

        if((offset-1)==Idz)   tIdz = 1;
        if((Nzp-offset)==Idz) tIdz = -1;

        *CondNote = ToLongId(Idx+tIdx,Idy+tIdy,Idz+tIdz,Nxp,Nyp,Nzp);
    }
    else
        *CondNote = -2;
}

void _searchdomain2(int * CondNote, int Idx, int Idy, int Idz, int Nxp, int Nyp, int Nzp, int offset)
{
    int State1 = 0, State2 = 0;
    if(offset < 0) offset = abs(offset);

    if(offset>Idx) SetBit(&State1,X);
    if(offset>Idy) SetBit(&State1,Y);
    if(offset>Idz) SetBit(&State1,Z);
    if(Idx>Nxp-1-offset) SetBit(&State1,DIM+X);
    if(Idy>Nyp-1-offset) SetBit(&State1,DIM+Y);
    if(Idz>Nzp-1-offset) SetBit(&State1,DIM+Z);

    if(Idx<(offset-1))   SetBit(&State2,X);
    if(Idx>(Nxp-offset)) SetBit(&State2,DIM+X);
    if(Idy<(offset-1))   SetBit(&State2,Y);
    if(Idy>(Nyp-offset)) SetBit(&State2,DIM+Y);
    if(Idz<(offset-1))   SetBit(&State2,Z);
    if(Idz>(Nzp-offset)) SetBit(&State2,DIM+Z);

    if(0==State1)
        *CondNote = -1;
    else if((0!=State1)&&(0==State2))
    {
        if((offset-1)==Idx)   SetBit(CondNote,X);
        if((Nxp-offset)==Idx) SetBit(CondNote,DIM + X);

        if((offset-1)==Idy)   SetBit(CondNote,Y);
        if((Nyp-offset)==Idy) SetBit(CondNote,DIM + Y);

        if((offset-1)==Idz)   SetBit(CondNote,Z);
        if((Nzp-offset)==Idz) SetBit(CondNote,DIM + Z);
    }
    else
        *CondNote = -2;
}


void _SearchDomainV(struct Mesh *_mesh,int offset,struct MeshInfo * _meshinfo)
{
    int * DomainMapV = (int *) malloc(_mesh->nv* sizeof(int));
    int * CondNote = (int *) malloc(_mesh->nv* sizeof(int));
    int Nxp = _mesh->nxp, Nyp=_mesh->nyp, Nzp = _mesh->nzp, Ncv=0, Ndv=0;

    for(int Idx=0;Idx<Nxp;Idx++)
    {
        for(int Idy=0;Idy<Nyp;Idy++)
        {
            for(int Idz=0;Idz<Nzp;Idz++)
            {
                int * VCondNote = CondNote + ToLongId(Idx,Idy,Idz,Nxp,Nyp,Nzp);
                _searchdomain2(VCondNote,Idx,Idy,Idz,Nxp,Nyp,Nzp,offset);

                if(0 <= *VCondNote)
                    Ncv++;
                else if(-1 == *VCondNote)
                    Ndv++;
            }
        }
    }
    struct ConditionV * CondV = (struct ConditionV *) malloc(Ncv * sizeof(struct ConditionV));
    struct Vertex **  DomainV = (struct Vertex **)    malloc(Ndv * sizeof(struct Veretx *));
    int kCondV = 0, kDomainV=0;
    for(int k=0;k<_mesh->nv;k++)
    {
        if(CondNote[k]>=0)
        {
            int Idx,Idy,Idz,tIdx=0,tIdy=0,tIdz=0;
            ToShortId(&Idx,&Idy,&Idz,k,Nxp,Nyp,Nzp);
            if(GetBit(CondNote[k],X))       tIdx = 1;
            if(GetBit(CondNote[k],X + DIM)) tIdx = -1;
            if(GetBit(CondNote[k],Y))       tIdy = 1;
            if(GetBit(CondNote[k],Y + DIM)) tIdy = -1;
            if(GetBit(CondNote[k],Z))       tIdz = 1;
            if(GetBit(CondNote[k],Z + DIM)) tIdz = -1;
            CondV[kCondV].dest = _mesh->Vertexes + k;

            for(int j=1;j<=NRCV;j++)
            {
                CondV[kCondV].refe[j-1] = _mesh->Vertexes + ToLongId(Idx+j*tIdx,Idy+j*tIdy,Idz+j*tIdz,Nxp,Nyp,Nzp);
            }
            kCondV ++;
        }
        else if(-1==CondNote[k])
        {
            DomainV[kDomainV] = _mesh->Vertexes + k;
            kDomainV ++;
        }
    }

    _mesh->CondV   = CondV;
    _mesh->ncv     = Ncv;
    _mesh->DomainV = DomainV;
    _mesh->ndv     = Ndv;
    free(DomainMapV);
    free(CondNote);
}

void FlushCondE(struct Mesh * _mesh)
{
    for(int k=0;k<_mesh->nce;k++)
    {
        struct ConditionE * ice = _mesh->CondE + k;
        ice->dest->Mass               = ice->refe->Mass;
        ice->dest->Pressure           = ice->refe->Pressure;
        ice->dest->ArtificialPressure = ice->refe->ArtificialPressure;
        ice->dest->Temperature        = ice->refe->Temperature;
        ice->dest->Density            = ice->refe->Density;

        for(int i=0;i<NMT;i++)
        {
            CopyEx(ice->refe->PhiL3[i],ice->dest->PhiL3[i],NPHI);
            CopyEx(ice->refe->sPhiL3[i],ice->dest->sPhiL3[i],NPHI);
        }
        CopyEx(ice->refe->PsiL3,ice->dest->PsiL3,NPSI);
        CopyEx(ice->refe->sPsiL3,ice->dest->sPsiL3,NPSI);

        // Velocity and momentum
        Copy(ice->refe->Momentum,ice->dest->Momentum);
    }
}

void FlushCondEL1(struct Mesh * _mesh)
{
    for(int k=0;k<_mesh->nce;k++)
    {
        struct ConditionE * ice = _mesh->CondE + k;
        CopyEx(ice->refe->VOF,ice->dest->VOF,NMT);
    }
}

void ApplyBC(int _type, double _tr,double _norm[], double _vd[], double _vr[])
{
    if(1==_type)
    {
        Copy(_vr,_vd);
        ScalerAddition(_vd,_norm,_tr* Dot(_vr,_norm));
    } else if(0==_type)
    {
        ScalerMove(_vd,_vr,-1.0);
    }
}

void FlushCondEL2(struct Mesh * _mesh)
{
    for(int k=0;k<_mesh->nce;k++)
    {
        struct ConditionE * ice = _mesh->CondE + k;
        for(int j=0;j<NMT;j++)
        {
            ApplyBC(ice->CEType,ice->Tr[0],ice->Wnorm[0],ice->dest->GradVOF[j],ice->refe->GradVOF[j]);
        }
        ConstructPlane(ice->dest);
    }
}

void FlushCondEL3(struct Mesh * _mesh)
{
    for(int k=0;k<_mesh->nce;k++)
    {
        struct ConditionE * ice = _mesh->CondE + k;
        ice->dest->Mass               = ice->refe->Mass;
        ice->dest->Csound             = ice->refe->Csound;
        if(1==ice->CEType)
        {
            Copy(ice->refe->Momentum,ice->dest->Momentum);
            ScalerAddition(ice->dest->Momentum,ice->Wnorm[0], ice->Tr[0]*Dot(ice->Wnorm[0],ice->refe->Momentum));
        } else if(2==ice->CEType)
        {
            ScalerMove(ice->dest->Momentum,ice->refe->Momentum,-1.0);
        } else if(0==ice->CEType)
        {
            Zero(ice->dest->Momentum);
        }
    }
}

void FlushCondEL4(struct Mesh * _mesh)
{
    for(int k=0;k<_mesh->nce;k++)
    {
        struct ConditionE * ice = _mesh->CondE + k;
        ice->dest->Pressure             = ice->refe->Pressure;
        ice->dest->ArtificialPressure   = ice->refe->ArtificialPressure;
        CopyEx(ice->refe->PsiL3,ice->dest->PsiL3,NPSI);
    }
}

void FlushCondEL5(struct Mesh * _mesh)
{
    for(int k=0;k<_mesh->nce;k++)
    {
        struct ConditionE * ice = _mesh->CondE + k;
        ice->dest->Courant      = ice->refe->Courant;
        for(int j=0;j<NMT;j++)
        {
            CopyEx(ice->refe->PhiL3[j],ice->dest->PhiL3[j],NPHI);
        }
        if(1==ice->CEType)
        {
            Copy(ice->refe->Momentum,ice->dest->Momentum);
            ScalerAddition(ice->dest->Momentum,ice->Wnorm[0], ice->Tr[0]*Dot(ice->Wnorm[0],ice->refe->Momentum));
        } else if(2==ice->CEType)
        {
            ScalerMove(ice->dest->Momentum,ice->refe->Momentum,-1.0);
        }
    }
}

typedef void (*FCELfunc)(struct Mesh *);
FCELfunc FlushCondEList[5] = {
        FlushCondEL1,
        FlushCondEL2,
        FlushCondEL3,
        FlushCondEL4,
        FlushCondEL5
};

void FlushCondEId(struct Mesh * _mesh, int _control)
{
    FlushCondEList[_control](_mesh);
}

void FlushCondV(struct Mesh * _mesh)
{

    for(int k=0;k<_mesh->ncv;k++)
    {
        struct ConditionV * icv = _mesh->CondV + k;
        struct Vertex * rv1 = icv->refe[0];
        struct Vertex * rv2 = icv->refe[1];
        struct Vertex * rv3 = icv->refe[2];
        struct Vertex * iv = icv->dest;
        // interpolation scheme 1
       /* iv->DistanceL1[Z] = rv1->DistanceL1[Z];
        iv->Velocity[Z] = rv1->Velocity[Z];*/
        iv->MassL1 = rv1->MassL1;
        Copy(rv1->Velocity,iv->Velocity);
        Copy(rv1->DistanceL1,iv->DistanceL1);
        Copy(rv1->VelocityL1,iv->VelocityL1);

       /* // interpolation scheme 2 3-order
        double tmpV[DIM],tmpD[DIM];
        LinearOp(tmpV,rv1->Velocity,rv2->Velocity,3.0,-3.0);
        ScalerAddition(tmpV,rv3->Velocity,1.0);
        LinearOp(tmpD,rv1->DistanceL1,rv2->DistanceL1,3.0,-3.0);
        ScalerAddition(tmpD,rv3->DistanceL1,1.0);
        Copy(tmpD,iv->DistanceL1);
        Copy(tmpV,iv->Velocity);*/

    }
}

void FlushVnode(struct Mesh * _mesh)
{
    for(int k=0;k<_mesh->nce;k++)
    {
        struct ConditionE * ice = _mesh->CondE + k;
        ConstructPlane(ice->dest);
    }
}

void FlushCondB(struct Mesh * _mesh)
{
    for(int k=0;k<_mesh->ncb;k++)
    {
        struct ConditionB * icb = _mesh->CondB + k;
        icb->dest->dMass = icb->refe->dMass;
        Copy(icb->refe->dMomentum,icb->dest->dMomentum);
        CopyEx(icb->refe->dPsi,icb->dest->dPsi,NPSI);
        for(int j=0;j<NMT;j++)
              CopyEx(icb->refe->dPhi[j],icb->dest->dPhi[j],NPHI);
    }
}

void SearchDomainE3(struct Mesh *_mesh, struct MeshInfo * _meshinfo)
{
    int * DomainMapE = (int *) malloc(_mesh->ne * sizeof(int));
    int * CondNote   = (int *) malloc(_mesh->ne * sizeof(int));
    int Nxp=_mesh->nxp-1, Nyp=_mesh->nyp-1, Nzp = _mesh->nzp-1, Nce=0, Nde=0;
    for(int Idx=0;Idx<Nxp;Idx++)
    {
        for(int Idy=0;Idy<Nyp;Idy++)
        {
            for(int Idz=0;Idz<Nzp;Idz++)
            {
                int * VState    = DomainMapE  + ToLongId(Idx,Idy,Idz,Nxp,Nyp,Nzp);
                int * VCondNote = CondNote + ToLongId(Idx,Idy,Idz,Nxp,Nyp,Nzp);
                _searchdomain(VCondNote,Idx,Idy,Idz,Nxp,Nyp,Nzp,2);
                if(*VCondNote >= 0)
                    Nce++;
                else if(*VCondNote == -1)
                    Nde++;
            }
        }
    }

    struct Element ** DomainE = (struct Element **) malloc(Nde* sizeof(struct Element *));
    struct ConditionE * CondE = (struct ConditionE *) malloc(Nce * sizeof(struct ConditionE));
    int kDomainE=0, kCondE=0;
    for(int k=0;k<_mesh->ne;k++)
    {
        if(CondNote[k] >= 0)
        {
            CondE[kCondE].dest = _mesh->Elements + k;
            CondE[kCondE].refe = _mesh->Elements + CondNote[k];
            kCondE++;
        }
        else if(-1==CondNote[k])
        {
            DomainE[kDomainE++] = _mesh->Elements + k;
        }
    }

    _mesh->CondE = CondE;
    _mesh->nce   = Nce;
    _mesh->DomainE = DomainE;
    _mesh->nde = Nde;


    _mesh->DomainB = (struct Boundary **) malloc((_mesh->nbx + _mesh->nby + _mesh->nbz)* sizeof(struct Boundary *));
    _mesh->ndb = 0;
    for(int i=0;i<_mesh->nbx;i++)
    {
        struct Boundary * ib = &_mesh->XBoundaries[i];
        if(NULL==ib->NeE[0] || NULL==ib->NeE[1]) continue;
        _mesh->DomainB[_mesh->ndb++] = ib;
    }
    for(int i=0;i<_mesh->nby;i++)
    {
        struct Boundary * ib = &_mesh->YBoundaries[i];
        if(NULL==ib->NeE[0] || NULL==ib->NeE[1]) continue;
        _mesh->DomainB[_mesh->ndb++] = ib;
    }
    for(int i=0;i<_mesh->nbz;i++)
    {
        struct Boundary * ib = &_mesh->ZBoundaries[i];
        if(NULL==ib->NeE[0] || NULL==ib->NeE[1]) continue;
        _mesh->DomainB[_mesh->ndb++] = ib;
    }


    int * DomainMapB = (int *) malloc(_mesh->ndb * sizeof(int));
    int Ndb = 0, Ncb=0;
    for(int k=0;k<_mesh->ndb;k++)
    {
        struct Boundary * ib = _mesh->DomainB[k];
        int kNe0 = (int) (ib->NeE[0] - _mesh->Elements);
        int kNe1 = (int) (ib->NeE[1] - _mesh->Elements);

        if((-1==CondNote[kNe0])&&(-1==CondNote[kNe1]))
        {
            Ndb++;
            DomainMapB[k] = -1;
        }
        else if((-1==CondNote[kNe0])&&(-1!=CondNote[kNe1]))
        {
            DomainMapB[k] = 0;
            Ncb++;
        }
        else if((-1!=CondNote[kNe0])&&(-1==CondNote[kNe1]))
        {
            DomainMapB[k] = 1;
            Ncb++;
        }
        else
            DomainMapB[k] = -2;

    }

    struct Boundary ** DomainB = (struct Boundary **) malloc(Ndb * sizeof(struct Boundary *));
    struct ConditionB * CondB = (struct ConditionB *) malloc(Ncb * sizeof(struct ConditionB));
    int kDomainB  = 0, kCondB = 0;
    for(int k=0;k<_mesh->ndb;k++)
    {
        if(-1 == DomainMapB[k])
            DomainB[kDomainB++] = _mesh->DomainB[k];
        else if(DomainMapB[k] >= 0)
        {
            struct Element * tmpE = _mesh->DomainB[k]->NeE[DomainMapB[k]];
            int Direction = -1;
            double DotRes = 0.0;
            for(int j=0;j<BPE;j++)
            {
                if(_mesh->DomainB[k] == tmpE->NeB[j]) continue;
                double jDotRes = fabs(Dot(_mesh->DomainB[k]->Normal,tmpE->NeB[j]->Normal));
                if(-1==Direction)
                {
                    Direction = j;
                    DotRes = jDotRes;
                }
                else if(jDotRes >= DotRes)
                {
                    Direction = j;
                    DotRes = jDotRes;
                }
            }
            CondB[kCondB].dest = _mesh->DomainB[k];
            CondB[kCondB].refe = tmpE->NeB[Direction];
            kCondB ++;
        }
    }

    free(_mesh->DomainB);
    _mesh->ndb = Ndb;
    _mesh->DomainB = DomainB;
    _mesh->ncb = Ncb;
    _mesh->CondB = CondB;
    free(CondNote);
    free(DomainMapE);
    free(DomainMapB);
}

void _SearchDomainE(struct Mesh *_mesh, int offset, struct MeshInfo * _meshinfo)
{
    int * DomainMapE = (int *) malloc(_mesh->ne * sizeof(int));
    int * CondNote   = (int *) malloc(_mesh->ne * sizeof(int));
    int Nxp=_mesh->nxp-1, Nyp=_mesh->nyp-1, Nzp = _mesh->nzp-1, Nce=0, Nde=0;
    for(int Idx=0;Idx<Nxp;Idx++)
    {
        for(int Idy=0;Idy<Nyp;Idy++)
        {
            for(int Idz=0;Idz<Nzp;Idz++)
            {
                int * VState    = DomainMapE  + ToLongId(Idx,Idy,Idz,Nxp,Nyp,Nzp);
                int * VCondNote = CondNote + ToLongId(Idx,Idy,Idz,Nxp,Nyp,Nzp);
                _searchdomain(VCondNote,Idx,Idy,Idz,Nxp,Nyp,Nzp,offset);
                if(*VCondNote >= 0)
                    Nce++;
                else if(*VCondNote == -1)
                    Nde++;
            }
        }
    }

    struct Element ** DomainE = (struct Element **) malloc(Nde* sizeof(struct Element *));
    struct ConditionE * CondE = (struct ConditionE *) malloc(Nce * sizeof(struct ConditionE));
    int kDomainE=0, kCondE=0;
    for(int k=0;k<_mesh->ne;k++)
    {
        if(CondNote[k] >= 0)
        {
            CondE[kCondE].dest = _mesh->Elements + k;
            CondE[kCondE].refe = _mesh->Elements + CondNote[k];
            kCondE++;
        }
        else if(-1==CondNote[k])
        {
            DomainE[kDomainE++] = _mesh->Elements + k;
        }
    }

    _mesh->CondE = CondE;
    _mesh->nce   = Nce;
    _mesh->DomainE = DomainE;
    _mesh->nde = Nde;

    _mesh->DomainB = (struct Boundary **) malloc((_mesh->nbx + _mesh->nby + _mesh->nbz)* sizeof(struct Boundary *));
    _mesh->ndb = 0;
    for(int i=0;i<_mesh->nbx;i++)
    {
        struct Boundary * ib = &_mesh->XBoundaries[i];
        if(NULL==ib->NeE[0] || NULL==ib->NeE[1]) continue;
        _mesh->DomainB[_mesh->ndb++] = ib;
    }
    for(int i=0;i<_mesh->nby;i++)
    {
        struct Boundary * ib = &_mesh->YBoundaries[i];
        if(NULL==ib->NeE[0] || NULL==ib->NeE[1]) continue;
        _mesh->DomainB[_mesh->ndb++] = ib;
    }
    for(int i=0;i<_mesh->nbz;i++)
    {
        struct Boundary * ib = &_mesh->ZBoundaries[i];
        if(NULL==ib->NeE[0] || NULL==ib->NeE[1]) continue;
        _mesh->DomainB[_mesh->ndb++] = ib;
    }


    int * DomainMapB = (int *) malloc(_mesh->ndb * sizeof(int));
    int Ndb = 0;
    for(int k=0;k<_mesh->ndb;k++)
    {
        struct Boundary * ib = _mesh->DomainB[k];
        int kNe0 = (int) (ib->NeE[0] - _mesh->Elements);
        int kNe1 = (int) (ib->NeE[1] - _mesh->Elements);

        if((-1==CondNote[kNe0])||(-1==CondNote[kNe1]))
        {
            Ndb++;
            DomainMapB[k] = 0;
            if(-1!=CondNote[kNe0])
                ib->UpE[1] = ib->NeE[1];
            else if(-1!=CondNote[kNe1])
                ib->UpE[0] = ib->NeE[0];
        }
        else
            DomainMapB[k] = 1;
    }

    struct Boundary ** DomainB = (struct Boundary **) malloc(Ndb * sizeof(struct Boundary *));
    int kDomainB  = 0;
    for(int k=0;k<_mesh->ndb;k++)
    {
        if(0 == DomainMapB[k])
            DomainB[kDomainB++] = _mesh->DomainB[k];
    }

    free(_mesh->DomainB);
    _mesh->ndb = Ndb;
    _mesh->DomainB = DomainB;
    free(CondNote);
    free(DomainMapE);
}

void SearchDomainE2(struct Mesh *_mesh, struct MeshInfo * _meshinfo)
{
    /*
     * different between SearchDomainE3:
     * remove the variable CondB and not necessary to call FlushCondB
     * A instance of SearchDomianE with offset=2
     */
    _SearchDomainE(_mesh,2,_meshinfo);
}

void SearchDomainV(struct Mesh *_mesh,struct MeshInfo * _meshinfo)
{
    _SearchDomainV(_mesh,3,_meshinfo);
}