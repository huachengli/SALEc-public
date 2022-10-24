//
// Created by huacheng on 2021/5/22.
//

#include "Parallelism.h"
#include "Error.h"

MPI_Comm SALEC_CART_COMM;


void SalecCartCommCreate(int dim[], int period[])
{
    MPI_Cart_create(MPI_COMM_WORLD,DIM,dim,period,1,&SALEC_CART_COMM);
}

int BoundaryDistance(int LongId,int nxp,int nyp,int nzp)
{
    int idx,idy,idz;
    ToShortId(&idx,&idy,&idz,LongId,nxp,nyp,nzp);
    int bDistance = (int) Min(Min(idx,nxp-1-idx),Min(Min(idy,nyp-1-idy),Min(idz,nzp-1-idz)));
    return bDistance;
}

void InitProcessInfo(struct Mesh * _mesh, struct ProcessInfo * _pinfo)
{
    int lrank,xrank,yrank,zrank,Npgx,Npgy,Npgz;
    Npgx = _pinfo->npgx;
    Npgy = _pinfo->npgy;
    Npgz = _pinfo->npgz;
    MPI_Comm_rank(SALEC_CART_COMM, &lrank);
    ToShortId(&xrank,&yrank,&zrank,lrank,Npgx,Npgy,Npgz);
    for(int ix=-1;ix<=1;ix++)
    {
        for(int iy=-1;iy<=1;iy++)
        {
            for(int iz=-1;iz<=1;iz++)
            {
                _pinfo->NeProc[ToLongId(ix+1,iy+1,iz+1,3,3,3)] =
                        ToLongId(xrank+ix,yrank+iy,zrank+iz,Npgx,Npgy,Npgz);
            }
        }
    }

    int Noffset,Nxp,Nyp,Nzp,Nsize;
    Noffset = _pinfo->noffset = 2;
    Nxp  = _pinfo->nwp[lrank][X] + 2*Noffset;
    Nyp  = _pinfo->nwp[lrank][Y] + 2*Noffset;
    Nzp  = _pinfo->nwp[lrank][Z] + 2*Noffset;
    Nsize = Nxp*Nyp*Nzp;
    struct RequireMap * tRMap = (struct RequireMap *) malloc(Nsize*sizeof(struct RequireMap));
    int bIdx,bIdy,bIdz;
    bIdx = _pinfo->BaseId[lrank][X];
    bIdy = _pinfo->BaseId[lrank][Y];
    bIdz = _pinfo->BaseId[lrank][Z];


    for(int k=0;k<Nsize;k++)
    {
        tRMap[k].npcs = 0;
    }

    for(int k=0;k<DIM3;k++)
    {
        int trank = _pinfo->NeProc[k];
        if(-1 == trank) continue;
        if(lrank == trank) continue;
        int tNxp,tNyp,tNzp,tIdx,tIdy,tIdz;
        tNxp = _pinfo->nwp[trank][X] + 2*Noffset;
        tNyp = _pinfo->nwp[trank][Y] + 2*Noffset;
        tNzp = _pinfo->nwp[trank][Z] + 2*Noffset;
        tIdx = _pinfo->BaseId[trank][X];
        tIdy = _pinfo->BaseId[trank][Y];
        tIdz = _pinfo->BaseId[trank][Z];

        for(int j=0;j<Nsize;j++)
        {
            int Idx,Idy,Idz;
            ToShortId(&Idx,&Idy,&Idz,j,Nxp,Nyp,Nzp);
            int destId = ToLongId(Idx + bIdx - tIdx,Idy + bIdy - tIdy,Idz + bIdz - tIdz,tNxp,tNyp,tNzp);
            if(-1 == destId) continue;

            if(BoundaryDistance(j,Nxp,Nyp,Nzp) >= Noffset)
            {
                /*
                 * In this element, information need to send their neighbour in other process.
                 */
                tRMap[j].prcs[tRMap[j].npcs] = trank;
                tRMap[j].dest[tRMap[j].npcs] = destId;
                tRMap[j].dist[tRMap[j].npcs] = BoundaryDistance(destId,tNxp,tNyp,tNzp);
                tRMap[j].npcs++;

            }
            else if(BoundaryDistance(destId,tNxp,tNyp,tNzp) >= Noffset)
            {
                /*
                 * Receive from other process, the Tag is local Id.
                 */
                tRMap[j].prcs[tRMap[j].npcs] = trank;
                tRMap[j].dest[tRMap[j].npcs] = destId;
                tRMap[j].dist[tRMap[j].npcs] = BoundaryDistance(j,Nxp,Nyp,Nzp);
//                tRMap[j].dest[tRMap[j].npcs] = j;
                tRMap[j].npcs++;
            }
        }
    }

    int Nstb=0,Nrtb=0;
    for(int k=0;k<Nsize;k++)
    {
        if(tRMap[k].npcs == 0) continue;
        if(BoundaryDistance(k,Nxp,Nyp,Nzp) < Noffset)
            Nrtb ++;
        else
            Nstb += tRMap[k].npcs;
    }

    struct RouteTableE * tSTB = (struct RouteTableE *) malloc(Nstb*sizeof(struct RouteTableE));
    struct RouteTableE * tRTB = (struct RouteTableE *) malloc(Nrtb*sizeof(struct RouteTableE));
    int kSTB=0,kRTB=0;
    for(int k=0;k<Nsize;k++)
    {
        if(tRMap[k].npcs == 0) continue;
        if(BoundaryDistance(k,Nxp,Nyp,Nzp) < Noffset)
        {
//            tRTB[kRTB].priority = Noffset - 1 - BoundaryDistance(k,Nxp,Nyp,Nzp);
            tRTB[kRTB].priority = Noffset - 1 - tRMap[k].dist[0];
            tRTB[kRTB].destProc = tRMap[k].prcs[0];
//            tRTB[kRTB].destTag  = tRMap[k].dest[0];
            /*
             * RecvTag is local Id --> iteration variable k.
             */
            tRTB[kRTB].destTag = k;
            tRTB[kRTB].local = _mesh->Elements + k;
            kRTB ++;
            if(tRMap[k].npcs > 1)
            {
                exit(0);
            }
        }
        else
        {
            for(int j=0;j<tRMap[k].npcs;j++)
            {
//                tSTB[kSTB].priority = BoundaryDistance(k,Nxp,Nyp,Nzp) - Noffset;
                tSTB[kSTB].priority = Noffset - 1 - tRMap[k].dist[j];
                tSTB[kSTB].local = _mesh->Elements + k;
                tSTB[kSTB].destTag  = tRMap[k].dest[j];
                tSTB[kSTB].destProc = tRMap[k].prcs[j];
                kSTB++;
            }
        }

    }

    _pinfo->RMAP = tRMap;
    _pinfo->nrmap = Nsize;

    _pinfo->RTB  = tRTB;
    _pinfo->STB  = tSTB;
    _pinfo->nrtb = Nrtb;
    _pinfo->nstb = Nstb;
    _pinfo->rank = lrank;
}


void ParaSearchDomainE(struct Mesh *_mesh, struct ProcessInfo *_pinfo)
{
    int offset = _pinfo->noffset;
    int * DomainMapE = (int *) malloc(_mesh->ne * sizeof(int));
    int * CondNote   = (int *) malloc(_mesh->ne * sizeof(int));
    int * CondNoteA   = (int *) malloc(_mesh->ne * sizeof(int));
    int Nxp=_mesh->nxp-1, Nyp=_mesh->nyp-1, Nzp = _mesh->nzp-1, Nce=0, Nde=0;
    int Nxp2 = _pinfo->Nxx,Nyp2 = _pinfo->Nyy, Nzp2 = _pinfo->Nzz, lrank= _pinfo->rank;
    int bIdx,bIdy,bIdz;
    bIdx = _pinfo->BaseId[lrank][X];
    bIdy = _pinfo->BaseId[lrank][Y];
    bIdz = _pinfo->BaseId[lrank][Z];

    for(int Idx=0;Idx<Nxp;Idx++)
    {
        for(int Idy=0;Idy<Nyp;Idy++)
        {
            for(int Idz=0;Idz<Nzp;Idz++)
            {
                int * VCondNoteA = CondNoteA + ToLongId(Idx,Idy,Idz,Nxp,Nyp,Nzp);
                int * VCondNote  = CondNote  + ToLongId(Idx,Idy,Idz,Nxp,Nyp,Nzp);
                _searchdomain(VCondNoteA,Idx+bIdx,Idy+bIdy,Idz+bIdz,Nxp2,Nyp2,Nzp2,offset);
                _searchdomain(VCondNote,Idx,Idy,Idz,Nxp,Nyp,Nzp,offset);
                if(*VCondNoteA >= 0 && *VCondNote >=0)
                {
                    Nce++;
                    int tNote = 0,tIdx=0,tIdy=0,tIdz=0;
                    _searchdomain2(&tNote,Idx+bIdx,Idy+bIdy,Idz+bIdz,Nxp2,Nyp2,Nzp2,offset);
                    if(GetBit(tNote,X))       tIdx = 1;
                    if(GetBit(tNote,X + DIM)) tIdx = -1;
                    if(GetBit(tNote,Y))       tIdy = 1;
                    if(GetBit(tNote,Y + DIM)) tIdy = -1;
                    if(GetBit(tNote,Z))       tIdz = 1;
                    if(GetBit(tNote,Z + DIM)) tIdz = -1;
                    *VCondNote = ToLongId(Idx+tIdx,Idy+tIdy,Idz+tIdz,Nxp,Nyp,Nzp);
                }

                if(*VCondNote == -1)
                    Nde++;
            }
        }
    }

    struct Element ** DomainE = (struct Element **) malloc(Nde* sizeof(struct Element *));
    struct ConditionE * CondE = (struct ConditionE *) malloc(Nce * sizeof(struct ConditionE));
    int kDomainE=0, kCondE=0;
    for(int k=0;k<_mesh->ne;k++)
    {
        if(CondNote[k] >= 0 && CondNoteA[k]>=0)
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
            if(-1!=CondNote[kNe1] && -1!=CondNoteA[kNe1])
                ib->UpE[1] = ib->NeE[1];
            else if(-1!=CondNote[kNe0] && -1!=CondNoteA[kNe0])
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


double * GetPhi(struct Element * ie)
{
    return &(ie->PhiL3[0][0]);
}
double * GetPsi(struct Element * ie)
{
    return &(ie->PsiL3[0]);
}
double * GetMass(struct Element * ie)
{
    return &(ie->Mass);
}
double * GetPres(struct Element * ie)
{
    return &(ie->Pressure);
}
double * GetPABV(struct Element * ie)
{
    return &(ie->ArtificialPressure);
}
double * GetCs(struct Element * ie)
{
    return &(ie->Csound);
}
double * GetVOF(struct Element * ie)
{
    return &(ie->VOF[0]);
}
double * GetMoment(struct Element * ie)
{
    return &(ie->Momentum[0]);
}
double * GetStress(struct Element * ie)
{
    return &(ie->PsiL3[1]);
}
double * GetCourant(struct Element * ie)
{
    return &(ie->Courant);
}
double * GetGradVOF(struct Element * ie)
{
    return &(ie->GradVOF[0][0]);
}

DataFunc IdToData[MaxField] = {
        GetPsi,
        GetPhi,
        GetMass,
        GetPres,
        GetPABV,
        GetCs,
        GetVOF,
        GetMoment,
        GetStress,
        GetCourant,
        GetGradVOF
};

int CacheLength[MaxField] = {
        NPSI, /*Psi*/
        NPHI*NMT, /*Phi*/
        1,    /*Mass*/
        1,    /*Pres*/
        1,    /*PresABV*/
        1,    /*Sound Spees*/
        NMT,  /*VOF*/
        DIM,  /*Moment*/
        6,    /*Stress*/
        1,    /*Courant*/
        DIM*NMT /*GradVOF*/
};

int SwitchR0Table[MaxField] = {
        0,    /* Psi */
        0,    /* Phi */
        0,    /* Mass */
        0,    /* Pres */
        0,    /* PresABV */
        0,    /* Sound Spees */
        1,    /* VOF */
        0,    /* Moment */
        0,    /* Stress */
        0,    /* Courant */
        0     /* GradVOF */
};

int SwitchR0(int _control,int Long,int Short)
{
    if(SwitchR0Table[_control])
        return Long;
    else
        return Short;
}


int rtTagCompare(const void *x, const void *y)
{
    struct RouteTableE * a = (struct RouteTableE *) x;
    struct RouteTableE * b = (struct RouteTableE *) y;
    if(a->destProc != b->destProc)
        return a->destProc - b->destProc;
    else if(a->priority != b->priority)
        return a->priority - b->priority;
    else
        return a->destTag - b->destTag;
}


int ParaPackRecvInit(struct ProcessInfo *_pinfo)
{
    qsort(_pinfo->RTB,_pinfo->nrtb, sizeof(struct RouteTableE), rtTagCompare);
    int RecvRoot[DIM3*2],Nrr=1,RecvRootR0[DIM3*2]={0};
    RecvRoot[0]   = 0;
    RecvRootR0[0] = 1;
    for(int k=1;k<_pinfo->nrtb;k++)
    {
        if(_pinfo->RTB[k-1].destProc!=_pinfo->RTB[k].destProc)
        {
            RecvRoot[Nrr++] = k;
        }
        if(_pinfo->RTB[k].priority==0)
        {
            RecvRootR0[Nrr-1]++;
        }
    }
    RecvRoot[Nrr] = _pinfo->nrtb;


    MPI_Request ** RecvReq = (MPI_Request **) malloc(MaxField*sizeof(MPI_Request *));
//    int ** PorityRecv = (int **) malloc(MaxField*sizeof(int *));
    for(int k=0;k<MaxField;k++)
    {
        RecvReq[k] = (MPI_Request *) malloc(sizeof(MPI_Request)*Nrr);
//        PorityRecv[k] = (int *) malloc(sizeof(int)*Nrr);
        int SegLen = CacheLength[k];
        for(int j=0;j<Nrr;j++)
        {
            int NumRecv  = SwitchR0(k,RecvRoot[j+1] - RecvRoot[j],RecvRootR0[j]);
            int DestProc = _pinfo->RTB[RecvRoot[j]].destProc;
            int DestTag  = k;
            MPI_Recv_init(_pinfo->PackRecvCache[k]+SegLen*RecvRoot[j],NumRecv*SegLen, MPI_DOUBLE,
                          DestProc, DestTag, SALEC_CART_COMM, RecvReq[k] + j);
        }
    }

    _pinfo->nrr = Nrr;
    _pinfo->RecvReq = RecvReq;
    return Nrr;
}

int ParaPackSendInit(struct ProcessInfo *_pinfo)
{
    qsort(_pinfo->STB,_pinfo->nstb, sizeof(struct RouteTableE), rtTagCompare);
    int SendRoot[DIM3*2], Nsr=1, SendRootR0[DIM3*2]={0};
    SendRoot[0]   = 0;
    SendRootR0[0] = 1;
    for(int k=1;k<_pinfo->nstb;k++)
    {
        if(_pinfo->STB[k-1].destProc!=_pinfo->STB[k].destProc)
        {
            SendRoot[Nsr++] = k;
        }
        if(_pinfo->STB[k].priority==0)
        {
            SendRootR0[Nsr-1] ++;
        }
    }
    SendRoot[Nsr] = _pinfo->nstb;


    MPI_Request ** SendReq = (MPI_Request **) malloc(MaxField*sizeof(MPI_Request *));
    for(int k=0;k<MaxField;k++)
    {
        SendReq[k] = (MPI_Request *) malloc(sizeof(MPI_Request)*Nsr);
        int SegLen = CacheLength[k];
        for(int j=0;j<Nsr;j++)
        {
            int NumSend = SwitchR0(k,SendRoot[j+1] - SendRoot[j],SendRootR0[j]);
            int DestProc = _pinfo->STB[SendRoot[j]].destProc;
            int DestTag  = k;
            MPI_Send_init(_pinfo->PackSendCache[k]+SendRoot[j]*SegLen,NumSend*SegLen, MPI_DOUBLE,
                          DestProc, DestTag, SALEC_CART_COMM, SendReq[k] + j);
        }
    }

    _pinfo->nsr     = Nsr;
    _pinfo->SendReq = SendReq;
    return Nsr;
}

int ParaPackInit(struct ProcessInfo *_pinfo)
{
    int SumLength = 0;
    for(int k=0;k<MaxField;k++) SumLength+= CacheLength[k];

    _pinfo->PackSendCache = (double **) malloc(MaxField*sizeof(double *));
    _pinfo->PackRecvCache = (double **) malloc(MaxField*sizeof(double *));
    double * TmpSendCache = (double *) malloc(SumLength*_pinfo->nstb*sizeof(double));
    double * TmpRecvCache = (double *) malloc(SumLength*_pinfo->nrtb*sizeof(double));
    for(int k=0;k<MaxField;k++)
    {
        _pinfo->PackSendCache[k] = TmpSendCache;
        _pinfo->PackRecvCache[k] = TmpRecvCache;
        TmpSendCache += _pinfo->nstb*CacheLength[k];
        TmpRecvCache += _pinfo->nrtb*CacheLength[k];
    }
}

int ParaPackId(struct ProcessInfo *_pinfo, int _control)
{
    int SegLen = CacheLength[_control];
    for(int k=0;k<_pinfo->nstb;k++)
    {
        if(_pinfo->STB[k].priority <= SwitchR0Table[_control])
        CopyEx(IdToData[_control](_pinfo->STB[k].local),_pinfo->PackSendCache[_control]+k*SegLen,SegLen);
    }
    return SegLen;
}

int ParaUnpackId(struct ProcessInfo *_pinfo, int _control)
{
    int SegLen = CacheLength[_control];
    for(int k=0;k<_pinfo->nrtb;k++)
    {
        if(_pinfo->RTB[k].priority <= SwitchR0Table[_control])
        CopyEx(_pinfo->PackRecvCache[_control]+k*SegLen,IdToData[_control](_pinfo->RTB[k].local),SegLen);
    }
    return SegLen;
}

int ParaPackBegin(struct ProcessInfo *_pinfo, int _control)
{
    ParaPackId(_pinfo,_control);
    MPI_Startall(_pinfo->nrr,_pinfo->RecvReq[_control]);
    MPI_Startall(_pinfo->nsr,_pinfo->SendReq[_control]);
}

int ParaPackEnd(struct ProcessInfo *_pinfo, int _control)
{
    MPI_Waitall(_pinfo->nrr,_pinfo->RecvReq[_control],MPI_STATUSES_IGNORE);
    MPI_Waitall(_pinfo->nsr,_pinfo->SendReq[_control],MPI_STATUSES_IGNORE);
    ParaUnpackId(_pinfo,_control);
}

int ParaPackEndAll(struct ProcessInfo *_pinfo)
{
    for(int k=0;k<MaxField;k++)
    {
        MPI_Waitall(_pinfo->nrr,_pinfo->RecvReq[k],MPI_STATUSES_IGNORE);
        MPI_Waitall(_pinfo->nsr,_pinfo->SendReq[k],MPI_STATUSES_IGNORE);
        ParaUnpackId(_pinfo,k);
    }
}


int ParaRecvSendInit(struct ProcessInfo *_pinfo)
{
    int Nrecv=_pinfo->nrtb,Nsend=_pinfo->nstb;
    MPI_Request * Send_PsiL3   = (MPI_Request *) malloc(Nsend * sizeof(MPI_Request));
    MPI_Request * Send_PhiL3   = (MPI_Request *) malloc(Nsend * sizeof(MPI_Request));
    MPI_Request * Send_Mass    = (MPI_Request *) malloc(Nsend * sizeof(MPI_Request));
    MPI_Request * Send_Pres    = (MPI_Request *) malloc(Nsend * sizeof(MPI_Request));
    MPI_Request * Send_PresABV = (MPI_Request *) malloc(Nsend * sizeof(MPI_Request));
    MPI_Request * Send_Csound  = (MPI_Request *) malloc(Nsend * sizeof(MPI_Request));
    MPI_Request * Send_VOF     = (MPI_Request *) malloc(Nsend * sizeof(MPI_Request));
    MPI_Request * Send_Momentum= (MPI_Request *) malloc(Nsend * sizeof(MPI_Request));
    MPI_Request * Send_Stress  = (MPI_Request *) malloc(Nsend * sizeof(MPI_Request));
    MPI_Request * Send_Courant = (MPI_Request *) malloc(Nsend * sizeof(MPI_Request));
    for(int k=0;k<Nsend;k++)
    {
        struct Element * LocalElement = _pinfo->STB[k].local;
        int BaseTag  = _pinfo->STB[k].destTag;
        int DestProc = _pinfo->STB[k].destProc;
        MPI_Send_init(LocalElement->PsiL3, NPSI, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField + IdPsiL3, SALEC_CART_COMM, Send_PsiL3 + k);
        MPI_Send_init(&(LocalElement->PhiL3[0][0]),NPHI*NMT, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField + IdPhiL3, SALEC_CART_COMM, Send_PhiL3 + k);
        MPI_Send_init(&(LocalElement->Mass), 1, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField + IdMass, SALEC_CART_COMM, Send_Mass + k);
        MPI_Send_init(&(LocalElement->Pressure), 1, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField + IdPres, SALEC_CART_COMM, Send_Pres + k);
        MPI_Send_init(&(LocalElement->ArtificialPressure), 1, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField + IdPresABV, SALEC_CART_COMM, Send_PresABV + k);
        MPI_Send_init(&(LocalElement->Csound), 1, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField + IdCsound, SALEC_CART_COMM, Send_Csound + k);
        MPI_Send_init(LocalElement->VOF, NMT, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField+IdVOF, SALEC_CART_COMM, Send_VOF + k);
        MPI_Send_init(LocalElement->Momentum, DIM, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField+IdMomentum, SALEC_CART_COMM, Send_Momentum + k);
        MPI_Send_init(LocalElement->PsiL3+1, 6, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField+IdStress, SALEC_CART_COMM, Send_Stress + k);
        MPI_Send_init(&(LocalElement->Courant), 1, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField+IdCourant, SALEC_CART_COMM, Send_Courant + k);
    }

    MPI_Request * Recv_PsiL3   = (MPI_Request *) malloc(Nrecv * sizeof(MPI_Request));
    MPI_Request * Recv_PhiL3   = (MPI_Request *) malloc(Nrecv * sizeof(MPI_Request));
    MPI_Request * Recv_Mass    = (MPI_Request *) malloc(Nrecv * sizeof(MPI_Request));
    MPI_Request * Recv_Pres    = (MPI_Request *) malloc(Nrecv * sizeof(MPI_Request));
    MPI_Request * Recv_PresABV = (MPI_Request *) malloc(Nrecv * sizeof(MPI_Request));
    MPI_Request * Recv_Csound  = (MPI_Request *) malloc(Nrecv * sizeof(MPI_Request));
    MPI_Request * Recv_VOF     = (MPI_Request *) malloc(Nrecv * sizeof(MPI_Request));
    MPI_Request * Recv_Momentum= (MPI_Request *) malloc(Nrecv * sizeof(MPI_Request));
    MPI_Request * Recv_Stress  = (MPI_Request *) malloc(Nrecv * sizeof(MPI_Request));
    MPI_Request * Recv_Courant = (MPI_Request *) malloc(Nrecv * sizeof(MPI_Request));
    for(int k=0;k<Nrecv;k++)
    {
        struct Element * LocalElement = _pinfo->RTB[k].local;
        int BaseTag  = _pinfo->RTB[k].destTag;
        int DestProc = _pinfo->RTB[k].destProc;
        MPI_Recv_init(LocalElement->PsiL3, NPSI, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField + IdPsiL3, SALEC_CART_COMM, Recv_PsiL3 + k);
        MPI_Recv_init(&(LocalElement->PhiL3[0][0]),NPHI*NMT, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField + IdPhiL3, SALEC_CART_COMM, Recv_PhiL3 + k);
        MPI_Recv_init(&(LocalElement->Mass), 1, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField + IdMass, SALEC_CART_COMM, Recv_Mass + k);
        MPI_Recv_init(&(LocalElement->Pressure), 1, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField + IdPres, SALEC_CART_COMM, Recv_Pres + k);
        MPI_Recv_init(&(LocalElement->ArtificialPressure), 1, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField + IdPresABV, SALEC_CART_COMM, Recv_PresABV + k);
        MPI_Recv_init(&(LocalElement->Csound), 1, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField + IdCsound, SALEC_CART_COMM, Recv_Csound + k);
        MPI_Recv_init(LocalElement->VOF, NMT, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField+IdVOF, SALEC_CART_COMM, Recv_VOF + k);
        MPI_Recv_init(LocalElement->Momentum, DIM, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField+IdMomentum, SALEC_CART_COMM, Recv_Momentum + k);
        MPI_Recv_init(LocalElement->PsiL3+1, 6, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField+IdStress, SALEC_CART_COMM, Recv_Stress + k);
        MPI_Recv_init(&(LocalElement->Courant), 1, MPI_DOUBLE, DestProc,
                      BaseTag*MaxField+IdCourant, SALEC_CART_COMM, Recv_Courant + k);
    }

    _pinfo->SendPsi     = Send_PsiL3;
    _pinfo->SendPhi     = Send_PhiL3;
    _pinfo->SendMass    = Send_Mass;
    _pinfo->SendPres    = Send_Pres;
    _pinfo->SendPABV    = Send_PresABV;
    _pinfo->SendCsound  = Send_Csound;
    _pinfo->SendVOF     = Send_VOF;
    _pinfo->SendMoment  = Send_Momentum;
    _pinfo->SendStress  = Send_Stress;
    _pinfo->SendCourant = Send_Courant;

    _pinfo->RecvPsi     = Recv_PsiL3;
    _pinfo->RecvPhi     = Recv_PhiL3;
    _pinfo->RecvMass    = Recv_Mass;
    _pinfo->RecvPres    = Recv_Pres;
    _pinfo->RecvPABV    = Recv_PresABV;
    _pinfo->RecvCsound  = Recv_Csound;
    _pinfo->RecvVOF     = Recv_VOF;
    _pinfo->RecvMoment  = Recv_Momentum;
    _pinfo->RecvStress  = Recv_Stress;
    _pinfo->RecvCourant = Recv_Courant;

    _pinfo->NumRecv = 0;
    _pinfo->NumSend = 0;
}

void IdToReq(struct ProcessInfo * _pinfo, int _control, MPI_Request ** NeedSend, MPI_Request ** NeedRecv)
{
    /*switch(_control)
    {
        case IdMass:
            NeedRecv = _pinfo->RecvMass;
            NeedSend = _pinfo->SendMass;
            break;
        case IdPres:
            NeedRecv = _pinfo->RecvPres;
            NeedSend = _pinfo->SendPres;
            break;
        case IdPresABV:
            NeedRecv = _pinfo->RecvPABV;
            NeedSend = _pinfo->SendPABV;
            break;
        case IdPsiL3:
            NeedRecv = _pinfo->RecvPsi;
            NeedSend = _pinfo->SendPsi;
            break;
        case IdPhiL3:
            NeedRecv = _pinfo->RecvPhi;
            NeedSend = _pinfo->SendPhi;
            break;
        case IdMomentum:
            NeedRecv = _pinfo->RecvMoment;
            NeedSend = _pinfo->SendMoment;
            break;
        case IdCsound:
            NeedRecv = _pinfo->RecvCsound;
            NeedSend = _pinfo->SendCsound;
            break;
        case IdVOF:
            NeedRecv = _pinfo->RecvVOF;
            NeedSend = _pinfo->SendVOF;
            break;
        case IdStress:
            NeedRecv = _pinfo->RecvStress;
            NeedSend = _pinfo->SendStress;
            break;
        default:
            NeedRecv = NULL;
            NeedSend = NULL;
            break;
    }*/

    switch(_control)
   {
       case IdMass:
           *NeedRecv = _pinfo->RecvMass;
           *NeedSend = _pinfo->SendMass;
           break;
       case IdPres:
           *NeedRecv = _pinfo->RecvPres;
           *NeedSend = _pinfo->SendPres;
           break;
       case IdPresABV:
           *NeedRecv = _pinfo->RecvPABV;
           *NeedSend = _pinfo->SendPABV;
           break;
       case IdPsiL3:
           *NeedRecv = _pinfo->RecvPsi;
           *NeedSend = _pinfo->SendPsi;
           break;
       case IdPhiL3:
           *NeedRecv = _pinfo->RecvPhi;
           *NeedSend = _pinfo->SendPhi;
           break;
       case IdMomentum:
           *NeedRecv = _pinfo->RecvMoment;
           *NeedSend = _pinfo->SendMoment;
           break;
       case IdCsound:
           *NeedRecv = _pinfo->RecvCsound;
           *NeedSend = _pinfo->SendCsound;
           break;
       case IdVOF:
           *NeedRecv = _pinfo->RecvVOF;
           *NeedSend = _pinfo->SendVOF;
           break;
       case IdStress:
           *NeedRecv = _pinfo->RecvStress;
           *NeedSend = _pinfo->SendStress;
           break;
       case IdCourant:
           *NeedRecv = _pinfo->RecvCourant;
           *NeedSend = _pinfo->SendCourant;
           break;
       default:
           *NeedRecv = NULL;
           *NeedSend = NULL;
           break;
   }

}


void ParaBegin(struct ProcessInfo *_pinfo, int _control)
{
    int Nrecv=_pinfo->nrtb,Nsend=_pinfo->nstb;
    MPI_Request * NeedSend;
    MPI_Request * NeedRecv;

    IdToReq(_pinfo,_control,&NeedSend,&NeedRecv);


    if((NULL == NeedSend)||(NULL==NeedRecv))
    {
        fprintf(stdout,"[ParaBegin] Error in NeedRecv/NeedSend with %d!\n",_control);
        exit(0);
    }

    _pinfo->OnRecv[_pinfo->NumRecv++] = NeedRecv;
    _pinfo->OnSend[_pinfo->NumSend++] = NeedSend;
#if COMM_ON
    MPI_Startall(Nrecv,NeedRecv);
    MPI_Startall(Nsend,NeedSend);
#endif
}

void ParaEndAll(struct ProcessInfo *_pinfo)
{
#if COMM_ON
    for(int k=0;k<_pinfo->NumRecv;k++)
        MPI_Waitall(_pinfo->nrtb,_pinfo->OnRecv[k],MPI_STATUSES_IGNORE);
    for(int k=0;k<_pinfo->NumSend;k++)
        MPI_Waitall(_pinfo->nstb,_pinfo->OnSend[k],MPI_STATUSES_IGNORE);
#endif
    _pinfo->NumSend = _pinfo->NumRecv = 0;
}

void ParaEnd(struct ProcessInfo *_pinfo, int _control)
{
    int Nrecv=_pinfo->nrtb,Nsend=_pinfo->nstb;
    MPI_Request * NeedSend;
    MPI_Request * NeedRecv;
    IdToReq(_pinfo,_control,&NeedSend,&NeedRecv);
    if((NULL == NeedSend)||(NULL==NeedRecv))
    {
        fprintf(stdout,"Error in NeedRecv/NeedSend!\n");
        exit(0);
    }
#if COMM_ON
    MPI_Waitall(Nrecv,NeedRecv,MPI_STATUSES_IGNORE);
    MPI_Waitall(Nsend,NeedSend,MPI_STATUSES_IGNORE);
#endif
}


void PIFClean(struct ProcessInfo *_pinfo)
{
    free(_pinfo->RMAP);
    for(int k=0;k<_pinfo->npgx*_pinfo->npgy*_pinfo->npgz;k++)
    {
        free(_pinfo->BaseId[k]);
        free(_pinfo->nwp[k]);
    }
    free(_pinfo->BaseId);
    free(_pinfo->nwp);
}

void ParaSearchDomainV(struct Mesh *_mesh, struct ProcessInfo *_pinfo)
{
    int * DistNote = (int *) malloc(_mesh->nv* sizeof(int));
    int Nxp = _mesh->nxp, Nyp=_mesh->nyp, Nzp = _mesh->nzp,Ndv1=0,Ndv2=0,Ndv3=0;
    for(int Idx=0;Idx<Nxp;Idx++)
    {
        for(int Idy=0;Idy<Nyp;Idy++)
        {
            for(int Idz=0;Idz<Nzp;Idz++)
            {
                int LongId = ToLongId(Idx,Idy,Idz,Nxp,Nyp,Nzp);
                int * DistNoteV = DistNote + LongId;
                *DistNoteV = BoundaryDistance(LongId,Nxp,Nyp,Nzp);

                if(*DistNoteV >= _pinfo->noffset + 1)
                    Ndv1 ++;
                else if(*DistNoteV == _pinfo->noffset)
                    Ndv2 ++;
                else
                    Ndv3 ++;
            }
        }
    }

    struct Vertex **  DomainV = (struct Vertex **)  malloc((Ndv1 + Ndv2) * sizeof(struct Veretx *));
    int kNdv1 = 0, kNdv2 = Ndv1;
    for(int k=0;k<_mesh->nv;k++)
    {
        if(DistNote[k] >= _pinfo->noffset + 1)
            DomainV[kNdv1++] = _mesh->Vertexes + k;
        else if(DistNote[k] == _pinfo->noffset)
            DomainV[kNdv2++] = _mesh->Vertexes + k;
    }

    _mesh->DomainV = DomainV;
    _mesh->ndv1 = Ndv1;
    _mesh->ndv  = Ndv1 + Ndv2;
    free(DistNote);
}

void ParaSortDomainE(struct Mesh *_mesh, struct ProcessInfo *_pinfo)
{
    int nde1 = _mesh->nde-1, nde0=0;
    int Nxp=_mesh->nxp-1, Nyp=_mesh->nyp-1, Nzp = _mesh->nzp-1;

    while(nde0 <= nde1)
    {
        int LongId = (int)(_mesh->DomainE[nde0] - _mesh->Elements);
        int bDistance = BoundaryDistance(LongId,Nxp,Nyp,Nzp);
        if(bDistance >= 2*_pinfo->noffset)
            nde0++;
        else if(bDistance <= _pinfo->noffset -1)
        {
            fprintf(stdout,"unexpected Boundary distance in DomainE!");
            exit(0);
        } else
        {
            struct Element * tmpE = _mesh->DomainE[nde0];
            _mesh->DomainE[nde0]  = _mesh->DomainE[nde1];
            _mesh->DomainE[nde1]  = tmpE;
            nde1--;
        }
    }
    _mesh->nde1 = nde1 + 1;
}

int ParaPackRecvSendInit(struct ProcessInfo *_pinfo)
{
    ParaPackInit(_pinfo);
    ParaPackSendInit(_pinfo);
    ParaPackRecvInit(_pinfo);
}

int PostPackEnd(struct ProcessInfo *_pinfo, int _control)
{
    if(IdGradVOF == _control)
    {
        for(int k=0;k<_pinfo->nrtb;k++)
        {
            ConstructPlane(_pinfo->RTB[k].local);
        }
    }
}

void ParaDebug()
{
    int debugvalue=1;
    while(debugvalue)
        sleep(5);
}

void InfoClean(struct MeshInfo * _minfo, struct ProcessInfo *_pinfo)
{
    for (int k = 0; k < _minfo->npz; k++)
    {
        for (int j = 0; j < _minfo->npy; j++)
        {
            for (int i = 0; i < _minfo->npx; i++)
            {
                free(_minfo->V[k][j][i]);
            }
            free(_minfo->V[k][j]);
        }
        free(_minfo->V[k]);
    }
    free(_minfo->V);

    int Npg = _pinfo->npgx * _pinfo->npgy * _pinfo->npgz;
    for(int k=0;k<Npg;k++)
    {
        free(_pinfo->nwp[k]);
        free(_pinfo->BaseId[k]);
    }
    free(_pinfo->nwp);
    free(_pinfo->BaseId);

    free(_pinfo->RTB);
    free(_pinfo->STB);
    free(_pinfo->RMAP);

    free(_pinfo->PackSendCache[0]);
    free(_pinfo->PackRecvCache[0]);

    // uncompleted
    // most memory will be free automatically ...
}