//
// Created by huacheng on 6/17/21.
//

#include "Reconstruction.h"
#include <math.h>

double CutH(double x)
{
    if(x > 0.0)
        return x*x*x;
    else
        return 0.0;
}

double CutVolume_(double d,double _p[])
{
    if(d<=0) return 0.0;
    double a = _p[0],b=_p[1],c=_p[2];
    double vstar = 0.0;
    if(b <= RecTOL*a)
    {
        if(d <= a)
            vstar = d;
        else
            vstar = a;
    } else if(c <= RecTOL*a)
    {
        if(d <= b)
            vstar = d*d;
        else if(d<=a)
            vstar = 2*b*d;
        else if(d <= a+b)
            vstar = 2*a*b - (d - a - b)*(d - a - b);
        else
            vstar = 2*a*b;
    } else
    {
        vstar = CutH(d) - CutH(d-a) - CutH(d-b) - CutH(d-c)
                + CutH(d-b-c) + CutH(d-a-c) + CutH(d - a -b) - CutH(d-a-b-c);
    }
    return vstar;
}

double CutVolume(double d,double _p[])
{
    double a = _p[0],b=_p[1],c=_p[2];
    double vstar = CutH(d) - CutH(d-a) - CutH(d-b) - CutH(d-c)
                   + CutH(d-b-c) + CutH(d-a-c) + CutH(d - a -b) - CutH(d-a-b-c);
    return vstar;
}

void Sort3(double _p[])
{
    for(int k=0;k<3;k++)
    {
        for(int j=k+1;j<3;j++)
        {
            if(_p[k] < _p[j])
            {
                double tmp = _p[k];
                _p[k] = _p[j];
                _p[j] = tmp;
            }
        }
    }
}

double CutLengthL(double _vf,double _p[])
{
    double a = _p[0],b=_p[1],c=_p[2];
    double dstar = 0.0;
    if(_vf <= 0.0)
        dstar = 0.0;
    else if(_vf <= 1.0)
        dstar = _vf*a;
    else
        dstar = a;
    return dstar + b + c;
}

double CutLengthS(double _vf,double _p[])
{
    double a = _p[0],b=_p[1],c=_p[2];
    double dstar;
    double crf1 = 0.5*b/a;
    if(_vf <= crf1)
        dstar = sqrt(2.0*a*b*_vf);
    else if(_vf <= 1.0 - crf1)
        dstar = 0.5*b + a*_vf;
    else if(_vf <= 1.0)
        dstar = a + b - sqrt(2.0*a*b*(1.0-_vf));
    else
        dstar = a + b;
    return dstar + c;
}

double CutLengthC(double _vf,double _p[])
{
    double a = _p[0],b=_p[1],c=_p[2];
    double dsL = 0.0,dsR=a+b+c,dsM=0.0,TargetV=6.0*a*b*c*_vf;
    for(int k=0;k<MaxCutCycle;k++)
    {
        dsM = 0.5*(dsL + dsR);
        double deltaV = dsM*dsM*dsM - CutH(dsM-a) - CutH(dsM-b) - CutH(dsM-c)
                        + CutH(dsM-b-c) + CutH(dsM-a-c) + CutH(dsM-a-b) - TargetV;
        if(deltaV > 0.0)
            dsR = dsM;
        else
            dsL = dsM;
    }
    return dsM;
}

double CutLength(double _vf,double _p[])
{
    double a = _p[0],b=_p[1],c=_p[2];
    if(b <= RecTOL*a)
    {
        return CutLengthL(_vf,_p);
    } else if(c <= RecTOL*b)
    {
        return CutLengthS(_vf,_p);
    } else
    {
        return CutLengthC(_vf,_p);
    }
}

double CutLength_1(double _vf,double _p[])
{
    /*
     * the error from CutLength_2 is much greater than CutLength !!!
     */
    double a = _p[0],b=_p[1],c=_p[2];
    double dstar = 0.0;
    if(b <= RecTOL*a)
    {
        dstar = _vf*a;
    } else if(c <= RecTOL*a)
    {
        double crf1 = 0.5*b/a;
        if(_vf <= crf1)
            dstar = sqrt(2.0*a*b*_vf);
        else if(_vf <= 1.0 - crf1)
            dstar = 0.5*b + a*_vf;
        else if(_vf <= 1.0)
            dstar = a + b - sqrt(2.0*a*b*(1.0-_vf));
        else
            dstar = a + b;
    } else
    {
        double dsL = 0.0,dsR=a+b+c,dsM=0.0,TargetV=6.0*a*b*c*_vf;
        for(int k=0;k<MaxCutCycle;k++)
        {
            dsM = 0.5*(dsL + dsR);
            double deltaV = dsM*dsM*dsM - CutH(dsM-a) - CutH(dsM-b) - CutH(dsM-c)
                    + CutH(dsM-b-c) + CutH(dsM-a-c) + CutH(dsM-a-b) - TargetV;
            if(fabs(deltaV) < 1e-4)
                break;
            else if(deltaV > 0.0)
                dsR = dsM;
            else
                dsL = dsM;
        }
        dstar = dsM;
    }
    return dstar;
}

double CutLength_2(double _vf,double _p[])
{
    double a = _p[0],b=_p[1],c=_p[2];
    double dstar = 0.0;
    if(a*b*c >= RecTOL)
    {
        double dsL = 0.0,dsR=a+b+c,dsM=0.0,TargetV = 6.0*a*b*c*_vf;
        for(int k=0;k<MaxCutCycle;k++)
        {
            dsM = 0.5*(dsL + dsR);
            double deltaV = dsM*dsM*dsM - CutH(dsM-a) - CutH(dsM-b) - CutH(dsM-c)
                            + CutH(dsM-b-c) + CutH(dsM-a-c) + CutH(dsM-a-b) - TargetV;
            if(fabs(deltaV) < 1e-4)
                break;
            else if(deltaV > 0.0)
                dsR = dsM;
            else
                dsL = dsM;
        }
        dstar = dsM;
    } else if(a*b >= RecTOL)
    {
        double crf1 = 0.5*b/a;
        if(_vf <= crf1)
            dstar = sqrt(2.0*a*b*_vf);
        else if(_vf <= 1.0 - crf1)
            dstar = 0.5*b + a*_vf;
        else
            dstar = a + b - sqrt(2.0*a*b*(1.0-_vf));
    } else
    {
        dstar = _vf*a;
    }
    return dstar;
}

double CutSf2(double dk,double dl,double dm,double dn)
{
    /*
     * dk-dn is sorted as
     * dk <= dl <= dm <= dn
     */
    double ResSf = 0.0;
    if(dk >= 0.0)
    {
        ResSf = 1.0;
    } else if(dn <= 0.0)
    {
        ResSf = 0.0;
    } else if(dl >= 0.0)
    {
        ResSf = 1.0 - 0.5*dk*dk/((dl-dk)*(dm-dk));
    } else if(dm <= 0.0)
    {
        ResSf = 0.5*dn*dn/((dl - dn)*(dm - dn));
    } else
    {
        ResSf = 0.5*(dm/(dm-dk) + dn/(dn-dl));
    }
    return ResSf;
}

double CutSf(double _p[])
{
    /*
     * dk-dn is sorted as
     * dk <= dl <= dm <= dn
     */
    double dk = _p[0];
    double dl = _p[1];
    double dm = _p[2];
    double dn = _p[3];
    double ResSf = 0.0;
    if(dk >= 0.0)
    {
        ResSf = 1.0;
    } else if(dn <= 0.0)
    {
        ResSf = 0.0;
    } else if(dl >= 0.0)
    {
        ResSf = 1.0 - 0.5*dk*dk/((dl-dk)*(dm-dk));
    } else if(dm <= 0.0)
    {
        ResSf = 0.5*dn*dn/((dl - dn)*(dm - dn));
    } else
    {
        ResSf = 0.5*(dm/(dm-dk) + dn/(dn-dl));
    }
    return ResSf;
}


const int bBID[BPE] = {
        24,  51, 153, 102,  15, 240
};

const int bVID[VPE] = {
        1,   2,   4,   8,
       16,  32,  64, 128
};

void InsertBoundaryBID(struct Boundary *ib)
{
    int b0=0,b1=0;

    for(int ibv=0;ibv<VPB;ibv++)
    {
        for(int jev=0;jev<VPE;jev++)
        {
            if(ib->NeE[0]->NeV[jev] == ib->NeV[ibv]) SetBit(&b0,jev);
            if(ib->NeE[1]->NeV[jev] == ib->NeV[ibv]) SetBit(&b1,jev);
        }
    }
    ib->bBID[0] = b0;
    ib->bBID[1] = b1;
    if((4!= SumBit(b0))||(4!= SumBit(b1)))
    {
        fprintf(stdout,"Error in BID\n");
        exit(0);
    }
}

int ndCompare(const void * a,const void * b)
{
    const struct Vnode * x = (const struct Vnode *) a;
    const struct Vnode * y = (const struct Vnode *) b;
    if(x->d < y->d)
        return -1;
    else
        return 1;
}

void ConstructPlane(struct Element * ie)
{
    if(ie->VOF[CutMaterial]>1.0-TOLVOF)
    {
        ie->CutTag = 1;
        ZeroEx(ie->CutRatio,VPE);
    }
    else if(ie->VOF[CutMaterial] < TOLVOF)
    {
        ie->CutTag = 0;
        CopyEx(ie->SubVolumeRatio,ie->CutRatio,VPE);
    }
    else
    {
        ie->CutTag = 2;
        double pCut[DIM] = {0.0};
        NormalizationEx(ie->GradVOF[CutMaterial],ie->Scale[Z]);
        for(int k=0;k<DIM;k++)
        {
            pCut[k] = fabs(ie->GradVOF[CutMaterial][k]*ie->Scale[k]);
        }
        Sort3(pCut);

        double dstar = CutLength(ie->VOF[CutMaterial],pCut);
        double dsum  = 0.0;
        for(int k=0;k<VPE;k++)
        {
            ie->dCut[k].vId = k;
            ie->dCut[k].d   = Dot(ie->GradVOF[CutMaterial],ie->NeV[k]->Position);
            dsum += ie->dCut[k].d;
        }
        qsort(ie->dCut,VPE,sizeof(struct Vnode),ndCompare);
        double dmax   = ie->dCut[VPE-1].d;
        dsum = 0.5*(pCut[0] + pCut[1] + pCut[2]) + dmax - dstar - 0.125*dsum;
        double negVol = 0.0;
        for(int k=0;k<VPE;k++)
        {
            ie->dCut[k].d = ie->dCut[k].d - dmax + dstar;
            unsigned short vId = ie->dCut[k].vId;
            if(ie->dCut[k].d < 0)
            {
                ie->CutRatio[vId] = CutVolume(dsum - ie->dCut[k].d,pCut);
                negVol += ie->CutRatio[vId];
            } else
            {
                ie->CutRatio[vId] = 0.0;
            }
            ie->dCut[k].vId = (1<<vId);
        }

        if(negVol > RecTOL)
        {
            for(int k=0;k<VPE;k++)
                ie->CutRatio[k] /= negVol;
        } else
        {
            CopyEx(ie->SubVolumeRatio,ie->CutRatio,VPE);
        }

    }
}

void SetVnodeCheck(struct Element * ie)
{
    if(ie->VOF[CutMaterial]>1.0-TOLVOF)
    {
        ie->CutTag = 1;
        ZeroEx(ie->CutRatio,VPE);
    }
    else if(ie->VOF[CutMaterial] < TOLVOF)
    {
        ie->CutTag = 0;
        CopyEx(ie->SubVolumeRatio,ie->CutRatio,VPE);
    }
    else
    {
        ie->CutTag = 2;
        double pCut[DIM] = {0.0};
        NormalizationEx(ie->GradVOF[CutMaterial],ie->Scale[Z]);
        for(int k=0;k<DIM;k++)
        {
            pCut[k] = fabs(ie->GradVOF[CutMaterial][k]*ie->Scale[k]);
        }
        Sort3(pCut);
        double dstar = CutLength(ie->VOF[CutMaterial],pCut);
        double dsum  = 0.0;
        CheckDouble(stdout,"CutMat",ie->VOF[CutMaterial]);
        CheckDoubleArray(stdout,"pCut array",pCut,3);
        CheckDoubleArrayE(stdout,"pCut array",pCut,3);
        CheckDouble(stdout,"DSTAR",dstar);

        double pCutAlt[3] = {565.68542,565.68542,0.0};
        CheckDouble(stdout,"DSTARAlt",CutLength(0.98,pCut));
        CheckDouble(stdout,"DSTARAlt2",CutLength(0.98,pCutAlt));
        for(int k=0;k<VPE;k++)
        {
            ie->dCut[k].vId = k;
            ie->dCut[k].d   = Dot(ie->GradVOF[CutMaterial],ie->NeV[k]->Position);
            dsum += ie->dCut[k].d;
            CheckDouble(stdout,"dk",ie->dCut[k].d);
        }
        qsort(ie->dCut,VPE,sizeof(struct Vnode),ndCompare);
        for(int k=0;k<VPE;k++)
        {
            CheckDouble(stdout,"dkSort",ie->dCut[k].d);
        }

        double dmax   = ie->dCut[VPE-1].d;
        dsum = 0.5*(pCut[0] + pCut[1] + pCut[2]) + dmax - dstar - 0.125*dsum;
        double negVol = 0.0;
        for(int k=0;k<VPE;k++)
        {
            ie->dCut[k].d = ie->dCut[k].d - dmax + dstar;
            unsigned short vId = ie->dCut[k].vId;
            if(ie->dCut[k].d < 0)
            {
                ie->CutRatio[vId] = CutVolume(dsum - ie->dCut[k].d,pCut);
                negVol += ie->CutRatio[vId];
            } else
            {
                ie->CutRatio[vId] = 0.0;
            }
            ie->dCut[k].vId = (1<<vId);
        }

        if(negVol > RecTOL)
        {
            for(int k=0;k<VPE;k++)
                ie->CutRatio[k] /= negVol;
        } else
        {
            CopyEx(ie->SubVolumeRatio,ie->CutRatio,VPE);
        }

    }
}

int dCompare(const void *a, const void *b)
{
    const double *aa = (const double *)a;
    const double *bb = (const double *)b;
    return (*aa > *bb) - (*aa < *bb);
}

double CutVOF(struct Boundary * ib, int tag)
{
    struct Element * ie = ib->NeE[tag];
    if(ie->CutTag == 0)
    {
        return 0.0;
    }
    else if(ie->CutTag == 1)
    {
        return 1.0;
    }
    else if(ie->CutTag == 2)
    {
        double tmpd[VPB] = {0.0};
        double tmpp[VPB] = {0.0};
        double tmpm[VPB] = {0.0};
        int ktmpd = 0;
        for(int k=0;k<VPE;k++)
        {
            if((ie->dCut[k].vId)&(ib->bBID[tag]))
            {
                double pDistance = Dot(ie->NeV[ReverseBit(ie->dCut[k].vId)]->DistanceL1,ie->GradVOF[CutMaterial]);
//                tmpp[ktmpd] = ie->dCut[k].d - pDistance;
                tmpp[ktmpd] = ie->dCut[k].d + pDistance;
                tmpm[ktmpd] = ie->dCut[k].d + 0.5*pDistance;
                tmpd[ktmpd++] = ie->dCut[k].d;
            }
        }
        qsort(tmpp,VPB,sizeof(double),dCompare);
        qsort(tmpm,VPB,sizeof(double),dCompare);
        double S1 = CutSf(tmpd);
        double S2 = CutSf(tmpp);
        double S3 = CutSf(tmpm);

        return S1;
//        return (S1 + S2 + S3)/3.0;
//        return (S1 + S2 + 2*S3 + sqrt(S1*S3) + sqrt(S2*S3))/6.0;
//        return (CutSf(tmpd) + CutSf(tmpp))*0.5;
//        return (S1 + S2 + sqrt(S1*S2))/3.0;
    }
}

void MapE2Vc(struct Vertex * iv)
{
    ZeroEx(iv->VOF,NMT);
    for(int j=0;j<EPV;j++)
    {
        ScalerAdditionEx(iv->VOF,iv->NeE[j]->VOF,iv->SubVolumeRatio[j],NMT);
    }
}

void CutGradVOF(struct Element * ie, int nm)
{
    for(int k=0;k<DIM;k++)
    {
        // calculate the direction of Volume function
        for (int iMat = 0; iMat < nm; iMat++)
        {
            ie->GradVOF[iMat][k] = 0.0;
            for (int i = 0; i < VPE; i++)
            {
                ie->GradVOF[iMat][k] += ie->NeV[i]->VOF[iMat] * ie->GradFactor[k][i];
            }
        }
    }
}

void CropVelocity(struct Element * ie, double MaxVel)
{
    if(MaxVel <= 0) return;
    double eVel = Length(ie->Momentum);
    if(eVel > MaxVel*ie->Mass)
    {
        for(int k=0;k<DIM;k++)
        {
            ie->Momentum[k] = ie->Momentum[k]/eVel * MaxVel*ie->Mass;
        }
    }
}





