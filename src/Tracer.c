//
// Created by huacheng on 9/14/21.
//

#include "Tracer.h"

int bSearch(double _target,const double * _list,int _len)
{
    int _lt = 0,_rt = _len-1;
    if(_target<_list[_lt])
        return -3;
    else if(_list[_rt]<_target)
        return -1;

    while(_lt+1<_rt)
    {
        int _mt = (_lt+_rt)/2;
        if(_target <= _list[_mt])
            _rt = _mt;
        else
            _lt = _mt;
    }
    return _lt;
}

void NeEDetect(Tracer * _self, TracerInfo * _tinfo,struct Mesh * _mesh)
{
    int ShortId[DIM];
    ShortId[X] = bSearch(_self->Position[X],_tinfo->Coordinate[X],_mesh->nxp);
    ShortId[Y] = bSearch(_self->Position[Y],_tinfo->Coordinate[Y],_mesh->nyp);
    ShortId[Z] = bSearch(_self->Position[Z],_tinfo->Coordinate[Z],_mesh->nzp);
    if(ShortId[X]<_mesh->noffset)
    {
        ShortId[X] = -3;
    } else if(ShortId[X]>_mesh->nxp-2-_mesh->noffset)
    {
        ShortId[X] = -1;
    }

    if(ShortId[Y]<_mesh->noffset)
    {
        ShortId[Y] = -3;
    } else if(ShortId[Y]>_mesh->nyp-2-_mesh->noffset)
    {
        ShortId[Y] = -1;
    }

    if(ShortId[Z]<_mesh->noffset)
    {
        ShortId[Z] = -3;
    } else if(ShortId[Z]>_mesh->nzp-2-_mesh->noffset)
    {
        ShortId[Z] = -1;
    }

    if(ShortId[X]>=0 && ShortId[Y]>=0 && ShortId[Z]>=0)
    {
        _self->NeE = ToLongId(ShortId[X],ShortId[Y],ShortId[Z],_mesh->nxp-1,_mesh->nyp-1,_mesh->nzp-1);
    } else
    {
        ShortId[X] = (ShortId[X]<0) ? (ShortId[X]+3) : 1;
        ShortId[Y] = (ShortId[Y]<0) ? (ShortId[Y]+3) : 1;
        ShortId[Z] = (ShortId[Z]<0) ? (ShortId[Z]+3) : 1;
        _self->NeE = ToLongId(ShortId[X],ShortId[Y],ShortId[Z],3,3,3) - DIM3;
    }
}


void TracerInfoExchangeBegin(TracerInfo * _tinfo, struct ProcessInfo * _pinfo, struct Mesh * _mesh)
{
    if(strcasecmp(_tinfo->Type,"NONE") == 0) return;
    for(int k=0;k<DIM3;k++) _tinfo->STI[k] = 0;

    Tracer * start = _tinfo->head;
    while(NULL!=start)
    {
        NeEDetect(start,_tinfo,_mesh);
        if(start->NeE < 0)
        {
            _tinfo->STI[DIM3+start->NeE] ++;
        }
        start = start->next;
    }

    _tinfo->nReq = 0;
    for(int i=0;i<DIM3;i++)
    {
        if(_pinfo->NeProc[i]<0 || _pinfo->rank == _pinfo->NeProc[i]) continue;
        /*
         * some tracers should be thrown out of domain !
         * now those will accumulate near the boundary.
         */
        MPI_Isend(_tinfo->STI+i,1,MPI_INT,_pinfo->NeProc[i],IdTracer,
                  SALEC_CART_COMM,_tinfo->SendRecvReq+_tinfo->nReq);
        MPI_Irecv(_tinfo->RTI+i,1,MPI_INT,_pinfo->NeProc[i],IdTracer,
                  SALEC_CART_COMM,_tinfo->SendRecvReq+_tinfo->nReq+1);
        _tinfo->nReq += 2;
    }

    for(int k=0;k<DIM3;k++)
    {
        if (_pinfo->NeProc[k] < 0 || _pinfo->rank == _pinfo->NeProc[k] || _tinfo->STI[k] <= 0)
        {
            _tinfo->SendTracerId[k] = NULL;
            _tinfo->SendTracer[k] = NULL;
            continue;
        }

        _tinfo->SendTracerId[k] = (int *) malloc(sizeof(int)*_tinfo->STI[k]*TRILEN);
        _tinfo->SendTracer[k]   = (double *) malloc(sizeof(double)*_tinfo->STI[k]*TRDLEN);
    }

    int kSTI[DIM3] = {0};
    start = _tinfo->head;
    Tracer * pstart = NULL;

/*
#ifdef TRACER_DEBUG
    static int cycle=0;
    cycle ++;
    int myrank;
    MPI_Comm_rank(SALEC_CART_COMM,&myrank);
    fprintf(stdout,"%d,%s,%d,%d\n",cycle,__func__,__LINE__,myrank);
    fflush(stdout);
#endif    
*/

    while(NULL!=start)
    {
        if(start->NeE < 0)
        {
            // send to its neighbour processor
            int DestProc = DIM3 +start->NeE;

            if(_pinfo->NeProc[DestProc]>=0)
            {
                int * kSendTracerId = _tinfo->SendTracerId[DestProc];
                double * kSendTracer= _tinfo->SendTracer[DestProc];
                kSendTracerId[kSTI[DestProc]*TRILEN + 0] = (int) start->Id;
                kSendTracerId[kSTI[DestProc]*TRILEN + 1] = (int) start->matId;
                Copy(start->Position,kSendTracer+kSTI[DestProc]*TRDLEN);
                kSendTracer[kSTI[DestProc]*TRDLEN + 3] = start->PeakPre;
                kSendTracer[kSTI[DestProc]*TRDLEN + 4] = start->PeakTem;
                kSendTracer[kSTI[DestProc]*TRDLEN + 5] = start->DumpTmp;
                kSTI[DestProc] ++;
                // Yue
            /*    TracerDel(_tinfo,pstart);
               if(NULL == pstart)
                   start = _tinfo->head;
               else
                   start = pstart->next;
            }*/
            }
            //  Yue
            TracerDel(_tinfo,pstart);
            if(NULL == pstart)
                start = _tinfo->head;
            else
            {
                start = pstart->next;
                //pstart = start;
                //start = start->next;
                // start = pstart->next;  
            }
        } else if(1 && (_mesh->Elements+start->NeE)->Density < TOLRHO)  // Yue
        {
            // remove tracers in the void cell
            TracerDel(_tinfo,pstart);
            if(NULL == pstart)
                start = _tinfo->head;
            else
            {
                //pstart = start;
                //start = start->next;
                start = pstart->next;
            }
        }
        else
        {
            pstart = start;
            start = start->next;
        }

    }
/*
#ifdef TRACER_DEBUG
    MPI_Comm_rank(SALEC_CART_COMM,&myrank);
    fprintf(stdout,"%d,%s,%d,%d\n",cycle,__func__,__LINE__,myrank);
    fflush(stdout);
#endif    
*/

}

void TracerInfoExchangeEnd(TracerInfo * _tinfo, struct ProcessInfo * _pinfo)
{
    if(strcasecmp(_tinfo->Type,"NONE") == 0) return;
    MPI_Waitall(_tinfo->nReq,_tinfo->SendRecvReq,MPI_STATUSES_IGNORE);
    for(int k=0;k<DIM3;k++)
    {
        if(_pinfo->NeProc[k] <0 || _pinfo->NeProc[k]==_pinfo->rank || _tinfo->RTI[k] <=0)
        {
            _tinfo->RecvTracer[k] = NULL;
            _tinfo->RecvTracerId[k] = NULL;
            continue;
        }

        _tinfo->RecvTracer[k] = (double *) malloc(sizeof(double)*_tinfo->RTI[k]*TRDLEN);
        _tinfo->RecvTracerId[k] = (int *) malloc(sizeof(int)*_tinfo->RTI[k]*TRILEN);
    }
}

void TracerExchangeBegin(TracerInfo * _tinfo, struct ProcessInfo * _pinfo, struct Mesh * _mesh)
{
    if(strcasecmp(_tinfo->Type,"NONE") == 0) return;
    int nReq = 0;
    for(int i=0;i<DIM3;i++)
    {
        if (_pinfo->NeProc[i] < 0 || _pinfo->rank == _pinfo->NeProc[i]) continue;
        if(_tinfo->STI[i]>=1)
        {
            MPI_Isend(_tinfo->SendTracer[i],_tinfo->STI[i]*TRDLEN,MPI_DOUBLE,
                      _pinfo->NeProc[i],IdTracer+2,SALEC_CART_COMM,_tinfo->SendRecvReq+nReq++);
            MPI_Isend(_tinfo->SendTracerId[i],_tinfo->STI[i]*TRILEN,MPI_INT,
                      _pinfo->NeProc[i],IdTracer+3,SALEC_CART_COMM,_tinfo->SendRecvReq+nReq++);
        }

        if(_tinfo->RTI[i]>=1)
        {
            MPI_Irecv(_tinfo->RecvTracer[i],_tinfo->RTI[i]*TRDLEN,MPI_DOUBLE,
                      _pinfo->NeProc[i],IdTracer+2,SALEC_CART_COMM,_tinfo->SendRecvReq+nReq++);
            MPI_Irecv(_tinfo->RecvTracerId[i],_tinfo->RTI[i]*TRILEN,MPI_INT,
                      _pinfo->NeProc[i],IdTracer+3,SALEC_CART_COMM,_tinfo->SendRecvReq+nReq++);
        }

    }
    _tinfo->nReq = nReq;
}

void TracerExchangeEnd(TracerInfo * _tinfo, struct ProcessInfo * _pinfo, struct Mesh * _mesh)
{
    if(strcasecmp(_tinfo->Type,"NONE") == 0) return;
    MPI_Waitall(_tinfo->nReq,_tinfo->SendRecvReq,MPI_STATUSES_IGNORE);
    for(int i=0;i<DIM3;i++)
    {
        if (_pinfo->NeProc[i] < 0 || _pinfo->rank == _pinfo->NeProc[i]) continue;
        if(_tinfo->STI[i]>=1)
        {
            free(_tinfo->SendTracerId[i]);
            free(_tinfo->SendTracer[i]);
            _tinfo->STI[i] = 0;
        }

        if(_tinfo->RTI[i]>=1)
        {
            for(int k=0;k<_tinfo->RTI[i];k++)
            {
                Tracer * tmpTracer = (Tracer *) malloc(sizeof(Tracer));
                tmpTracer->Id = _tinfo->RecvTracerId[i][k*TRILEN + 0];
                tmpTracer->matId = _tinfo->RecvTracerId[i][k*TRILEN + 1];
                Copy(_tinfo->RecvTracer[i]+k*TRDLEN,tmpTracer->Position);
                tmpTracer->PeakPre = _tinfo->RecvTracer[i][k*TRDLEN+3];
                tmpTracer->PeakTem = _tinfo->RecvTracer[i][k*TRDLEN+4];
                tmpTracer->DumpTmp = _tinfo->RecvTracer[i][k*TRDLEN+5];
                NeEDetect(tmpTracer,_tinfo,_mesh);
                TailInsert(_tinfo,tmpTracer);
            }
            free(_tinfo->RecvTracerId[i]);
            free(_tinfo->RecvTracer[i]);
            _tinfo->RTI[i] = 0;
        }
    }
}

void SingleTracerMove(Tracer * _self, struct Mesh * _mesh,double dt)
{

    struct Element * ie = _mesh->Elements + _self->NeE;

#ifdef TRACER_DETAIL
    int myrank;
    MPI_Comm_rank(SALEC_CART_COMM,&myrank);
    if(myrank==0)
    {
        fprintf(stdout,"c-t distance of %d:%e,vel:%e:(mass)%e:(dt)%e:(density)%e:(vol)%e\n",
                _self->NeE,
                Distance(ie->Center,_self->Position),
                Length(ie->Momentum)
                ,ie->Mass,dt,
                ie->Density,ie->Volume);
    }
#endif

    if(ie->Density < TOLRHO)
    {
        // in void element, just move according to gravity/body force (.)
        // just delete those tracer
        //
    }
    else
    {
        ScalerAddition(_self->Position,ie->Momentum,dt/ie->Mass);
        if(ie->Pressure > _self->PeakPre) _self->PeakPre = ie->Pressure;   // Yue
        if(ie->Temperature > _self->PeakTem) _self->PeakTem = ie->Temperature;   // Yue
    }

}



void TracerMove(TracerInfo * _tinfo, struct Mesh * _mesh,double dt)
{
    // move tracer position with velocity of cell center
    //
    Tracer * start = _tinfo->head;
    while(NULL!=start)
    {
        if(start->NeE >= 0)
        {
            SingleTracerMove(start,_mesh,dt);
        }
        start = start->next;
    }
}

void TailInsert(TracerInfo * _tinfo,Tracer * _x)
{
    if(NULL == _tinfo->tail)
    {
        _tinfo->tail = _x;
        _tinfo->head = _x;
    }
    else
    {
        _tinfo->tail->next = _x;
        _tinfo->tail = _x;
    }

    _x->next = NULL;
    _tinfo->TracerNum ++;
}

void HeadInsert(TracerInfo * _tinfo,Tracer * _x)
{
    _x->next = _tinfo->head;
    _tinfo->head = _x;
    _tinfo->TracerNum ++;
}

void TracerDel(TracerInfo * _tinfo, Tracer * _x)
{
    if(_x==NULL)
    {
        // _x == NULL is invalid calling of Del
        Tracer * tmpDel = _tinfo->head;
        if(_tinfo->head->next==NULL) _tinfo->tail=NULL;
        _tinfo->head = tmpDel->next;
        //if(_tinfo->TracerNum==1) _tinfo->tail=NULL;
        free(tmpDel);
        _tinfo->TracerNum --;
    } else
    {
        Tracer * tmpDel = _x->next;
        if(tmpDel != NULL)
        {
           if(_x->next->next == NULL) _tinfo->tail = _x;
           _x->next = tmpDel->next;
            // Yue
           free(tmpDel);
           _tinfo->TracerNum --;
        }
    }

}

void TailInsertList(TracerInfo * _tinfo,Tracer * _x,unsigned int _Len)
{
    for(int k=0;k<_Len-1;k++)
    {
        _x[k].next = &(_x[k+1]);
    }
    _x[_Len-1].next = NULL;
    _tinfo->tail->next = &(_x[0]);
}

void InitTracer(TracerInfo * _tinfo, InputFile * ifp,struct Mesh * _mesh)
{
    char TracerOpt[MaxStrLen];
    GetValueS(ifp,"tracer.type",TracerOpt,"NONE");
    GetValueS(ifp,"tracer.prefix",_tinfo->Prefix,"default");
    strcpy(_tinfo->Type,TracerOpt);

    int myrank;
    MPI_Comm_rank(SALEC_CART_COMM,&myrank);
    if(0==strcasecmp(TracerOpt,"NONE"))
    {
        if(0==myrank)
            fprintf(stdout,"tracer is trun off\n");
    } else if(0==strcasecmp(TracerOpt,"BenchmarkSP"))
    {
        double PointA[DIM],PointB[DIM];
        for(int k=0;k<DIM;k++)
        {
            PointA[k] = GetValueDk(ifp,"tracer.A",k,"0.0");
            PointB[k] = GetValueDk(ifp,"tracer.B",k,"0.0");
        }

        int NumTracer = GetValueI(ifp,"tracer.ntracer","2");

        GetValueS(ifp,"tracer.prefix",_tinfo->Prefix,"TestTracer");

        if(NumTracer>=2) {
            for (int k = 0; k < NumTracer; ++k) {
                Tracer *tmpTracer = (Tracer *) malloc(sizeof(Tracer));
                tmpTracer->Id = k;
                double _rNum = 1.0/(NumTracer -1.);
                LinearOp(tmpTracer->Position, PointB, PointA, (NumTracer - k - 1) *_rNum,
                         k*_rNum);
                NeEDetect(tmpTracer, _tinfo, _mesh);
                if (tmpTracer->NeE >= 0) {
                    TailInsert(_tinfo, tmpTracer);
                } else {
                    free(tmpTracer);
                }
            }
        } else if(NumTracer == 1)
        {
            Tracer *tmpTracer = (Tracer *) malloc(sizeof(Tracer));
            tmpTracer->Id = 0;
            Copy(PointA,tmpTracer->Position);
            NeEDetect(tmpTracer, _tinfo, _mesh);
            if (tmpTracer->NeE >= 0) {
                TailInsert(_tinfo, tmpTracer);
            } else {
                free(tmpTracer);
            }
        } else if(NumTracer <= 0)
        {
            fprintf(stdout,"NumTracer = %d!\n",NumTracer);
            exit(0);
        }
    } else if(0==strcasecmp(TracerOpt,"Field"))
    {
        if(0==myrank)
            fprintf(stdout,"Tracer type is %s\n",TracerOpt);
        for(unsigned int k=0;k<_mesh->nde;++k)
        {
            // skip some void element
            struct Element * _ie = _mesh->DomainE[k];
            if(_ie->Density < TOLRHO) continue;
            Tracer *tmpTracer = (Tracer *) malloc(sizeof(Tracer));
            Copy(_ie->Center,tmpTracer->Position);
            tmpTracer->PeakPre = _ie->Pressure;  // Yue
            tmpTracer->PeakTem = _ie->Temperature;  //Yue
            tmpTracer->DumpTmp = 0.0;
            tmpTracer->matId = 0;
            for(int iMat=1;iMat<NMT;++iMat)
            {
                if(_ie->VOF[iMat] > _ie->VOF[tmpTracer->matId]) tmpTracer->matId = iMat;
            }

            tmpTracer->Id = (int)(_ie - _mesh->Elements) + _mesh->ne*myrank;

            NeEDetect(tmpTracer,_tinfo,_mesh);
            if (tmpTracer->NeE >= 0)
            {
                // inside of this progress cells
                TailInsert(_tinfo, tmpTracer);
            } else
            {
                // outside, delete this tracer
                free(tmpTracer);
            }
        }
    } else
    {
        // unknown tracer type
        if(0==myrank)
            fprintf(stdout,">[ ] unknown Tracer Type [%s]\n",TracerOpt);
    }


    /*
     * init the ejecta post processing
     */
    char EjectaOpt[MaxStrLen];

    double EjectaThreshold = GetValueD(ifp,"ejecta.threshold","2.0");
    double ProjRadius = GetValueD(ifp,"projectile.radiu","7000.0");
    _tinfo->Hthreshold = EjectaThreshold*ProjRadius;
    GetValueS(ifp,"ejecta.type",EjectaOpt,"NONE");
    strcpy(_tinfo->EjectaModel,EjectaOpt);
    if( (0 == strcasecmp(EjectaOpt,"Flat")) || (0 == strcasecmp(EjectaOpt,"Plane")))
    {
        _tinfo->EjectaFunc = FlatEjectaDetect;

        if(0 == myrank)
        {
            fprintf(stdout,"%s:ejecta model is %s\n",__func__, EjectaOpt);
            fprintf(stdout,"%s:threshold is %e\n",__func__,_tinfo->Hthreshold);
        }
    }
    else if(0 == strcasecmp(EjectaOpt,"Sphere"))
    {
        _tinfo->Pcenter[X] = 0;
        _tinfo->Pcenter[Y] = 0;
        _tinfo->Pcenter[Z] = -1.0*GetValueD(ifp,"target.r_planet","1738.0e3");
        _tinfo->Hthreshold += fabs(_tinfo->Pcenter[Z]);
        _tinfo->EjectaFunc = PlanetEjectaDetect;

        if(0 == myrank)
        {
            fprintf(stdout,"%s:ejecta model is %s\n",__func__, EjectaOpt);
            fprintf(stdout,"%s:threshold is %e\n",__func__,_tinfo->Hthreshold);
        }
    }
    else if(0 == strcasecmp(EjectaOpt,"NONE"))
    {
        _tinfo->EjectaFunc = NULL;
        Zero(_tinfo->Pcenter);
        if(0 == myrank)
        {
            fprintf(stdout,"%s:ejecta model is OFF\n",__func__);
        }
    }
    else
    {
        _tinfo->EjectaFunc = NULL;
        fprintf(stdout,"%s:Unimplemeted ejecta option:%s\n",__func__,EjectaOpt);
        exit(0);
    }
    _tinfo->efp = NULL;

    /*
     * output some info about tracer
     */
#ifdef TRACER_DETAIL
    fprintf(stdout,">[%d] Tracer[%s] is activated; \n"
                   ">[%d] %d tracer initialized;\n",myrank,TracerOpt,myrank,_tinfo->TracerNum);
#endif
}

int PlanetEjectaDetect(Tracer * _t, TracerInfo * _tinfo, struct Element * _mesh)
{
    if((_t->matId >= 0) && (Distance(_t->Position,_tinfo->Pcenter) > _tinfo->Hthreshold))
    {
        const struct Element * _ie = _mesh + _t->NeE;
        double nSuf[DIM] ={0.0};
        LinearOp(nSuf,_t->Position,_tinfo->Pcenter,1.0,-1.0);
        double rVel = Dot(nSuf,_ie->Momentum);
        if(rVel > 0.0)
            return 1;
        else
            return 0;
    }
    else
    {
        return 0;
    }
}


int FlatEjectaDetect(Tracer * _t, TracerInfo * _tinfo, struct Element * _mesh)
{
    if((_t->matId >= 0)&&(_t->Position[Z] > _tinfo->Hthreshold))
    {
        const struct Element * _ie = _mesh + _t->NeE;
        double nSuf[DIM] ={0.0, 0.0, 1.0};
        double rVel = Dot(nSuf,_ie->Momentum);
        if(rVel > 0.0)
            return 1;
        else
            return 0;
    }
    else
    {
        return 0;
    }
}


void OutputTracer(TracerInfo * _tinfo,struct Mesh * _mesh,int step,int rank)
{
    /*
     * write all tracer info into txt
     */
    #if WRITETRACERTXT==0
    return;
    #endif

    if(strcasecmp(_tinfo->Type,"NONE") == 0) return;
    if(_tinfo->EjectaFunc == NULL) return;
    char _tname[MaxStrLen*2];
    sprintf(_tname,"txt/%s.proc%d.%04d.tracer",_tinfo->Prefix,rank,step);
    FILE * fp = fopen(_tname,"w");
    fprintf(fp,"    %d\n",_tinfo->TracerNum);
    Tracer * start = _tinfo->head;
    while(NULL!=start)
    {
        if(start->NeE < 0)
        {
            fprintf(fp,"%15d, %15d, %15.5f, %15.5f, %15.5f, %15.5e, %15.5f\n",
                    start->Id,start->NeE,
                    start->Position[X],start->Position[Y],start->Position[Z],
                    0.0,
                    0.0);
        } else
        {
            fprintf(fp,"%15d, %15d, %15.5f, %15.5f, %15.5f, %15.5e, %15.5f\n",
                    start->Id,start->NeE,
                    start->Position[X],start->Position[Y],start->Position[Z],
                    _mesh->Elements[start->NeE].Pressure,
                    _mesh->Elements[start->NeE].Density);
        }
        start = start->next;
    }
    fclose(fp);
}

void OutputTracerVtp(TracerInfo * _tinfo,struct Mesh * _mesh,int step,int rank)
{
    if(strcasecmp(_tinfo->Type,"NONE") == 0) return;
    // collect info from tracers
    const unsigned int tr_len = _tinfo->TracerNum + 1;
    float * tr_pos = malloc(sizeof(float)*tr_len*DIM);
    float * tr_vel = malloc(sizeof(float)*tr_len*DIM);
    float * tr_mid  = malloc(sizeof(float)*tr_len);
    unsigned int * tr_gid  = malloc(sizeof(unsigned int)*tr_len);
    float * tr_mpre = malloc(sizeof(float)*tr_len);
    float * tr_mtem = malloc(sizeof(float)*tr_len);
    float * tr_epre = malloc(sizeof(float)*tr_len);
    float * tr_etem = malloc(sizeof(float)*tr_len);
    float * tr_dump = malloc(sizeof(float)*tr_len);
    float * tr_eden = malloc(sizeof(float)*tr_len);

    // write some dump info
    tr_pos[X] = tr_pos[Y] = tr_pos[Z] = 0.0;
    tr_vel[X] = tr_vel[Y] = tr_vel[Z] = 0.0;
    tr_mid[0] = 0.0;
    tr_gid[0] = 0.0;
    tr_mpre[0] = 0.0;
    tr_mtem[0] = 0.0;
    tr_epre[0] = 0.0;
    tr_etem[0] = 0.0;
    tr_dump[0] = 0.0;
    tr_eden[0] = 0.0;

    Tracer * start = _tinfo->head;
    unsigned int tr_index = 0;
    while(NULL!=start)
    {
        if(start->NeE < 0)
        {
            // skip, this tracer will be ignored
        } else
        {
            if(tr_index > _tinfo->TracerNum)
            {
                fprintf(stdout,"%s:tr_index(%d) > TracerNum(%d)\n",__func__,tr_index,_tinfo->TracerNum);
                exit(0);
            }
            // info in tracer_t
            tr_pos[tr_index*DIM + X] = start->Position[X];
            tr_pos[tr_index*DIM + Y] = start->Position[Y];
            tr_pos[tr_index*DIM + Z] = start->Position[Z];

            tr_mid[tr_index] = 1.0*start->matId;
            tr_gid[tr_index] = start->Id;
            tr_mpre[tr_index] = start->PeakPre;
            tr_mtem[tr_index] = start->PeakTem;
            tr_dump[tr_index] = start->DumpTmp;

            // info in element
            struct Element * _tie = _mesh->Elements + start->NeE;
            tr_vel[tr_index*DIM + X] = _tie->Momentum[X]/_tie->Mass;
            tr_vel[tr_index*DIM + Y] = _tie->Momentum[Y]/_tie->Mass;
            tr_vel[tr_index*DIM + Z] = _tie->Momentum[Z]/_tie->Mass;
            tr_epre[tr_index] = _tie->Pressure;
            tr_etem[tr_index] = _tie->Temperature;
            tr_eden[tr_index] = _tie->Density;
            tr_index++;
        }
        start = start->next;
    }
    if(tr_index >= 1)
    {
        char _tname[MaxStrLen*2];
        sprintf(_tname,"vtp/%s.proc%d.%d.vtp",_tinfo->Prefix,rank,step);
        const char * vtp_data_format = "binary";
        FILE * fp = fopen(_tname,"w");
        vtp_file_header(fp,tr_index);
        vtk_point_data_header(NULL,fp);
        vtk_dataarray_vec_f(fp,"matid",vtp_data_format,tr_mid,tr_index,1);
        vtk_dataarray_vec_u32(fp,"id",vtp_data_format,tr_gid,tr_index,1);
        vtk_dataarray_vec_f(fp,"pre",vtp_data_format,tr_epre,tr_index,1);
        vtk_dataarray_vec_f(fp,"mpre",vtp_data_format,tr_mpre,tr_index,1);
        vtk_dataarray_vec_f(fp,"tem",vtp_data_format,tr_etem,tr_index,1);
        vtk_dataarray_vec_f(fp,"mtem",vtp_data_format,tr_mtem,tr_index,1);
        vtk_dataarray_vec_f(fp,"vel",vtp_data_format,tr_vel,tr_index,3);
        vtk_dataarray_vec_f(fp,"dump",vtp_data_format,tr_dump,tr_index,1);
        vtk_dataarray_vec_f(fp,"den",vtp_data_format,tr_eden,tr_index,1);
        vtk_point_data_trailer(NULL,fp);
        vtk_point_header(fp);
        vtk_dataarray_vec_f(fp,"coordinate",vtp_data_format,tr_pos,tr_index,3);
        vtk_point_trailer(fp);
        vtp_file_trailer(fp);
        fclose(fp);
    }
    else
    {
        char _tname[MaxStrLen*2];
        sprintf(_tname,"vtp/%s.proc%d.%d.vtp",_tinfo->Prefix,rank,step);
        const char * vtp_data_format = "ascii";
        FILE * fp = fopen(_tname,"w");
        vtp_file_header(fp,0);
        vtk_point_data_header(NULL,fp);
        vtk_dataarray_vec_f(fp,"matid",vtp_data_format,tr_mid,tr_index,1);
        vtk_dataarray_vec_u32(fp,"id",vtp_data_format,tr_gid,tr_index,1);
        vtk_dataarray_vec_f(fp,"pre",vtp_data_format,tr_epre,tr_index,1);
        vtk_dataarray_vec_f(fp,"mpre",vtp_data_format,tr_mpre,tr_index,1);
        vtk_dataarray_vec_f(fp,"tem",vtp_data_format,tr_etem,tr_index,1);
        vtk_dataarray_vec_f(fp,"mtem",vtp_data_format,tr_mtem,tr_index,1);
        vtk_dataarray_vec_f(fp,"vel",vtp_data_format,tr_vel,tr_index,3);
        vtk_dataarray_vec_f(fp,"dump",vtp_data_format,tr_dump,tr_index,1);
        vtk_dataarray_vec_f(fp,"den",vtp_data_format,tr_eden,tr_index,1);
        vtk_point_data_trailer(NULL,fp);
        vtk_point_header(fp);
        vtk_dataarray_vec_f(fp,"coordinate",vtp_data_format,tr_pos,tr_index,3);
        vtk_point_trailer(fp);
        vtp_file_trailer(fp);
        fclose(fp);
    }

    free(tr_pos);
    free(tr_vel);
    free(tr_mid);
    free(tr_gid);
    free(tr_mpre);
    free(tr_mtem);
    free(tr_epre);
    free(tr_etem);
    free(tr_dump);
    free(tr_eden);
}

void OutputEjecta(TracerInfo * _tinfo,struct Mesh * _mesh,double time,int rank)
{
    // if(rank == 0)
    // {
    //     fprintf(stdout,"%s: %f,%d\n",__func__ ,time,rank);
    // }
    if(strcasecmp(_tinfo->Type,"NONE") == 0) return;
    if(_tinfo->EjectaFunc == NULL) return;
    if(NULL == _tinfo->efp)
    {
        fprintf(stdout,"%s:%d try write ejecta to null file\n",__func__,rank);
        exit(0);
    }
    Tracer * start = _tinfo->head;
    Tracer * pstart = NULL;
    unsigned int ej_no=0;
    while(NULL!=start)
    {
        if((_tinfo->EjectaFunc)(start,_tinfo,_mesh->Elements))
        {
            ej_no++;
        }
        start = start->next;
    }
    // only output the file when the ejecta number is non-zero
    if(ej_no > 0)
    {
        Tracer * start = _tinfo->head;
        Tracer * pstart = NULL;

        while(NULL!=start)
        {
            if((_tinfo->EjectaFunc)(start,_tinfo,_mesh->Elements))
            {
                start->matId *= -1;
                struct Element * ie = _mesh->Elements + start->NeE;
                double vp[DIM] = {0.0};
                for(int k=0;k<DIM;k++)  vp[k] = ie->Momentum[k]/ie->Mass;

                if(start->NeE < 0)
                {
                    fprintf(_tinfo->efp,"%d, %d, %d, %d, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f\n",
                            start->Id,rank,0,start->NeE,time,
                            start->Position[X],start->Position[Y],start->Position[Z],
                            start->PeakPre,start->PeakTem,
                            0.0,0.0,0.0,
                            -1.0,-1.0);
                } else
                {
                    fprintf(_tinfo->efp,"%d, %d, %d, %d,  %15.5f, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f, %15.5f\n",
                            start->Id,rank,start->matId,start->NeE,time,
                            start->Position[X],start->Position[Y],start->Position[Z],
                            start->PeakPre,start->PeakTem,
                            vp[X],vp[Y],vp[Z],
                            _mesh->Elements[start->NeE].Pressure,
                            _mesh->Elements[start->NeE].Temperature);
                }
                // keep ejecta tracer in mem
                /*TracerDel(_tinfo, pstart);
                if(NULL == pstart)
                {
                    start = _tinfo->head;
                }else
                {
                    start = pstart->next;
                }*/
                pstart = start;
                start = start->next;

            }else
            {
                pstart = start;
                start = start->next;
            }
        }
    }
}

void FlushEjecta(TracerInfo * _tinfo, int step, int rank)
{
    if(strcasecmp(_tinfo->Type,"NONE") == 0) return;
    if(_tinfo->EjectaFunc == NULL) return;
    char _tname[MaxStrLen*2];
    sprintf(_tname,"ejecta/%s.proc%d.%04d.ejecta",_tinfo->Prefix,rank,step);
    if(NULL != _tinfo->efp) fclose(_tinfo->efp);
    _tinfo->efp = fopen(_tname,"w");
    // fprintf(stdout,"make %s in %d/%d step\n",_tname,step,rank);
    if(NULL == _tinfo->efp)
    {
        fprintf(stdout,"%s(%d):cannot open %s\n",__func__,rank,_tname);
        exit(0);
    }
}

void TracerClean(TracerInfo * _tinfo)
{
    Tracer * start = _tinfo->head;
    Tracer * pstart;
    while (NULL!=start)
    {
        pstart = start;
        start = start->next;
        free(pstart);
    }

    for(int k=0;k<DIM;++k)
    {
        free(_tinfo->Coordinate[k]);
    }
    free(_tinfo->Coordinate);
}

void InitTracerInfo(TracerInfo * _tinfo, struct MeshInfo * _minfo)
{
    // set coordinate line
    _tinfo->Coordinate = (double **) malloc(DIM* sizeof(double *));
    _tinfo->Coordinate[X] = (double *) malloc(sizeof(double)*_minfo->npx);
    for(int i=0;i<_minfo->npx;++i)
    {
        _tinfo->Coordinate[X][i] = _minfo->V[1][1][i][X];
    }
    _tinfo->Coordinate[Y] = (double *) malloc(sizeof(double)*_minfo->npy);
    for(int j=0;j<_minfo->npy;++j)
    {
        _tinfo->Coordinate[Y][j] = _minfo->V[1][j][1][Y];
    }
    _tinfo->Coordinate[Z] = (double *) malloc(sizeof(double)*_minfo->npz);
    for(int k=0;k<_minfo->npz;++k)
    {
        _tinfo->Coordinate[Z][k] = _minfo->V[k][1][1][Z];
    }

    _tinfo->TracerNum = 0;
    _tinfo->tail = NULL;
    _tinfo->head = NULL;

    for(int k=0;k<DIM3;++k)
    {
        _tinfo->RTI[k] = 0;
        _tinfo->STI[k] = 0;
    }
}

