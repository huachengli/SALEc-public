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

        _tinfo->SendTracerId[k] = (int *) malloc(sizeof(int)*_tinfo->STI[k]);
        _tinfo->SendTracer[k]   = (double *) malloc(sizeof(double)*_tinfo->STI[k]*DIM);
    }

    int kSTI[DIM3] = {0};
    start = _tinfo->head;
    Tracer * pstart = NULL;


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
                kSendTracerId[kSTI[DestProc]] = (int) start->Id;
                Copy(start->Position,kSendTracer+kSTI[DestProc]*DIM);
                kSTI[DestProc] ++;
            }

            TracerDel(_tinfo,pstart);
            if(NULL == pstart)
                start = _tinfo->head;
            else
                start = pstart->next;
        } else if((_mesh->Elements+start->NeE)->Density < TOLRHO)
        {
            // remove tracers in the void cell
            TracerDel(_tinfo,pstart);
            if(NULL == pstart)
                start = _tinfo->head;
            else
                start = pstart->next;
        }
        else
        {
            pstart = start;
            start = start->next;
        }

    }

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

        _tinfo->RecvTracer[k] = (double *) malloc(sizeof(double)*_tinfo->RTI[k]*DIM);
        _tinfo->RecvTracerId[k] = (int *) malloc(sizeof(int)*_tinfo->RTI[k]);
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
            MPI_Isend(_tinfo->SendTracer[i],_tinfo->STI[i]*DIM,MPI_DOUBLE,
                      _pinfo->NeProc[i],IdTracer+2,SALEC_CART_COMM,_tinfo->SendRecvReq+nReq++);
            MPI_Isend(_tinfo->SendTracerId[i],_tinfo->STI[i],MPI_INT,
                      _pinfo->NeProc[i],IdTracer+3,SALEC_CART_COMM,_tinfo->SendRecvReq+nReq++);
        }

        if(_tinfo->RTI[i]>=1)
        {
            MPI_Irecv(_tinfo->RecvTracer[i],_tinfo->RTI[i]*DIM,MPI_DOUBLE,
                      _pinfo->NeProc[i],IdTracer+2,SALEC_CART_COMM,_tinfo->SendRecvReq+nReq++);
            MPI_Irecv(_tinfo->RecvTracerId[i],_tinfo->RTI[i],MPI_INT,
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
                Copy(_tinfo->RecvTracer[i]+k*DIM,tmpTracer->Position);
                tmpTracer->Id = _tinfo->RecvTracerId[i][k];
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

    if(ie->Density < TOLRHO) return;
    ScalerAddition(_self->Position,ie->Momentum,dt/ie->Mass);
}

void TracerMove(TracerInfo * _tinfo, struct Mesh * _mesh,double dt)
{
    // move tracer position with velocity of cell center
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
        _tinfo->head = tmpDel->next;
        free(tmpDel);
    } else
    {
        Tracer * tmpDel = _x->next;
        _x->next = tmpDel->next;
        free(tmpDel);
    }

    _tinfo->TracerNum --;

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
     * output some info about tracer
     */
#ifdef TRACER_DETAIL
    fprintf(stdout,">[%d] Tracer[%s] is activated; \n"
                   ">[%d] %d tracer initialized;\n",myrank,TracerOpt,myrank,_tinfo->TracerNum);
#endif
}

void OutputTracer(TracerInfo * _tinfo,struct Mesh * _mesh,int step,int rank)
{
    if(strcasecmp(_tinfo->Type,"NONE") == 0) return;
    char _tname[MaxStrLen*2];
    sprintf(_tname,"%s.proc%d.%04d.tracer",_tinfo->Prefix,rank,step);
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

