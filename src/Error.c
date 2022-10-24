//
// Created by huacheng on 2020/12/22.
//

//#include <arpa/nameser.h>
#include "Error.h"
void NullPointer(FILE *fp, const char *_message) {
    fprintf(fp, "%s", _message);
    exit(0);
}

void CheckInt(FILE *fp, const char *_message, int x)
{
    fprintf(fp, "%s", _message);
    fprintf(fp, "Value=%d\n", x);
}

void CheckDouble(FILE *fp, const char *_message, double x)
{
    fprintf(fp, "%s", _message);
    fprintf(fp, "Value=%f\n", x);
}

void RunPoint(FILE *fp, const char *_message)
{
#if __RUNPOINT__
    static int Truns = 0;
    fprintf(fp,"# check trun #%d\n",Truns++);
    fprintf(fp,_message);
#endif
}

void RunPoint_s(FILE *fp, const char *_message)
{
    static int Truns = 0;
    fprintf(fp,"# check turn #%d\n",Truns++);
    fprintf(fp,_message);
}

void TestMeshInfo_air(struct MeshInfo *_info){

    _info->npx = 400;
    _info->npy = 10;
    _info->npz = 10;
    _info->ARTVIS   = 2.0;
    _info->CISSAMK1 = 1.0;
    _info->MaxVelocity = 10.0*20000.0;

    double dx = 400.0/(_info->npx-1);
    double dy = 20.0/(_info->npy-1);
    double dz = 20.0/(_info->npz-1);
    double ddxyz = 0;

    double pX0 = -200.0;
    double pY0 = -10.0;
    double pZ0 = -10.0;

    _info->V = (double ****) malloc(_info->npz * sizeof(double ***));
    for (int k = 0; k < _info->npz; k++)
    {
        _info->V[k] = (double ***) malloc(_info->npy * sizeof(double **));
        for (int j = 0; j < _info->npy; j++)
        {
            _info->V[k][j] = (double **) malloc(_info->npx * sizeof(double *));
            for (int i = 0; i < _info->npx; i++)
            {
                _info->V[k][j][i] = (double *) malloc(DIM * sizeof(double));
                _info->V[k][j][i][X] = pX0 + i * dx + ddxyz*RandomDouble() * dx;
                _info->V[k][j][i][Y] = pY0 + j * dy + ddxyz*RandomDouble() * dy;
                _info->V[k][j][i][Z] = pZ0 + k * dz + ddxyz*RandomDouble() * dz;
            }
        }
    }

    _info->MaterialNum = 2;
    char MaterialName[3][30] = {
            "../eos/basalt_.aneos",
            "../eos/AirEos2.txt"
    };

    for(int i=1;i<_info->MaterialNum;i++)
    {
        _info->Materialfp[i] = fopen(MaterialName[i],"r");
        if(NULL == _info->Materialfp[i])
        {
            fprintf(stdout,"Cannot open %s\n",MaterialName[i]);
            exit(0);
        }
    }

}

double StopWatch(int _control, int _k) {
    static struct timeval _start;
    static double _sequence[500];
    static int k = 0;

    if (0 == _control)
    {
        if (0 == k)
        {
            k++;
            gettimeofday(&_start, NULL);
            return 0.0;
        }
        else{
            struct timeval _end;
            gettimeofday(&_end, NULL);
            _sequence[k++] = (_end.tv_sec - _start.tv_sec) + (_end.tv_usec - _start.tv_usec) / 1000000.0;
            return _sequence[k-1];
        }
    }
    else if (1 == _control)
    {
        return _sequence[_k];
    }
    else if(2 == _control)
    {
        if(_k > 0)
            return _sequence[_k] - _sequence[_k-1];
        else
            return _sequence[_k];
    }
    else if(3 == _control)
    {
        return _sequence[k-1];
    }
}

void CheckDoubleArray(FILE *fp, const char *_message, const double *x, int length) {
    fprintf(fp, "%s\n", _message);
    for (int k = 0; k < length; k++)
    {
        fprintf(fp, "%10.5f, ", x[k]);
        if(k%15 == 14) fprintf(fp,"\n");
    }

    fprintf(fp, "#\n");
}

void CheckDoubleArrayE(FILE *fp, const char *_message, const double *x, int length) {
    fprintf(fp, "%s\n", _message);
    for (int k = 0; k < length; k++)
        fprintf(fp, "%10.5e, ", x[k]);
    fprintf(fp, "#\n");
}

void CheckElement(const struct Element *_element) {
    CheckDoubleArray(stdout, "The SubVolumRatio:", _element->SubVolumeRatio, NSEPE);
}

double RandomDouble()
{
    static int k = 0;
    if (0 == k) {
        k++;
        srand((unsigned) time(NULL));
        return (rand() % 10000) / 10000.0;
    } else
        return (rand() % 10000) / 10000.0;
}

void Shuffle(int a[], int len)
{
    /*
     * Fisher and Yate Algorithm
     */
    static int State = 0;
    if(0==State)
    {
        State++;srand((unsigned) time(NULL));
    }
    for(int k=0;k<len;k++) a[k] = k;

    for(int k=len-1;k>=0;k--)
    {
        int j = rand()%(k+1);
        if(j!=k)
        {
            int tmp = a[j];
            a[j] = a[k];
            a[k] = tmp;
        }
    }
}


void CheckANEOSTable(struct ANEOSTable * _t)
{
    /*
     *  Here We only check the size and the first 10
     *  data in Table
     */

    fprintf(stdout,"#head of CheckANEOSTable\n");
    fprintf(stdout,"The size of this aneos table is (%d,%d),  ",_t->nTem,_t->nDen);
    fprintf(stdout,"and the first 10 data in this Tem is:\n");
    for(int k=0;k<5;k++)
    {
        fprintf(stdout,"%1.10e,%1.10e\n",_t->xTem[2*k],_t->xTem[2*k+1]);
    }
    fprintf(stdout,"the first 10 data in this Table is:\n");
    for(int k=0;k<5;k++)
    {
        fprintf(stdout,"%1.12e,%1.12e,%1.12e\n",_t->Data[0][k][0],_t->Data[0][k][1],_t->Data[0][k][2]);
    }
    fprintf(stdout,"The last 10 data in Den is:\n");
    for(int k=0;k<5;k++)
    {
        fprintf(stdout,"%1.10e,%1.10e\n",_t->yDen[_t->nDen -10 + 2*k],_t->yDen[_t->nDen - 10 + 2*k+1]);
    }

    fprintf(stdout,"#end of CheckANEOSTable\n");

}

void CheckEformat()
{
    char TestNum[] = "0.451759751638-11";
    double TestN = 0.0;
    sscanf(TestNum,"%le",&TestN);
    fprintf(stdout,"%le\n",TestN);
//    sscanf(TestNum,"%le",&TestN);
    TestN = strtod(TestNum,NULL);
    fprintf(stdout,"%le\n",TestN);
}

void CheckStateReference(struct StateReference * _s, struct ANEOSTable * _t)
{
    fprintf(stdout,"#head of CheckStateReference\n");
    fprintf(stdout,"The reference in this SRF is:\n");
    fprintf(stdout,"norm(P,T,Rho) = (%lf,%lf,%lf)\n",_t->Pnorm,_t->Tnorm,_t->Dnorm);
    fprintf(stdout,"#end of CheckStateReference\n");
}

void TestVelStavlity(struct Vertex * iv)
{
    /*
     * Set dx_min which related the stablity
     */
    int First = 1.0;
    for(int i=0;i<VPE;i++)
    {
        if(NULL == iv->NeE[i])
            continue;
        else if(First)
        {
            First = 0.0;
            iv->dx_min = 1.0/iv->NeE[i]->rdx_min;
            if(iv->dx_min < 0)
                exit(0);
        }
        else
            iv->dx_min = Min(1.0/iv->NeE[i]->rdx_min,iv->dx_min);
    }
    iv->dx_min *= COURANT0;
}

void TestMesh(struct Mesh * E, const char * VtsName)
{
    FILE * vtsfp = fopen(VtsName,"w");
    vtk_output(E,vtsfp);
    fclose(vtsfp);
}

void DataClean(struct ANEOSTable * ATB, struct Mesh * _mesh,int nm)
{
    free(_mesh->Vertexes);
    free(_mesh->Elements);
    free(_mesh->XBoundaries);
    free(_mesh->YBoundaries);
    free(_mesh->ZBoundaries);

    free(_mesh->DomainB);
    free(_mesh->DomainE);
    free(_mesh->DomainV);

    free(_mesh->CondE);

    for(int k=1;k<nm;k++)
        UnAllocateANEOS(ATB + k);
}