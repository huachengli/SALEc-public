//
// Created by huacheng on 6/23/21.
//

#include "InputParser.h"
char * rtrim(char *str)
{
    if (str == NULL || *str == '\0')
    {
        return str;
    }

    int len = strlen(str);
    char *p = str + len - 1;
    while (p >= str  && isspace(*p))
    {
        *p = '\0';
        --p;
    }
    return str;
}

char * ltrim(char *str)
{
    if (str == NULL || *str == '\0')
    {
        return str;
    }

    int len = 0;
    char *p = str;
    while (*p != '\0' && (isspace(*p) || (*p=='.')))
    {
        ++p;
        ++len;
    }
    memmove(str, p, strlen(str) - len + 1);
    return str;
}

char *trim(char *str)
{
    str = rtrim(str);
    str = ltrim(str);
    return str;
}

int IsComment(const char c[])
{
    if(strlen(c)==0)
        return 1;
    char CommentHead[5]   = "\"#/-*";
    int IsComm = 0;
    for(int k=0;k<5;k++)
    {
        if(c[0] == CommentHead[k])
        {
            IsComm = 1;
            break;
        }
    }
    return IsComm;
}

int InField(const char c[])
{
    const char * p = c;
    for(;*p!='\0';p++)
    {
        if(*p=='.') return 0;
    }
    return 1;
}

int InSubset(char _c,const char _set[])
{
    const char * p = _set;
    for(;*p!='\0';p++)
    {
        if(*p==_c)
            return 1;
    }
    return 0;
}

int Strok(const char _str[],const char _delim[], char value[])
{
    const char * p = _str;
    int k=0,i=0;
    for(;*p!='\0';p++)
    {
        if(InSubset(*p,_delim))
        {
            if(k==0)
                continue;
            else
                break;
        } else
        {
            value[k++] = *p;
        }
        i++;
    }
    for(;*p!='\0';p++)
    {
        if(!InSubset(*p,_delim)) break;
        i++;
    }
    value[k] = '\0';
    return i;
}

InputFile * OpenInputFile(const char fname[])
{
    FILE * fp = fopen(fname,"r");
    InputFile * ifp = (InputFile *) malloc(sizeof(InputFile));
    if(NULL == ifp)
    {
        fprintf(stdout,"Cannot open input file * %s * \n",fname);
        exit(0);
    }
    ifp->Len = 0;

    char LineBuffer[300];
    char MainDelimiter[] = "= ";
    char SubDelimiter[]  = "/#*\"";
    char Field[MaxStrLen] = "mesh";

    while(fscanf(fp,"%[^\n]",LineBuffer)!=EOF)
    {
        fgetc(fp);
        trim(LineBuffer);
        int LineBufferLen = strlen(LineBuffer);
        if((LineBufferLen==0) || IsComment(LineBuffer)) continue;
        if((LineBuffer[0]=='[') && (LineBuffer[LineBufferLen-1]==']'))
        {
            LineBuffer[LineBufferLen-1] = '\0';
            strcpy(Field,LineBuffer+1);
            continue;
        }

        char tKey[MaxStrLen], tValue[MaxStrLen];
        int r = Strok(LineBuffer,MainDelimiter,tKey);
        Strok(LineBuffer+r,SubDelimiter,tValue);
        if((0== strlen(tKey)) || (0== strlen(tValue)))
        {
            LineBuffer[0] = '\0';
            continue;
        }

        if(InField(tKey))
        {
            sprintf(ifp->Key[ifp->Len],"%s.%s",Field,tKey);
            strcpy(ifp->Value[ifp->Len], tValue);
        } else
        {
            strcpy(ifp->Key[ifp->Len], tKey);
            strcpy(ifp->Value[ifp->Len], tValue);
        }

        trim(ifp->Key[ifp->Len]);
        trim(ifp->Value[ifp->Len]);
        ifp->Len++;
        LineBuffer[0] = '\0';
    }
    fclose(fp);
    SortInputFile(ifp,0,ifp->Len-1);
    return ifp;
}

InputFile * ParseInputFile(FILE * fp)
{
    InputFile * ifp = (InputFile *) malloc(sizeof(InputFile));
    if(NULL == ifp)
    {
        fprintf(stdout,"Cannot open input file with fp \n");
        exit(0);
    }
    ifp->Len = 0;

    char LineBuffer[300];
    char MainDelimiter[] = "= ";
    char SubDelimiter[]  = "/#*\"";
    char Field[MaxStrLen] = "mesh";

    while(fscanf(fp,"%[^\n]",LineBuffer)!=EOF)
    {
        fgetc(fp);
        trim(LineBuffer);
        int LineBufferLen = strlen(LineBuffer);
        if((LineBufferLen==0) || IsComment(LineBuffer)) continue;
        if((LineBuffer[0]=='[') && (LineBuffer[LineBufferLen-1]==']'))
        {
            LineBuffer[LineBufferLen-1] = '\0';
            strcpy(Field,LineBuffer+1);
            continue;
        }

        char tKey[MaxStrLen], tValue[MaxStrLen];
        int r = Strok(LineBuffer,MainDelimiter,tKey);
        Strok(LineBuffer+r,SubDelimiter,tValue);
        if((0== strlen(tKey)) || (0== strlen(tValue)))
        {
            LineBuffer[0] = '\0';
            continue;
        }

        if(InField(tKey))
        {
            sprintf(ifp->Key[ifp->Len],"%s.%s",Field,tKey);
            strcpy(ifp->Value[ifp->Len], tValue);
        } else
        {
            strcpy(ifp->Key[ifp->Len], tKey);
            strcpy(ifp->Value[ifp->Len], tValue);
        }

        trim(ifp->Key[ifp->Len]);
        trim(ifp->Value[ifp->Len]);
        ifp->Len++;
        LineBuffer[0] = '\0';
    }
    fclose(fp);
    SortInputFile(ifp,0,ifp->Len-1);
    return ifp;
}


void * CloseInputFile(InputFile * ifp)
{
    free(ifp);
}

void CheckInput(FILE * fp, const InputFile * ifp)
{
    for(int k=0;k<ifp->Len;k++)
    {
        fprintf(fp,"<%3d>:<%3d>:<%s>: = \"%s\"\n",k,strcmp(ifp->Key[k],ifp->Key[k-1>=0?k-1:0]),ifp->Key[k],ifp->Value[k]);
    }
}


void Merge(InputFile * ifp, int start, int mid, int end)
{
    char (*tKey)[MaxStrLen]   = (char (*)[MaxStrLen]) malloc((end+1-start)*sizeof(char [MaxStrLen]));
    char (*tValue)[MaxStrLen] = (char (*)[MaxStrLen]) malloc((end+1-start)*sizeof(char [MaxStrLen]));

    int i = start;
    int j = mid + 1;
    int k = 0;

    while (i <= mid && j <= end)
    {
        if(strcasecmp(ifp->Key[i],ifp->Key[j]) <= 0)
        {
            strcpy(tKey[k],ifp->Key[i]);
            strcpy(tValue[k],ifp->Value[i]);
            k++;i++;
        } else
        {
            strcpy(tKey[k],ifp->Key[j]);
            strcpy(tValue[k],ifp->Value[j]);
            k++;j++;
        }
    }

    while (i <= mid)
    {
        strcpy(tKey[k],ifp->Key[i]);
        strcpy(tValue[k],ifp->Value[i]);
        k++;i++;
    }
    while (j <= end)
    {
        strcpy(tKey[k],ifp->Key[j]);
        strcpy(tValue[k],ifp->Value[j]);
        k++;j++;
    }

    for (i = 0; i < k; i++)
    {
        strcpy(ifp->Key[i+start],tKey[i]);
        strcpy(ifp->Value[i+start],tValue[i]);
    }

    free(tKey);
    free(tValue);
}

void SortInputFile(InputFile * ifp, int start, int end)
{
    if(start >= end)
        return;
    else
    {
        int mid = start + (end - start) / 2;
        SortInputFile(ifp,start,mid);
        SortInputFile(ifp,mid+1,end);
        Merge(ifp,start,mid,end);
    }
}

void Swap(char a[],char b[])
{
    char t[MaxStrLen];
    strcpy(t,a);
    strcpy(a,b);
    strcpy(b,t);
}

void SimpleSort(InputFile * ifp, int start, int end)
{
    for(int i=start;i<=end;i++)
    {
        for(int j=i+1;j<=end;j++)
        {
            if(strcasecmp(ifp->Key[i],ifp->Key[j]) >= 0)
            {
                Swap(ifp->Key[i],ifp->Key[j]);
                Swap(ifp->Value[i],ifp->Value[j]);
            }
        }
    }
}

int SearchInput(InputFile * ifp,const char key[])
{
    int il=0,ir=ifp->Len-1;
    if(strcasecmp(ifp->Key[il],key)>0 || strcasecmp(ifp->Key[ir],key) <0)
        return -1;
    while(il + 1 < ir)
    {
        int mid = (il+ir)/2;
        int caseresult = strcasecmp(ifp->Key[mid],key);
        if(caseresult > 0)
        {
            ir = mid;
        } else if(caseresult < 0)
        {
            il = mid;
        } else
        {
            return mid;
        }
    }
    if(0==strcasecmp(ifp->Key[il],key))
        return il;
    else if(0==strcasecmp(ifp->Key[ir],key))
        return ir;
    else
        return -1;
}

void SearchTest(InputFile * ifp, const char key[])
{
    int r = SearchInput(ifp,key);
    fprintf(stdout,"Expected:%-25s;Searched:%3d|%-25s\n",key,r,r>0?ifp->Key[r]:"NONE");
}

void GetTest(InputFile * ifp)
{
    char TestD[] = "cycle.dt0";
    fprintf(stdout,"Test[%s]=%e\n",TestD, GetValueD(ifp,TestD,"1.0e-100"));
    char TestI[] = "output.ostep";
    fprintf(stdout,"Test[%s]=%d\n",TestI, GetValueI(ifp,TestI,"-100"));
    char TestDk[] = "material.Cvib";
    for(int k=0;k<3;k++)
    {
        fprintf(stdout,"Test[%s|%d]=%e\n",TestDk,k,GetValueDk(ifp,TestDk,k,"1.0e-100"));
    }
    char TestSk[] = "model.layer";
    char tmp[MaxStrLen];
    GetValueSk(ifp,TestSk,tmp,0,"1.0e-100");
    fprintf(stdout,"Test[%s|%d]=%s|\n",TestSk,0,tmp);
    char TestU[] = "cycle.dt000";
    fprintf(stdout,"Test[%s]=%e\n",TestU, GetValueD(ifp,TestU,"1.0e-100"));
}

int GetValueI(InputFile * ifp,const char key[], char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r < 0)
        return strtol(dvalue,NULL,10);
    else
        return strtol(ifp->Value[r],NULL,10);
}

double GetValueD(InputFile * ifp,const char key[], char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r < 0)
        return strtod(dvalue,NULL);
    else
        return strtod(ifp->Value[r],NULL);
}

int GetValueS(InputFile * ifp,const char key[], char value[], char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r < 0)
        return strlen(strcpy(value,dvalue));
    else
        return strlen(strcpy(value,ifp->Value[r]));
}

int GetStrk(char strlist[],int pos, char value[])
{
    int len = strlen(strlist);
    if(('['!=strlist[0]) || (']'!=strlist[len-1]))
    {
        return -1;
    }
    int lpos = 0,rpos=0,k=0;
    for(int i=1;i<len-1;i++)
    {
        if(','==strlist[i])
        {
            lpos++;
            continue;
        }
        if(lpos>pos) break;
        if(lpos==pos)
        {
            value[k++] = strlist[i];
        }
    }
    value[k]='\0';
    return k;
}

double GetValueDk(InputFile * ifp,const char key[], int k,char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r < 0)
        return strtod(dvalue,NULL);
    else
    {
        char tmp[MaxStrLen];
        if(GetStrk(ifp->Value[r],k,tmp))
        {
            return strtod(tmp,NULL);
        } else
        {
            return strtod(dvalue,NULL);
        }
    }
}

int GetValueSk(InputFile * ifp,const char key[], char value[],int k,char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r<0)
        return strlen(strcpy(value,dvalue));
    else
    {
        char tmp[MaxStrLen];
        GetStrk(ifp->Value[r],k,tmp);
        trim(tmp);
        if(strlen(tmp))
        {
            return strlen(strcpy(value,tmp));
        } else
        {
            return strlen(strcpy(value,dvalue));
        }
    }
}

int GetValueIk(InputFile * ifp,const char key[], int k,char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r < 0)
        return strtod(dvalue,NULL);
    else
    {
        char tmp[MaxStrLen];
        if(GetStrk(ifp->Value[r],k,tmp))
        {
            return strtol(tmp,NULL,10);
        } else
        {
            return strtol(dvalue,NULL,10);
        }
    }
}


void SetMeshInfo(InputFile * ifp, struct MeshInfo * _minfo, struct ProcessInfo *_pinfo, struct Mesh * _mesh)
{
   /*
    int myrank;
    MPI_Comm_rank(SALEC_CART_COMM,&myrank);
    if(0==myrank)
        CheckInput(stdout,ifp);
   */

    _pinfo->npgx = GetValueI(ifp,"processor.npgx","2");
    _pinfo->npgy = GetValueI(ifp,"processor.npgy","2");
    _pinfo->npgz = GetValueI(ifp,"processor.npgz","2");
    _pinfo->npgs = _pinfo->npgx*_pinfo->npgy*_pinfo->npgz;
    _pinfo->noffset = GetValueI(ifp,"processor.noffset","2");
    _mesh->noffset = _pinfo->noffset;

    int CommDim[DIM], CommPeriod[DIM];
    CommDim[X] = _pinfo->npgx; CommPeriod[X] = 0;
    CommDim[Y] = _pinfo->npgy; CommPeriod[Y] = 0;
    CommDim[Z] = _pinfo->npgz; CommPeriod[Z] = 0;
    SalecCartCommCreate(CommDim, CommPeriod);

    int mesh_npx = GetValueI(ifp,"mesh.npx","32");
    int mesh_npy = GetValueI(ifp,"mesh.npy","32");
    int mesh_npz = GetValueI(ifp,"mesh.npz","32");

    int Npg = _pinfo->npgx * _pinfo->npgy * _pinfo->npgz;
    _pinfo->nwp    = (int **) malloc(sizeof(int*)*Npg);
    _pinfo->BaseId = (int **) malloc(sizeof(int*)*Npg);
    for(int k=0;k<Npg;k++)
    {
        _pinfo->nwp[k]    = (int *) malloc(sizeof(int)*DIM);
        _pinfo->BaseId[k] = (int *) malloc(sizeof(int)*DIM);
    }

    _pinfo->Nxx = mesh_npx*_pinfo->npgx + 2*_pinfo->noffset;
    _pinfo->Nyy = mesh_npy*_pinfo->npgy + 2*_pinfo->noffset;
    _pinfo->Nzz = mesh_npz*_pinfo->npgz + 2*_pinfo->noffset;

    int nwpx = (_pinfo->Nxx - 2*_pinfo->noffset) / _pinfo->npgx;
    int nwpy = (_pinfo->Nyy - 2*_pinfo->noffset) / _pinfo->npgy;
    int nwpz = (_pinfo->Nzz - 2*_pinfo->noffset) / _pinfo->npgz;

    for(int igx=0;igx<_pinfo->npgx;igx++)
    {
        for(int igy=0;igy<_pinfo->npgy;igy++)
        {
            for(int igz=0;igz<_pinfo->npgz;igz++)
            {
                int tId = ToLongId(igx,igy,igz,_pinfo->npgx,_pinfo->npgy,_pinfo->npgz);
                _pinfo->nwp[tId][X]    = nwpx;
                _pinfo->nwp[tId][Y]    = nwpy;
                _pinfo->nwp[tId][Z]    = nwpz;

                _pinfo->BaseId[tId][X] = nwpx*igx;
                _pinfo->BaseId[tId][Y] = nwpy*igy;
                _pinfo->BaseId[tId][Z] = nwpz*igz;
            }
        }
    }

    int lrank,bIdx,bIdy,bIdz;
    MPI_Comm_rank(SALEC_CART_COMM,&lrank);
    _pinfo->rank = lrank;

    char t_pname[100];
    int t_length;
    MPI_Get_processor_name(_pinfo->hostname,&t_length);

    _minfo->npx = _pinfo->nwp[lrank][X] + 2*_pinfo->noffset + 1;
    _minfo->npy = _pinfo->nwp[lrank][Y] + 2*_pinfo->noffset + 1;
    _minfo->npz = _pinfo->nwp[lrank][Z] + 2*_pinfo->noffset + 1;
    bIdx = _pinfo->BaseId[lrank][X];
    bIdy = _pinfo->BaseId[lrank][Y];
    bIdz = _pinfo->BaseId[lrank][Z];

    _minfo->ARTVIS   = GetValueD(ifp,"numerical.ARTVIS","0.5");
    _minfo->ARTVIS2  = GetValueD(ifp,"numerical.ARTVIS2","0.5");
    _minfo->CISSAMK1 = GetValueD(ifp,"numerical.CISSAMK1","0.5");
    _minfo->PartPressure = GetValueI(ifp,"numerical.PARTPRES","0");
    _minfo->TensileFailure = GetValueI(ifp,"numerical.TENSILE","0");

    double dx = GetValueD(ifp,"mesh.dx","800.0");
    double dy = GetValueD(ifp,"mesh.dy","800.0");
    double dz = GetValueD(ifp,"mesh.dz","800.0");

    char GridExtOpt[MaxStrLen];
    GetValueS(ifp,"mesh.ext",GridExtOpt,"off");
    if(strcasecmp(GridExtOpt,"off")==0)
    {
        double Lx = mesh_npx * _pinfo->npgx * dx;
        double Ly = mesh_npy * _pinfo->npgy * dy;
        double Lz = mesh_npz * _pinfo->npgz * dz;

        double pX0 = -1.0*Lx* GetValueDk(ifp,"mesh.O",0,"0.5") - _pinfo->noffset*dx;
        double pY0 = -1.0*Ly* GetValueDk(ifp,"mesh.O",1,"0.5") - _pinfo->noffset*dy;
        double pZ0 = -1.0*Lz* GetValueDk(ifp,"mesh.O",2,"0.5") - _pinfo->noffset*dz;

        _minfo->V = (double ****) malloc(_minfo->npz * sizeof(double ***));
        for (int k = 0; k < _minfo->npz; k++)
        {
            _minfo->V[k] = (double ***) malloc(_minfo->npy * sizeof(double **));
            for (int j = 0; j < _minfo->npy; j++)
            {
                _minfo->V[k][j] = (double **) malloc(_minfo->npx * sizeof(double *));
                for (int i = 0; i < _minfo->npx; i++)
                {
                    _minfo->V[k][j][i] = (double *) malloc(DIM * sizeof(double));
                    _minfo->V[k][j][i][X] = pX0 + (i+bIdx) * dx;
                    _minfo->V[k][j][i][Y] = pY0 + (j+bIdy) * dy;
                    _minfo->V[k][j][i][Z] = pZ0 + (k+bIdz) * dz;
                }
            }
        }

        char DampOpt[MaxStrLen];
        GetValueS(ifp,"numerical.damping",DampOpt,"off");
        if(strcasecmp(DampOpt,"off")==0)
        {
            _minfo->DampLambda = -1.0;
        } else
        {
            int DampZk = GetValueIk(ifp,"numerical.LAMBDA",0,"10");
            _minfo->DampLambda = GetValueDk(ifp,"numerical.LAMBDA",1,"0.5");
            _minfo->DampZ0 = _pinfo->noffset*dz + pZ0;
            _minfo->DampZ1 = (_pinfo->noffset + DampZk)*dz + pZ0;
        }

    } else if(strcasecmp(GridExtOpt,"on")==0)
    {
        double Xext = GetValueDk(ifp,"mesh.ex",0,"1.05");
        int    XeL  = GetValueIk(ifp,"mesh.ex",1,"0");
        int    XeR  = GetValueIk(ifp,"mesh.ex",2,"0");

        double Yext = GetValueDk(ifp,"mesh.ey",0,"1.05");
        int    YeL  = GetValueIk(ifp,"mesh.ey",1,"0");
        int    YeR  = GetValueIk(ifp,"mesh.ey",2,"0");

        double Zext = GetValueDk(ifp,"mesh.ez",0,"1.05");
        int    ZeL  = GetValueIk(ifp,"mesh.ez",1,"0");
        int    ZeR  = GetValueIk(ifp,"mesh.ez",2,"0");

//        double Lx = Spacing(_pinfo->Nxx-1,Xext,XeL,XeR,_pinfo->Nxx)*dx;
//        double Ly = Spacing(_pinfo->Nyy-1,Yext,YeL,YeR,_pinfo->Nyy)*dy;
//        double Lz = Spacing(_pinfo->Nzz-1,Zext,ZeL,ZeR,_pinfo->Nzz)*dz;

        double Lx = Spacing(_pinfo->Nxx-1-_pinfo->noffset,Xext,XeL,XeR,_pinfo->Nxx)*dx
                - Spacing(_pinfo->noffset,Xext,XeL,XeR,_pinfo->Nxx)*dx;
        double Ly = Spacing(_pinfo->Nyy-1-_pinfo->noffset,Yext,YeL,YeR,_pinfo->Nyy)*dy
                - Spacing(_pinfo->noffset,Yext,YeL,YeR,_pinfo->Nyy)*dy;
        double Lz = Spacing(_pinfo->Nzz-1-_pinfo->noffset,Zext,ZeL,ZeR,_pinfo->Nzz)*dz
                - Spacing(_pinfo->noffset,Zext,ZeL,ZeR,_pinfo->Nzz)*dz;

        double pX0 = -1.0*Lx* GetValueDk(ifp,"mesh.O",0,"0.5");
        double pY0 = -1.0*Ly* GetValueDk(ifp,"mesh.O",1,"0.5");
        double pZ0 = -1.0*Lz* GetValueDk(ifp,"mesh.O",2,"0.5");

        pX0 -= Spacing(_pinfo->noffset,Xext,XeL,XeR,_pinfo->Nxx)*dx;
        pY0 -= Spacing(_pinfo->noffset,Yext,YeL,YeR,_pinfo->Nyy)*dy;
        pZ0 -= Spacing(_pinfo->noffset,Zext,ZeL,ZeR,_pinfo->Nzz)*dz;

        _minfo->V = (double ****) malloc(_minfo->npz * sizeof(double ***));
        for (int k = 0; k < _minfo->npz; k++)
        {
            _minfo->V[k] = (double ***) malloc(_minfo->npy * sizeof(double **));
            for (int j = 0; j < _minfo->npy; j++)
            {
                _minfo->V[k][j] = (double **) malloc(_minfo->npx * sizeof(double *));
                for (int i = 0; i < _minfo->npx; i++)
                {
                    _minfo->V[k][j][i] = (double *) malloc(DIM * sizeof(double));
                    _minfo->V[k][j][i][X] = pX0 + Spacing(i+bIdx,Xext,XeL,XeR,_pinfo->Nxx) * dx;
                    _minfo->V[k][j][i][Y] = pY0 + Spacing(j+bIdy,Yext,YeL,YeR,_pinfo->Nyy) * dy;
                    _minfo->V[k][j][i][Z] = pZ0 + Spacing(k+bIdz,Zext,ZeL,ZeR,_pinfo->Nzz) * dz;
                }
            }
        }

        char DampOpt[MaxStrLen];
        GetValueS(ifp,"numerical.damping",DampOpt,"off");
        if(strcasecmp(DampOpt,"off")==0)
        {
            _minfo->DampLambda = -1.0;
        } else
        {
            int DampZk = GetValueIk(ifp,"numerical.LAMBDA",0,"10");
            _minfo->DampLambda = GetValueDk(ifp,"numerical.LAMBDA",0,"0.5");
            _minfo->DampZ0 = pZ0 + Spacing(_pinfo->noffset,Zext,ZeL,ZeR,_pinfo->Nzz) * dz;
            _minfo->DampZ1 = pZ0 + Spacing(_pinfo->noffset+DampZk,Zext,ZeL,ZeR,_pinfo->Nzz) * dz;
        }

    } else
    {
        fprintf(stdout,"Undefine Value in mesh.ext!\n");
        exit(0);
    }


    _minfo->MaterialNum = GetValueI(ifp,"material.nm","3");
    char MaterialName[NMT][30];
    strcpy(MaterialName[0],"../eos/basalt_.aneos");

    for(int iMat=1;iMat<_minfo->MaterialNum;iMat++)
    {
        char postfix[MaxStrLen];
        char matname[MaxStrLen];
        GetValueSk(ifp,"material.name",matname,iMat,"basalt_");
        GetValueSk(ifp,"material.postfix",postfix,iMat,"aneos");
        sprintf(MaterialName[iMat],"../eos/%s.%s",matname,postfix);
        strcpy(_minfo->EosType[iMat],postfix);
    }

    for(int i=1;i<_minfo->MaterialNum;i++)
    {
        _minfo->Materialfp[i] = fopen(MaterialName[i],"r");
        if(NULL == _minfo->Materialfp[i])
        {
            fprintf(stdout,"Cannot open %s\n",MaterialName[i]);
            exit(0);
        }
    }

    GetValueS(ifp,"output.format",_minfo->output_format,"binary");
    _minfo->MaxVelocity = GetValueD(ifp,"numerical.MaxVelocity","2.0e4");

}


void SRFInit(InputFile * ifp,struct StateReference * _s, int iMat)
{
    // (ROCK)
    _s->Yint0    = GetValueDk(ifp,"material.yint0",iMat,"1.0");
    _s->Yfricint = GetValueDk(ifp,"material.yintfri",iMat,"1.0");
    _s->Ylimint  = GetValueDk(ifp,"material.yintlim",iMat,"2.50e9");
    _s->Ydam0    = GetValueDk(ifp,"material.ydam0",iMat,"1.0");
    _s->Yfricdam = GetValueDk(ifp,"material.ydamfri",iMat,"1.0");
    _s->Ylimdam  = GetValueDk(ifp,"material.ydamlim",iMat,"2.50e9");
    char strengthmodel[MaxStrLen];
    GetValueSk(ifp,"material.yshear",strengthmodel,iMat,"SimpleRock");
    _s->Yfunc = SelectStrength(strengthmodel);

    // (IVANOV)
    _s->IvanA = GetValueDk(ifp,"material.IvanA",iMat,"1.0e-4");
    _s->IvanB = GetValueDk(ifp,"material.IvanB",iMat,"1.0e-11");
    _s->IvanC = GetValueDk(ifp,"material.IvanC",iMat,"3.0e8");
    GetValueSk(ifp,"material.damage",strengthmodel,iMat,"SimpleShear");
    _s->dSDf = SelectDamage(strengthmodel);

    // (LINEAR) tensile

    _s->Yten0 = GetValueDk(ifp,"material.yten0",iMat,"3.0e8");
    GetValueSk(ifp,"material.ytens",strengthmodel,iMat,"SimpleTensile");
    _s->Tfunc = SelectStrength(strengthmodel);

    // (Other)
    double Poisson = GetValueDk(ifp,"material.Poisson",iMat,"0.3");
    _s->GaCoef = 1.5 * (1.0 - 2 * Poisson) / (1.0 + Poisson);

    // set the parameters related to pressure_min : ROCK
    _s->minPint = _s->Yint0/_s->Yfricint * (_s->Yint0/_s->Ylimint - 1.0);
    _s->minPdam = -1.0*_s->Ydam0/_s->Yfricdam;

    // set acfl parameters
    double BlockSize = GetValueDk(ifp,"material.BlockSize",iMat,"8000.0");
    _s->GammaEta  = GetValueDk(ifp, "material.GammaEta", iMat, "8.0e-3");
    _s->GammaBeta = GetValueDk(ifp, "material.GammaBeta", iMat, "1.15e2");
    _s->Toff     = GetValueDk(ifp,"material.Toff",iMat,"16.0");
    _s->Cvib     = GetValueDk(ifp,"material.Cvib",iMat,"0.1");
    _s->VibMax   = GetValueDk(ifp,"material.VibMax",iMat,"200.0");
    _s->Pvlim    = GetValueDk(ifp,"material.Pvlim",iMat,"2.50e10");
    _s->Acvis    = _s->GammaEta * BlockSize * 5000.0;
    _s->Tdamp    = _s->GammaBeta * BlockSize / 5000.0;

    // set Johson-Cook strength parameters
    _s->JCA      = GetValueDk(ifp,"material.JcA",iMat,"4.90e7");
    _s->JCB      = GetValueDk(ifp,"material.JcB",iMat,"1.57e8");
    _s->JCN      = GetValueDk(ifp,"material.JcN",iMat,"1.67e-1");
    _s->JCC      = GetValueDk(ifp,"material.JcC",iMat,"1.60e-2");
    _s->JCM      = GetValueDk(ifp,"material.JcM",iMat,"1.70");
    _s->JCTREF   = GetValueDk(ifp,"material.JcTref",iMat,"8.00e2");

    // set thermal temperature parameters
    _s->Asimon   = 1.0/GetValueDk(ifp,"material.SimonA",iMat,"6.00e9");
    _s->Csimon   = 1.0/GetValueDk(ifp,"material.SimonC",iMat,"3.00");
    _s->Tmelt0   = GetValueDk(ifp,"material.SimonT0",iMat,"933.0");
    _s->Tfrac    = GetValueDk(ifp,"material.OhnakaXi",iMat,"0.5");

    char StrengthModel[100];
    GetValueSk(ifp,"material.yshear",StrengthModel,iMat,"SimpleRock");
    if((0==strcasecmp(StrengthModel,"JohnsonCook2"))||(0==strcasecmp(StrengthModel,"JohnsonCook1")))
    {
        _s->minPdam = _s->minPint = GetValueDk(ifp,"material.JcminP",iMat,"-2.44e9");
        _s->JCMaxPlastic = GetValueDk(ifp,"material.JcMaxPlastic",iMat,"10.0");
    }
}

int SelectBC(char s[], double * tr, int *type)
{
    if(strcasecmp(s,"outflow")==0)
    {
        *type = 1;
        *tr   = 0.0;
    } else if(strcasecmp(s,"freeslip")==0)
    {
        *type = 1;
        *tr   = -2.0;
    } else if(strcasecmp(s,"fixed")==0)
    {
        *type = 0;
        *tr   = 0.0;
    } else
    {
        *type = 0;
        *tr   = 0.0;
        fprintf(stdout,"undefined boundary type: %s\n",s);
        exit(0);
    }
    return *type;
}

const char BoundaryName[BPE][MaxStrLen] = {
        "boundary.left",
        "boundary.right",
        "boundary.front",
        "boundary.back",
        "boundary.bottom",
        "boundary.top"
};

const double BoundaryNorm[BPE][DIM] = {
        { 0.0, 1.0, 0.0},
        { 0.0,-1.0, 0.0},
        {-1.0, 0.0, 0.0},
        { 1.0, 0.0, 0.0},
        { 0.0, 0.0, 1.0},
        { 0.0, 0.0,-1.0}
};

void SetBC(InputFile * ifp, struct Mesh * _mesh)
{
    double tmpTr[BPE];
    int tmpType[BPE];
    for(int k=0;k<BPE;k++)
    {
        char tBC[MaxStrLen];
        GetValueS(ifp,BoundaryName[k],tBC,"fixed");
        SelectBC(tBC,tmpTr+k,tmpType+k);
    }

    for(int k=0;k<_mesh->nce;k++)
    {
        struct ConditionE * ice = _mesh->CondE + k;
        struct Element * dest = ice->dest;
        struct Element * refe = ice->refe;
        double vLine[DIM] = {0.0}, mPoint[DIM]={0.0};
        LinearOp(vLine,dest->Center,refe->Center,-1.0,1.0);
        LinearOp(mPoint,dest->Center,refe->Center,0.5,0.5);
        Normalization(vLine);
        Copy(vLine,ice->Wnorm[0]);

        int ib = 0;
        for(int j=1;j<BPE;j++)
        {
            if(Dot(BoundaryNorm[j],vLine) > Dot(BoundaryNorm[ib],vLine))
                ib = j;
        }

        ice->Tr[0] = tmpTr[ib];
        ice->CEType = tmpType[ib];
    }
}


void SetDamping(struct MeshInfo * _minfo, struct Mesh * _mesh)
{
    for(int k=0;k<_mesh->ndv;k++)
    {
        struct Vertex * iv = _mesh->DomainV[k];
        double ivZ = iv->Position[Z];
        if(_minfo->DampLambda>0.0 && ivZ >= _minfo->DampZ0 && ivZ <=_minfo->DampZ1)
        {
            double tmp = (_minfo->DampZ1 - ivZ)/(_minfo->DampZ1 - _minfo->DampZ0);
            iv->Damping = _minfo->DampLambda * pow(tmp,4.0);
        } else
        {
            iv->Damping = 0.0;
        }
    }
}


void InitMeshANEOS(InputFile * ifp, struct Mesh * _mesh, struct ANEOSTable * ATB)
{
    char Material[NMT][MaxStrLen];
    int NumMaterial = GetValueI(ifp,"material.nm","1.0");
    for(int iMat=0;iMat<NumMaterial;iMat++)
    {
        GetValueSk(ifp,"material.name",Material[iMat],iMat,"dunite_");
    }
    // process the projectile
    double speed = GetValueDk(ifp,"projectile.velocity",0,"1.0e4");
    double theta = M_PI/180.0*GetValueDk(ifp,"projectile.velocity",1,"0.0");
    double phi   = M_PI/180.0*GetValueDk(ifp,"projectile.velocity",2,"0.0");
    double Radiu = GetValueD(ifp,"projectile.radiu","7000.0");
    double ImpactVelocity[DIM],Center[DIM];
    ImpactVelocity[X] = speed * sin(theta) * cos(phi);
    ImpactVelocity[Y] = speed * sin(theta) * sin(phi);
    ImpactVelocity[Z] = speed * cos(theta);
    Center[X] = GetValueDk(ifp,"projectile.center",X,"0.0");
    Center[Y] = GetValueDk(ifp,"projectile.center",Y,"0.0");
    Center[Z] = GetValueDk(ifp,"projectile.center",Z,"0.0");

    char iName[MaxStrLen],pDamage[MaxStrLen],pTem[MaxStrLen],pPre[MaxStrLen],pDen[MaxStrLen];
    double iDamage=0.0,iTem=0.0,iPre=0.0,iDen=0.0,iEng=0.0;
    GetValueS(ifp,"projectile.material",iName,"dunite_");
    GetValueSk(ifp,"projectile.damage",pDamage,0,"unknown");
    GetValueSk(ifp,"projectile.temperature",pTem,0,"unknown");
    GetValueSk(ifp,"projectile.pressure",pPre,0,"unknown");
    GetValueSk(ifp,"projectile.density",pDen,0,"unknown");

    int IdImpact = 1;
    for(int iMat=1;iMat<NumMaterial;iMat++)
    {
        if(0==strcasecmp(iName,Material[iMat])) IdImpact = iMat;
    }
    if(0==strcasecmp(pDamage,"const"))
        iDamage = GetValueDk(ifp,"projectile.damge",1,"1.0");
    if((0==strcasecmp(pTem,"const")) && (0==strcasecmp(pPre,"const")) && (0==strcasecmp(pDen,"derived")))
    {
        iPre = GetValueDk(ifp,"projectile.pressure",1,"1.0");
        iTem = GetValueDk(ifp,"projectile.temperature",1,"1.0");
        iDen = ANEOSInterpolateTP(ATB + IdImpact, iTem, iPre, -1);
        iEng = ANEOSInterpolateTP(ATB + IdImpact, iTem, iPre, 0);
    }

    int NumTarget = GetValueI(ifp,"target.number","1");
    double dz = 0.5*GetValueD(ifp,"mesh.dz","1.0e4");
    double depth[NMT],SumDepth;
    int IdTarget[NMT];
    for(int iTar=0;iTar<NumTarget;iTar++)
    {
        char TarName[MaxStrLen];
        GetValueSk(ifp,"target.material",TarName,iTar,"unknown");
        IdTarget[iTar] = -1;
        for(int iMat=1;iMat<NumMaterial;iMat++)
        {
            if(0==strcasecmp(TarName,Material[iMat]))
                IdTarget[iTar] = iMat;
        }
        if(-1==IdTarget[iTar])
        {
            fprintf(stdout,"unknown material in target %s!\n",TarName);
            exit(0);
        }
        depth[iTar] = GetValueDk(ifp,"target.depth",iTar,"0.0");
        SumDepth += depth[iTar];
    }

    int ProfileLen = (int)(SumDepth/dz) + 1;
    double * tGrav = (double *) malloc(ProfileLen* sizeof(double ));
    double * tTemp = (double *) malloc(ProfileLen* sizeof(double ));
    double * tPres = (double *) malloc(ProfileLen* sizeof(double ));
    double * tDens = (double *) malloc(ProfileLen* sizeof(double ));
    double * tEng  = (double *) malloc(ProfileLen* sizeof(double ));
    double * tDam  = (double *) malloc(ProfileLen* sizeof(double ));
    int * profId   = (int *) malloc(ProfileLen* sizeof(int ));

    char TarGrav[MaxStrLen],TarTemp[MaxStrLen],TarPres[MaxStrLen],TarDam[MaxStrLen];
    double gravity = 0.0,temperature=0.0;
    GetValueSk(ifp,"condition.gravity",TarGrav,0,"unknown");
    GetValueSk(ifp,"target.temperature",TarTemp,0,"unknown");
    GetValueSk(ifp,"target.damage",TarDam,0,"unknown");
    GetValueSk(ifp,"target.pressure",TarPres,0,"unknown");

    if(0==strcasecmp(TarGrav,"const"))
    {
        gravity = GetValueDk(ifp,"condition.gravity",1,"-9.8");
        for(int k=0;k<ProfileLen;k++) tGrav[k] = gravity;
    }
    if(0==strcasecmp(TarTemp,"const"))
    {
        temperature = GetValueDk(ifp,"target.temperature",1,"295.0");
        for(int k=0;k<ProfileLen;k++) tTemp[k] = temperature;
    }
    if(0==strcasecmp(TarDam,"const"))
    {
        double Damage = GetValueDk(ifp,"target.damage",1,"0.0");
        for(int k=0;k<ProfileLen;k++) tDam[k] = Damage;
    }

    if(0==strcasecmp(TarPres,"derived"))
    {
        tPres[0] = 1.0;
        int lnd = 0;
        profId[0] = IdTarget[0];
        for(int iTar=0;iTar<NumTarget;iTar++)
        {
            int iInterval = (int)(depth[iTar]/dz+1);
            ANEOSPresProfRK3(ATB + IdTarget[iTar], tPres + lnd, tGrav + lnd, tTemp + lnd, iInterval, dz);
            for(int j=1;j<=iInterval-1;j++)
            {
                profId[j + lnd] = IdTarget[iTar];
                tDens[j+lnd] = ANEOSInterpolateTP(ATB + IdTarget[iTar], tTemp[j + lnd], tPres[j + lnd], -1);
                tEng[j+lnd] = ANEOSInterpolateTP(ATB + IdTarget[iTar], tTemp[j + lnd], tPres[j + lnd], 0);
            }
            lnd += iInterval-1;
        }
    }

    double toplevel = GetValueD(ifp,"target.toplevel","-1000.0");

    for(int i=0;i<_mesh->ne;i++)
    {
        struct Element * ie = _mesh->Elements + i;
        double * VibVelocity = ie->PsiL3 + 8;
        ZeroEx(ie->VOF,NMT);
        ZeroEx(ie->sPsiL3,NPSI);
        ZeroEx(ie->PsiL3,NPSI);
        ie->VibPressure = *VibVelocity = 0.0;
        for(int k=0;k<NMT;k++)
        {
            ZeroEx(ie->PhiL3[k],NPHI);
            ZeroEx(ie->sPhiL3[k],NPHI);
            ie->PhiL3[k][0] = 1.0;
        }
        if(Distance(Center,ie->Center) <= Radiu)
        {
            ie->Temperature = iTem;
            ie->Pressure = iPre;
            ie->Density  = iDen;

            ie->VOF[IdImpact] = 1.0;
            ie->PhiL3[IdImpact][1] = iDen;
            ie->PhiL3[IdImpact][2] = iEng;

            ScalerMoveEx(ie->sPhiL3[IdImpact],ie->PhiL3[IdImpact],1.0,NPHI);
            ie->Damage = ie->PsiL3[0] = iDamage;
            ie->Mass = ie->Density * ie->Volume;
            ie->sPsiL3[0] = ie->Mass*iDamage;
            ScalerMove(ie->Momentum,ImpactVelocity,ie->Mass);
        } else if(ie->Center[Z] < toplevel)
        {
            int dpz = (int) Wind((toplevel - ie->Center[Z])/dz,0.0,ProfileLen-1);
            ie->Temperature = tTemp[dpz];
            ie->Pressure    = tPres[dpz];
            ie->Density     = tDens[dpz];

            ie->VOF[profId[dpz]] = 1.0;
            ie->PhiL3[profId[dpz]][1] = tDens[dpz];
            ie->PhiL3[profId[dpz]][2] = tEng[dpz];

            ScalerMoveEx(ie->sPhiL3[profId[dpz]],ie->PhiL3[profId[dpz]],1.0,NPHI);
            ie->Damage = ie->PsiL3[0] = tDam[dpz];
            ie->Mass = ie->Density * ie->Volume;
            ie->sPsiL3[0] = ie->Mass*ie->Damage;
            Zero(ie->Momentum);
        } else
        {
            ie->Temperature = 0.0;
            ie->Pressure    = 0.0;
            ie->Density     = 0.0;
            ie->VOF[0] = 1.0;
            ie->sPhiL3[0][0] = 1.0;
            ie->Mass = 0.0;
            Zero(ie->Momentum);
        }
    }

    for(int k=0;k<_mesh->nv;k++)
    {
        struct Vertex * iv = _mesh->Vertexes + k;
        Zero(iv->BodyF);
        iv->BodyF[Z] = gravity;
    }

    for(int k=0;k<_mesh->ndb;k++)
    {
        struct Boundary * ib = _mesh->DomainB[k];
        ib->dMass = 0.0;
        Zero(ib->dMomentum);
        ZeroEx(ib->dPsi,NPSI);
        for(int j=0;j<NMT;j++)
        {
            ZeroEx(ib->dPhi[j],NPHI);
        }
    }

    free(tGrav);
    free(tTemp);
    free(tPres);
    free(tDens);
    free(tEng);
    free(tDam);
    free(profId);
}


void InitMesh(InputFile * ifp, struct Mesh * _mesh, struct EosTable * _etb)
{
    char Material[NMT][MaxStrLen];
    int NumMaterial = GetValueI(ifp,"material.nm","1.0");
    for(int iMat=0;iMat<NumMaterial;iMat++)
    {
        GetValueSk(ifp,"material.name",Material[iMat],iMat,"dunite_");
    }
    // process the projectile
    double speed = GetValueDk(ifp,"projectile.velocity",0,"1.0e4");
    double theta = M_PI/180.0*GetValueDk(ifp,"projectile.velocity",1,"0.0");
    double phi   = M_PI/180.0*GetValueDk(ifp,"projectile.velocity",2,"0.0");
    double Radiu = GetValueD(ifp,"projectile.radiu","7000.0");
    double ImpactVelocity[DIM],Center[DIM];
    ImpactVelocity[X] = speed * sin(theta) * cos(phi);
    ImpactVelocity[Y] = speed * sin(theta) * sin(phi);
    ImpactVelocity[Z] = speed * cos(theta);
    Center[X] = GetValueDk(ifp,"projectile.center",X,"0.0");
    Center[Y] = GetValueDk(ifp,"projectile.center",Y,"0.0");
    Center[Z] = GetValueDk(ifp,"projectile.center",Z,"0.0");

    char iName[MaxStrLen],pDamage[MaxStrLen],pTem[MaxStrLen],pPre[MaxStrLen],pDen[MaxStrLen];
    double iDamage=0.0,iTem=0.0,iPre=0.0,iDen=0.0,iEng=0.0;
    GetValueS(ifp,"projectile.material",iName,"dunite_");
    GetValueSk(ifp,"projectile.damage",pDamage,0,"unknown");
    GetValueSk(ifp,"projectile.temperature",pTem,0,"unknown");
    GetValueSk(ifp,"projectile.pressure",pPre,0,"unknown");
    GetValueSk(ifp,"projectile.density",pDen,0,"unknown");

    int IdImpact = 1;
    for(int iMat=1;iMat<NumMaterial;iMat++)
    {
        if(0==strcasecmp(iName,Material[iMat])) IdImpact = iMat;
    }
    if(0==strcasecmp(pDamage,"const"))
        iDamage = GetValueDk(ifp,"projectile.damge",1,"1.0");
    if((0==strcasecmp(pTem,"const")) && (0==strcasecmp(pPre,"const")) && (0==strcasecmp(pDen,"derived")))
    {
        iPre = GetValueDk(ifp,"projectile.pressure",1,"1.0");
        iTem = GetValueDk(ifp,"projectile.temperature",1,"1.0");
        iDen = _etb[IdImpact].InterpolateTP(_etb + IdImpact, iTem, iPre, -1);
        iEng = _etb[IdImpact].InterpolateTP(_etb + IdImpact, iTem, iPre, 0);
    }

    int NumTarget = GetValueI(ifp,"target.number","1");
    double dz = 0.5*GetValueD(ifp,"mesh.dz","1.0e4");
    double depth[NMT],SumDepth;
    int IdTarget[NMT];
    for(int iTar=0;iTar<NumTarget;iTar++)
    {
        char TarName[MaxStrLen];
        GetValueSk(ifp,"target.material",TarName,iTar,"unknown");
        IdTarget[iTar] = -1;
        for(int iMat=1;iMat<NumMaterial;iMat++)
        {
            if(0==strcasecmp(TarName,Material[iMat]))
                IdTarget[iTar] = iMat;
        }
        if(-1==IdTarget[iTar])
        {
            fprintf(stdout,"unknown material in target %s!\n",TarName);
            exit(0);
        }
        depth[iTar] = GetValueDk(ifp,"target.depth",iTar,"0.0");
        SumDepth += depth[iTar];
    }

    int ProfileLen = (int)(SumDepth/dz) + 1;
    double * tGrav = (double *) malloc(ProfileLen* sizeof(double ));
    double * tTemp = (double *) malloc(ProfileLen* sizeof(double ));
    double * tPres = (double *) malloc(ProfileLen* sizeof(double ));
    double * tDens = (double *) malloc(ProfileLen* sizeof(double ));
    double * tEng  = (double *) malloc(ProfileLen* sizeof(double ));
    double * tDam  = (double *) malloc(ProfileLen* sizeof(double ));
    int * profId   = (int *) malloc(ProfileLen* sizeof(int ));

    char TarGrav[MaxStrLen],TarTemp[MaxStrLen],TarPres[MaxStrLen],TarDam[MaxStrLen];
    double gravity = 0.0,temperature=0.0;
    GetValueSk(ifp,"condition.gravity",TarGrav,0,"unknown");
    GetValueSk(ifp,"target.temperature",TarTemp,0,"unknown");
    GetValueSk(ifp,"target.damage",TarDam,0,"unknown");
    GetValueSk(ifp,"target.pressure",TarPres,0,"unknown");

    if(0==strcasecmp(TarGrav,"const"))
    {
        gravity = GetValueDk(ifp,"condition.gravity",1,"-9.8");
        for(int k=0;k<ProfileLen;k++) tGrav[k] = gravity;
    } else
    {
        fprintf(stdout,"error gravity type:%s\n",TarGrav);
    }

    if(0==strcasecmp(TarTemp,"const"))
    {
        temperature = GetValueDk(ifp,"target.temperature",1,"295.0");
        for(int k=0;k<ProfileLen;k++) tTemp[k] = temperature;
    }
    if(0==strcasecmp(TarDam,"const"))
    {
        double Damage = GetValueDk(ifp,"target.damage",1,"0.0");
        for(int k=0;k<ProfileLen;k++) tDam[k] = Damage;
    }

    if(0==strcasecmp(TarPres,"derived"))
    {
        tPres[0] = 1.0;
        int lnd = 0;
        profId[0] = IdTarget[0];
        for(int iTar=0;iTar<NumTarget;iTar++)
        {
            int iInterval = (int)(depth[iTar]/dz+1);
            PresProfRK3(_etb + IdTarget[iTar], tPres + lnd, tGrav + lnd, tTemp + lnd, iInterval, dz);
            for(int j=1;j<=iInterval-1;j++)
            {
                profId[j + lnd] = IdTarget[iTar];
                tDens[j+lnd] = _etb[IdTarget[iTar]].InterpolateTP(_etb + IdTarget[iTar], tTemp[j + lnd], tPres[j + lnd], -1);
                tEng[j+lnd] = _etb[IdTarget[iTar]].InterpolateTP(_etb + IdTarget[iTar], tTemp[j + lnd], tPres[j + lnd], 0);
            }
            lnd += iInterval-1;
        }

    }

    double toplevel = GetValueD(ifp,"target.toplevel","-1000.0");
    double topnorm[3] = {0.,0.,1.};
    char TopTypeOpt[MaxStrLen];
    GetValueSk(ifp,"target.toptype",TopTypeOpt,0,"plane");

    if(0== strcasecmp(TopTypeOpt,"plane"))
    {
        /*
         * do nothing, as special case of slope
         */

        for(int i=0;i<_mesh->ne;i++)
        {
            struct Element * ie = _mesh->Elements + i;
            double * VibVelocity = ie->PsiL3 + 8;
            ie->Viscosity = 0.0;
            ZeroEx(ie->VOF,NMT);
            ZeroEx(ie->sPsiL3,NPSI);
            ZeroEx(ie->PsiL3,NPSI);
            ie->VibPressure = *VibVelocity = 0.0;
            for(int k=0;k<NMT;k++)
            {
                ZeroEx(ie->PhiL3[k],NPHI);
                ZeroEx(ie->sPhiL3[k],NPHI);
                ie->PhiL3[k][0] = 1.0;
            }
            if(Distance(Center,ie->Center) <= Radiu)
            {
                ie->Temperature = iTem;
                ie->Pressure = iPre;
                ie->Density  = iDen;

                ie->VOF[IdImpact] = 1.0;
                ie->PhiL3[IdImpact][1] = iDen;
                ie->PhiL3[IdImpact][2] = iEng;

                ScalerMoveEx(ie->sPhiL3[IdImpact],ie->PhiL3[IdImpact],1.0,NPHI);
                ie->Damage = ie->PsiL3[0] = iDamage;
                ie->Mass = ie->Density * ie->Volume;
                ie->sPsiL3[0] = ie->Mass*iDamage;
                ScalerMove(ie->Momentum,ImpactVelocity,ie->Mass);
            } else if(ie->Center[Z] < toplevel)
            {
                int dpz = (int) (Wind((toplevel - ie->Center[Z])/dz,0,ProfileLen-1));
                ie->Temperature = tTemp[dpz];
                ie->Pressure    = tPres[dpz];
                ie->Density     = tDens[dpz];

                ie->VOF[profId[dpz]] = 1.0;
                ie->PhiL3[profId[dpz]][1] = tDens[dpz];
                ie->PhiL3[profId[dpz]][2] = tEng[dpz];

                ScalerMoveEx(ie->sPhiL3[profId[dpz]],ie->PhiL3[profId[dpz]],1.0,NPHI);
                ie->Damage = ie->PsiL3[0] = tDam[dpz];
                ie->Mass = ie->Density * ie->Volume;
                ie->sPsiL3[0] = ie->Mass*ie->Damage;
                Zero(ie->Momentum);
            } else
            {
                ie->Temperature = 0.0;
                ie->Pressure    = 0.0;
                ie->Density     = 0.0;
                ie->VOF[0] = 1.0;
                ie->sPhiL3[0][0] = 1.0;
                ie->Mass = 0.0;
                Zero(ie->Momentum);
            }
        }

    } else if(0 == strcasecmp(TopTypeOpt,"slope"))
    {
        /*
         * theta,phi (degree)
         */
        double SlopeTheta = M_PI/180.0*GetValueDk(ifp,"target.toptype",1,"0.0");
        double SlopePhi   = M_PI/180.0*GetValueDk(ifp,"target.toptype",2,"0.0");
        topnorm[0] = sin(SlopeTheta)* cos(SlopePhi);
        topnorm[1] = sin(SlopeTheta)* sin(SlopePhi);
        topnorm[2] = cos(SlopeTheta);

        for(int i=0;i<_mesh->ne;i++)
        {
            struct Element * ie = _mesh->Elements + i;
            double * VibVelocity = ie->PsiL3 + 8;
            double topDistance = Dot(topnorm,ie->Center);
            ie->Viscosity = 0.0;
            ZeroEx(ie->VOF,NMT);
            ZeroEx(ie->sPsiL3,NPSI);
            ZeroEx(ie->PsiL3,NPSI);
            ie->VibPressure = *VibVelocity = 0.0;
            for(int k=0;k<NMT;k++)
            {
                ZeroEx(ie->PhiL3[k],NPHI);
                ZeroEx(ie->sPhiL3[k],NPHI);
                ie->PhiL3[k][0] = 1.0;
            }
            if(Distance(Center,ie->Center) <= Radiu)
            {
                ie->Temperature = iTem;
                ie->Pressure = iPre;
                ie->Density  = iDen;

                ie->VOF[IdImpact] = 1.0;
                ie->PhiL3[IdImpact][1] = iDen;
                ie->PhiL3[IdImpact][2] = iEng;

                ScalerMoveEx(ie->sPhiL3[IdImpact],ie->PhiL3[IdImpact],1.0,NPHI);
                ie->Damage = ie->PsiL3[0] = iDamage;
                ie->Mass = ie->Density * ie->Volume;
                ie->sPsiL3[0] = ie->Mass*iDamage;
                ScalerMove(ie->Momentum,ImpactVelocity,ie->Mass);
            } else if(topDistance < toplevel)
            {
                int dpz = (int) (Wind((toplevel - topDistance)/dz,0,ProfileLen-1));
                ie->Temperature = tTemp[dpz];
                ie->Pressure    = tPres[dpz];
                ie->Density     = tDens[dpz];

                ie->VOF[profId[dpz]] = 1.0;
                ie->PhiL3[profId[dpz]][1] = tDens[dpz];
                ie->PhiL3[profId[dpz]][2] = tEng[dpz];

                ScalerMoveEx(ie->sPhiL3[profId[dpz]],ie->PhiL3[profId[dpz]],1.0,NPHI);
                ie->Damage = ie->PsiL3[0] = tDam[dpz];
                ie->Mass = ie->Density * ie->Volume;
                ie->sPsiL3[0] = ie->Mass*ie->Damage;
                Zero(ie->Momentum);
            } else
            {
                ie->Temperature = 0.0;
                ie->Pressure    = 0.0;
                ie->Density     = 0.0;
                ie->VOF[0] = 1.0;
                ie->sPhiL3[0][0] = 1.0;
                ie->Mass = 0.0;
                Zero(ie->Momentum);
            }
        }

    } else
    {
        fprintf(stdout,"unknown slopetype [%s]\n",TopTypeOpt);
        exit(0);
    }


    // initialization of vertex
    for(int k=0;k<_mesh->nv;k++)
    {
        struct Vertex * iv = _mesh->Vertexes + k;
        Zero(iv->BodyF);
        iv->BodyF[Z] = gravity;
    }

    for(int k=0;k<_mesh->ndb;k++)
    {
        struct Boundary * ib = _mesh->DomainB[k];
        ib->dMass = 0.0;
        Zero(ib->dMomentum);
        ZeroEx(ib->dPsi,NPSI);
        for(int j=0;j<NMT;j++)
        {
            ZeroEx(ib->dPhi[j],NPHI);
        }
    }

    free(tGrav);
    free(tTemp);
    free(tPres);
    free(tDens);
    free(tEng);
    free(tDam);
    free(profId);
}

void InitChk(InputFile * ifp, struct Mesh * _mesh, struct ProcessInfo * _pinfo, struct MeshInfo * _minfo)
{
    char ChkLoadOpt[MaxStrLen];
    GetValueS(ifp,"cycle.chkload",ChkLoadOpt,"off");
    if(strcasecmp(ChkLoadOpt,"on")==0)
    {
        char ChkLoadPrefix[MaxStrLen];
        GetValueS(ifp,"cycle.chklpre", ChkLoadPrefix,"null");
        ChkLoad(_mesh,_pinfo->rank,ChkLoadPrefix);
        if(0 == _pinfo->rank)
        {
            fprintf(stdout,"ChkLoadOpt is ON, load checkpoint from %s*\n",ChkLoadPrefix);
        }
    }
    else
    {
        if(_pinfo->rank==0)
        {
            fprintf(stdout,"ChkLoadOpt is OFF.\n");
        }
    }

    char ChkStoreOpt[MaxStrLen];
    GetValueS(ifp,"cycle.chkstore",ChkStoreOpt,"off");
    if(strcasecmp(ChkStoreOpt,"on")==0)
    {
        _minfo->chkStep = GetValueI(ifp,"cycle.chkstep","5000");
        GetValueS(ifp,"cycle.chkspre",_minfo->chkPrefix,"tmpchk");
        if(0 == _pinfo->rank)
        {
            fprintf(stdout,"ChkStoreOpt is ON, check point will be written to %s* every %d cycle.\n",_minfo->chkPrefix,_minfo->chkStep);
        }
    }
    else
    {
        _minfo->chkStep = -1;
        if(_pinfo->rank==0)
        {
            fprintf(stdout,"ChkStoreOpt is OFF\n");
        }
    }
}

void UpdateChk(struct Mesh * _mesh, struct MeshInfo * _minfo, struct ProcessInfo * _pinfo, int _step)
{
    char tmp[200];
    sprintf(tmp,"%s.%04d",_minfo->chkPrefix,_step);
    ChkWrite(_mesh,_pinfo->rank,tmp);
}

void CompareChk(struct Mesh * _mesh, struct MeshInfo * _minfo, struct ProcessInfo * _pinfo, int _step)
{
    // compare data in checkpoint with data in mesh
    struct Mesh test_mesh;
    test_mesh.Elements = (struct Element *) malloc(_mesh->ne * sizeof(struct Element));
    test_mesh.ne = _mesh->ne;

    char tmp_prefix[200];
    sprintf(tmp_prefix,"%s.%04d",_minfo->chkPrefix,_step);

    ChkLoad(&test_mesh,_pinfo->rank,tmp_prefix);

    double max_diff = 0.,avg_diff = 0.;
#define DIFFPERCENT(a,b) fabs((a)-(b))/(fabs((b)) + 1.0e-6)
    for(int k=0;k<_mesh->ne;++k)
    {
        struct Element * ie = _mesh->Elements + k;
        struct Element * te = test_mesh.Elements + k;

        double sum_diff = 0.0;
        for(int j=0;j<NPSI;++j)
        {
            sum_diff += DIFFPERCENT(te->PsiL3[j],ie->PsiL3[j]);
        }
        for(int j=0;j<NMT;++j)
        {
            sum_diff += DIFFPERCENT(te->VOF[j],ie->VOF[j]);
        }
        avg_diff += sum_diff;
        max_diff = Max(max_diff,sum_diff);
    }
#undef DIFFPERCENT

    fprintf(stdout,"the max relative difference is %%%f, average difference is %%%f\n",max_diff*100.,avg_diff*100.);
    // clean test_mesh
    if(NULL != test_mesh.Elements) free(test_mesh.Elements);
}


double Sk(int k,double ext)
{
    return ext*(pow(ext,k) - 1.0)/(ext - 1.0);
}

double Spacing(int k,double ext, int eL, int eR, int N)
{
    double Xk = 0.0;
    if(k<eL)
        Xk = Sk(eL,ext) - Sk(eL-k,ext);
    else if(k<N-eR)
        Xk = Sk(eL,ext) + (k - eL);
    else
        Xk = Sk(eL,ext) + (N - eL - eR) + Sk(k-N+eR,ext);
    return Xk;
}

double LoadTillEOS(struct TillotsonTable * _t, FILE * fp)
{
    InputFile * ifp = ParseInputFile(fp);
    _t->TLRho0 = GetValueD(ifp,"Tillostson.Rho0","2.7e3");
    _t->TLCv   = GetValueD(ifp,"ColdEnergy.Cv","8.96e2");
    _t->TLA    = GetValueD(ifp,"Tillostson.BulkA","7.52e10");
    _t->TLB    = GetValueD(ifp,"Tillostson.BulkB","6.50e10");
    _t->TLE0   = GetValueD(ifp,"Tillostson.E0","5.0e6");
    _t->TLa    = GetValueD(ifp,"Tillostson.Ta","0.50");
    _t->TLb    = GetValueD(ifp,"Tillostson.Tb","1.63");
    _t->TLAlpha = GetValueD(ifp,"Tillostson.Alpha","5.0");
    _t->TLBeta = GetValueD(ifp,"Tillostson.Beta","5.0");
    _t->TLEiv = GetValueD(ifp,"Tillostson.Eiv","3.00e6");
    _t->TLEcv = GetValueD(ifp,"Tillostson.Ecv","1.39e7");
    _t->TLTref = GetValueD(ifp,"ColdEnergy.Tref","293.0");
    _t->nDen = GetValueI(ifp,"ColdEnergy.RhoStep","1501");
    CloseInputFile(ifp);
    TillColdEnergy(_t);
}

void EosInit(struct MeshInfo * _minfo, struct EosTable * _etb, struct StateReference * _srf)
{
    _etb[0]._atb = NULL;
    _etb[0]._ttb = NULL;
    _etb[0]._gtb = NULL;
    _etb[0].CalL3t = NULL;
    _etb[0].InterpolateTP = NULL;

    for(int k=1;k<_minfo->MaterialNum;k++)
    {
        struct EosTable * ketb = _etb + k;
        struct StateReference * ksrf = _srf + k;
        if(0==strcasecmp("aneos",_minfo->EosType[k]))
        {
            ketb->_atb = (struct ANEOSTable *) malloc(sizeof(struct ANEOSTable));
            ketb->_ttb = NULL;
            ketb->_gtb = NULL;
            LoadANEOS(ketb->_atb,_minfo->Materialfp[k]);
            ANEOSInitStateRef(ketb->_atb,ksrf);
            ketb->CalL3t = ANEOSCalL3tCompose;
            ketb->InterpolateTP = ANEOSInterpolateTPCompose;
        } else if(0== strcasecmp("tillotson",_minfo->EosType[k]))
        {
            ketb->_atb = NULL;
            ketb->_ttb = (struct TillotsonTable *) malloc(sizeof(struct TillotsonTable));
            ketb->_gtb = NULL;
            LoadTillEOS(ketb->_ttb,_minfo->Materialfp[k]);
            TillInitStateRef(ketb->_ttb,ksrf);
            ketb->CalL3t = TillEOSCalL3tCompose;
            ketb->InterpolateTP = TillEOSInterpolateTPCompose;
        } else
        {
            fprintf(stdout,"unknown type of Eos\n");
            exit(0);
        }
    }
}

void LoadStateRef(InputFile * ifp, struct MeshInfo * _minfo,struct StateReference * _srf)
{
    for(int k=1;k<_minfo->MaterialNum;k++)
    {
        SRFInit(ifp,_srf+k,k);
    }
}