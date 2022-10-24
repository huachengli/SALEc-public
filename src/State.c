//
// Created by huacheng on 2021/1/22.
//

#include <malloc.h>
#include "State.h"


void AllocateANEOS(struct ANEOSTable *_t)
{
    _t->yDen = (double *) malloc(_t->nDen * sizeof(double));
    _t->xTem = (double *) malloc(_t->nTem * sizeof(double));
    _t->Data = (double ***) malloc(_t->nTem * sizeof(double **));

    for (int k = 0; k < _t->nTem; k++)
    {
        _t->Data[k] = (double **) malloc(_t->nDen * sizeof(double *));
        for (int i = 0; i < _t->nDen; i++)
        {
            _t->Data[k][i] = (double *) malloc(3 * sizeof(double));
        }
    }
}

void UnAllocateANEOS(struct ANEOSTable *_t)
{

    free(_t->xTem);
    free(_t->yDen);
    for (int k = 0; k < 3; k++) // need to be fixed
    {
        for (int i = 0; i < _t->nTem; i++)
        {
            free(_t->Data[k][i]);
        }
        free(_t->Data[k]);
    }
    free(_t->Data);
}

void LoadANEOS(struct ANEOSTable *_t, FILE *fp)
{
    // the comment line should be processed before this function
    // process them here is not a right solution
    Over(fp,3);
    if(2!=fscanf(fp, "%d %d", &_t->nTem, &_t->nDen))
    {
        fprintf(stdout,"Can not get number of indexes for density or temperature!\n");
        exit(0);
    }
    AllocateANEOS(_t);

    if(3!=fscanf(fp, "%le %le %le", &_t->Tnorm, &_t->Dnorm, &_t->Pnorm))
    {
        fprintf(stdout, "Can not get the standary state from EOS file!\n");
        exit(0);
    }


    for (int i = 0; i < _t->nTem; i++)
    {
        if(1!=fscanf(fp, "%le", &_t->xTem[i]))
        {
            fprintf(stdout,"Can not get table value!!\n");
            exit(0);
        }
    }

    for (int i = 0; i < _t->nTem; i++)
    {
        for (int j = 0; j < _t->nDen; j++)
        {
            if(3!=fscanf(fp, "%le %le %le\n", &_t->Data[i][j][0],&_t->Data[i][j][1],&_t->Data[i][j][2]))
            {
                fprintf(stdout, "Can not get table value from EOS table!\n");
                exit(0);
            }
        }
    }

    for (int i = 0; i < _t->nDen; i++)
    {
        if(1!=fscanf(fp, "%le", &_t->yDen[i]))
            fprintf(stdout, "Can not get index values of density!\n");
    }

    // in the raw table , the first column is energy
    // but in calculation energy*density is used more usually
    // this will make a difference to the linear interpolation of energy in
    // (temperature,density). (a linear interpolation of energy*density)
    for (int i = 0; i < _t->nTem; i++)
    {
        for (int j = 0; j < _t->nDen; j++)
        {
            _t->Data[i][j][0] *= _t->yDen[j];
        }
    }

    fclose(fp);
}

void ANEOSCalL3t(struct Element *_e, const struct ANEOSTable *_t, int _k)
{
    // update the ThetaL3 in element for _k-th material
    // ThetaL3 is (Temperature,Pressure,Csound)
    int DI = 0;
    int DJ = _t->nDen - 1;
    while (DI + 1 < DJ)
    {
        int DM = (DI + DJ) / 2;
        if (_t->yDen[DM] > _e->PhiL3[_k][1])
            DJ = DM;
        else
            DI = DM;
    }

    // the local coordinate for Density
    double Dxr = (_e->PhiL3[_k][1] - _t->yDen[DI]) / (_t->yDen[DJ] - _t->yDen[DI]);
    Dxr = Wind(Dxr,0.0,1.0);
    double Dxl = 1.0 - Dxr;

//    double Energy = _e->PsiL3[7] * _e->PhiL3[_k][1];
    double Energy = _e->PhiL3[_k][2];
    int TI = 0;
    double IRE = _t->Data[TI][DI][0] * Dxl + _t->Data[TI][DJ][0] * Dxr;
    int TJ = _t->nTem - 1;
    double JRE = _t->Data[TJ][DI][0] * Dxl + _t->Data[TJ][DJ][0] * Dxr;
    while (TI + 1 < TJ)
    {
        int TM = (TI + TJ) / 2;
        double MRE = _t->Data[TM][DI][0] * Dxl + _t->Data[TM][DJ][0] * Dxr;
        if (MRE > Energy)
        {
            TJ = TM;
            JRE = MRE;
        }
        else
        {
            TI = TM;
            IRE = MRE;
        }
    }

    // the local coordinate for temperature
    double Txr = (Energy - IRE) / (JRE - IRE);
    Txr = Wind(Txr,0.0,1.0);
    double Txl = 1.0 - Txr;
    _e->ThetaL3[_k][0] = _t->xTem[TI] * Txl + _t->xTem[TJ] * Txr;

    // pressure and speed of sound
    for (int m = 1; m < 3; m++)
    {
        _e->ThetaL3[_k][m] = (_t->Data[TI][DI][m] * Dxl + _t->Data[TI][DJ][m] * Dxr) * Txl
                             + (_t->Data[TJ][DI][m] * Dxl + _t->Data[TJ][DJ][m] * Dxr) * Txr;
    }

#if __RUNPOINT__
    if(_e->ThetaL3[_k][2] < 0)
    {
        fprintf(stdout,"pre=%f,L=%f,R=%f\n",_e->ThetaL3[_k][1],_t->Data[TI][DI][1],_t->Data[TJ][DJ][1]);
        fprintf(stdout,"%d,(%d,%d) of data\n",_k,TI,DI);
        fprintf(stdout,"%d,(%1.10e,%1.10e) of data\n",_k,_t->xTem[TI],_t->yDen[DI]);
        for(int i=0;i<_t->xTem;i++)
        {
            fprintf(stdout,"%1.10e,",_t->Data[i][DI][1]);
            if(0==i%6)
                fprintf(stdout,"\n");
        }
        exit(12);
    }
#endif
    /*
     * This function can be rewritten using Ip2d*
     * but it's not a nice choice when performance is considered.
     * A potential improved method/implement is to interpolate this result on a
     * regular grid. This implement need rearrange the ANEOS data table.
     */
}

void ANEOSInitStateRef(struct ANEOSTable *_t, struct StateReference *_s)
{
    /*
     * change the Dnorm accord to EOS
     * make (Tnorm,Dnorm,Pnorm) is consistent with ANEOSTable
     * this segment is copied from ANEOSCalL3t
     */
    int TI = 0;
    int TJ = _t->nTem;
    while (TI + 1 < TJ)
    {
        int TM = (TI + TJ) / 2;
        if (_t->xTem[TM] > _t->Tnorm)
            TJ = TM;
        else
            TI = TM;
    }

    double Txr = (_t->Tnorm - _t->xTem[TI]) / (_t->xTem[TJ] - _t->xTem[TI]);
    double Txl = 1.0 - Txr;

    int DI = 0;
    int DJ = _t->nDen - 1;
    double IRP = _t->Data[TI][DI][1] * Txl + _t->Data[TJ][DI][1] * Txr;
    double JRP = _t->Data[TI][DJ][1] * Txl + _t->Data[TJ][DJ][1] * Txr;
    while (DI + 1 < DJ)
    {
        int DM = (DI + DJ) / 2;
        double MRP = _t->Data[TI][DM][1] * Txl + _t->Data[TJ][DM][1] * Txr;
        if (MRP > _t->Pnorm)
        {
            DJ = DM;
            JRP = MRP;
        }
        else
        {
            DI = DM;
            IRP = MRP;
        }
    }

    double Dxr = (_t->Pnorm - IRP) / (JRP - IRP);
    double Dxl = 1.0 - Dxr;
    _t->Dnorm = _t->yDen[DI] * Dxl + _t->yDen[DJ] * Dxr;

    // Get the density of vapor and melt accord to the Dnorm
    _s->MeltDen =  0.85 * _t->Dnorm;
    _s->VaporDen = 0.04 * _t->Dnorm;
    _s->Viscosity = 0.0;
}

int GetState(struct Element *_e, const struct StateReference _s[], int nm)
{
    /*
     * (Elastic    = [state] 0) -- SOLID
     * (Invicosity = [state] 1) -- MELT
     * (Gas        = [state] 3) -- VAPOR
     * (hasVaccum  = [state] 2) -- VOID
     */

    if(_e->VOF[VACUUM] > TOLVOF)
    {
        /*
         * a void element;
         * stress = 0.0
         */
        _e->State = VoidElement;
    }
    else
    {
        double AverageMeltDen  = 0.0;
        double AverageVaporDen = 0.0;
        _e->Viscosity = 0.0;
        for (int i = 1; i < nm; i++)
        {
            if (_e->VOF[i] > TOLVOF)
            {
                AverageMeltDen  += _e->VOF[i] * _s[i].MeltDen;
                AverageVaporDen += _e->VOF[i] * _s[i].VaporDen;
                _e->Viscosity   += _e->VOF[i] * _s[i].Viscosity;
            }
        }

        if (_e->Density > AverageMeltDen)
        {
            _e->State = SolidElement;
        }
        else if (_e->Density > AverageVaporDen)
        {
            _e->State = FluidElement;
        }
        else
        {
            _e->State = VaporElement;
        }
    }
    return _e->State;
}

double SimpleRockStrength(struct Element *_e, const struct StateReference *_s)
{
    /*
     * this is a example for the strength function
     * This code encourage user to write the Strength function
     * accord to different materials.
     * The parameters required by this function for different materials can be added
     * in the struct StateReference.
     */

    // Intact strength uses the lundborg approximation
    double Yint = Ylundborg(_e->Pressure,_s->Yint0,_s->Yfricint,_s->Ylimint);
    double Ydam = Ydrucker(_e->Pressure,_s->Ydam0,_s->Yfricdam,_s->Ylimdam);
    double Yac = BlockVibStrength(_e,_s);
    Ydam = Min(Yint,Ydam);
    return Min(Ydam * _e->Damage + Yint * (1.0 - _e->Damage),Yac);
}

double SimpleRockNoacfl(struct Element *_e, const struct StateReference *_s)
{
    /*
    * a special case for SimpleRockStrength with (GammaBeta,GammaEta=0)
     */

    // Intact strength uses the lundborg approximation
    double Yint = Ylundborg(_e->Pressure,_s->Yint0,_s->Yfricint,_s->Ylimint);
    double Ydam = Ydrucker(_e->Pressure,_s->Ydam0,_s->Yfricdam,_s->Ylimdam);

    Ydam = Min(Yint,Ydam);
    return Ydam * _e->Damage + Yint * (1.0 - _e->Damage);
}

double JohnsonCook2(struct Element *_e, const struct StateReference *_s)
{
    double PlasticStrain    = Wind(_e->PsiL3[7],0.0,_s->JCMaxPlastic);
    double JcY1 = _s->JCA + _s->JCB * pow(PlasticStrain,_s->JCN);
    double DotPlasticStrain = _e->sInvSta2;
    double JcY2 = 1.0;
    if(DotPlasticStrain > 1.0 && _s->JCC > 0.0)
    {
        JcY2 = 1.0 + _s->JCC * log(DotPlasticStrain);
    }

    return JcY1*JcY2*JNCKSoft(_e,_s)*BETA;
}

double JohnsonCook1(struct Element *_e, const struct StateReference *_s)
{
    double PlasticStrain    = Wind(_e->PsiL3[7],0.0,_s->JCMaxPlastic);
    double JcY1 = _s->JCA + _s->JCB * pow(PlasticStrain,_s->JCN);
    double DotPlasticStrain = _e->sInvSta2;
    double JcY2 = 1.0;
    if(DotPlasticStrain > 1.0 && _s->JCC > 0.0)
    {
        JcY2 = 1.0 + _s->JCC * log(DotPlasticStrain);
    }
    return JcY1*JcY2*BETA;
}


double SimpleTensile(struct Element *_e, const struct StateReference *_s)
{
    return _s->Yten0 * (1.0 - _e->Damage);
}

double SimpleSheardam(struct Element *_e, const struct StateReference *_s)
{
    /*
     * this function return the increment of damage caused by shear failure
     * This function return a reverse of the counterpart of iSALE2d
     * Example for IVANOV damage model
     */
    return 1.0/Max(_s->IvanA, _s->IvanB * (_e->Pressure - _s->IvanC));
}

double SimpleSoft(struct Element *_e, const struct StateReference *_s)
{
    double Tmelt = SimonMelt(_e,_s);
    return OnnakaSoft(_e,_s,Tmelt);
}

double JNCKSoft(struct Element *_e, const struct StateReference *_s)
{
    double Tmelt = SimonMelt(_e,_s);
    double MeltFrac = (_e->Temperature - _s->JCTREF)/(Tmelt - _s->JCTREF);
    MeltFrac = Wind(MeltFrac,0.0,1.0);
    return Wind(1.0 - pow(MeltFrac,_s->JCM),0.0,1.0);
}

double SimonMelt(struct Element *_e, const struct StateReference *_s)
{
    return _s->Tmelt0 * pow(_s->Asimon*_e->Pressure + 1.0,_s->Csimon);
}

double OnnakaSoft(struct Element *_e, const struct StateReference *_s, double Tmelt)
{
    return tanh(Tmelt/_e->Temperature - 1.0)*_s->Tfrac;
}

void Over(FILE * fp,int n)
{
    /*
     * To process the head of ANEOS table
     * The first version is only able to process some comment line
     * in the head.
     */
    char buf[1024];
    for(int k=0;k<n;k++)
    {
        if(NULL==fgets(buf,sizeof(buf),fp))
            fprintf(stdout,"error in fgets(). skip over some lines is not valid.\n");
    }
}

double InterpolateTD(struct ANEOSTable *_t,double tTem, double tDen, int DataId)
{
    /*
     * Interpolate the data using (Tem,Den)
     */

    int TI = 0;
    int TJ = _t->nTem - 1;
    while (TI + 1 < TJ)
    {
        int TM = (TI + TJ) / 2;
        if (_t->xTem[TM] > tTem)
            TJ = TM;
        else
            TI = TM;
    }

    double Txr = (tTem - _t->xTem[TI]) / (_t->xTem[TJ] - _t->xTem[TI]);
    Txr = Wind(Txr,0.0,1.0);
    double Txl = 1.0 - Txr;

    int DI = 0;
    int DJ = _t->nDen - 1;
    while (DI + 1 < DJ)
    {
        int DM = (DI + DJ) / 2;
        if (_t->yDen[DM] > tDen)
            DJ = DM;
        else
            DI = DM;
    }

    double Dxr = (tDen - _t->yDen[DI])/(_t->yDen[DJ] - _t->yDen[DI]);
    Dxr = Wind(Dxr,0.0,1.0);
    double Dxl = 1.0 - Dxr;

    return (_t->Data[TI][DI][DataId] * Dxl + _t->Data[TI][DJ][DataId] * Dxr) * Txl
           + (_t->Data[TJ][DI][DataId] * Dxl + _t->Data[TJ][DJ][DataId] * Dxr) * Txr;
}

double ANEOSInterpolateTP(struct ANEOSTable *_t, double tTem, double tPre, int DataId)
{
    /*
     * get data using (Ten,Pre)
     * copy from ANEOSInitStateRef
     */

    int TI = 0;
    int TJ = _t->nTem;
    while (TI + 1 < TJ)
    {
        int TM = (TI + TJ) / 2;
        if (_t->xTem[TM] > tTem)
            TJ = TM;
        else
            TI = TM;
    }

    double Txr = (tTem - _t->xTem[TI]) / (_t->xTem[TJ] - _t->xTem[TI]);
    double Txl = 1.0 - Txr;

    int DI = 0;
    int DJ = _t->nDen - 1;
    double IRP = _t->Data[TI][DI][1] * Txl + _t->Data[TJ][DI][1] * Txr;
    double JRP = _t->Data[TI][DJ][1] * Txl + _t->Data[TJ][DJ][1] * Txr;
    while (DI + 1 < DJ)
    {
        int DM = (DI + DJ) / 2;
        double MRP = _t->Data[TI][DM][1] * Txl + _t->Data[TJ][DM][1] * Txr;
        if (MRP > tPre)
        {
            DJ = DM;
            JRP = MRP;
        }
        else
        {
            DI = DM;
            IRP = MRP;
        }
    }

    double Dxr = (tPre - IRP) / (JRP - IRP);
    Dxr = Wind(Dxr,0.0,1.0);
    double Dxl = 1.0 - Dxr;
    if(DataId == -1)
        return _t->yDen[DI] * Dxl + _t->yDen[DJ] * Dxr;
    else
        return (_t->Data[TI][DI][DataId] * Dxl + _t->Data[TI][DJ][DataId] * Dxr) * Txl
           + (_t->Data[TJ][DI][DataId] * Dxl + _t->Data[TJ][DJ][DataId] * Dxr) * Txr;

}

double BilinearInverseM[4][5] = {
        { 1.0, 1.0, 1.0, 1.0,-4.0},
        { 1.0,-1.0,-1.0, 1.0, 0.0},
        { 1.0, 1.0,-1.0,-1.0, 0.0},
        { 1.0,-1.0, 1.0,-1.0, 0.0}};
int BilinearInverse(const double xi[],const double yi[], double xl[])
{
    double a[4] = {0.0, 0.0, 0.0, 0.0};
    double b[4] = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            a[i] += 0.25*BilinearInverseM[i][j]*xi[j];
            b[i] += 0.25*BilinearInverseM[i][j]*yi[j];
        }
    }

    double tA = a[1]*b[3] - a[3]*b[1];
    double tB = a[1]*b[2] - a[2]*b[1] + a[0]*b[3] - a[3]*b[0];
    double tC = a[0]*b[2] - a[2]*b[0];

    double Delta = tB*tB - 4.0*tA*tC;
    if(Delta < 0.0)
        return 0;
    else
    {
        xl[0] = (-tB + sqrt(Delta))/(2.0*tA);
        if((xl[0] > 1.0 + 1e-5)|| (xl[0] < -1.0 -1.0e-5))
            xl[0] = (-tB - sqrt(Delta))/(2.0*tA);
        xl[1] = -1.0*(a[0] + a[1]*xl[0])/(a[2] + a[3]*xl[0]);
    }

    if((xl[0] > 1.0 + 1e-5)|| (xl[0] < -1.0-1.0e-5) || (xl[1] > 1.0 + 1e-5) || (xl[1] < -1.0-1.0e-5))
        return 0;
    else
        return 1;
}


void AllocateAirEOS(struct AirEOSTable * _air)
{
    _air->xEng = (double *) malloc(sizeof(double) * _air->nEng);
    _air->yDen = (double *) malloc(sizeof(double) * _air->nDen);
    _air->Data = (double ***) malloc(sizeof(double **) * _air->nEng);
    for(int k=0;k<_air->nEng;k++)
    {
        _air->Data[k] = (double **) malloc(sizeof(double *)*_air->nDen);
        for(int i=0;i<_air->nDen;i++)
            _air->Data[k][i] = (double *) malloc(sizeof(double)*3);
    }
}

void UnAllocateAirEOS(struct AirEOSTable * _air)
{
    free(_air->xEng);
    free(_air->yDen);
    for(int k=0;k<_air->nEng;k++)
    {
        for(int i=0;i<_air->nDen;i++)
        {
            free(_air->Data[k][i]);
        }
        free(_air->Data[k]);
    }
    free(_air->Data);
}

void LoadAirEOS(struct AirEOSTable * _t, FILE *fp)
{
    // the comment line is the first 3 line
    Over(fp,3);
    // a line for nTem, nPre
    if(2!=fscanf(fp, "%d %d", &_t->nEng, &_t->nDen))
        fprintf(stdout,"Can not get number of indexes for temperature or pressure!\n");
    AllocateAirEOS(_t);
    // the index for temperature
    for(int i=0;i<_t->nEng;i++)
    {
        if(1!=fscanf(fp, "%le", &_t->xEng[i]))
            fprintf(stdout,"Can not read index for temperature!\n");
    }

    // the index for Pressure
    for(int i=0;i<_t->nDen;i++)
    {
        if(1!=fscanf(fp, "%le", &_t->yDen[i]))
            fprintf(stdout,"Can not read index for pressure!\n");
    }
        ;
    // the data array for (Energy, Density, Speed of sound)
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<_t->nDen;j++)
            for(int k=0;k<_t->nEng;k++)
            {
                if(1!=fscanf(fp, "%le",&_t->Data[k][j][i]))
                    fprintf(stdout,"Can not read data from aireos file!\n");
            }
    }
    fclose(fp);
}

int Cmpfunc(const void * a, const void * b)
{
    const double * tl = (double *) b;
    const double * tr = tl + 1;
    const double * t  = (double *) a;

    if(*t < *tl)
        return -1;
    else if(*t > *tr)
        return 1;
    else
        return 0;
}

int Bsearch(double _key, double * _base,int _len)
{
    int TI = 0;
    int TJ = _len - 1;
    while(TI + 1 < TJ)
    {
        int TM = (TI + TJ)/2;
        if(_base[TM] > _key)
            TJ = TM;
        else
            TI = TM;
    }
    return TI;
}


void CalculateL3atm(struct Element *_e, const struct AirEOSTable *_gas,int _k)
{
    double logRo = log10(_e->PhiL3[_k][1]);
    double logEn = log10(_e->PhiL3[_k][2]) - logRo;

    int TIE = Bsearch(logEn,_gas->xEng,_gas->nEng);
    int TIR = Bsearch(logRo,_gas->yDen,_gas->nDen);

    double Exr = (logEn - _gas->xEng[TIE])/(_gas->xEng[TIE+1] - _gas->xEng[TIE]);
    double Rxr = (logRo - _gas->yDen[TIR])/(_gas->yDen[TIR+1] - _gas->yDen[TIR]);

    if(Exr > 1.0) Exr = 1.0;
    if(Rxr > 1.0) Rxr = 1.0;

    double logPre = ((1.0 - Exr)*_gas->Data[TIE][TIR][0] + Exr*_gas->Data[TIE+1][TIR][0])*(1.0-Rxr)
                    + ((1.0 - Exr)*_gas->Data[TIE][TIR+1][0] + Exr*_gas->Data[TIE+1][TIR+1][0])*Rxr;
    double logTem = ((1.0 - Exr)*_gas->Data[TIE][TIR][1] + Exr*_gas->Data[TIE+1][TIR][1])*(1.0-Rxr)
                    + ((1.0 - Exr)*_gas->Data[TIE][TIR+1][1] + Exr*_gas->Data[TIE+1][TIR+1][1])*Rxr;
    double logCs  = ((1.0 - Exr)*_gas->Data[TIE][TIR][2] + Exr*_gas->Data[TIE+1][TIR][2])*(1.0-Rxr)
                    + ((1.0 - Exr)*_gas->Data[TIE][TIR+1][2] + Exr*_gas->Data[TIE+1][TIR+1][2])*Rxr;
    _e->ThetaL3[_k][0] = pow(10.0,logTem);
    _e->ThetaL3[_k][1] = pow(10.0,logPre);
    _e->ThetaL3[_k][2] = pow(10.0,logCs );

    if(_e->ThetaL3[_k][0] < 500.0)
    {
        _e->ThetaL3[_k][0] = _e->PhiL3[_k][2] / _e->PhiL3[_k][1] / (2.52 * 287.12835);
        _e->ThetaL3[_k][1] = _e->PhiL3[_k][2] / 2.52;
    }
}


double ANEOSRhoGrav(double z, double Pres, struct ANEOSTable *_t, const double Grav[], const double Temp[])
{
    double lGrav = (1-z)*Grav[0] + z*Grav[1];
    double lTemp = (1-z)*Temp[0] + z*Temp[1];
    lGrav = fabs(lGrav);
    double lRho  = ANEOSInterpolateTP(_t, lTemp, Pres, -1);
    return lRho*lGrav;
}


double ANEOSPresProfRK3(struct ANEOSTable *_t, double Pres[], double Grav[], double Temp[], int steps, double dh)
{
    double K1,K2,K3;
    for(int k=1;k<steps;k++)
    {
        K1 = ANEOSRhoGrav(0.0 , Pres[k - 1]             , _t, Grav + k - 1, Temp + k - 1);
        K2 = ANEOSRhoGrav(0.50, Pres[k - 1] + 0.50 * K1 * dh, _t, Grav + k - 1, Temp + k - 1);
        K3 = ANEOSRhoGrav(0.75, Pres[k - 1] + 0.75 * K2 * dh, _t, Grav + k - 1, Temp + k - 1);
        Pres[k] = Pres[k-1] + (2.0*K1 + 3.0*K2 + 4.0*K3)*dh/9.0;
    }
    return Pres[steps-1];
}

double TillEOSRhoGrav(double z, double Pres, struct TillotsonTable *_t, const double Grav[], const double Temp[])
{
    double lGrav = (1-z)*Grav[0] + z*Grav[1];
    double lTemp = (1-z)*Temp[0] + z*Temp[1];
    lGrav = fabs(lGrav);
    double lRho  = TillEOSInterpolateTP(_t, lTemp, Pres, -1);
    return lRho*lGrav;
}

double TillEOSPresProfRK3(struct TillotsonTable *_t, double Pres[], double Grav[], double Temp[], int steps, double dh)
{
    double K1,K2,K3;
    for(int k=1;k<steps;k++)
    {
        K1 = TillEOSRhoGrav(0.0 , Pres[k - 1]             , _t, Grav + k - 1, Temp + k - 1);
        K2 = TillEOSRhoGrav(0.50, Pres[k - 1] + 0.50 * K1 * dh, _t, Grav + k - 1, Temp + k - 1);
        K3 = TillEOSRhoGrav(0.75, Pres[k - 1] + 0.75 * K2 * dh, _t, Grav + k - 1, Temp + k - 1);
        Pres[k] = Pres[k-1] + (2.0*K1 + 3.0*K2 + 4.0*K3)*dh/9.0;
    }
    return Pres[steps-1];
}

double RhoGrav(double z, double Pres, struct EosTable *_t, const double Grav[], const double Temp[])
{
    double lGrav = (1-z)*Grav[0] + z*Grav[1];
    double lTemp = (1-z)*Temp[0] + z*Temp[1];
    lGrav = fabs(lGrav);
    double lRho  = _t->InterpolateTP(_t,lTemp, Pres, -1);
    return lRho*lGrav;
}

double PresProfRK3(struct EosTable *_t, double Pres[], double Grav[], double Temp[], int steps, double dh)
{
    double K1,K2,K3;
    for(int k=1;k<steps;k++)
    {
        K1 = RhoGrav(0.0 , Pres[k - 1]             , _t, Grav + k - 1, Temp + k - 1);
        K2 = RhoGrav(0.50, Pres[k - 1] + 0.50 * K1 * dh, _t, Grav + k - 1, Temp + k - 1);
        K3 = RhoGrav(0.75, Pres[k - 1] + 0.75 * K2 * dh, _t, Grav + k - 1, Temp + k - 1);
        Pres[k] = Pres[k-1] + (2.0*K1 + 3.0*K2 + 4.0*K3)*dh/9.0;
    }
    return Pres[steps-1];
}


/*
double PresProfBS23(struct ANEOSTable *_t,double Pres[], double Grav[], double Temp[],int steps, double dh,double ptol)
{
    return 0.0;
    // uncompleted
    double K1 = Grav[0]*ANEOSInterpolateTP(_t,Temp[0],Pres[0],-1),K2,K3,K4;
    double lPres=Pres[0],rPres,lz=0.0,rz,dh1=dh;
    int k = 1;

    while(1)
    {
        K2 = ANEOSRhoGrav( 0.50*dh1/dh, lPres + 0.50*K1*dh1);
        K3 = ANEOSRhoGrav( 0.75*dh1/dh, lPres + 0.75*K2*dh1);
        rPres = lPres + (2.0*K1 + 3.0*K2 + 4.0*K3)/9.0*dh1;
        K4 = ANEOSRhoGrav( 1.00*dh1/dh, rPres);
        double Err = (-5.0*K1 + 6.0*K2 + 8.0*K3 - 9.0*K4)*dh1/72.0;
        if(Err <= ptol)
        {
            lPres = rPres;
            dh1  = 0.0;
        } else
        {
            dh1 = 0.8*pow(Err/ptol,1.0/3.0)*dh1;
        }
    }
}
*/

int BlockVibration(struct Element * ie, double dt,struct StateReference * _s, int nm,double time)
{
    /*
     * PsiL3[8] is the square of vib-velocity(vibration energy density)
     */
    double * VibVelocity = ie->PsiL3 + 8;

    if(ie->CutTag || (ie->Damage < 0.5))
    {
        ie->VibPressure = 0.0;
        *VibVelocity = 0.0;
        return 0;
    }

    double Tdeacy = 0.0, Cvib=0.0, Toff=0.0, VibMax=0.0, Pvlim=0.0;
    for(int k=1;k<nm;k++)
    {
        if(ie->VOF[k] < TOLVOF) continue;
        Tdeacy = Max(_s[k].Tdamp,Tdeacy);
        Cvib   = Max(_s[k].Cvib , Cvib );
        Toff   = Max(_s[k].Toff , Toff);
        VibMax = Max(_s[k].VibMax, VibMax);
        Pvlim  = Max(_s[k].Pvlim, Pvlim);
    }

    if((Tdeacy <= 0.0) || (Toff <= 0.0))
    {
        ie->VibPressure = 0.0;
        *VibVelocity = 0.0;
        return 0;
    }

    double newVibVelocity = Cvib * Length(ie->Momentum) / ie->Mass;
    double newVibVelocity2 = newVibVelocity*newVibVelocity;

    if( (time <= Toff) && (newVibVelocity2 > *VibVelocity) && (newVibVelocity > 1.0) && (ie->ArtificialPressure > 0.0))
    {
        *VibVelocity = Min(newVibVelocity2, VibMax*VibMax);
    } else
    {
        *VibVelocity = Max(*VibVelocity*(1.0 - 2.0*dt/Tdeacy + 2.0*(dt/Tdeacy)*(dt/Tdeacy)),0.0);
    }

    if( (ie->Pressure < 0.0) || (ie->Pressure > Pvlim) || (*VibVelocity < 1.0))
    {
        ie->VibPressure = 0.0;
    } else
    {
        ie->VibPressure = Min(ie->Pressure,ie->Density * ie->Csound*sqrt(*VibVelocity));
    }
    return 1;
}

double BlockVibStrength(struct Element *_e, const struct StateReference *_s)
{
    double ploc = Max(_e->Pressure - _e->VibPressure*_e->Damage*_e->Damage,0.0);
    double Yint = Ylundborg(ploc,_s->Yint0,_s->Yfricint,_s->Ylimint);
    double Ydam = Ydrucker(ploc,_s->Ydam0,_s->Yfricdam,_s->Ylimdam);
    Ydam = Min(Ydam,Yint);
    double Yvib = Yint*(1.0 - _e->Damage) + Ydam*_e->Damage + _e->Density*_s->Acvis*_e->sInvSta2;
    return Yvib;
}

double Ylundborg(double p,double y0,double fric,double ylim)
{
    /*
     * Smooth lundborg function from y0 to ylim
     */
    return y0 + p*fric/(1 + p*fric/(ylim-y0));
}

double Ydrucker(double p,double y0,double fric,double ylim)
{
    /*
     * Linear with pressure from y0 up to ylim
     */
    return Min(y0 + fric*p,ylim);
}

double LowdensitySoft(double Density, double MeltDensity)
{
    return pow(Density/MeltDensity,4.0);
}

double TillCold(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs)
{
    /*
     * Tillotson equation of state
     * cold compressed states, Eta >= 1 or EngIn <= Eiv
     */
    double TillEta = RhoIn/_t->TLRho0;
    double TillMu  = TillEta - 1.0;
    double TillPhi0 = EngIn/(_t->TLE0*TillEta*TillEta);
    double TillPhi1 = 1.0/(TillPhi0 + 1.0);
    *Pres   = (_t->TLa + _t->TLb*TillPhi1)*EngIn*RhoIn + _t->TLA*TillMu + _t->TLB*TillMu*TillMu;

    double CsSquare1 = (1.0/_t->TLRho0)*(_t->TLA + 2.0*_t->TLB*TillMu) + EngIn*(_t->TLa + _t->TLb*TillPhi1*(3.0 - 2.0*TillPhi1));
    double CsSquare2 = (_t->TLa + _t->TLb*TillPhi1*TillPhi1)*(*Pres)/RhoIn;
    *Cs = Sqrt_trunc(CsSquare1 + CsSquare2);

    return *Cs;
}

double TillHot(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs)
{
    /*
     * Tillotson equation of state
     * expanded states, Eta <= 1 && EngIn >= Ecv
     */
    double TillEta = RhoIn/_t->TLRho0;
    double TillMu  = TillEta - 1.0;
    double TillXi  = _t->TLRho0/RhoIn - 1.0;
    double TillPhi0 = EngIn/(_t->TLE0*TillEta*TillEta);
    double TillPhi1 = 1.0/(TillPhi0 + 1.0);

    double Exp_A2   = 0.0;
    if(_t->TLAlpha*TillXi*TillXi < 100.0)
        Exp_A2 = exp(-1.0*_t->TLAlpha*TillXi*TillXi);
    double Exp_B1  = 0.0;
    if(_t->TLBeta*TillXi < 100.0)
        Exp_B1 = exp(-1.0*_t->TLBeta*TillXi);
    *Pres = _t->TLa*RhoIn*EngIn + (_t->TLb*EngIn*RhoIn*TillPhi1 + _t->TLA*TillMu*Exp_B1)*Exp_A2;

    double CsSquare1 = _t->TLa*EngIn + _t->TLb*EngIn*TillPhi1*(3.0 + 2.0*_t->TLAlpha*TillXi/TillEta - 2.0*TillPhi1)*Exp_A2
            + _t->TLA/_t->TLRho0*(1.0 - (TillXi/TillEta)*(_t->TLBeta+ 2.0*_t->TLAlpha*TillXi))*Exp_A2*Exp_B1;
    double CsSquare2 = (_t->TLa + _t->TLb*TillPhi1*TillPhi1*Exp_A2)*(*Pres)/RhoIn;
    *Cs = Sqrt_trunc(CsSquare1 + CsSquare2);

    return *Cs;
}

double TillTran(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs)
{
    /*
     * Tillotson equation of state
     * transition between Cold and Hot
     */

    double TillETran = _t->TLEcv - _t->TLEiv;
    double TilldEr   = EngIn - _t->TLEiv;
    double TilldEl   = _t->TLEcv - EngIn;

    double TillEta = RhoIn/_t->TLRho0;
    double TillMu  = TillEta - 1.0;
    double TillPhi0 = EngIn/(_t->TLE0*TillEta*TillEta);
    double TillPhi1 = 1.0/(TillPhi0 + 1.0);
    double TillXi  = _t->TLRho0/RhoIn - 1.0;

    double PresCold = (_t->TLa + _t->TLb*TillPhi1)*EngIn*RhoIn + _t->TLA*TillMu + _t->TLB*TillMu*TillMu;

    double Exp_A2   = 0.0;
    if(_t->TLAlpha*TillXi*TillXi < 100.0)
    {
        Exp_A2 = exp(-1.0*_t->TLAlpha*TillXi*TillXi);
    }
    double Exp_B1  = 0.0;
    if(_t->TLBeta*TillXi < 100.0)
    {
        Exp_B1 = exp(-1.0*_t->TLBeta*TillXi);
    }

    double PresHot = _t->TLa*RhoIn*EngIn + (_t->TLb*EngIn*RhoIn*TillPhi1 + _t->TLA*TillMu*Exp_B1)*Exp_A2;

    *Pres = (TilldEl*PresCold + TilldEr*PresHot)/TillETran;

    double CsColdSquare1 = (1.0/_t->TLRho0)*(_t->TLA + 2.0*_t->TLB*TillMu) + EngIn*(_t->TLa + _t->TLb*TillPhi1*(3.0 - 2.0*TillPhi1));
    double CsColdSquare2 = (_t->TLa + _t->TLb*TillPhi1*TillPhi1)*(*Pres)/RhoIn;

    double CsHotSquare1 = _t->TLa*EngIn + _t->TLb*EngIn*TillPhi1*(3.0 + 2.0*_t->TLAlpha*TillXi/TillEta - 2.0*TillPhi1)*Exp_A2
            + _t->TLA/_t->TLRho0*(1.0 - (TillXi/TillEta)*(_t->TLBeta+ 2.0*_t->TLAlpha*TillXi))*Exp_A2*Exp_B1;
    double CsHotSquare2 = (_t->TLa + _t->TLb*TillPhi1*TillPhi1*Exp_A2)*(*Pres)/RhoIn;

    double Cs2Cold = CsColdSquare1 + CsColdSquare2;
    double Cs2Hot  = CsHotSquare1  + CsHotSquare2 ;

    double CsSquare1 = (TilldEl*Cs2Cold + TilldEr*Cs2Hot)/TillETran;
    double CsSquare2 = (PresHot - PresCold)/TillETran * (*Pres)/(RhoIn * RhoIn);
    *Cs = Sqrt_trunc(CsSquare1 + CsSquare2);

    return *Cs;
}


double TillPres(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs)
{

    if((RhoIn >= _t->TLRho0) || (EngIn < _t->TLEiv))
        TillCold(_t,EngIn,RhoIn,Pres,Cs);
    else if(EngIn > _t->TLEcv)
        TillHot(_t,EngIn,RhoIn,Pres,Cs);
    else
        TillTran(_t,EngIn,RhoIn,Pres,Cs);

    if(*Pres <= 1.0e-2)
        *Cs = Sqrt_trunc(_t->TLA/_t->TLRho0);
    return *Cs;
}

double TilldEdr(const struct TillotsonTable * _t, double EngIn, double RhoIn)
{
    double tCs,tPres;
    TillPres(_t,EngIn,RhoIn,&tPres,&tCs);
    return tPres/(RhoIn*RhoIn);
}

void TillColdEnergy(struct TillotsonTable * _t)
{
    double RefDen = _t->TLRho0;
    double MaxDen = _t->TLRho0*15.0;
    int nDen = _t->nDen;
    double StepDen = MaxDen/(nDen-1);
    double * xDen = (double *) malloc(sizeof(double)*nDen);
    double * yEng = (double *) malloc(sizeof(double)*nDen);
    for(int k=0;k<nDen;k++)
    {
        xDen[k] = k*StepDen;
        yEng[k] = 0.0;
    }

    int IndexM = (int)(floor(RefDen/MaxDen*nDen)+1);
    double K1,K2,K3;
    for(int k=IndexM+1;k<nDen;k++)
    {
        K1 = TilldEdr(_t,yEng[k-1]                  ,xDen[k-1]);
        K2 = TilldEdr(_t,yEng[k-1] + 0.50*K1*StepDen,xDen[k-1]+0.50*StepDen);
        K3 = TilldEdr(_t,yEng[k-1] + 0.75*K2*StepDen,xDen[k-1]+0.75*StepDen);
        yEng[k] = yEng[k-1] + (2.0*K1 + 3.0*K2 + 4.0*K3)*StepDen/9.0;
    }

    StepDen *= -1.0;
    for(int k=IndexM-1;k>=1;k--)
    {
        K1 = TilldEdr(_t,yEng[k+1],xDen[k+1]);
        K2 = TilldEdr(_t,yEng[k+1] + 0.50*K1*StepDen,xDen[k+1]+0.50*StepDen);
        K3 = TilldEdr(_t,yEng[k+1] + 0.75*K2*StepDen,xDen[k+1]+0.75*StepDen);
        yEng[k] = yEng[k+1] + (2.0*K1 + 3.0*K2 + 4.0*K3)*StepDen/9.0;
    }

    yEng[0] = yEng[1];

    _t->xDen = xDen;
    _t->yEng = yEng;
    _t->StepDen = MaxDen/(nDen-1);
}

double TillTemp(const struct TillotsonTable * _t,double EngIn, double RhoIn)
{
    int DL = (int)(Min(floor(RhoIn/_t->StepDen),_t->nDen-2));
    double Dxl = DL + 1.0 - RhoIn/_t->StepDen;
    double Dxr = 1.0 - Dxl;
    double rEngIn = Dxl*_t->yEng[DL] + Dxr*_t->yEng[DL+1];
    return Max((EngIn - rEngIn)/_t->TLCv,0.0) + _t->TLTref;
}


void TillEOSCalL3t(struct Element *_e, const struct TillotsonTable * _t, int _k)
{
    double EngIn = _e->PhiL3[_k][2]/_e->PhiL3[_k][1];
    double RhoIn = _e->PhiL3[_k][1];
    _e->ThetaL3[_k][0] = TillTemp(_t,EngIn,RhoIn);
    TillPres(_t,EngIn,RhoIn,&(_e->ThetaL3[_k][1]),&(_e->ThetaL3[_k][2]));
}

void TillEOSCalL3tCompose(struct Element *_e, const struct EosTable * _t, int _k)
{
    TillEOSCalL3t(_e,_t->_ttb,_k);
}

void ANEOSCalL3tCompose(struct Element *_e, const struct EosTable * _t, int _k)
{
    ANEOSCalL3t(_e,_t->_atb,_k);
}

void TillInitStateRef(struct TillotsonTable *_t, struct StateReference *_s)
{
    _s->MeltDen   = 0.85*_t->TLRho0;
    double GasMin = 1.0e-4*_t->TLRho0;
    _s->VaporDen = GasMin*10.0;
    _s->Viscosity = 0.0;
}


StrengthFunc SelectStrength(const char s[])
{
    if(strcasecmp(s,"SimpleRock")==0)
    {
        return SimpleRockStrength;
    }
    if(strcasecmp(s,"SimpleRockNoacfl")==0)
    {
        return SimpleRockNoacfl;
    }
    else if(strcasecmp(s,"SimpleTensile")==0)
    {
        return SimpleTensile;
    }
    else if(strcasecmp(s,"JohnsonCook1")==0)
    {
        return JohnsonCook1;
    }
    else if(strcasecmp(s,"JohnsonCook2")==0)
    {
        return JohnsonCook2;
    }
    else
        return NULL;
}

DamageFunc SelectDamage(const char s[])
{
    if(strcasecmp(s,"SimpleShear")==0)
    {
        return SimpleSheardam;
    } else
        return NULL;
}

double TillEOSInterpolateTP(struct TillotsonTable *_t, double tTem, double tPre, int DataId)
{
    double ThermalE = (tTem - _t->TLTref)*_t->TLCv;
    int DI = 0;
    int DJ = _t->nDen - 1;
    while(DI + 1 < DJ)
    {
        int DM = (DI + DJ)/2;
        double DEk = _t->yEng[DM] + ThermalE;
        double DPk,Csk;
        TillPres(_t,DEk,_t->xDen[DM],&DPk,&Csk);
        if(DPk > tPre)
            DJ = DM;
        else
            DI = DM;
    }

    double LIPre,LCs,RJPre,RCs;
    TillPres(_t,_t->yEng[DI]+ThermalE,_t->xDen[DI],&LIPre,&LCs);
    TillPres(_t,_t->yEng[DJ]+ThermalE,_t->xDen[DJ],&RJPre,&RCs);

    double LIDen = _t->xDen[DI], RJDen = _t->xDen[DJ], LMDen, LMEng, LMPre, LMCs;
    for(int k=0;k<30;k++)
    {
        double DenError = (RJDen - LIDen)/RJDen;
        LMDen = (LIDen + RJDen)/2.0;
        if(DenError < 1e-4) break;
        double Dxl = (_t->xDen[DJ] - LMDen)/_t->StepDen;
        double Dxr = 1.0 - Dxl;
        LMEng = Dxl*_t->yEng[DI] + Dxr*_t->yEng[DJ] + ThermalE;
        TillPres(_t,LMEng,LMDen,&LMPre,&LMCs);
        if(LMPre > tPre)
            RJDen = LMDen;
        else
            LIDen = LMDen;
    }

    if(-1==DataId)
        return LMDen;
    else if(0 == DataId)
        return LMEng*LMDen;
    else if(1 == DataId)
        return tPre;
    else if(2 == DataId)
        return LMCs;
    else
    {
        fprintf(stdout,"unsupported data id in InterpolateTP\n");
        exit(0);
    }
}

double TillEOSInterpolateTPCompose(struct EosTable *_t, double tTem, double tPre, int DataId)
{
    return TillEOSInterpolateTP(_t->_ttb,tTem,tPre,DataId);
}

double ANEOSInterpolateTPCompose(struct EosTable *_t, double tTem, double tPre, int DataId)
{
    return ANEOSInterpolateTP(_t->_atb,tTem,tPre,DataId);
}

