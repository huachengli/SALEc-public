//
// Created by huacheng on 11/8/22.
//

#include "Planet.h"

void AllocatePlanetProfile(planet_profile * _prof)
{
    if(_prof->n > 0)
    {
        _prof->pre = (double *) malloc(sizeof(double)*_prof->n);
        _prof->rad = (double *) malloc(sizeof(double)*_prof->n);
        _prof->den = (double *) malloc(sizeof(double)*_prof->n);
        _prof->tem = (double *) malloc(sizeof(double)*_prof->n);
        _prof->gra = (double *) malloc(sizeof(double)*_prof->n);
        _prof->dam = (double *) malloc(sizeof(double)*_prof->n);
        _prof->eng = (double *) malloc(sizeof(double)*_prof->n);
        _prof->mat = (int *) malloc(sizeof(int)*_prof->n);
    }
}

void UnAllocatePlanetProfile(planet_profile * _prof)
{
    free(_prof->pre);
    free(_prof->rad);
    free(_prof->den);
    free(_prof->tem);
    free(_prof->mat);
    free(_prof->gra);
    free(_prof->dam);
    free(_prof->eng);
}

void SetPlanetProfile_Tem(planet_profile * _prof, state_ref * _s, planet_definition * _pdef)
{
    for(int k=0;k<_prof->n;++k)
    {
        _prof->den[k] = (_s + _prof->mat[k])->den0;
        _prof->rad[k] = _pdef->r_planet*k/(_prof->n - 1);
    }

    FILE * file;
    double num1,num2;
    int kl=0;

    switch (_pdef->tem_type) {
        case condconv_profile:
        case convective_profile:
            for(int k=_prof->n-1;k>=0;--k)
            {
                double k_depth = _pdef->r_planet - _prof->rad[k];
                if(_prof->n== (k+1))
                {
                    _prof->tem[k] = _pdef->surf_tem;
                    continue;
                }

                double dtdz_cond = _pdef->surf_dtdz * (1.0 + k_depth/_pdef->r_planet);
                double conv_coef = _s[_prof->mat[k]].thermexp * _pdef->surf_gravity / _s[_prof->mat[k]].c_heat;
                double dtdz_conv = _pdef->t_lith * conv_coef * exp(conv_coef * (k_depth - _pdef->d_lith));
                double conduction = 1.0 - (_prof->tem[k+1] - 0.9*_pdef->t_lith)/(0.11*_pdef->t_lith);
                conduction = Wind(conduction,0.0,1.0);
                _prof->tem[k] = _prof->tem[k+1] + (dtdz_cond*conduction + dtdz_conv*(1.0 - conduction))*(_prof->rad[k+1] - _prof->rad[k]);
            }
            break;
        case user_profile:
            file = fopen("../data.txt","r");
            if (file == NULL)
            {
                printf("error in read the user input temperature profile");
                exit(0);
            }

            while(fscanf(file,"%lf %lf", &num1, &num2) == 2)
            {
                _prof->tem[kl] = num2;
                kl = kl + 1;
            }
            if (kl != _prof->n)
            {
                printf("The user input temperature profile is not correct, length should be %d\n",_prof->n);
                exit(0);
            }
            fclose(file);
            break;
        default:
            fprintf(stdout,"%s:UnImplemented type for tem_type!\n",__func__);
            exit(0);
    }

}

void SetPlanetProfile_Mat(planet_profile * _prof, double * depth, int * IdTarget, int NumTarget)
{
    int cur_cell = _prof->n - 1;
    for(int iTar=0;iTar<NumTarget;iTar++)
    {
        int TarIntervals = (int)floor(depth[iTar]/_prof->dz + 1);
        for(int j=0;j<TarIntervals && cur_cell >= 0;++j)
        {
            _prof->mat[cur_cell--] = IdTarget[iTar];
        }
    }

    for(;cur_cell>=0;--cur_cell)
    {
        _prof->mat[cur_cell] = IdTarget[NumTarget-1];
    }
}

void CapPlanetProfile(planet_profile * _prof, eos_table * _t, state_ref * _s)
{
    for(int k=0;k<_prof->n;++k)
    {
        // cap the temperature
        state_ref * local_sref = _s + _prof->mat[k];
        double melt_tem = SolidusTem(local_sref, _prof->pre[k]);
        double local_tem = Min(_prof->tem[k],melt_tem);
        _prof->tem[k] = local_tem;
        eos_table * local_eos = _t + _prof->mat[k];
        _prof->eng[k] = local_eos->InterpolateTP(local_eos,local_tem,_prof->pre[k],EOSENG);
    }
}


void ConsistentLith(planet_definition * _pdef)
{
    if(_pdef->t_lith >= 0 && _pdef->d_lith < 0)
    {
        // calculated d_lith automatically
        double Delta = 1.0 - 2.0*(_pdef->t_lith - _pdef->surf_tem)/(_pdef->r_planet*_pdef->surf_dtdz);
        Delta = Max(0.0,Delta);
        _pdef->d_lith = _pdef->r_planet * (1.0 - sqrt(Delta));
    } else if(_pdef->t_lith < 0 && _pdef->d_lith >= 0)
    {
        _pdef->t_lith = _pdef->surf_tem + _pdef->surf_dtdz*_pdef->d_lith*(1.0 - 0.5*_pdef->d_lith/_pdef->r_planet);
    } else
    {
        fprintf(stdout,"cannot set d_lith and t_lith at same time!\n");
        exit(0);
    }
}

void CorrectDensity(planet_profile * _prof, eos_table * _t, state_ref * _s)
{
    for(int k=0;k<_prof->n;++k)
    {
        // cap the temperature
        state_ref * local_sref = _s + _prof->mat[k];
        double melt_tem = SolidusTem(local_sref, _prof->pre[k]);
        double local_tem = Min(_prof->tem[k],melt_tem);
        // interpolate the density through eos
        eos_table * local_eos = _t + _prof->mat[k];
        _prof->den[k] = local_eos->InterpolateTP(local_eos,local_tem,_prof->pre[k],EOSDEN);
    }
}

void CorrectionDensity2(planet_profile * _prof, planet_definition * _pdef, eos_table * _t, state_ref * _s)
{
    // conduction lithosphere + convection mantle
    for(int k=_prof->n -1 ;k >= 0; k--)
    {
        // fix the surface temperature; just skip temperature updating of [n-1] index
        double k_depth = _prof->rad[_prof->n-1] - _prof->rad[k];
        // recalculate temperature according to the gravity/thermexp, using local gravity not the surf_gravity
        double dtdz_cond = _pdef->surf_dtdz * (1.0 + k_depth/_pdef->r_planet);
        _prof->tem[k] = _pdef->surf_tem;
        if( (_prof->n - 1) != k)
        {
            double conv_coef = _s[_prof->mat[k]].thermexp * fabs(_prof->gra[k]) / _s[_prof->mat[k]].c_heat;
            double dtdz_conv = _prof->tem[k+1] * conv_coef;
            // double dtdz_conv = _pdef->t_lith * conv_coef * exp(conv_coef * (k_depth - _pdef->d_lith));
            // double conduction = 1.0 - (_prof->tem[k+1] - 0.9*_pdef->t_lith)/(0.10*_pdef->t_lith);
            double conduction = 1.0 - (k_depth - 0.9*_pdef->d_lith)/(0.11*_pdef->d_lith);
            // double conduction = 1.0 - (_prof->tem[k+1] - 0.9*_pdef->t_lith)/(0.11*_pdef->t_lith);
            conduction = Wind(conduction,0.0,1.0);
            _prof->tem[k] = _prof->tem[k+1] + (dtdz_cond*conduction + dtdz_conv*(1.0 - conduction))*(_prof->rad[k+1] - _prof->rad[k]);
        }

        // cap the temperature greater than melt temperature
        state_ref * local_sref = _s + _prof->mat[k];
        // double melt_tem = local_sref->Tmelt0 * pow(local_sref->Asimon*_prof->pre[k] + 1.0,local_sref->Csimon);
        // replace melt_tem as solidus; if polysol is defined
        double melt_tem = SolidusTem(local_sref,_prof->pre[k]);
        double local_tem = _prof->tem[k];

        local_tem = Min(_prof->tem[k],melt_tem);

        _prof->tem[k] = local_tem;
        eos_table * local_eos = _t + _prof->mat[k];
        _prof->eng[k] = local_eos->InterpolateTP(local_eos,local_tem,_prof->pre[k],EOSENG);
        _prof->den[k] = local_eos->InterpolateTP(local_eos,local_tem,_prof->pre[k],EOSDEN);
    }
}

void WriteTargetProfile_All(FILE * fp,planet_profile * _prof, state_ref * _s, planet_definition * _pdef)
{
    fputs("      Depth, Matid,    Density,    Gravity, Pressure,Temperature,solidus\n",fp);
    for(int k= _prof->n-1; k>= 0; --k)
    {
        state_ref * local_sref = _s + _prof->mat[k];
        double melt_tem = SolidusTem(local_sref,_prof->pre[k]);
        fprintf(fp,"%10.5f,   %3d, %10.5f, %10.5f,%10.5e,%10.5e, %10.5e\n",_prof->rad[k], _prof->mat[k],_prof->den[k],
                _prof->gra[k] ,_prof->pre[k],_prof->tem[k],
                melt_tem);
    }
}

void IntergateGravity(planet_profile * _prof)
{
    /// typo error, do not change it
    _prof->gra[0] = 0.0;
    for(int k=1;k<_prof->n;++k)
    {
        _prof->gra[k] = _prof->gra[k-1] - 4.0*M_PI*6.674e-11*(_prof->den[k]+_prof->den[k-1])
                                          *(_prof->rad[k] + _prof->rad[k-1])*(_prof->rad[k] + _prof->rad[k-1])/8.0
                                          *(_prof->rad[k] - _prof->rad[k-1]);
    }
    for(int k=1;k<_prof->n;++k)
    {
        _prof->gra[k] = _prof->gra[k]/(_prof->rad[k]*_prof->rad[k]);
    }

    _prof->pre[_prof->n-1] = 0.0;
    for(int k=_prof->n-2;k>=0;--k)
    {
        _prof->pre[k] = _prof->pre[k+1] + (_prof->den[k] + _prof->den[k+1]) * (_prof->gra[k] + _prof->gra[k+1])
                                          * (_prof->rad[k] - _prof->rad[k+1])/4.0;
    }
}

void IntegrateFlatGravity(planet_profile * _prof)
{
    /// only update the pressure, keep gravity constant
    _prof->pre[_prof->n-1] = 0.0;
    for(int k=_prof->n-2;k>=0;--k)
    {
        _prof->pre[k] = _prof->pre[k+1] + (_prof->den[k] + _prof->den[k+1]) * (_prof->gra[k] + _prof->gra[k+1])
                                          * (_prof->rad[k] - _prof->rad[k+1])/4.0;
    }
}


void LoadPlanetDef(InputFile * ifp, planet_definition * _pdef)
{
    char TemProfOpt[MaxStrLen];
    GetValueSk(ifp,"target.temperature",TemProfOpt,0,"const");
    if(strcasecmp(TemProfOpt,"convective") == 0)
    {
        _pdef->tem_type = convective_profile;
    } else if(strcasecmp(TemProfOpt,"user") == 0)
    {
        _pdef->tem_type = user_profile;
    } else if(strcasecmp(TemProfOpt,"const") == 0)
    {
        _pdef->tem_type = const_profile;
    }
    else
    {
        _pdef->tem_type = unknown_profile;
        fprintf(stdout,"%s: UnImplemented tem_type:%s\n",__func__,TemProfOpt);
    }

    _pdef->surf_gravity = GetValueD(ifp,"target.surf_gravity","-1.63");
    _pdef->surf_dtdz = GetValueD(ifp,"target.surf_dtdz","290.0");
    _pdef->t_lith = GetValueD(ifp,"target.t_lith","-1.0");
    _pdef->d_lith = GetValueD(ifp,"target.d_lith","-1.0");
    _pdef->r_planet = GetValueD(ifp,"target.r_planet","1738.0e3");
    _pdef->surf_tem = GetValueD(ifp,"target.surf_temperature","293.0");
    _pdef->center[Z] = - _pdef->r_planet;
    _pdef->center[X] = _pdef->center[Y] = 0.0;
    ConsistentLith(_pdef);
}

planet_profile * BuildTargetProf(InputFile * ifp, struct EosTable * _etb, struct StateReference * _sref)
{
    /// build thermal profile for planet model, not used !
    char Material[NMT][MaxStrLen];
    int NumMaterial = GetValueI(ifp,"material.nm","1.0");
    for(int iMat=0;iMat<NumMaterial;iMat++)
    {
        GetValueSk(ifp,"material.name",Material[iMat],iMat,"dunite_");
    }

    int NumTarget = GetValueI(ifp,"target.number","1");
    double dz = 0.25*GetValueD(ifp,"mesh.dz","1.0e4");
    double depth[NMT*4],SumDepth;
    int IdTarget[NMT*4];
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

    struct PlanetDefinition pdef;
    LoadPlanetDef(ifp,&pdef);

    struct PlanetProfile tprofile;
    tprofile.n = (int)(pdef.r_planet/dz) + 1;
    tprofile.dz = dz;
    AllocatePlanetProfile(&tprofile);

    SetPlanetProfile_Mat(&tprofile, depth, IdTarget, NumTarget);
    SetPlanetProfile_Tem(&tprofile,_sref,&pdef);

    char TarDamOpt[MaxStrLen];
    GetValueSk(ifp,"target.damage",TarDamOpt,0,"unknown");

    if(0==strcasecmp(TarDamOpt,"const"))
    {
        double Damage = GetValueDk(ifp,"target.damage",1,"0.0");
        for(int k=0;k<tprofile.n;++k)
        {
            tprofile.dam[k] = Damage;
        }
    } else
    {
        fprintf(stdout,"undefined damage profile option: %s\n",TarDamOpt);
        exit(0);
    }


    for(int cycle=0;cycle<100;++cycle)
    {
        IntergateGravity(&tprofile);

        if(pdef.tem_type == convective_profile || pdef.tem_type == condconv_profile)
        {
            CorrectionDensity2(&tprofile,&pdef,_etb,_sref);
        }
        else
            CorrectDensity(&tprofile,_etb,_sref);
    }

    CapPlanetProfile(&tprofile, _etb, _sref);
    struct PlanetProfile * _prof = malloc(sizeof(struct PlanetProfile));
    *_prof = tprofile;
    return _prof;
}

planet_profile * BuildFlatTargetProf(InputFile * ifp, struct EosTable * _etb, struct StateReference * _sref)
{
    /// copy from BuildTargetProf
    /// change the gravity to constant
    char Material[NMT][MaxStrLen];
    int NumMaterial = GetValueI(ifp,"material.nm","1.0");
    for(int iMat=0;iMat<NumMaterial;iMat++)
    {
        GetValueSk(ifp,"material.name",Material[iMat],iMat,"dunite_");
    }

    int NumTarget = GetValueI(ifp,"target.number","1");
    double dz = 0.25*GetValueD(ifp,"mesh.dz","1.0e4");
    double depth[NMT*4],SumDepth;
    int IdTarget[NMT*4];
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

    struct PlanetDefinition pdef;
    LoadPlanetDef(ifp,&pdef);

    struct PlanetProfile tprofile;
    tprofile.n = (int)(pdef.r_planet/dz) + 1;
    tprofile.dz = dz;
    AllocatePlanetProfile(&tprofile);

    SetPlanetProfile_Mat(&tprofile, depth, IdTarget, NumTarget);
    SetPlanetProfile_Tem(&tprofile,_sref,&pdef);

    char TarDamOpt[MaxStrLen];
    GetValueSk(ifp,"target.damage",TarDamOpt,0,"unknown");

    if(0==strcasecmp(TarDamOpt,"const"))
    {
        double Damage = GetValueDk(ifp,"target.damage",1,"0.0");
        for(int k=0;k<tprofile.n;++k)
        {
            tprofile.dam[k] = Damage;
        }
    } else
    {
        fprintf(stdout,"undefined damage profile option: %s\n",TarDamOpt);
        exit(0);
    }

    char TarGrav[MaxStrLen];
    double gravity = 0.0;
    GetValueSk(ifp,"condition.gravity",TarGrav,0,"unknown");
    if(0==strcasecmp(TarGrav,"const"))
    {
        gravity = GetValueDk(ifp, "condition.gravity", 1, "-9.8");
        for(int k = 0; k < tprofile.n; ++k)
        {
            tprofile.gra[k] = gravity;
        }
    } else
    {
        fprintf(stdout,"%s:error gravity type:%s\n",__func__,TarGrav);
        exit(0);
    }

    for(int cycle=0;cycle<100;++cycle)
    {
        //IntergateGravity(&tprofile);
        IntegrateFlatGravity(&tprofile);

        if(pdef.tem_type == convective_profile || pdef.tem_type == condconv_profile)
        {
            CorrectionDensity2(&tprofile,&pdef,_etb,_sref);
        }
        else
            CorrectDensity(&tprofile,_etb,_sref);
    }

    CapPlanetProfile(&tprofile, _etb, _sref);
    struct PlanetProfile * _prof = malloc(sizeof(struct PlanetProfile));
    *_prof = tprofile;
    return _prof;
}


int CopyProfToE(planet_profile * _prof, int dpz, struct Element * ie)
{
    assert(_prof != NULL);
    if(dpz >= 0 && dpz <= _prof->n-1)
    {
        // set as full materials
        ie->Temperature = _prof->tem[dpz];
        ie->Pressure    = _prof->pre[dpz];
        ie->Density     = _prof->den[dpz];

        int local_mat = _prof->mat[dpz];
        ie->VOF[local_mat] = 1.0;
        ie->PhiL3[local_mat][1] = _prof->den[dpz];
        ie->PhiL3[local_mat][2] = _prof->eng[dpz];

        ScalerMoveEx(ie->sPhiL3[local_mat],ie->PhiL3[local_mat],1.0,NPHI);
        ie->Damage = ie->PsiL3[0] = _prof->dam[dpz];
        ie->Mass = ie->Density * ie->Volume;
        ie->sPsiL3[0] = ie->Mass*ie->Damage;
        return local_mat;
    }
    else
    {
        // set to vacuum
        ie->Temperature = 0.0;
        ie->Pressure    = 0.0;
        ie->Density     = 0.0;
        ie->VOF[0] = 1.0;
        ie->sPhiL3[0][0] = 1.0;
        ie->Mass = 0.0;
        return 0;
    }
}

void InitPlanetMesh(InputFile * ifp, struct Mesh * _mesh, struct EosTable * _etb, struct StateReference * _sref)
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
    double depth[NMT*4],SumDepth;
    int IdTarget[NMT*4];
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

    struct PlanetDefinition pdef;
    LoadPlanetDef(ifp,&pdef);
    struct PlanetProfile tprofile;
    tprofile.n = (int)(pdef.r_planet/dz) + 1;
    tprofile.dz = dz;
    AllocatePlanetProfile(&tprofile);

    SetPlanetProfile_Mat(&tprofile, depth, IdTarget, NumTarget);
    SetPlanetProfile_Tem(&tprofile,_sref,&pdef);

    char TarDamOpt[MaxStrLen];
    GetValueSk(ifp,"target.damage",TarDamOpt,0,"unknown");

    if(0==strcasecmp(TarDamOpt,"const"))
    {
        double Damage = GetValueDk(ifp,"target.damage",1,"0.0");
        for(int k=0;k<tprofile.n;++k)
        {
            tprofile.dam[k] = Damage;
        }
    } else
    {
        fprintf(stdout,"undefined damage profile option: %s\n",TarDamOpt);
        exit(0);
    }

    for(int cycle=0;cycle<100;++cycle)
    {
        IntergateGravity(&tprofile);

        if(pdef.tem_type == convective_profile || pdef.tem_type == condconv_profile)
        {
            CorrectionDensity2(&tprofile,&pdef,_etb,_sref);
        }
        else
            CorrectDensity(&tprofile,_etb,_sref);
    }

    CapPlanetProfile(&tprofile, _etb, _sref);
    // initialization of element
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
        } else if(Distance(pdef.center,ie->Center) <= pdef.r_planet)
        {
            int dpz = (int) (Wind(Distance(pdef.center,ie->Center)/tprofile.dz,0,tprofile.n-1));
            ie->Temperature = tprofile.tem[dpz];
            ie->Pressure    = tprofile.pre[dpz];
            ie->Density     = tprofile.den[dpz];

            int local_mat = tprofile.mat[dpz];
            ie->VOF[local_mat] = 1.0;
            ie->PhiL3[local_mat][1] = tprofile.den[dpz];
            ie->PhiL3[local_mat][2] = tprofile.eng[dpz];

            ScalerMoveEx(ie->sPhiL3[local_mat],ie->PhiL3[local_mat],1.0,NPHI);
            ie->Damage = ie->PsiL3[0] = tprofile.dam[dpz];
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

    // initialization of vertex
    for(int k=0;k<_mesh->nv;k++)
    {
        struct Vertex * iv = _mesh->Vertexes + k;
        double local_gravity[DIM];
        LinearOp(local_gravity,iv->Position,pdef.center,1.0,-1.0);
        double r = Length(local_gravity);
        double r0 = pdef.r_planet;
        NormalizationEx(local_gravity,1.0);
        int dpz = (int) (Wind(Distance(pdef.center,iv->Position)/tprofile.dz,0,tprofile.n-1));
	if( r > r0){
        Scaler(local_gravity,tprofile.gra[dpz]*(r0/r)*(r0/r));
	}
	else
	{
        Scaler(local_gravity,tprofile.gra[dpz]);
	}
        Copy(local_gravity,iv->BodyF);
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

    int lrank;
    MPI_Comm_rank(SALEC_CART_COMM, &lrank);
    if(lrank == 1)
    {
        FILE * fp = fopen("txt/profile.txt","w");
        WriteTargetProfile_All(fp,&tprofile,_sref,&pdef);
        fclose(fp);
    }

    UnAllocatePlanetProfile(&tprofile);
}

void LoadCycleControl(InputFile * ifp, cycle_control * _ccl)
{
    GetValueS(ifp, "output.prefix", _ccl->prefix, "ParaTest");
    _ccl->step_time = GetValueD(ifp,"cycle.time0","0.0");
    _ccl->start_time = _ccl->step_time;
    _ccl->end_time = GetValueD(ifp,"cycle.time1","5.0");
    _ccl->dt = GetValueD(ifp,"cycle.dt0","1.0e-3");
    _ccl->max_dt = GetValueD(ifp,"cycle.maxdt","1.0e-2");
    _ccl->out_time = GetValueD(ifp,"output.time","1.0");
    _ccl->max_step = GetValueI(ifp,"cycle.maxstep","-1");
    if(_ccl->max_step < 0)
    {
        _ccl->max_step = INT_MAX;
    }
    _ccl->k_out_step = 0;
    _ccl->init_step = GetValueI(ifp,"output.init","0");
    _ccl->i_step = 0;

    _ccl->out_tracer = 0;
    _ccl->cycle_continue = 1;
    _ccl->out_log = 0;
    _ccl->out_vts = 0;
    _ccl->out_init = 0;
}

void CheckCycleControl(cycle_control * _ccl)
{
    if(0 == _ccl->i_step%100)
        _ccl->out_log = 1;
    else
        _ccl->out_log = 0;

    if(fabs(_ccl->step_time - _ccl->k_out_step * _ccl->out_time - _ccl->start_time) <= 0.5*_ccl->max_dt)
        _ccl->out_vts = 1;
    else
        _ccl->out_vts = 0;

    if(_ccl->i_step <= _ccl->init_step)
        _ccl->out_init = 1;
    else
        _ccl->out_init = 0;
}

void UpdateCycleControl(cycle_control * _ccl)
{
    _ccl->i_step ++;
    _ccl->step_time += _ccl->dt;

    if(_ccl->i_step < _ccl->max_step)
        _ccl->cycle_continue = 1;
    else
        _ccl->cycle_continue = 0;

    if(_ccl->step_time > _ccl->end_time + _ccl->max_dt)
        _ccl->cycle_continue = 0;
}
