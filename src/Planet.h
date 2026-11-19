//
// Created by huacheng on 11/8/22.
//
// although named after planet, this file includes subroutines that initialize thermal prof for flat target.

#ifndef SALEC_PLANET_H
#define SALEC_PLANET_H
#include "State.h"
#include "InputParser.h"

typedef enum TemProfType
{
    conductive_profile,
    convective_profile,
    condconv_profile,
    const_profile,
    lunar_profile,
    user_profile,
    unknown_profile
} tem_prof_type;

typedef struct PlanetDefinition
{
    double t_lith;
    double d_lith;
    double surf_gravity;
    double surf_tem;
    double surf_dtdz;
    double r_planet;
    double center[DIM];
    tem_prof_type tem_type;
} planet_definition;

void LoadPlanetDef(InputFile * ifp,planet_definition * _pdef);

typedef struct PlanetProfile
{
    double * pre;
    double * rad;
    double * den;
    double * tem;
    double * gra;
    double * dam;
    double * eng;
    int * mat;
    int n;
    double dz;
} planet_profile;

planet_profile * BuildTargetProf(InputFile * ifp, struct EosTable * _etb, struct StateReference * _sref);
planet_profile * BuildFlatTargetProf(InputFile * ifp, struct EosTable * _etb, struct StateReference * _sref);
void WriteTargetProfile_All(FILE * fp,planet_profile * _prof, state_ref * _s, planet_definition * _pdef);
int CopyProfToE(planet_profile * _prof, int dpz, struct Element * ie);
// void BuildPlanetProf(InputFile * ifp, planet_profile * _prof, struct EosTable * _etb, struct StateReference * _sref);
void AllocatePlanetProfile(planet_profile * _prof);
void UnAllocatePlanetProfile(planet_profile * _prof);
void ConsistentLith(planet_definition * _pdef);
void SetPlanetProfile_Tem(planet_profile * _prof, state_ref * _s, planet_definition * _pdef);
void SetPlanetProfile_Mat(planet_profile * _prof, double * depth, int * IdTarget, int NumTarget);
void CapPlanetProfile(planet_profile * _prof, eos_table * _t, state_ref * _s);
void CorrectDensity(planet_profile * _prof, eos_table * _t, state_ref * _s);
void CorrectionDensity2(planet_profile * _prof, planet_definition * _pdef, eos_table * _t, state_ref * _s);
void IntergateGravity(planet_profile * _prof);
void IntegrateFlatGravity(planet_profile * _prof);
void InitPlanetMesh(InputFile * ifp, struct Mesh * _mesh, struct EosTable * _etb, struct StateReference * _sref);

typedef struct CycleControl
{
    double start_time;
    double end_time;
    int max_step;
    double step_time;
    double dt;
    double max_dt;
    double out_time;
    char prefix[MaxStrLen];

    int init_step;
    int k_out_step;
    int i_step;

    int out_tracer;
    int out_vts;
    int out_log;
    int out_init;
    int cycle_continue;
} cycle_control;

void LoadCycleControl(InputFile * ifp, cycle_control * _ccl);
void CheckCycleControl(cycle_control * _ccl);
void UpdateCycleControl(cycle_control * _ccl);
#endif //SALEC_PLANET_H
