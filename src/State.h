//
// Created by huacheng on 2021/1/22.
//

#ifndef SALEC_STATE_H
#define SALEC_STATE_H
/*
 * updating the pressure / eos
 * this first favor type of eos is ANEOS which is inherited from sale2d
 * only the Temperature-Density type is implemented
 */

#include "Variable.h"

#define SolidElement 0
#define FluidElement 1
#define VaporElement 3
#define VoidElement 2


struct ANEOSTable
{
    /*
     * structured table, interpolate the value linearly
     * between data points
     */
    int nTem;
    int nDen;
    double Pnorm;
    double Dnorm;
    double Tnorm;
    double *yDen;
    double *xTem;
    double *** Data;
    // data is (density*energy,pressure,csound)
    // not same with iSALE2d
};

struct StateReference;
typedef double (*StrengthFunc)(struct Element *, const struct StateReference *);
typedef double (*DamageFunc)(struct Element *, const struct StateReference *);
typedef double (*SoftFunc)(struct Element *, const struct StateReference *);
struct StateReference
{
    /*
     * The parameters shared in one materials
     * which can be used in EOS or strength model.
     */
    int MaterialId;
    double VaporDen;
    double MeltDen;
    double GaCoef; // Constant to convert bulk modulus to shear modulus

    // parameters for simple ROCK strength model
    double Yint0;
    double Yfricint;
    double Ylimint;
    double Ydam0;
    double Yfricdam;
    double Ylimdam;
    double Yten0;

    double IvanA;
    double IvanB;
    double IvanC;
    StrengthFunc Yfunc; // function to calculte the yield strength(shear)
    StrengthFunc Tfunc; //          for tensile strength
    DamageFunc   dSDf;  // the increment of damage (shear failure)

    double minPint;
    double minPdam;


    // parameters for Johnson-Cook strength model
    double JCA;
    double JCB;
    double JCN;
    double JCC;
    double JCM;
    double JCTREF;
    double JCMaxPlastic;

    // soft parameters
    double Asimon;
    double Csimon;
    double Tmelt0;
    double Tfrac;
    double Viscosity;

    // Acoustic Fluidization parameters
    double Toff;
    double Cvib;
    double VibMax;
    double Tdamp;
    double Pvlim;
    double Acvis;
    double GammaEta;
    double GammaBeta;
};

void AllocateANEOS(struct ANEOSTable *_t);
void UnAllocateANEOS(struct ANEOSTable *_t);
void LoadANEOS(struct ANEOSTable *_t, FILE *fp);
void ANEOSCalL3t(struct Element *_e, const struct ANEOSTable *_t, int _k);
void ANEOSInitStateRef(struct ANEOSTable *_t, struct StateReference *_s);
int GetState(struct Element *_e, const struct StateReference _s[], int nm);
void Over(FILE * fp,int n);
double InterpolateTD(struct ANEOSTable *_t,double tTem, double tDen, int DataId);
double ANEOSInterpolateTP(struct ANEOSTable *_t, double tTem, double tPre, int DataId);

double SimpleTensile(struct Element *_e, const struct StateReference *_s);
double SimpleRockStrength(struct Element *_e, const struct StateReference *_s);
double SimpleRockNoacfl(struct Element *_e, const struct StateReference *_s);
double SimpleSheardam(struct Element *_e, const struct StateReference *_s);
double SimpleSoft(struct Element *_e, const struct StateReference *_s);
double JNCKSoft(struct Element *_e, const struct StateReference *_s);
double SimonMelt(struct Element *_e, const struct StateReference *_s);
double OnnakaSoft(struct Element *_e, const struct StateReference *_s, double Tmelt);
double Ylundborg(double p,double y0,double fric,double ylim);
double Ydrucker(double p,double y0,double fric,double ylim);
double LowdensitySoft(double Density, double MeltDensity);
double BlockVibStrength(struct Element *_e, const struct StateReference *_s);
int BlockVibration(struct Element * ie, double dt,struct StateReference * _s, int nm,double time);

struct AirEOSTable
{
    int nEng;
    int nDen;
    double * xEng;
    double * yDen;
    double *** Data;
    // (All the variables is log10. form)
    // Data = log.10 (Pressure,Temperature,Speed of sound)
};

int BilinearInverse(const double xi[],const double yi[], double xl[]);
void AllocateAirEOS(struct AirEOSTable * _air);
void UnAllocateAirEOS(struct AirEOSTable * _air);
void LoadAirEOS(struct AirEOSTable * _t, FILE *fp);
void CalculateL3atm(struct Element *_e, const struct AirEOSTable *_gas,int _k);
double ANEOSPresProfRK3(struct ANEOSTable *_t, double Pres[], double Grav[], double Temp[], int steps, double dh);

double JohnsonCook2(struct Element *_e, const struct StateReference *_s);
double JohnsonCook1(struct Element *_e, const struct StateReference *_s);

struct TillotsonTable
{
    double TLRho0;
    double TLCv;
    double TLA;
    double TLB;
    double TLE0;
    double TLa;
    double TLb;
    double TLAlpha;
    double TLBeta;
    double TLEiv;
    double TLEcv;
    double TLTref;

    double * xDen;
    double * yEng;
    double StepDen;
    int nDen;
};

double TillPres(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs);
double TillTemp(const struct TillotsonTable * _t,double EngIn, double RhoIn);
void TillEOSCalL3t(struct Element *_e, const struct TillotsonTable * _t, int _k);

struct EosTable;
typedef void (*CalL3tFunc)(struct Element *_e, const struct EosTable * _t, int _k);
typedef double (*InterpolateFunc)(struct EosTable *_t, double tTem, double tPre, int DataId);
struct EosTable
{
    struct ANEOSTable     * _atb;
    struct TillotsonTable * _ttb;
    struct AirEOSTable    * _gtb;
    CalL3tFunc CalL3t;
    InterpolateFunc InterpolateTP;
};

void TillEOSCalL3tCompose(struct Element *_e, const struct EosTable * _t, int _k);
void ANEOSCalL3tCompose(struct Element *_e, const struct EosTable * _t, int _k);
void TillColdEnergy(struct TillotsonTable * _t);
void TillInitStateRef(struct TillotsonTable *_t, struct StateReference *_s);
double TillEOSInterpolateTP(struct TillotsonTable *_t, double tTem, double tPre, int DataId);
double TillEOSInterpolateTPCompose(struct EosTable *_t, double tTem, double tPre, int DataId);
double ANEOSInterpolateTPCompose(struct EosTable *_t, double tTem, double tPre, int DataId);
double TillEOSPresProfRK3(struct TillotsonTable *_t, double Pres[], double Grav[], double Temp[], int steps, double dh);
double PresProfRK3(struct EosTable *_t, double Pres[], double Grav[], double Temp[], int steps, double dh);

StrengthFunc SelectStrength(const char s[]);
DamageFunc SelectDamage(const char s[]);
#endif //SALEC_STATE_H
