//
// Created by huacheng on 2020/12/22.
//

#ifndef SALEC_VARIABLE_H
#define SALEC_VARIABLE_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
// The basic dimension information related to the structured points
#define DIM  3  // mesh dimension
#define DIM3 27
#define VPB  4  // nodes/vertexes per boundary
#define EPB  2  // neighbor elements per boundary
#define EPV  8  // neighbor elements per node/vertex
#define VPE  8  // neighbor nodes per element
#define BPE  6  // neighbor boundaries per element
#define TOL    1e-7
#define TOLVOF 1e-5
#define TOLRHO 1e-3
#define TOLSTR 1e-6
#define TOLMASS  1e-6
#define TOLVEL   1e-5
#define COURANT0 0.5
#define COURANT1 0.5
#define COURANT2 2.0
#define NMT    4   // the max number of material types in model
#define NPHI   3   // the number of conservation variables in L3
#define NTHETA 3   // the number of variables determined by EOS
#define NPSI   10  // the variable effected by advection
// Index for coordinate line
#define X    0
#define Y    1
#define Z    2
#define BETA 0.57735026918962576451
#define NTHREAD 2
#define OABF 3    // order of Adams-Bashforth formula

#define TDSIZE  100
// Index for stress tensor component
#define XX 0
#define YY 1
#define ZZ 2
#define XY 3
#define YX 3
#define XZ 4
#define ZX 4
#define ZY 5
#define YZ 5

// Index for position of cell boundaries
#define Left                0
#define Right               1
#define Front               2
#define Back                3
#define Bottom              4
#define Top                 5

// Abbreviation : Index for position of cell boundaries
#define L                   0  // (Left)
#define R                   1  //
#define F                   2
#define K                   3  // (Back)
#define B                   4  // (Bottom)
#define T                   5

// Index for position of vertexs of a cell
#define BFR    0
#define BKR    1
#define BKL    2
#define BFL    3
#define TFR    4
#define TKR    5
#define TKL    6
#define TFL    7

// other values to define the element
#define NSBPB 4
#define NSEPE 8

// some fixed material Id
#define VACUUM 0
struct Vnode{
    int  vId;
    double d;
};
// An incomplete statement used in struct Element,
// which will be complete stated following struct Element.
struct Vertex;
struct Boundary;

#define ELEMENT_DEBUG 0
struct Element
{
    double Center[DIM];
    double Scale[DIM];
    double Volume;
    double rVolume;
    double rdx_min;

    double Mass;
    double Density;
    double Damage;
    double Momentum[DIM];

    struct Vertex   *NeV[VPE];
    struct Boundary *NeB[BPE];
    struct Vnode    dCut[VPE];
    double CutRatio[VPE];

    unsigned short CutTag;
    double SubFace[VPE][DIM];
    double ABVArea;

    double GradFactor[DIM][VPE];
    double GradVel[DIM][DIM];
    double ArtificialPressure;
    double Strain[DIM][DIM];
    double DivVel;
    double InvSta2;
    double InvSta3;
    double sInvSte2; // sqrt of 2rd invariant of stress
    double sInvSte2_old; // used in plastic .
    double sInvSta2; // sqrt of 2rd invariant of strain
    // the invariants of strain

    // used in VOF related functions
    double GradVOF[NMT][DIM];
    double VOF[NMT];
    double Courant;

    /*
     * the variables used for L3 calculation
     * they conservation in calculation domain
     */

    double PhiL3[NMT][NPHI];
    double sPhiL3[NMT][NPHI];
    double PsiL3[NPSI];
    double sPsiL3[NPSI];

    // ThetaL3 is some variables updated by EOS
    double ThetaL3[NMT][NTHETA];
    double Temperature;
    double Pressure;
    double Csound;
    double ShearModule;

    double SubVolumeRatio[NSEPE];

    unsigned short State;
    double ShearStrength;
    double TensileStrength;
    unsigned short FailureState;

    double Viscosity;

    double VibPressure;

#if ELEMENT_DEBUG
    double Debug;
#endif
};

struct Vertex
{
    unsigned int State;
    double dx_min;
    double Position[DIM];
    struct Element *NeE[EPV];

    double Velocity[DIM];

    double VelocityL1[DIM];
    double MtVelocityL1[OABF-1][DIM];
    double DistanceL1[DIM];
    double BodyF[DIM];

    double VOF[NMT];
    double MassL1;
    double CsoundL1;
    double Damping;

    double Volume;
    double rVolume;
    double SubVolumeRatio[NSEPE];
};

struct Boundary
{
    struct Vertex  *NeV[VPB];
    struct Element *UpE[EPB];
    struct Element *NeE[EPB];

    int bBID[EPB];
    double Normal[DIM];
    double FluxFactor[VPB][DIM];

    // calculated  CalculateBFV per cycle
    double FluxVolume;

    // the conservation flows get through this boundary
    double dPhi[NMT][NPHI];
    double dMomentum[DIM];
    double dMass;
    double dPsi[NPSI];
};

struct ConditionE
{/*
 * Now only
 * [CEType] = 0 -- Copy
 * [CEType] = 1 -- Iverse in Wnorm
 * [CEType] = 2 -- Fixed
 */
    int CEType;
    double Tr[DIM];
    double Wnorm[DIM][DIM];
    struct Element * dest;
    struct Element * refe;
//    struct Boundary * NeB[2];
};

#define NRCV 3
struct ConditionV
{
    struct Vertex * dest;
    struct Vertex * refe[NRCV];
};

struct ConditionB
{
    struct Boundary * dest;
    struct Boundary * refe;
};


struct Mesh
{
    // some global variables that can be adjusted by setup information
    int nxp; // number of points in dimension X
    int nyp; // number of points in dimension Y
    int nzp; // number of points in dimension Z
    int ne; // number of Elements
    int nv; // number of vertexes
    int nbx; // number of bounderies in direction X
    int nby; // ..                   in direction Y
    int nbz; // ..                   in direction Z
    struct Element  *Elements;
    struct Vertex   *Vertexes;
    struct Boundary *XBoundaries;
    struct Boundary *YBoundaries;
    struct Boundary *ZBoundaries;

    struct Vertex ** DomainV;
    int ndv;
    int ndv1;
    struct Boundary ** DomainB;
    int ndb;
    struct Element ** DomainE;
    struct ConditionE * CondE;
    struct ConditionV * CondV;
    struct ConditionB * CondB;
    int nde;
    int nde1;
    int nce;
    int ncv;
    int ncb;

    char vtk_format[100];
    int nvof;
    int noffset;
};

struct MeshInfo
{
    int npx;
    int npy;
    int npz;
    double ARTVIS;
    double ARTVIS2;
    double CISSAMK1;
    double MaxVelocity;
    double DampZ0;
    double DampZ1;
    double DampLambda;
    double MinPres;

    int PartPressure;
    int TensileFailure;
    int MaterialNum;
    FILE * Materialfp[NMT];
    char EosType[NMT][100];

    double **** V; // the position of points
    char output_format[100];
    char chkPrefix[100];
    int chkStep;
};

#define NGI2d 4
#define NGI3d 8
// #define NSPE   9    // number of sub-element in every big element
// value of 1/sqrt(3) used in position of gaussian integral
//extern const double BETA;
// points of Gaussian integral for 2*2 in 2d
extern const double GIPS2d[NGI2d][2];
// points of Gaussian integral for 2*2*2 in 3d
extern const double GIPS3d[NGI3d][3];
// weights of Gaussian integral for 2*2 in 2d
extern const double GIWS2d[NGI2d];
// weights of Gaussian integral for 2*2*2 in 3d
extern const double GIWS3d[NGI3d];

// definition for the interpolation function of boundary/volume
// interpolation function for boundary within 4 points
typedef double (*Interpolation)(const double *_x);

#define NIpB 4
#define NIpV 8

extern Interpolation IpB[NIpB];
extern Interpolation IpB_X[NIpB][2];
// interpolation function for Volume with in 8 points
extern Interpolation IpV[NIpV];
extern Interpolation IpV_X[NIpV][3];

struct Vertex *GetVertex(struct Mesh *_mesh, int Idx, int Idy, int Idz);
struct Element *GetElement(struct Mesh *_mesh, int Idx, int Idy, int Idz);
struct Boundary *GetXBoundaries(struct Mesh *_mesh, int Idx, int Idy, int Idz);
struct Boundary *GetYBoundaries(struct Mesh *_mesh, int Idx, int Idy, int Idz);
struct Boundary *GetZBoundaries(struct Mesh *_mesh, int Idx, int Idy, int Idz);

// sub-element in every element : local coordinate
extern double SEPE[NSEPE][VPE][DIM];
// points of every element : local coordinate
extern const double ELC[VPE][DIM];
extern const int PCVE[VPE];

// sub-boundary in every boundary : local coordinate
extern double SBPB[NSBPB][VPB][DIM - 1];
// points of every boundary : local coordinate
extern const double BLC[VPB][DIM - 1];
// interal surface in every element : local coordinate
extern double ISFE[DIM][VPB][DIM];
void DeriveSEPE();
void DeriveSBPB();
extern const int ToOuterFace[BPE][VPB];

// functions for linear operation of vector
double Sqrt_trunc(double x);
double Max(double a, double b);
double Min(double a, double b);
void Copy(const double _a[], double _b[]);
void Addition(double *_a, const double *_b);
void Scaler(double _a[], double _p);
void Normalization(double p[]);
void NormalizationEx(double p[], double Scale);
double Length(const double p[]);
void Zero(double _a[]);
double Dot(const double *_a, const double *_b);
double Contraction(const double _a[][DIM], const double _b[][DIM]);
void ScalerAddition(double *_a, const double *_b, double _p);
void LinearOp(double *_a, const double * _b, const double * _c, double _pb, double _pc);

void ScalerMove(double *_a, const double *_b, double _p);
double MatrixTime(double _a[], const double _m[][DIM], const double _b[]);
double Determinant(double tM[][DIM]);
double SecondInvariant(double tM[][DIM]);
double MaxScaler(const double Xi[], double p);
double MaxComponent(const double Xi[], int d);
double Wind(double x, double limitL, double limitR);

// a set of extended linear operation
void ScalerAdditionEx(double *_a, const double *_b, double _p, int d);
void ScalerAdditionEx1(double *_a, const double *_b, double _p);
void ScalerAdditionEx2(double *_a, const double *_b, double _p);
void ZeroEx(double _a[], int d);
void CopyEx(const double _a[], double _b[], int _d);
void WindEx(double limitL, double limitR, double _a[], int _d);
void ScalerMoveEx(double *_a, const double *_b, double _p, int d);
void ScalerEx(double _a[], double _p, int _d);

// functions to calculated the parameters related to the geometry of Mesh
void DeriveSubVolume(double Xi[][DIM], double Vol[]);
void DeriveInterFace(double Xi[][DIM], double InterFace[][NSBPB][DIM]);
void DeriveSubFace(double InterFace[][NSBPB][DIM], double SubFace[][DIM]);
void DeriveGradFactor(double Xi[][DIM], double GradFactor[][VPE]);
void DeriveOuterFace(double Xi[][DIM], double InterFace[][DIM]);
void DeriveFlowFactor(double Xi[][DIM], double FlowFactor[][DIM]);
void DeriveArea(double Xi[][DIM], double a[]);
double DeriveCenter(double Xi[][DIM], double _c[]);
double DeriveScale(double Xi[][DIM], double _c[]);

int GetBit(int x,int pos);
void SetBit(int * x,int pos);
int SumBit(int n);
int ReverseBit(int n);
double MinMod(double a, double b);
double Distance(const double x1[], const double x2[]);
double Sign(double x);

// function for Symmetrical matrix
int SymIndex(int i,int j);
double MatrixTimeSym(double _a[], const double _m[], const double _b[]);
double ContractionSym(const double _a[],double _b[][DIM]);
double SecondInvariantSym(const double tM[]);
double DeterminantSym(const double tM[]);

int ToLongId(int Idx, int Idy, int Idz, int nxp, int nyp, int nzp);
int ToShortId(int * Idx, int * Idy, int * Idz,int LongId, int nxp, int nyp, int nzp);

#endif //SALEC_VARIABLE_H
