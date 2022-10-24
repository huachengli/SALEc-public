//
// Created by huacheng on 2020/12/22.
//
#include "Variable.h"

struct Vertex *GetVertex(struct Mesh *_mesh, int Idx, int Idy, int Idz)
{
    if (Idx >= _mesh->nxp || Idy >= _mesh->nyp || Idz >= _mesh->nzp)
        return NULL;
    if (Idx < 0 || Idy < 0 || Idz < 0)
        return NULL;

    int LongId = Idx + Idy * _mesh->nxp + Idz * _mesh->nxp * _mesh->nyp;
    return &(_mesh->Vertexes[LongId]);
}

struct Element *GetElement(struct Mesh *_mesh, int Idx, int Idy, int Idz)
{
    if (Idx >= (_mesh->nxp - 1) || Idy >= (_mesh->nyp - 1) || Idz >= (_mesh->nzp - 1))
        return NULL;
    if (Idx < 0 || Idy < 0 || Idz < 0)
        return NULL;

    int LongId = Idx + Idy * (_mesh->nxp - 1) + Idz * (_mesh->nxp - 1) * (_mesh->nyp - 1);
    return &(_mesh->Elements[LongId]);
}

struct Boundary *GetXBoundaries(struct Mesh *_mesh, int Idx, int Idy, int Idz)
{
    if (Idx >= _mesh->nxp || Idy >= (_mesh->nyp - 1) || Idz >= (_mesh->nzp - 1))
        return NULL;
    if (Idx < 0 || Idy < 0 || Idz < 0)
        return NULL;

    int LongId = Idx + Idy * _mesh->nxp + Idz * _mesh->nxp * (_mesh->nyp - 1);
    return &(_mesh->XBoundaries[LongId]);
}

struct Boundary *GetYBoundaries(struct Mesh *_mesh, int Idx, int Idy, int Idz)
{
    if (Idx >= (_mesh->nxp - 1) || Idy >= _mesh->nyp || Idz >= (_mesh->nzp - 1))
        return NULL;
    if (Idx < 0 || Idy < 0 || Idz < 0)
        return NULL;
    int LongId = Idx + Idy * (_mesh->nxp - 1) + Idz * (_mesh->nxp - 1) * _mesh->nyp;
    return &(_mesh->YBoundaries[LongId]);
}

struct Boundary *GetZBoundaries(struct Mesh *_mesh, int Idx, int Idy, int Idz)
{
    if (Idx >= (_mesh->nxp - 1) || Idy >= (_mesh->nyp - 1) || Idz >= _mesh->nzp)
        return NULL;
    if (Idx < 0 || Idy < 0 || Idz < 0)
        return NULL;
    int LongId = Idx + Idy * (_mesh->nxp - 1) + Idz * (_mesh->nxp - 1) * (_mesh->nyp - 1);
    return &(_mesh->ZBoundaries[LongId]);
}

//const double BETA = 0.57735026918962576451;
const double GIPS2d[NGI2d][2] = {
        {-BETA, -BETA},
        {-BETA, BETA},
        {BETA,  -BETA},
        {BETA,  BETA}};
const double GIPS3d[NGI3d][3] = {
        {-BETA, -BETA, -BETA},
        {-BETA, -BETA, BETA},
        {-BETA, BETA,  -BETA},
        {-BETA, BETA,  BETA},
        {BETA,  -BETA, -BETA},
        {BETA,  -BETA, BETA},
        {BETA,  BETA,  -BETA},
        {BETA,  BETA,  BETA}};
const double GIWS2d[NGI2d] = {1.0, 1.0, 1.0, 1.0};
const double GIWS3d[NGI3d] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

// functions used in 4 points interpolate
double Ip2d0(const double *_x)
{
    return 0.25 * (1 + _x[0]) * (1 + _x[1]);
}

double Ip2d0_X1(const double *_x)
{
    return 0.25 * (1 + _x[1]);
}

double Ip2d0_X2(const double *_x)
{
    return 0.25 * (1 + _x[0]);
}

double Ip2d1(const double *_x)
{
    return 0.25 * (1 - _x[0]) * (1 + _x[1]);
}

double Ip2d1_X1(const double *_x)
{
    return -0.25 * (1 + _x[1]);
}

double Ip2d1_X2(const double *_x)
{
    // correct by huacheng
    return 0.25 * (1 - _x[0]);
}

double Ip2d2(const double *_x)
{
    return 0.25 * (1 - _x[0]) * (1 - _x[1]);
}

double Ip2d2_X1(const double *_x)
{
    return -0.25 * (1 - _x[1]);
}

double Ip2d2_X2(const double *_x)
{
    return -0.25 * (1 - _x[0]);
}


// functions used in 4 points interpolate
double Ip2d3(const double *_x)
{
    return 0.25 * (1 + _x[0]) * (1 - _x[1]);
}

double Ip2d3_X1(const double *_x)
{
    return 0.25 * (1 - _x[1]);
}

double Ip2d3_X2(const double *_x)
{
    return -0.25 * (1 + _x[0]);
}
// the functions used in 4 points for a boundary
// see the array IpB and IpB_X

// functions used in 8 points interpolate
// see the function array IpV and IpV_B
double Ip3d0(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 + _x[1]) * (1 - _x[2]);
}

double Ip3d1(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 + _x[1]) * (1 - _x[2]);
}

double Ip3d2(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 - _x[1]) * (1 - _x[2]);
}

double Ip3d3(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 - _x[1]) * (1 - _x[2]);
}

double Ip3d4(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 + _x[1]) * (1 + _x[2]);
}

double Ip3d5(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 + _x[1]) * (1 + _x[2]);
}

double Ip3d6(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 - _x[1]) * (1 + _x[2]);
}

double Ip3d7(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 - _x[1]) * (1 + _x[2]);
}

double Ip3d0_X1(const double *_x)
{
    return 0.125 * (1 + _x[1]) * (1 - _x[2]);
}

double Ip3d1_X1(const double *_x)
{
    return -0.125 * (1 + _x[1]) * (1 - _x[2]);
}

double Ip3d2_X1(const double *_x)
{
    return -0.125 * (1 - _x[1]) * (1 - _x[2]);
}

double Ip3d3_X1(const double *_x)
{
    return 0.125 * (1 - _x[1]) * (1 - _x[2]);
}

double Ip3d4_X1(const double *_x)
{
    return 0.125 * (1 + _x[1]) * (1 + _x[2]);
}

double Ip3d5_X1(const double *_x)
{
    return -0.125 * (1 + _x[1]) * (1 + _x[2]);
}

double Ip3d6_X1(const double *_x)
{
    return -0.125 * (1 - _x[1]) * (1 + _x[2]);
}

double Ip3d7_X1(const double *_x)
{
    return 0.125 * (1 - _x[1]) * (1 + _x[2]);
}

double Ip3d0_X2(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 - _x[2]);
}

double Ip3d1_X2(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 - _x[2]);
}

double Ip3d2_X2(const double *_x)
{
    return -0.125 * (1 - _x[0]) * (1 - _x[2]);
}

double Ip3d3_X2(const double *_x)
{
    return -0.125 * (1 + _x[0]) * (1 - _x[2]);
}

double Ip3d4_X2(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 + _x[2]);
}

double Ip3d5_X2(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 + _x[2]);
}

double Ip3d6_X2(const double *_x)
{
    return -0.125 * (1 - _x[0]) * (1 + _x[2]);
}

double Ip3d7_X2(const double *_x)
{
    return -0.125 * (1 + _x[0]) * (1 + _x[2]);
}

double Ip3d0_X3(const double *_x)
{
    return -0.125 * (1 + _x[0]) * (1 + _x[1]);
}

double Ip3d1_X3(const double *_x)
{
    return -0.125 * (1 - _x[0]) * (1 + _x[1]);
}

double Ip3d2_X3(const double *_x)
{
    return -0.125 * (1 - _x[0]) * (1 - _x[1]);
}

double Ip3d3_X3(const double *_x)
{
    return -0.125 * (1 + _x[0]) * (1 - _x[1]);
}

double Ip3d4_X3(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 + _x[1]);
}

double Ip3d5_X3(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 + _x[1]);
}

double Ip3d6_X3(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 - _x[1]);
}

double Ip3d7_X3(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 - _x[1]);
}

Interpolation IpB[NIpB] = {Ip2d0, Ip2d1, Ip2d2, Ip2d3};
Interpolation IpB_X[NIpB][2] = {
        {Ip2d0_X1, Ip2d0_X2},
        {Ip2d1_X1, Ip2d1_X2},
        {Ip2d2_X1, Ip2d2_X2},
        {Ip2d3_X1, Ip2d3_X2},
};
Interpolation IpV[NIpV] = {Ip3d0, Ip3d1, Ip3d2, Ip3d3, Ip3d4, Ip3d5, Ip3d6, Ip3d7};
Interpolation IpV_X[NIpV][3] = {
        {Ip3d0_X1, Ip3d0_X2, Ip3d0_X3},
        {Ip3d1_X1, Ip3d1_X2, Ip3d1_X3},
        {Ip3d2_X1, Ip3d2_X2, Ip3d2_X3},
        {Ip3d3_X1, Ip3d3_X2, Ip3d3_X3},
        {Ip3d4_X1, Ip3d4_X2, Ip3d4_X3},
        {Ip3d5_X1, Ip3d5_X2, Ip3d5_X3},
        {Ip3d6_X1, Ip3d6_X2, Ip3d6_X3},
        {Ip3d7_X1, Ip3d7_X2, Ip3d7_X3},
};

double Max(double a, double b)
{
    if (a > b)
        return a;
    else
        return b;
}

double Min(double a, double b)
{
    if (a < b)
        return a;
    else
        return b;
}

double Sqrt_trunc(double x)
{
    if(x < 0.) x = 0.;
    return sqrt(x);
}

void Cross(const double *_a, const double *_b, double *_c)
{
    // _c = _a x _b
    _c[X] = _a[Y] * _b[Z] - _a[Z] * _b[Y];
    _c[Y] = _a[Z] * _b[X] - _a[X] * _b[Z];
    _c[Z] = _a[X] * _b[Y] - _a[Y] * _b[X];
}

void CrossAddition(const double *_a, const double *_b, double *_c)
{
    // _c += _a x _b
    _c[X] += _a[Y] * _b[Z] - _a[Z] * _b[Y];
    _c[Y] += _a[Z] * _b[X] - _a[X] * _b[Z];
    _c[Z] += _a[X] * _b[Y] - _a[Y] * _b[X];
}

double Dot(const double *_a, const double *_b)
{
    return _a[X] * _b[X] + _a[Y] * _b[Y] + _a[Z] * _b[Z];
}

void Addition(double *_a, const double *_b)
{
    // _a = _a + _b
    for (int k = 0; k < DIM; k++)
        _a[k] = _a[k] + _b[k];
}

//void ScalerAddition(double *_a, const double *_b, double _p)
//{
//    // _a = _a + p * _b
//    for (int k = 0; k < DIM; k++)
//        _a[k] += _p * _b[k];
//}
void ScalerMove(double *_a, const double *_b, double _p)
{
    // _a = p * _b
    for (int k = 0; k < DIM; k++)
        _a[k] = _p * _b[k];
}


void Copy(const double _a[], double _b[])
{
    // _b = _a
    for (int k = 0; k < DIM; k++)
        _b[k] = _a[k];
}

void Scaler(double _a[], double _p)
{
    for (int k = 0; k < DIM; k++)
        _a[k] *= _p;
}


void Zero(double _a[])
{
    for (int k = 0; k < DIM; k++)
        _a[k] = 0.0;
}

double Determinant(double tM[][DIM])
{
    // operator det of 3*3 matrix
    double tmp = tM[X][X] * (tM[Y][Y] * tM[Z][Z] - tM[Z][Y] * tM[Y][Z])
                 - tM[Y][X] * (tM[X][Y] * tM[Z][Z] - tM[X][Z] * tM[Z][Y])
                 + tM[Z][X] * (tM[X][Y] * tM[Y][Z] - tM[X][Z] * tM[Y][Y]);
    return tmp;
}

double MatrixTime(double _a[], const double _m[][DIM], const double _b[])
{
    // _a = _m * _b
    for (int k = 0; k < DIM; k++)
    {
        _a[k] = Dot(_m[k], _b);
    }
}



void Normalization(double p[])
{
    double tmp = 0;
    for (int k = 0; k < DIM; k++)
    {
        tmp += p[k] * p[k];
    }
    tmp = sqrt(tmp);

    for (int k = 0; k < DIM; k++)
    {
        p[k] = p[k] / tmp;
    }
}

void NormalizationEx(double p[], double Scale)
{
    double tmp = 0;
    for (int k = 0; k < DIM; k++)
    {
        p[k] *= Scale;
        tmp += p[k] * p[k];
    }
    tmp = sqrt(tmp);

    if(tmp > 1e-4)
    {
        for (int k = 0; k < DIM; k++)
        {
            p[k] = p[k] / tmp;
        }
    } else
    {
        p[Z] = 1.0;
        p[Y] = 0.0;
        p[X] = 0.0;
    }

}


double Length(const double p[])
{
    double tmp = 0.0;
    for (int k = 0; k < DIM; k++)
    {
        tmp += p[k] * p[k];
    }
    return sqrt(tmp);
}

double Contraction(const double _a[][DIM], const double _b[][DIM])
{
    double tmp = 0.0;
    for (int k = 0; k < DIM; k++)
        for (int j = 0; j < DIM; j++)
        {
            tmp += _a[k][j] * _b[k][j];
        }
    return tmp;
}

void ScalerAddition(double *_a, const double *_b, double _p)
{
    // _a = _a + p * _b
    for (int k = 0; k < DIM; k++)
        _a[k] += _p * _b[k];
}

void LinearOp(double *_a, const double * _b, const double * _c, double _pb, double _pc)
{
    for(int k=0;k<DIM;k++)
    {
        _a[k] = _pb*_b[k] + _pc*_c[k];
    }
}



double SecondInvariant(double tM[][DIM])
{
    return tM[X][Y] * tM[X][Y] + tM[Y][Z] * tM[Y][Z] + tM[Z][X] * tM[Z][X]
           - tM[X][X] * tM[Y][Y] - tM[Y][Y] * tM[Z][Z] - tM[Z][Z] * tM[X][X];
}

void ScalerAdditionEx(double *_a, const double *_b, double _p, int d)
{
    // _a = _a + p * _b
    // this is a extended version of ScalerAddition
    // ScalerAddition(a,b,p) = ScalerAdditionEx(a,b,p,DIM)
    for (int k = 0; k < d; k++)
        _a[k] += _p * _b[k];
}

void ScalerAdditionEx1(double *_a, const double *_b, double _p)
{
    for(int k=0;k<NPHI;k++)
        _a[k] += _p * _b[k];
}

void ScalerAdditionEx2(double *_a, const double *_b, double _p)
{
    for(int k=0;k<NPSI;k++)
        _a[k] += _p * _b[k];
}

void ScalerMoveEx(double *_a, const double *_b, double _p, int d)
{
    // _a = p * _b
    for (int k = 0; k < d; k++)
        _a[k] = _p * _b[k];
}

void ScalerEx(double _a[], double _p, int _d)
{
    for (int k = 0; k < _d; k++)
        _a[k] *= _p;
}

void ZeroEx(double _a[], int d)
{
    // a extended version of Zero
    for (int k = 0; k < d; k++)
        _a[k] = 0.0;
}

void CopyEx(const double _a[], double _b[], int _d)
{
    // _b = _a
    for (int k = 0; k < _d; k++)
        _b[k] = _a[k];
}

void WindEx(double limitL, double limitR, double _a[], int _d)
{
    for (int i = 0; i < _d; i++)
    {
        if (_a[i] < limitL)
            _a[i] = limitL;
        else if (_a[i] > limitR)
            _a[i] = limitR;
    }
}

double Wind(double x, double limitL, double limitR)
{
    if(x < limitL)
        return limitL;
    else if(x > limitR)
        return limitR;
    else
        return x;
}

//void DeriveArea(struct Boundary * _b)
//{
//    double tArea[3] = {0.0,0.0,0.0};
//    for(int k=0;k<NGI2d;k++)
//    {
//        double tA[3] = {0.0,0.0,0.0};
//        double tX1[3] = {0.0,0.0,0.0};
//        double tX2[3] = {0.0,0.0,0.0};
//        for(int i=0;i<VPB;i++)
//        {
//            ScalerAddition(tX1,_b->NeV[i]->Position,IpB_X[i][0](GIPS2d[k]));
//            ScalerAddition(tX2,_b->NeV[i]->Position,IpB_X[i][1](GIPS2d[k]));
//        }
//        Cross(tX1,tX2,tA);
//        ScalerAddition(tArea,tA,GIWS2d[k]);
//    }
//    Copy(tArea,_b->Area);
//}

double SEPE[NSEPE][VPE][DIM];

void ConstructCube(const double a[], const double b[], double c[][DIM])
{
    // construct a cube with point a and b
    double cMin[DIM];
    double cMax[DIM];
    for (int k = 0; k < DIM; k++)
    {
        if (a[k] >= b[k])
        {
            cMax[k] = a[k];
            cMin[k] = b[k];
        }
        else
        {
            cMax[k] = b[k];
            cMin[k] = a[k];
        }
    }
    c[BFR][X] = c[BFL][X] = c[TFR][X] = c[TFL][X] = cMax[X];
    c[BKR][X] = c[BKL][X] = c[TKR][X] = c[TKL][X] = cMin[X];

    c[BFR][Y] = c[BKR][Y] = c[TFR][Y] = c[TKR][Y] = cMax[Y];
    c[BFL][Y] = c[BKL][Y] = c[TFL][Y] = c[TKL][Y] = cMin[Y];

    c[BFR][Z] = c[BKR][Z] = c[BKL][Z] = c[BFL][Z] = cMin[Z];
    c[TFR][Z] = c[TKR][Z] = c[TKL][Z] = c[TFL][Z] = cMax[Z];
}

const double ELC[VPE][DIM] = { // the local coordinates of vertexes in an element (ELC)
        {1.0,   1.0,  -1.0},
        {-1.0,  1.0,  -1.0},
        {-1.0, -1.0,  -1.0},
        {1.0,  -1.0,  -1.0}, // BFR, BKR, BKL, BFL
        {1.0,   1.0,   1.0},
        {-1.0,  1.0,   1.0},
        {-1.0, -1.0,   1.0},
        {1.0,  -1.0,   1.0}  // TFR, TKR, TKL, TFL
};

const int PCVE[VPE] ={
        /*
         * Position commutation between Vertex and Element
         * the relative position for in indexes is different in Vertex.NeE and Element.NeV
         * PCVE bind Vertex.NeE to Element.NeV
         */
    TKL, TFL, TFR, TKR,
    BKL, BFL, BFR, BKR
};

void DeriveSEPE()
{
    double tb[DIM] = {0.0, 0.0, 0.0};
    for (int k = 0; k < NSEPE; k++)
    {
        ConstructCube(tb, ELC[k], SEPE[k]);
    }
//    ConstructCube(ELC[TFR],ELC[BKL],SEPE[VPE]);
}


double SBPB[NSBPB][VPB][DIM - 1];

void ConstructRect(const double a[], const double b[], double c[][DIM - 1])
{
    // construct a rectangle with point a and b
    double rMin[DIM];
    double rMax[DIM];
    for (int k = 0; k < DIM - 1; k++)
    {
        if (a[k] > b[k])
        {
            rMax[k] = a[k];
            rMin[k] = b[k];
        }
        else
        {
            rMax[k] = b[k];
            rMin[k] = a[k];
        }
    }
    c[0][X] = c[3][X] = rMin[X];
    c[1][X] = c[2][X] = rMax[X];

    c[0][Y] = c[1][Y] = rMin[Y];
    c[2][Y] = c[3][Y] = rMax[Y];
}

const double BLC[VPB][DIM - 1] = {
        {-1, -1},
        { 1, -1},
        { 1,  1},
        {-1,  1}
};

void DeriveSBPB()
{
    double tb[2] = {0.0, 0.0};
    for (int k = 0; k < VPB; k++)
    {
        ConstructRect(tb, BLC[k], SBPB[k]);
    }
}

void XgIpV(double Xg[], double Xi[][DIM], const double xl[])
{
    // Xg is the coordinate global
    // Xi is the Vertexes of Volume
    // xl is the local coordinate
    double tXg[DIM] = {0.0, 0.0, 0.0};
    for (int i = 0; i < NIpV; i++)
    {
        ScalerAddition(tXg, Xi[i], IpV[i](xl));
    }
    Copy(tXg, Xg);
}

void XgIpV_X(double Xg[], double Xi[][DIM], const double xl[], const int _d)
{
    // _d is in {X,Y,Z}
    double tXg[DIM] = {0.0, 0.0, 0.0};
    for (int i = 0; i < NIpV; i++)
    {
        ScalerAddition(tXg, Xi[i], IpV_X[i][_d](xl));
    }
    Copy(tXg, Xg);
}

double detJV(double Xi[][DIM], const double xl[])
{
    // the Determinant of jacobi matrix
    // Xi is the Vertexes of Volume
    // xl is the local coordinate
    // tM is the transpose of the jacobi matrix
    double tM[DIM][DIM] = {0.0};
    for (int k = 0; k < DIM; k++)
    {
        XgIpV_X(tM[k], Xi, xl, k);
    }
    return Determinant(tM);
}

double DeriveVolume(double Xi[][DIM])
{
    double tV = 0.0;
    for (int k = 0; k < NGI3d; k++)
    {
        tV += GIWS3d[k] * detJV(Xi, GIPS3d[k]);
    }
    return tV;
}

void XgIpB(double Xg[], double Xi[][DIM], const double xl[])
{
    // Xg is the coordinate global
    // Xi is the Vertexes of Volume
    // xl is the local coordinate
    Zero(Xg);
    for (int i = 0; i < NIpB; i++)
    {
        ScalerAddition(Xg, Xi[i], IpB[i](xl));
    }
}

void XgIpB_X(double Xg[], double Xi[][DIM], const double xl[], const int _d)
{
    Zero(Xg);
    for (int i = 0; i < NIpB; i++)
    {
        ScalerAddition(Xg, Xi[i], IpB_X[i][_d](xl));
    }
}

//void DriveArea(const double p[][DIM], double a[])
//{
//    double tArea[DIM] = {0.0,0.0,0.0};
//    for(int k=0;k<NGI2d;k++)
//    {
//        double tA[3] = {0.0,0.0,0.0};
//        double tX1[3] = {0.0,0.0,0.0};
//        double tX2[3] = {0.0,0.0,0.0};
//        for(int i=0;i<NIpB;i++)
//        {
//            ScalerAddition(tX1,p[i],IpB_X[i][0](GIPS2d[k]));
//            ScalerAddition(tX2,p[i],IpB_X[i][1](GIPS2d[k]));
//        }
//        Cross(tX1,tX2,tA);
//        ScalerAddition(tArea,tA,GIWS2d[k]);
//    }
//    Copy(tArea,a);
//}

void DeriveArea(double Xi[][DIM], double a[])
{
    // calculate the norm vector of a face
    double tArea[DIM] = {0.0, 0.0, 0.0};
    for (int k = 0; k < NGI2d; k++)
    {
        double tA[3] = {0.0, 0.0, 0.0};
        double tX1[3] = {0.0, 0.0, 0.0};
        double tX2[3] = {0.0, 0.0, 0.0};

        XgIpB_X(tX1, Xi, GIPS2d[k], X);
        XgIpB_X(tX2, Xi, GIPS2d[k], Y);

        Cross(tX1, tX2, tA);
        ScalerAddition(tArea, tA, GIWS2d[k]);
    }
    Copy(tArea, a);
}

void XgSubElement(double Xg[][DIM], double Xi[][DIM], double xl[][DIM])
{
    for (int k = 0; k < VPE; k++)
        XgIpV(Xg[k], Xi, xl[k]);
}

void XgSubBoundary(double Xg[][DIM], double Xi[][DIM], double xl[][DIM - 1])
{
    for (int k = 0; k < NSBPB; k++)
        XgIpB(Xg[k], Xi, xl[k]);
}

void XgInterFace(double Xg[][DIM], double Xi[][DIM], double xl[][DIM])
{
    for (int k = 0; k < VPB; k++)
        XgIpV(Xg[k], Xi, xl[k]);
}

void DeriveSubVolume(double Xi[][DIM], double Vol[])
{
    double tXg[VPE][DIM] = {0.0};
    for (int k = 0; k < NSEPE; k++)
    {
        XgSubElement(tXg, Xi, SEPE[k]);
        Vol[k] = DeriveVolume(tXg);
    }
}

void DeriveSubBoundary(double Xi[][DIM], double SubBound[][DIM])
{
    double tXg[VPB][DIM] = {0.0};
    for (int k = 0; k < NSBPB; k++)
    {
        XgSubBoundary(tXg, Xi, SBPB[k]);
        DeriveArea(tXg, SubBound[k]);
    }
}

double ISFE[DIM][VPB][DIM] = {
        {{0.0,  -1.0, -1.0}, {0.0,  1.0,  -1.0}, {0.0, 1.0, 1.0}, {0.0,  -1.0, 1.0}},
        {{-1.0, 0.0,  -1.0}, {-1.0, 0.0,  1.0},  {1.0, 0.0, 1.0}, {1.0,  0.0,  -1.0}},
        {{-1.0, -1.0, 0.0},  {1.0,  -1.0, 0.0},  {1.0, 1.0, 0.0}, {-1.0, 1.0,  0.0}}
};

void DeriveInterFace(double Xi[][DIM], double InterFace[][NSBPB][DIM])
{
    for (int k = 0; k < DIM; k++)
    {
        double tXg[VPB][DIM];
        XgInterFace(tXg, Xi, ISFE[k]);
        DeriveSubBoundary(tXg, InterFace[k]);
    }
}

const int ToSubFace0[VPE][DIM] = {
        {1, 3, 2},
        {1, 0, 3},
        {0, 0, 0},
        {0, 3, 1},
        {2, 2, 2},
        {2, 1, 3},
        {3, 1, 0},
        {3, 2, 1}};

const double ToSubFace1[VPE][DIM] = {
        {1.0,  1.0,  -1.0},
        {-1.0, 1.0,  -1.0},
        {-1.0, -1.0, -1.0},
        {1.0,  -1.0, -1.0},
        {1.0,  1.0,  1.0},
        {-1.0, 1.0,  1.0},
        {-1.0, -1.0, 1.0},
        {1.0,  -1.0, 1.0}};

/*
 * All the coefficients in ToSubFace need to be opposite to the order
 */
void DeriveSubFace(double InterFace[][NSBPB][DIM], double SubFace[][DIM])
{
    for (int i = 0; i < VPE; i++)
    {
        double tA[DIM] = {0.0, 0.0, 0.0};
        for (int j = 0; j < DIM; j++)
        {
            int idsub = ToSubFace0[i][j];
            ScalerAddition(tA, InterFace[j][idsub], -1.0 * ToSubFace1[i][j]);
            // adjust according the definition of area vector
            // old version is " ScalerAddition(tA, InterFace[j][idsub],ToSubFace1[i][j]);
        }
        Copy(tA, SubFace[i]);
    }
}

void XiCopy(double Xi[][DIM], double tXi[][DIM], int Length)
{
    for (int i = 0; i < Length; i++)
    {
        Copy(Xi[i], tXi[i]);
    }
}

void SetEi(double Xi[][DIM], int _d, int _gV)
{
    for (int i = 0; i < VPE; i++)
    {
        Xi[i][_d] = 0.0;
    }
    Xi[_gV][_d] = 1.0;
}

void DeriveGradFactor(double Xi[][DIM], double GradFactor[][VPE])
{
    for (int k = 0; k < DIM; k++)
    {
        double tXi[VPE][DIM];
        XiCopy(Xi, tXi, VPE);
        for (int j = 0; j < VPE; j++)
        {
            SetEi(tXi, k, j);
            GradFactor[k][j] = DeriveVolume(tXi);
        }
    }
}

const int ToOuterFace[BPE][VPB] = {
        {2, 3, 7, 6}, /* Left */
        {0, 1, 5, 4}, /* R    */
        {0, 4, 7, 3}, /* F    */
        {1, 2, 6, 5}, /* Back */
        {0, 3, 2, 1}, /* Bottom */
        {4, 5, 6, 7}  /* T */
};

void XgOuterFace(double Xg[][DIM], double Xi[][DIM], const int ids[])
{
    for (int k = 0; k < VPB; k++)
    {
        int Id = ids[k];
        Copy(Xi[Id], Xg[k]);
    }
}

void DeriveOuterFace(double Xi[][DIM], double InterFace[][DIM])
{
    // calculate the normal of every face of a Hexahedron
    for (int k = 0; k < BPE; k++)
    {
        double tXg[VPB][DIM];
        XgOuterFace(tXg, Xi, ToOuterFace[k]);
        DeriveArea(tXg, InterFace[k]);
    }
}

double DeriveArea2(double Si[][DIM], double Xi[][DIM])
{
    double tFlux = 0.0;
    for (int k = 0; k < NGI2d; k++)
    {
        double tA[3] = {0.0, 0.0, 0.0};
        double tX1[3] = {0.0, 0.0, 0.0};
        double tX2[3] = {0.0, 0.0, 0.0};
        double tS[3] = {0.0, 0.0, 0.0};

        XgIpB_X(tX1, Xi, GIPS2d[k], X);
        XgIpB_X(tX2, Xi, GIPS2d[k], Y);
        XgIpB(tS, Si, GIPS2d[k]);
        Cross(tX1, tX2, tA);
        tFlux += GIWS2d[k] * Dot(tS, tA);
    }
    return tFlux;
}

void SetSi(double Si[][DIM], int _a, int _b)
{
    for (int k = 0; k < DIM; k++)
        for (int j = 0; j < VPB; j++)
        {
            Si[j][k] = 0.0;
        }
    Si[_a][_b] = 1.0;
}

void DeriveFlowFactor(double Xi[][DIM], double FlowFactor[][DIM])
{
    double tSi[VPB][DIM];
    for (int k = 0; k < DIM; k++)
    {
        for (int j = 0; j < VPB; j++)
        {
            SetSi(tSi, j, k);
            FlowFactor[j][k] = DeriveArea2(tSi, Xi);
        }
    }
}

double DeriveCenter(double Xi[][DIM], double _c[])
{
    Zero(_c);
    for(int i=0;i<VPE;i++)
        ScalerAddition(_c,Xi[i],0.125);

    double rdx_min = fabs(Xi[0][Z] - _c[Z]);
    for(int i=0;i<VPE;i++)
    {
        double t1 = Min(Min(fabs(Xi[i][X] - _c[X]),fabs(Xi[i][Y] - _c[Y])),fabs(Xi[i][Z] - _c[Z]));
        rdx_min = Min(t1,rdx_min);
    }
    if(rdx_min > 0.0)
        return 0.5/rdx_min;
    else
    {
        fprintf(stdout,"Coordinate of Element Center is not set correctly!\n");
        exit(0);
    }
}

double DeriveScale(double Xi[][DIM], double _c[])
{
    Zero(_c);
    for(int k=0;k<VPE;k++)
    {
        for(int j=k+1;j<VPE;j++)
        {
            for(int i=0;i<DIM;i++)
            {
                _c[i] = Max(_c[i], fabs(Xi[j][i] - Xi[k][i]));
            }
        }
    }
}

double MaxScaler(const double Xi[], double p)
{
    double tx1 = fabs(Xi[0] * p);
    double tx2 = fabs(Xi[1] * p);
    double tx3 = fabs(Xi[2] * p);
    if(tx1 > tx2 && tx1 > tx3)
    {
        return tx1;
    }
    else if(tx2 > tx3)
    {
        return tx2;
    }
    else
        return tx3;
}

double MaxComponent(const double Xi[], int d)
{
    double ResVal = fabs(Xi[0]);
    for(int k=1;k<d;k++)
    {
        ResVal = Max(ResVal,fabs(Xi[k]));
    }
    return ResVal;
}

int GetBit(int x,int pos)
{
    return (x >> pos) & 1;
}

void SetBit(int * x,int pos)
{
    *x |= (1<<pos);
}

double MinMod(double a, double b)
{
    if(a*b < 0.0)
        return 0.0;
    else if(a > 0.0)
        return Min(fabs(a),fabs(b));
    else
        return -1.0*Min(fabs(a),fabs(b));
}

double Distance(const double x1[], const double x2[])
{
    double tmp[DIM];
    for(int k=0;k<DIM;k++)
        tmp[k] = x1[k] - x2[k];
    return Length(tmp);
}

double Sign(double x)
{
    if(x >= 0)
        return 1.0;
    else
        return -1.0;
}

int sIndex[9] = {XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ};

int SymIndex(int i,int j)
{
    return sIndex[i + 3*j];
}

double MatrixTimeSym(double _a[], const double _m[], const double _b[])
{
    // _a = _m * _b
    // _m is a Symmetrical 3*3 matrix

    _a[X] = _m[XX]*_b[X] + _m[XY]*_b[Y] + _m[XZ]*_b[Z];
    _a[Y] = _m[YX]*_b[X] + _m[YY]*_b[Y] + _m[YZ]*_b[Z];
    _a[Z] = _m[ZX]*_b[X] + _m[ZY]*_b[Y] + _m[ZZ]*_b[Z];
}

double ContractionSym(const double _a[], double _b[][DIM])
{
    double tmp = 0.0;
    for (int k = 0; k < DIM; k++)
        for (int j = 0; j < DIM; j++)
        {
            tmp += _a[SymIndex(k,j)] * _b[k][j];
        }
    return tmp;
}

double SecondInvariantSym(const double tM[])
{
    return tM[XY]*tM[XY] + tM[YZ] * tM[YZ] + tM[ZX] * tM[ZX]
           - tM[XX] * tM[YY] - tM[YY] * tM[ZZ] - tM[ZZ] * tM[XX];
}

double DeterminantSym(const double tM[])
{
    // operator det of 3*3 matrix
    double tmp =   tM[XX] * (tM[YY] * tM[ZZ] - tM[ZY] * tM[YZ])
                 - tM[YX] * (tM[XY] * tM[ZZ] - tM[XZ] * tM[ZY])
                 + tM[ZX] * (tM[XY] * tM[YZ] - tM[XZ] * tM[YY]);
    return tmp;
}

int SumBit( int n )
{
    n = (n&0x55555555) + ((n>>1)&0x55555555);
    n = (n&0x33333333) + ((n>>2)&0x33333333);
    n = (n&0x0f0f0f0f) + ((n>>4)&0x0f0f0f0f);
    n = (n&0x00ff00ff) + ((n>>8)&0x00ff00ff);
    n = (n&0x0000ffff) + ((n>>16)&0x0000ffff);
    return n;
}

int ReverseBit(int n)
{
    int k = n;
    int r = 0;
    while( n > 1)
    {
        n = n/2;
        r++;
    }
    return r;
}


int ToLongId(int Idx, int Idy, int Idz, int nxp, int nyp, int nzp)
{
    if(Idx >= nxp || Idy >= nyp || Idz >= nzp)
        return -1;
    if (Idx < 0 || Idy < 0 || Idz < 0)
        return -1;

    int LongId = Idx + Idy * nxp + Idz * nxp * nyp;
    return LongId;
}

int ToShortId(int * Idx, int * Idy, int * Idz,int LongId, int nxp, int nyp, int nzp)
{
    *Idx = LongId%nxp;
    *Idy = (LongId/nxp)%nyp;
    *Idz = LongId/(nxp*nyp);

    if(*Idx >= nxp || *Idy >= nyp || *Idz >= nzp)
        return -1;
    else if (*Idx < 0 || *Idy < 0 || *Idz < 0)
        return -1;
    else
        return 0;
}

