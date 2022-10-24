//
// Created by huacheng on 2020/12/22.
//

#include "Setup.h"

int GetMeshGeometry(struct Mesh *_mesh, const struct MeshInfo *_info)
{
    if (NULL == _mesh)
    {
        NullPointer(stdout, "# the struct Mesh is NULL!\n");
        return 1;
    }
    // get/set (npx,npy,npz) of the mesh
    // the allocate the memory for mesh

    _mesh->nxp = _info->npx;
    _mesh->nyp = _info->npy;
    _mesh->nzp = _info->npz;
    _mesh->nvof = _info->MaterialNum;

    _mesh->ne = (_mesh->nxp - 1) * (_mesh->nyp - 1) * (_mesh->nzp - 1);
    _mesh->nv = _mesh->nxp * _mesh->nyp * _mesh->nzp;
    _mesh->nbx = _mesh->nxp * (_mesh->nyp - 1) * (_mesh->nzp - 1);
    _mesh->nby = (_mesh->nxp - 1) * _mesh->nyp * (_mesh->nzp - 1);
    _mesh->nbz = (_mesh->nxp - 1) * (_mesh->nyp - 1) * _mesh->nzp;

    _mesh->Vertexes = (struct Vertex *) malloc(_mesh->nv * sizeof(struct Vertex));
    _mesh->Elements = (struct Element *) malloc(_mesh->ne * sizeof(struct Element));
    _mesh->XBoundaries = (struct Boundary *) malloc(_mesh->nbx * sizeof(struct Boundary));
    _mesh->YBoundaries = (struct Boundary *) malloc(_mesh->nby * sizeof(struct Boundary));
    _mesh->ZBoundaries = (struct Boundary *) malloc(_mesh->nbz * sizeof(struct Boundary));

    if (NULL == _mesh->Vertexes || NULL == _mesh->Elements || NULL == _mesh->XBoundaries ||
        NULL == _mesh->YBoundaries || NULL == _mesh->ZBoundaries)
    {
        NullPointer(stdout, "# Fail to allocate the array of Vertex or Element\n");
    }

    if((0==strcmp("binary",_info->output_format))||(0==strcmp("ascii",_info->output_format)))
        strcpy(_mesh->vtk_format,_info->output_format);
    else
        strcpy(_mesh->vtk_format,"ascii");

    return 0;
}

int BuildMeshGeometry(struct Mesh *_mesh, const struct MeshInfo *_info)
{
    #pragma omp parallel for num_threads(NTHREAD)
    for (int i = 0; i < _mesh->nxp; i++)
        for (int j = 0; j < _mesh->nyp; j++)
            for (int k = 0; k < _mesh->nzp; k++)
            {
                // Search the neighbors of every boundary in X direction
                struct Boundary *ibx = GetXBoundaries(_mesh, i, j, k);
                if (NULL != ibx)
                {
                    ibx->NeE[0] = GetElement(_mesh, i - 1, j, k);
                    ibx->NeE[1] = GetElement(_mesh, i    , j, k);

                    ibx->NeV[0] = GetVertex(_mesh, i, j    , k    );
                    ibx->NeV[1] = GetVertex(_mesh, i, j + 1, k    );
                    ibx->NeV[2] = GetVertex(_mesh, i, j + 1, k + 1);
                    ibx->NeV[3] = GetVertex(_mesh, i, j    , k + 1);

                    /*
                     * Search the upwind element required by VOF algorithm
                     * Those pointer can also be used to determined the element that on
                     * boundaries.
                     */

                    ibx->UpE[0] = GetElement(_mesh, i - 2, j, k);
                    ibx->UpE[1] = GetElement(_mesh, i + 1, j, k);
                }

//                RunPoint(stdout, "Finish Search the neighbors of every boundary in X direction\n");

                // Search the neighbor of every boundary in Y direction
                struct Boundary *iby = GetYBoundaries(_mesh, i, j, k);
                if (NULL != iby)
                {
                    iby->NeE[0] = GetElement(_mesh, i, j - 1, k);
                    iby->NeE[1] = GetElement(_mesh, i, j    , k);

                    iby->NeV[0] = GetVertex(_mesh, i    , j, k);
                    iby->NeV[1] = GetVertex(_mesh, i    , j, k + 1);
                    iby->NeV[2] = GetVertex(_mesh, i + 1, j, k + 1);
                    iby->NeV[3] = GetVertex(_mesh, i + 1, j, k);

                    iby->UpE[0] = GetElement(_mesh, i, j - 2, k);
                    iby->UpE[1] = GetElement(_mesh, i, j + 1, k);
                }
//                RunPoint(stdout, "Finish Search the neighbors of every boundary in Y direction\n");

                // Search the neighbor of every boundary in Z direction
                struct Boundary *ibz = GetZBoundaries(_mesh, i, j, k);
                if (NULL != ibz)
                {
                    ibz->NeE[0] = GetElement(_mesh, i, j, k - 1);
                    ibz->NeE[1] = GetElement(_mesh, i, j, k);

                    ibz->NeV[0] = GetVertex(_mesh, i    , j    , k);
                    ibz->NeV[1] = GetVertex(_mesh, i + 1, j    , k);
                    ibz->NeV[2] = GetVertex(_mesh, i + 1, j + 1, k);
                    ibz->NeV[3] = GetVertex(_mesh, i    , j + 1, k);

                    ibz->UpE[0] = GetElement(_mesh, i, j, k - 2);
                    ibz->UpE[1] = GetElement(_mesh, i, j, k + 1);
                }
//                RunPoint(stdout, "Finish Search the neighbors of every boundary in Z direction\n");

                // Search the neighbor of every element
                struct Element *ie = GetElement(_mesh, i, j, k);
                if (NULL != ie)
                {
                    ie->NeV[BFR] = GetVertex(_mesh, i + 1, j + 1, k);
                    ie->NeV[BKR] = GetVertex(_mesh, i    , j + 1, k);
                    ie->NeV[BKL] = GetVertex(_mesh, i    , j    , k);
                    ie->NeV[BFL] = GetVertex(_mesh, i + 1, j    , k);
                    ie->NeV[TFR] = GetVertex(_mesh, i + 1, j + 1, k + 1);
                    ie->NeV[TKR] = GetVertex(_mesh, i    , j + 1, k + 1);
                    ie->NeV[TKL] = GetVertex(_mesh, i    , j    , k + 1);
                    ie->NeV[TFL] = GetVertex(_mesh, i + 1, j    , k + 1);

                    ie->NeB[Top]    = GetZBoundaries(_mesh, i    , j    , k + 1);
                    ie->NeB[Bottom] = GetZBoundaries(_mesh, i    , j    , k);
                    ie->NeB[Left]   = GetYBoundaries(_mesh, i    , j    , k);
                    ie->NeB[Right]  = GetYBoundaries(_mesh, i    , j + 1, k);
                    ie->NeB[Back]   = GetXBoundaries(_mesh, i    , j    , k);
                    ie->NeB[Front]  = GetXBoundaries(_mesh, i + 1, j    , k);

                }
//                RunPoint(stdout, "Finish Search the neighbor of every element\n");

                // Search the neighbor of every vertex
                struct Vertex *iv = GetVertex(_mesh, i, j, k);
                if (NULL != iv)
                {
                    Copy(_info->V[k][j][i], iv->Position);

                    iv->NeE[4] = GetElement(_mesh, i    , j    , k);
                    iv->NeE[5] = GetElement(_mesh, i - 1, j    , k);
                    iv->NeE[6] = GetElement(_mesh, i - 1, j - 1, k);
                    iv->NeE[7] = GetElement(_mesh, i    , j - 1, k);

                    iv->NeE[0] = GetElement(_mesh, i    , j    , k - 1);
                    iv->NeE[1] = GetElement(_mesh, i - 1, j    , k - 1);
                    iv->NeE[2] = GetElement(_mesh, i - 1, j - 1, k - 1);
                    iv->NeE[3] = GetElement(_mesh, i    , j - 1, k - 1);
                }
//                RunPoint(stdout, "Finish Search the neighbor of every vertex\n");
            }
    return 0;
}

int Calculate_SV(struct Mesh *_mesh)
{
    /*
     * Calculate the sub-volume in any element
     * Those parameter will be used in Map between Vertex and Element.
     * The factors of volume that gets through the outerfaces for any given displacement of vertex
     * are also calculated.
     */

    DeriveSEPE();
    DeriveSBPB();

    #pragma omp parallel for num_threads(NTHREAD)
    for(int i=0;i<_mesh->ne;i++)
    {
        struct Element *ie = &_mesh->Elements[i];
        if(NULL==ie) continue;

        double tXi[VPE][DIM];
        for(int m = 0; m < VPE; m++)
        {
            // no need for check NeV == NULL
            Copy(ie->NeV[m]->Position, tXi[m]);
        }

        ie->rdx_min = DeriveCenter(tXi,ie->Center);
        DeriveScale(tXi,ie->Scale);
        double TmpSubVolume[NSEPE];
        DeriveSubVolume(tXi, TmpSubVolume);

        double TmpInterSubFace[DIM][NSBPB][DIM];
        // ie->InterSubFace is removed. a Tmp is set here
        DeriveInterFace(tXi, TmpInterSubFace);

        /*
         * SubFace is the area vector for every NeV
         * Used in update the velocity
         */
        DeriveSubFace(TmpInterSubFace, ie->SubFace);

        /*
         * Other geometry parameters determined by sub-volume
         */
        ie->Volume = 0.0;
        for (int m = 0; m < NSEPE; m++)
        {
            ie->Volume += TmpSubVolume[m];
        }
        ie->rVolume = 1.0 / ie->Volume;

        for (int m = 0; m < NSEPE; m++)
        {
            ie->SubVolumeRatio[m] = TmpSubVolume[m] * ie->rVolume;
        }

        /*
         * NOTICE: DeriveGradFactor give the factor integral with in the element
         * the grad should be Integral(Grad(V),Element)/Volume(Element)
         */
        DeriveGradFactor(tXi, ie->GradFactor);
        for(int fi=0;fi<DIM;fi++)
            for(int fj=0;fj<VPE;fj++)
            {
                ie->GradFactor[fi][fj] *= ie->rVolume;
            }

        // calculated the vector vec(S) on the outer face
        // ie->OuterFace is removed and instead by TmpOuterFace
        double TmpOuterFace[BPE][DIM];
        DeriveOuterFace(tXi, TmpOuterFace);
        ie->ABVArea = 0.0;
        for (int fk = 0; fk < DIM; fk++)
        {
            for (int fi = 0; fi < BPE; fi++)
                ie->ABVArea += fabs(TmpOuterFace[fi][fk]);
        }
        /*
         *  the norm of the vec(S) is also used in the inner boundary reconstruction
         *  so we calculated and store in the element
         *  the value of vec(S) is calculated repeated 2 times in the setup subroutine !!!
         */

    }

    // The factors of volume that gets through the outer-faces
    #pragma omp parallel for num_threads(NTHREAD)
    for(int i=0;i<_mesh->nbx;i++)
    {
        struct Boundary *ibx = &_mesh->XBoundaries[i];
        if (NULL != ibx) Calculate_SVboundary(ibx);
    }

    #pragma omp parallel for num_threads(NTHREAD)
    for(int i=0;i<_mesh->nby;i++)
    {
        struct Boundary *iby = &_mesh->YBoundaries[i];
        if (NULL != iby) Calculate_SVboundary(iby);
    }

    #pragma omp parallel for num_threads(NTHREAD)
    for(int i=0;i<_mesh->nbz;i++)
    {
        struct Boundary *ibz = &_mesh->ZBoundaries[i];
        if (NULL != ibz) Calculate_SVboundary(ibz);
    }

    /*
     * Derive the geometry parameters with in Vertex
     */
    #pragma omp parallel for num_threads(NTHREAD)
    for(int i=0;i<_mesh->nv;i++)
    {
        struct Vertex * iv = &_mesh->Vertexes[i];
        iv->Volume = 0.0;
        for(int j=0;j<NSEPE;j++)
        {
            if(NULL == iv->NeE[j])
            {
                iv->SubVolumeRatio[j] = 0.0;
            }
            else
            {
                iv->Volume += iv->NeE[j]->SubVolumeRatio[PCVE[j]] * iv->NeE[j]->Volume;
                iv->SubVolumeRatio[j] = iv->NeE[j]->SubVolumeRatio[PCVE[j]] * iv->NeE[j]->Volume;
            }
        }
        iv->rVolume = 1.0/iv->Volume;
        for(int j=0;j<NSEPE;j++)
            iv->SubVolumeRatio[j] *= iv->rVolume;
    }
    return 0;
}

void Calculate_SVboundary(struct Boundary *_boundary)
{
    double tXi[VPB][DIM];
    for (int k = 0; k < VPB; k++)
    {
        Copy(_boundary->NeV[k]->Position, tXi[k]);
    }
    DeriveArea(tXi, _boundary->Normal);
    Normalization(_boundary->Normal);
    DeriveFlowFactor(tXi, _boundary->FluxFactor);
}

void SetInitVertex_gas(struct Mesh *_mesh, InitConditionV _icf, int _nm)
{
    for(int k=0;k<_mesh->nv;k++)
    {
        struct Vertex * iv = &_mesh->Vertexes[k];
        if(NULL != iv) _icf(iv);
    }
}

void SetInitElement(struct Mesh * _mesh, InitConditionE _icf)
{
    #pragma omp parallel for num_threads(NTHREAD)
    for(int i=0;i<_mesh->ne;i++)
    {
        struct Element * ie = &_mesh->Elements[i];
        if(NULL == ie) continue;
        _icf(ie);
    }
}

void SetVelStavlity(struct Mesh *_mesh, InitConditionV _icf)
{
    for(int k=0;k<_mesh->ndv;k++)
    {
        struct Vertex * iv = _mesh->DomainV[k];
        if(NULL == iv) continue;
        _icf(iv);
    }
}


void SearchDomain(struct Mesh *_mesh, void (*_cef)(int, struct ConditionE *))
{
    /*
     * Search the calculate domain of vertexes
     */

    _mesh->DomainV = (struct Vertex **) malloc(_mesh->nv * sizeof(struct Vertex *));
    _mesh->ndv = 0;

    /*
     * Don't make parallel, the ndv is important for the domain
     * make it critical may improve performance.
     */
    for(int i=0;i<_mesh->nv;i++)
    {

        int InDomain = 0x0000;
        for(int j=0;j<EPV;j++)
        {
            if(NULL != _mesh->Vertexes[i].NeE[j])
            {
                InDomain |= 1<<j;
            }
        }
        if(InDomain == 0x00ff)
            _mesh->DomainV[_mesh->ndv++] = &_mesh->Vertexes[i];
        else if(0x0033 == InDomain)
        {
            // the left-boundary
        }
        else if(0x00cc == InDomain)
        {
            // the right-boundary
        }
        else if(0x000f == InDomain)
        {
            // the top-boundary
        }
        else if(0x00f0 == InDomain)
        {
            // the bottom-boundary
        }
        else if(0x0066 == InDomain)
        {
            // the back-boundary
        }
        else if(0x0099 == InDomain)
        {
            // the front-boundary
        }
    }

    /*
     * Search the domain of boundary[X,Y,Z]
     */
    int dbSize = _mesh->nbx + _mesh->nby + _mesh->nbz;
    _mesh->DomainB = (struct Boundary **) malloc(dbSize* sizeof(struct Boundary *));
    _mesh->ndb = 0;
    for(int i=0;i<_mesh->nbx;i++)
    {
        struct Boundary * ib = &_mesh->XBoundaries[i];
        if(NULL==ib->NeE[0] || NULL==ib->NeE[1]) continue;
        _mesh->DomainB[_mesh->ndb++] = ib;
    }
    for(int i=0;i<_mesh->nby;i++)
    {
        struct Boundary * ib = &_mesh->YBoundaries[i];
        if(NULL==ib->NeE[0] || NULL==ib->NeE[1]) continue;
        _mesh->DomainB[_mesh->ndb++] = ib;
    }
    for(int i=0;i<_mesh->nbz;i++)
    {
        struct Boundary * ib = &_mesh->ZBoundaries[i];
        if(NULL==ib->NeE[0] || NULL==ib->NeE[1]) continue;
        _mesh->DomainB[_mesh->ndb++] = ib;
    }

    /*
     * adjust the upwind elements of boundaries that near the boundary-condition;
     * this operation can be moved to BuildGeometry if necessary.
     */
    for(int i=0;i<_mesh->ndb;i++)
    {
        struct Boundary * ib = _mesh->DomainB[i];
        for(int j=0;j<2;j++)
        {
            if(NULL==ib->UpE[j])
                ib->UpE[j] = ib->NeE[j];
        }
    }

    /*
     * Add DomainE and ghost cell in _mesh struct
     */
    _mesh->DomainE = (struct Element **) malloc(_mesh->ne * sizeof(struct Element *));
    _mesh->CondE   = (struct ConditionE *) malloc(2*(_mesh->nxp*_mesh->nyp + _mesh->nyp*_mesh->nzp + _mesh->nzp*_mesh->nxp)*
                                                    sizeof(struct ConditionE));
    _mesh->nde = _mesh->nce = 0;
    for(int i=0;i<_mesh->ne;i++)
    {
        struct Element * ie = &_mesh->Elements[i];
        int InDomain = 0;
        for(int j=0;j<BPE;j++)
        {
            if(NULL != ie->NeB[j]->NeE[0] && NULL != ie->NeB[j]->NeE[1])
                InDomain |= 1<<j;
        }
        if(63 == InDomain)
        {
            _mesh->DomainE[_mesh->nde++] = ie;
        }
        else
        {
            struct ConditionE * tCE = &_mesh->CondE[_mesh->nce];
            tCE->dest   = ie;
            tCE->refe   = ie;
            tCE->CEType = 0;

            if(0 == GetBit(InDomain,Left))
            {
                _cef(Left,tCE);
                struct Element * tmpE = tCE->refe;
                tCE->refe = tmpE->NeB[Right]->NeE[1];
                tCE->CEType ++;
            }
            else if(0 == GetBit(InDomain,Right))
            {
                _cef(Right,tCE);
                struct Element * tmpE = tCE->refe;
                tCE->refe = tmpE->NeB[Left]->NeE[0];
                tCE->CEType ++;
            }

            if(0 == GetBit(InDomain,Top))
            {
                _cef(Top,tCE);
                struct Element * tmpE = tCE->refe;
                tCE->refe = tmpE->NeB[Bottom]->NeE[0];
                tCE->CEType ++;
            }
            else if(0 == GetBit(InDomain,Bottom))
            {
                _cef(Bottom,tCE);
                struct Element * tmpE = tCE->refe;
                tCE->refe = tmpE->NeB[Top]->NeE[1];
                tCE->CEType ++;
            }

            if(0 == GetBit(InDomain,Front))
            {
                _cef(Front,tCE);
                struct Element * tmpE = tCE->refe;
                tCE->refe = tmpE->NeB[Back]->NeE[0];
                tCE->CEType ++;
            }
            else if(0 == GetBit(InDomain,Back))
            {
                _cef(Back,tCE);
                struct Element * tmpE = tCE->refe;
                tCE->refe = tmpE->NeB[Front]->NeE[1];
                tCE->CEType ++;
            }

            if(tCE->CEType > 0)
                _mesh->nce ++;
            else
            {
                fprintf(stdout,"Check for tce");
                exit(0);
            }
        }
    }

    _mesh->CondV = (struct ConditionV *) malloc(2*((_mesh->nxp+1)*(_mesh->nyp + 1)
            + (_mesh->nyp + 1)*(_mesh->nzp + 1) + (_mesh->nzp+1)*(_mesh->nxp+1))*sizeof(struct ConditionV));
    _mesh->ncv = 0;

}
