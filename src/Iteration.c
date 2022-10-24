//
// Created by huacheng on 2020/12/25.
//

#include "Iteration.h"

void CalculateGradM(struct Element *_element, int nm)
{
    if(_element->CutTag)
    {
        for(int i=0;i<DIM;i++)
        {
            for(int j=0;j<DIM;j++)
            {
                _element->GradVel[i][j] = 0.0;
                _element->Strain[i][j]  = 0.0;
            }
        }
        _element->DivVel   = 0.0;
        _element->InvSta2  = 0.0;
        _element->InvSta3  = 0.0;
        _element->sInvSta2 = 0.0;
    }
    else
    {
        for (int k = 0; k < DIM; k++)
        {
            for (int j = 0; j < DIM; j++)
            {
                _element->GradVel[k][j] = 0.0;
                for (int i = 0; i < VPE; i++)
                {
                    _element->GradVel[k][j] += _element->NeV[i]->Velocity[k] * _element->GradFactor[j][i];
                }
            }
        }
        _element->DivVel = _element->GradVel[X][X] + _element->GradVel[Y][Y] + _element->GradVel[Z][Z];
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                _element->Strain[i][j] = (_element->GradVel[i][j] + _element->GradVel[j][i]) / 2.0;
            }
            _element->Strain[i][i] -= 1.0 / 3.0  * _element->DivVel;
        }

        // the second invariant
        _element->InvSta2 = SecondInvariant(_element->Strain);
        _element->InvSta3 = Determinant(_element->Strain);
        _element->sInvSta2 = sqrt(_element->InvSta2);
    }
}

void NoStressCell(struct  Element *_element)
{
    ZeroEx(_element->PsiL3+1,6);
    _element->VibPressure = 0.0;
    _element->PsiL3[8] = 0.0;
    _element->InvSta2  = 0.0;
}



void RheologyUpdate(struct Element *_element, double dt, const struct StateReference _s[], int nm)
{
    /*
     * calculate the stress using rheology model.
     * the first version only include 2 types of rheology.
     * Elastic Model and Invicosity Fluid
     * (Elastic    = [state] 0)
     * (Invicosity = [state] 1)
     * (Nostress   = [state] 2)
     * (Newton Fluid = [ state] 3)
     */


    // step 0: determine the state of this element
    GetState(_element, _s, nm);
    double * Stress = _element->PsiL3 + 1;
    double * VibVelocity = _element->PsiL3 + 8;
    // step 1: Get the increment of stress (Elastic)
    // accord to the Jaumann Equation
    if (SolidElement == _element->State)
    {

        _element->sInvSte2_old = sqrt(SecondInvariantSym(Stress));
        double Rotation[DIM][DIM];
        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                Rotation[i][j] = 0.5 * (_element->GradVel[i][j] - _element->GradVel[j][i]);

        // calculation Rotation * Stress(old)
        double tStrainEla[DIM][DIM] = {0.0};
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                for (int k = 0; k < DIM; k++)
                    tStrainEla[i][j] += Rotation[i][k] * Stress[SymIndex(k,j)];
            }
        }

        // += 2*dt*G*Strain + dt*(transpose<R*Stress> + R*Stress)
        double tmpDt = 2.0 * dt * _element->ShearModule;
        for (int i = 0; i < DIM; i++)
        {
            for (int j = i; j < DIM; j++)
            {
                Stress[SymIndex(i,j)] += tmpDt * _element->Strain[i][j] + dt * (tStrainEla[i][j] + tStrainEla[j][i]);
            }
        }
    }
    else if (FluidElement == _element->State)
    {
        // viscosity Newton fluid
        for(int i=0;i<DIM;i++)
        {
            for(int j=i;j<DIM;j++)
            {
                Stress[SymIndex(i,j)] = 2.0*_element->Viscosity * _element->Strain[i][j];
            }
        }

        *VibVelocity = 0.0;
        _element->VibPressure = 0.0;
        _element->PsiL3[0] = 1.0;   //Damage
    }
    else if (VoidElement == _element->State)
    {
        /*
         * no stress in this element
         */
        ZeroEx(_element->PsiL3+1,6);
        _element->VibPressure = 0.0;
        _element->PsiL3[8] = 0.0;
        _element->InvSta2  = 0.0;
    }
    else if(VaporElement == _element->State)
    {
        /*
         * Vapor Element
         * zero the stress and plastic strain
         */

        _element->sInvSta2 = 0.0;
        ZeroEx(_element->PsiL3,9);
        _element->VibPressure = 0.0;
    }
}

void FailureUpdate(struct Element *_element, double dt, const struct StateReference _s[], struct MeshInfo *_minfo)
{
    /*
     * Update the yield strength of element and
     * cut the stress tensor
     * When sqrt(J2(Stress)) > Yield Strength ==> Shear Failure
     *                         [FailureState = 1]
     * When abd(Sigma1 - Sigma3) > Tensile Strength ==> Tensile Failure
     *                         [FailureState = 2]
     * When Shear Failure + Tensile Failure appear together
     *                         [FailureState = 3]
     */

    int nm = _minfo->MaterialNum;
    if (SolidElement != _element->State)
    {
        _element->ShearStrength   = 0.0;
        _element->TensileStrength = 0.0;
        _element->PsiL3[7] = Wind(_element->PsiL3[7]+dt*_element->sInvSta2,0.0,10.0);
    }
    else
    {
        double * Stress = _element->PsiL3 + 1;
        double InvSte2 = SecondInvariantSym(Stress);
        _element->sInvSte2        = sqrt(InvSte2);
        _element->ShearStrength   = 0.0;
        _element->TensileStrength = 0.0;
        _element->FailureState    = 0;


        for(int i = 1; i < nm; i++)
        {
            // the strength contributed by vacuum is ignored.
            // mix-strength model
            _element->ShearStrength   += _element->VOF[i] * _s[i].Yfunc(_element, &_s[i]);
            _element->TensileStrength += _element->VOF[i] * _s[i].Tfunc(_element, &_s[i]);
//            _element->ShearStrength   = Min(_s[i].Yfunc(_element, &_s[i]),_element->ShearStrength);
//            _element->TensileStrength = Min(_s[i].Tfunc(_element, &_s[i]),_element->TensileStrength);
        }

        if(_element->sInvSte2 > _element->ShearStrength)
        {
            _element->FailureState += 1;
            double tmpS = _element->ShearStrength / _element->sInvSte2;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = i; j < DIM; j++)
                    Stress[SymIndex(i,j)] *= tmpS;
            }

            double _sInvSta2_2d = _element->sInvSta2/Max(0.5,(_element->sInvSte2 - _element->sInvSte2_old)/(_element->sInvSte2 - _element->ShearStrength));
            _element->sInvSta2 = Min((_element->sInvSte2 - _element->ShearStrength)/(2.0*_element->ShearModule*dt),_element->sInvSta2);

            /*
             * The general damage increase model
             * LINEAR　is default
             * The exponential model will be added in future code
             */
//            double sInvSta2dt = _sInvSta2_2d * dt;
            double sInvSta2dt = _element->sInvSta2*dt;
            double dSDamage   = 0.0;
            for (int i = 1; i < nm; i++)
            {
                dSDamage += _element->VOF[i] * _s[i].dSDf(_element, &_s[i]);
            }
            _element->Damage = Min(1.0, _element->Damage +  sInvSta2dt*dSDamage);
        }


        if(1 == _minfo->TensileFailure)
        {
            // calculate the tensile failure
            double tmpx = _element->Pressure + _element->TensileStrength;
            double Invste3 = DeterminantSym(Stress);
            if ((tmpx > 0.0) && (tmpx * tmpx * tmpx - InvSte2 * tmpx - Invste3 > 0.0) && (3.0 * tmpx * tmpx - InvSte2 > 0.0))
            {
                /*
                 * no tensile failure;
                 * do nothing.
                 */
            }
            else
            {
                _element->FailureState |= 2;

                /*
                 * it is very strange to update the damage increment here
                 * just copy it from iSALE2d
                 */

                double DamageOld = _element->Damage;
                double tD = pow(DamageOld, 1.0 / 3.0);
                tD += 0.4 * _element->Csound * dt * _element->rdx_min;
                _element->Damage = Min(1.0, tD * tD * tD);
                double stressfrac = 0.0;
                if (DamageOld < 0.98)
                    stressfrac = (1.0 - _element->Damage) / (1.0 - DamageOld);

                for (int i = 0; i < DIM; i++)
                {
                    for (int j = i; j < DIM; j++)
                        Stress[SymIndex(i,j)] *= stressfrac;
                }
            }
        }



        /*
         * NOTICE:
         * 3 function pointers: StateReference->Tfunc, StateReference->Yfunc, and StateReference->dSDf
         * can be defined by user. See State.h
         */

        if (0 != _element->FailureState)
        {
            /*
             *  Update the total PLASTIC strain in the cell. For cells that are
             * yielding, the plastic strain increment for the current timestep is
             * the total strain rate in the timestep (erate) times the timestep
             * size (dt). We limit this to a strain of 10--strains larger than this
             * are somewhat meaningless. . .
             */
            _element->PsiL3[7] = Wind(_element->PsiL3[7]+dt*_element->sInvSta2,0.0,10.0);
//            _element->PsiL3[7] = Wind(_element->PsiL3[7]+_element->sInvSta2,0.0,10.0);
        }
    }

    _element->PsiL3[0]  = _element->Damage;
    ScalerMoveEx(_element->sPsiL3,_element->PsiL3,_element->Mass,NPSI);
}

void SyncForward(struct Element *_element)
{
    // synchronize the data between Damage and (Psi[0],sPsi[0])
    // ...                          Stress and (Psi[1:6],sPsi[1:6])
    _element->PsiL3[0] = _element->Damage;
    /*
    _element->PsiL3[1] = _element->Stress[0][0];
    _element->PsiL3[2] = _element->Stress[1][1];
    _element->PsiL3[3] = _element->Stress[2][2];
    _element->PsiL3[4] = _element->Stress[0][1];
    _element->PsiL3[5] = _element->Stress[0][2];
    _element->PsiL3[6] = _element->Stress[1][2];
    */
    for(int i=0;i<=6;i++)
        _element->sPsiL3[i] = _element->PsiL3[i] * _element->Mass;

    // synchronize data between PhiL3 and sPhiL3
    for(int i=0;i<NMT;i++)
        ScalerMoveEx(_element->sPhiL3[i],_element->PhiL3[i],_element->VOF[i],NPHI);
//        _element->sPhiL3[i][2] = _element->PhiL3[i][2] * _element->VOF[i];
}

void SyncBackward(struct Element * _e)
{
    // the average density is sum by weight of materials
    // Sync-Backward
    // Copy Psi[0] to damage
    _e->Damage = _e->PsiL3[0];
    // stress dev, is stored in sPsiL3[1:6]
    /*
    _e->Stress[0][0] = _e->PsiL3[1];
    _e->Stress[1][1] = _e->PsiL3[2];
    _e->Stress[2][2] = _e->PsiL3[3];
    _e->Stress[0][1] = _e->Stress[1][0] = _e->PsiL3[4];
    _e->Stress[0][2] = _e->Stress[2][0] = _e->PsiL3[5];
    _e->Stress[1][2] = _e->Stress[2][1] = _e->PsiL3[6];
    */
    for(int k=0;k<NMT;k++)

        _e->VOF[k] = _e->sPhiL3[k][0];
}

void CalculateABV(struct Element *_element, double _artvis)
{
    // get the artificial bulk viscosity
    // it is a [ pressure ]
    // _artvis is the lambda in sale3d manual and should be < 0.25 to avoid excessive damping.
    // 0.1 has been satisfactory for many application
    if(2 == _element->State || _element->DivVel > 0.0)
    {
        _element->ArtificialPressure = 0.0;
    }
    else
    {
        double Area = _element->ABVArea;
        _element->ArtificialPressure = _artvis * _element->DivVel * _element->DivVel * _element->Density * Area / 6.0;
        double AAP = 0.1;
        if(fabs(_element->ArtificialPressure) > AAP*fabs(_element->Pressure))
            _element->ArtificialPressure *= (fabs(_element->Pressure) + 1e-6)/(fabs(_element->ArtificialPressure) + 1e-6) * AAP;
    }
}

void CalculateABV2(struct Element *_element, double _artvis, double _artvis2)
{
    // get the artificial bulk viscosity
    // it is a [ pressure ]
    // _artvis is the lambda in sale3d manual and should be < 0.25 to avoid excessive damping.
    // 0.1 has been satisfactory for many application
    if(VoidElement == _element->State || _element->DivVel > 0.0)
    {
        _element->ArtificialPressure = 0.0;
    }
    else
    {
        double Area = _element->ABVArea/6.0;
        double Qline = _artvis2 * sqrt(fabs(Area)) * _element->Csound;
        double Qquad = _artvis * Area * _element->DivVel;
        _element->ArtificialPressure = (Qquad - Qline) * _element->DivVel * _element->Density;
    }
}

void CalculateBFV(struct Boundary *_boundary)
{
    /*
     * just calculate the flux volume that across the bounary .
     * before Courant number of donor cell
     */
    double _boundaryFluxVolume = 0.0;
    for (int k = 0; k < VPB; k++)
    {
        struct Vertex *kNeV = _boundary->NeV[k];
        _boundaryFluxVolume += -1.0 * Dot(kNeV->DistanceL1, _boundary->FluxFactor[k]);
    }

    _boundary->FluxVolume = _boundaryFluxVolume;
}

void CalculateCourant(struct Element * doe)
{
    doe->Courant = 0.0;
    for(int k=0;k<BPE;k++)
        doe->Courant += Max(0.0,doe->NeB[k]->FluxVolume*FlowWeight[k]);
    doe->Courant *= doe->rVolume;
}

void CalculateVOFM(struct Boundary *_boundary, double CISSAMK1, int nm)
{
    // calculate the volume across this face
    // notice the minus
    // vg = -vp ==> distance is - Distance L1

    double _boundaryFluxVolume = _boundary->FluxVolume;

    // nm is is the max index of materials needed to be updated
    double VOFU, VOFD, VOFA; // they are temporary variables
    struct Element * upe;
    struct Element * doe;
    struct Element * ade;
    double CutMatVOF;

    if(_boundaryFluxVolume < 0)
    {
        upe = _boundary->UpE[0];
        doe = _boundary->NeE[0];
        ade = _boundary->NeE[1];
        CutMatVOF = CutVOF(_boundary, 0);
    }
    else
    {
        upe = _boundary->UpE[1];
        doe = _boundary->NeE[1];
        ade = _boundary->NeE[0];
        CutMatVOF = CutVOF(_boundary, 1);
    }


    double CourantD,Courantf;
    Courantf = fabs(_boundaryFluxVolume * doe->rVolume);
    CourantD = doe->Courant;

    double _boundaryVOF[NMT] = {0.0};
    if(CutMatVOF > 1.0 - TOLVOF)
    {
        _boundaryVOF[CutMaterial] = 1.0;
    }
    else
    {
        for (int iMat = 0; iMat < nm; iMat++)
        {
            if (iMat == CutMaterial) continue;
            VOFU = upe->VOF[iMat];
            VOFD = doe->VOF[iMat];
            VOFA = ade->VOF[iMat];

            double tmp = VOFA - VOFU;
            double nVOFD = (VOFD - VOFU) / tmp;

            if (fabs(tmp) < TOLVOF)
            {
                _boundaryVOF[iMat] = VOFD;
            } else if (nVOFD < TOLVOF || nVOFD > 1 - TOLVOF)
            {
                _boundaryVOF[iMat] = VOFD;
            } else
            {
                double CosTf;
                CosTf = Dot(doe->GradVOF[iMat], _boundary->Normal);
                CosTf = CosTf * CosTf /
                        (doe->GradVOF[iMat][X] * doe->GradVOF[iMat][X] + doe->GradVOF[iMat][Y] * doe->GradVOF[iMat][Y]
                         + doe->GradVOF[iMat][Z] * doe->GradVOF[iMat][Z]);
                double Gf, nVOFf;
                Gf = Min(CISSAMK1 * CosTf * CosTf, 1.0);
                nVOFf = vofMSTACS(nVOFD, Gf, CourantD);
                _boundaryVOF[iMat] = VOFU + nVOFf * (VOFA - VOFU);
            }
        }

        double Zeta  = 0.0;
        double ZetaD = 0.0;
        _boundaryVOF[CutMaterial] = CutMatVOF;
        for(int iMat=0;iMat<nm;iMat++)
        {
            if(iMat != CutMaterial)
            {
                Zeta  += _boundaryVOF[iMat];
                ZetaD += doe->VOF[iMat];
            }
        }
        if(Zeta > TOLVOF)
        {
            for(int iMat=0;iMat<nm;iMat++)
            {
                if(iMat != CutMaterial)
                    _boundaryVOF[iMat] =  _boundaryVOF[iMat]/Zeta * (1.0 - CutMatVOF);
            }
        } else if(ZetaD > TOLVOF)
        {
            for(int iMat=0;iMat<nm;iMat++)
            {
                if(iMat != CutMaterial)
                    _boundaryVOF[iMat] =  doe->VOF[iMat]/ZetaD * (1.0 - CutMatVOF);
            }
        }
        else
        {
            ZeroEx(_boundaryVOF,nm);
            _boundaryVOF[CutMaterial] = 1.0;
        }
    }


    /*
     * Following code is merged from CalculateL3a:
     * Calculating the materials that get through this boundary.
     */
    _boundary->dMass = 0.0;
    for(int k = 0; k < nm; k++)
    {
        double tmpVol  = _boundaryVOF[k] * _boundaryFluxVolume;
//        if(fabs(tmpVol) > doe->VOF[k]*doe->Volume)
//            tmpVol = Sign(_boundaryFluxVolume) * doe->VOF[k]*doe->Volume;
        ScalerMoveEx(_boundary->dPhi[k],doe->PhiL3[k],tmpVol,NPHI);
        _boundary->dMass += tmpVol*doe->PhiL3[k][1];
    }

    if(doe->Mass > doe->Volume*TOLRHO)
    {
        ScalerMove(_boundary->dMomentum,doe->Momentum,_boundary->dMass/doe->Mass);
        ScalerMoveEx(_boundary->dPsi,doe->PsiL3,_boundary->dMass,NPSI);
    }
    else
    {
        Zero(_boundary->dMomentum);
        ZeroEx(_boundary->dPsi,NPSI);
    }
}

void CalculateVOFTVD(struct Boundary *_boundary, double CISSAMK1, int nm)
{
    // this is another version of CalculateVOF to
    // estimate the value of VOF for multi-materials
    // Calculate the gradient direction of the inner boundary
    // this boundary is the division between the vacuum and the fluid occupied field

    // calculate the volume across this face
    // notice the minus
    // vg = -vp ==> distance is - Distance L1

    double _boundaryFluxVolume = _boundary->FluxVolume;

    // nm is is the max index of materials needed to be updated
    double VOFU, VOFD, VOFA; // they are temporary variables
    struct Element * upe;
    struct Element * doe;
    struct Element * ade;

    if(_boundaryFluxVolume < 0)
    {
        upe = _boundary->UpE[0];
        doe = _boundary->NeE[0];
        ade = _boundary->NeE[1];
    }
    else
    {
        upe = _boundary->UpE[1];
        doe = _boundary->NeE[1];
        ade = _boundary->NeE[0];
    }

    doe->Courant = 0.0;
    for(int k=0;k<BPE;k++)
        doe->Courant += Max(0.0,doe->NeB[k]->FluxVolume*FlowWeight[k]);
    doe->Courant *= doe->rVolume;
    double CourantD,Courantf;
    CourantD = doe->Courant;
    Courantf = fabs(_boundaryFluxVolume * doe->rVolume);

    double _boundaryVOF[NMT];
    for (int iMat = 0; iMat < nm; iMat++)
    {
        VOFU = upe->VOF[iMat];
        VOFD = doe->VOF[iMat];
        VOFA = ade->VOF[iMat];

        double tmp = VOFA - VOFU;
        double nVOFD = (VOFD - VOFU) / tmp;

        if (fabs(tmp) < TOLVOF)
        {
            _boundaryVOF[iMat] = VOFD;
        }
        else if(nVOFD < TOLVOF || nVOFD > 1 - TOLVOF)
        {
            _boundaryVOF[iMat] = VOFD;
        }
        else
        {
            double CosTf;
            CosTf = Dot(doe->GradVOF[iMat],_boundary->Normal);
            CosTf = CosTf*CosTf/(doe->GradVOF[iMat][X]*doe->GradVOF[iMat][X] + doe->GradVOF[iMat][Y]*doe->GradVOF[iMat][Y]
                                 + doe->GradVOF[iMat][Z]*doe->GradVOF[iMat][Z]);

            double Gf = Min(CISSAMK1 * CosTf, 1.0);
            Gf = Max(0.0,Gf);
            double nVOFCBC, nVOFUQ;
            nVOFCBC = Min(nVOFD / CourantD, 1.0);
            nVOFUQ  = Min(CourantD * nVOFD + (1.0 - CourantD) * (6.0 * nVOFD + 3.0) / 8.0, nVOFCBC);
            double nVOFf = Gf * nVOFCBC + (1.0 - Gf) * nVOFUQ;
            _boundaryVOF[iMat] = VOFD + (nVOFf - nVOFD) * (VOFA - VOFU);

            /*
            if(fabs(_boundaryVOF[iMat]) < TOLVOF)
                _boundaryVOF[iMat] = 0.0;
            */
            /*
             * Crop on _boundaryVOF is removed, when the dt is very small, all of
             * the _boundaryVOF[iMat] is much less than TOLVOF and all of the flux diminish
             * System is stuck in a trap arisen by numerical.
            */
        }
    }

    /*
     * Check the validity of boundaryVOF
     */
    double SumVOF = 0.0;
    for(int iMat = 0;iMat< nm;iMat++)
    {
        SumVOF += _boundaryVOF[iMat];
    }

    /*
     * refactor interface value by tvd
     */
    double tPhiL3[NMT][NPHI];
    double tMass = 0.0;
    double tMomentum[DIM];
    Refactor(upe->Momentum,doe->Momentum,ade->Momentum,tMomentum,DIM);
    for(int i=0;i<nm;i++)
    {
        Refactor(upe->PhiL3[i],doe->PhiL3[i],ade->PhiL3[i],tPhiL3[i],NPHI);
        tMass += tPhiL3[i][1]*tPhiL3[i][0]*doe->Volume;
    }


    /*
     * Following code is merged from CalculateL3a:
     * Calculating the materials that get through this boundary.
     */
    _boundary->dMass = 0.0;
    for(int k = 0; k < nm; k++)
    {
        double tmpVol  = _boundaryVOF[k] * _boundaryFluxVolume / SumVOF;
        if(fabs(tmpVol) > doe->VOF[k]*doe->Volume)
            tmpVol = Sign(_boundaryFluxVolume) * doe->VOF[k]*doe->Volume;

        ScalerMoveEx(_boundary->dPhi[k],tPhiL3[k],tmpVol,NPHI);
        _boundary->dMass += tmpVol*tPhiL3[k][1];
    }

    if(tMass > doe->Volume*TOLVOF*TOLRHO)
        ScalerMove(_boundary->dMomentum,tMomentum,_boundary->dMass/tMass);
    else
        Zero(_boundary->dMomentum);
    ScalerMoveEx(_boundary->dPsi,doe->PsiL3,_boundary->dMass,NPSI);
}

void Refactor(const double * Uu, const double * Du, const double * Au, double * fu, int _len)
{
    for(int i=0;i<_len;i++)
    {
          double ul = Du[i] - Uu[i];
          double ur = Au[i] - Du[i];
          double su = (2.0*ul*ur + 1.0e-6)/(ul*ul + ur*ur + 1.0e-6);
          fu[i]  = Du[i] + 0.25*su*((1.0 - su/3.0)*ul + (1.0 + su/3.0)*ur);
    }

}


void EnergyUpdate(struct Element *_element, const struct MeshInfo * _minfo , double dt)
{
    /*
     * _name changed: CalculateL2() -> EnergyUpdate()
     * Update Kinetic Energy.
     * Ignoring the empty element
     */

    int nm = _minfo->MaterialNum;

    if(_element->VOF[0]<1.0-TOLVOF)
    {
        double EVelocityL1[DIM] = {0.0}, tmpTotalMass = 0.0;
        for(int i=0;i<VPE;i++)
        {
            tmpTotalMass += _element->CutRatio[i];
            ScalerAddition(EVelocityL1,_element->NeV[i]->VelocityL1,_element->CutRatio[i]);
        }
        ScalerAddition(_element->Momentum,EVelocityL1,dt*_element->Mass);
    }

    /*
     * Update Internal energy.
     * Ignoring the empty and void element
     */
    if(_element->VOF[VACUUM] > TOLVOF || _element->Density < TOLRHO) return;

    double dE = -1.0*(_element->Pressure + _element->ArtificialPressure)* _element->DivVel
                + ContractionSym(_element->PsiL3+1,_element->GradVel);
    double ldE = dE * dt;

    for(int i=1;i<nm;i++)
    {
        if(_element->sPhiL3[i][0] < TOLVOF) continue;

        _element->PhiL3[i][2]  = Max(0.0,_element->PhiL3[i][2] + ldE);
        _element->sPhiL3[i][2] = _element->PhiL3[i][2] * _element->sPhiL3[i][0];
    }
}


// the coefficient of the flow: represent the boundary direction
// the order is Left/Right/Front/Back/Bottom/Top
const double FlowWeight[BPE] = {-1.0, 1.0, 1.0, -1.0, -1.0, 1.0};
void CalculateL3b(struct Element *_e, const struct MeshInfo * _minfo)
{

    /*
     * The information of momentum or velocity in elements are
     * stored in boundary.
     * here, _e->Momentum is cleared and its definition is changed]
     * to the increment of momentum because of advection.
     */

    int nm = _minfo->MaterialNum;
    for(int i=0;i<BPE;i++)
    {
        _e->Mass += FlowWeight[i]*_e->NeB[i]->dMass;
        ScalerAddition(_e->Momentum,_e->NeB[i]->dMomentum,FlowWeight[i]);
        ScalerAdditionEx(_e->sPsiL3,_e->NeB[i]->dPsi,FlowWeight[i],NPSI);
    }

#ifdef SALEC_ADD_TENSILE_VOID
    /*
     * add void space according to Senft and Stewart(2007), “Modeling Impact Cratering in Layered Surfaces.”
     */
    if(_e->FailureState & 2)
    {
        _e->sPhiL3[VACUUM][0] += _e->Damage*_e->Courant*_e->Volume;
    }
#endif


    double Zeta = 0.0;
    for(int i = 0; i < nm; i++)
    {
        for (int j = 0; j < BPE; j++)
        {
            ScalerAdditionEx(_e->sPhiL3[i], _e->NeB[j]->dPhi[i], FlowWeight[j] * _e->rVolume, NPHI);
        }
        _e->sPhiL3[i][0] = Max(_e->sPhiL3[i][0],0.0);
        if(i != 0)
            Zeta += _e->sPhiL3[i][0];
    }

    /* this is similar with sale2d
    Compress or expand each material proportionally with
    volume fraction of the material.  In essence, this assumes
    the compressibility of each material is the same. . . */

    double ZetaLim = 32.0;
    if((Zeta > 1.0-TOLVOF) || (_e->sPhiL3[0][0] < ZetaLim*TOLVOF && Zeta > 0.0))
    {
        _e->sPhiL3[0][0] = 0.0;
        for (int i = 1; i < nm; i++)
            _e->sPhiL3[i][0] = _e->sPhiL3[i][0] / Zeta;
    }
    /*else if((1==_minfo->PartPressure)&&(Zeta<=1.0-_e->sPhiL3[0][0]-TOLVOF))
    {
        Zeta += _e->sPhiL3[0][0];
        for (int i = 0; i < nm; i++)
        {
            _e->sPhiL3[i][0] = _e->sPhiL3[i][0] / Zeta;
        }
    }*/
    else
    {
        _e->sPhiL3[0][0] = 1.0 - Zeta;
    }



    /*
     * for some extreme cases, this method is not valid
     * i.e, the volume Phi[k][0] limits to 0 when residual mass Phi[k][1] doesn't limit to 0.
     * the density will be infinity.
     */
    double tMass = 0.0;
    for(int i = 1; i < nm; i++)
    {
        if (_e->sPhiL3[i][0] > TOLVOF)
        {
            _e->PhiL3[i][0] = 1.0;
            _e->PhiL3[i][1] = _e->sPhiL3[i][1] / _e->sPhiL3[i][0];
            _e->PhiL3[i][2] = _e->sPhiL3[i][2] / _e->sPhiL3[i][0];
        }
        else
        {
            _e->PhiL3[i][0] = 1.0;
            _e->PhiL3[i][1] = 0.0;
            _e->PhiL3[i][2] = 0.0;
            ZeroEx(_e->sPhiL3[i],NPHI);
        }
        _e->VOF[i]  = _e->sPhiL3[i][0];
        tMass      += _e->sPhiL3[i][1];
    }
    tMass *= _e->Volume;

    if((tMass < TOLRHO*_e->Volume) || (_e->sPhiL3[0][0] > 1.0 - TOLVOF))
    {
        for(int k=0;k<nm;k++)
        {
            ZeroEx(_e->sPhiL3[k],NPHI);
            ZeroEx(_e->PhiL3[k],NPHI);
            _e->PhiL3[k][0] = 1.0;
        }
        Zero(_e->Momentum);
        ZeroEx(_e->sPsiL3,NPSI);
        ZeroEx(_e->PsiL3,NPSI);
        ZeroEx(_e->VOF,NMT);
        _e->sPhiL3[0][0] = 1.0;
        _e->VOF[0]       = 1.0;
        _e->Mass         = 0.0;
        // need improve !!
    }
    else
    {
        _e->Mass = tMass;
        _e->VOF[0] = _e->sPhiL3[0][0];
        ScalerMoveEx(_e->PsiL3,_e->sPsiL3,1.0/_e->Mass,NPSI);
    }


    // the average density is sum by weight of materials
    // Sync-Backward
    _e->Density = _e->Mass * _e->rVolume;
    // Copy Psi[0] to damage
    _e->PsiL3[0]  = Wind(_e->PsiL3[0],0.0,1.0);
    _e->sPsiL3[0] = _e->PsiL3[0] * _e->Mass;
    _e->Damage    = _e->PsiL3[0];

    // stress dev, is stored in sPsiL3[1:6]

}


void CalculateL3bOut_(struct Element * _e, int nm)
{
    if(_e->Density < TOLRHO) return;

    double tMass        = 0.0;
    int DirtyBit        = 0;
    int OutTable[BPE]   = {0};

    for(int j=0;j<nm;j++)
    {
        if(_e->sPhiL3[j][0] < TOLVOF) continue;
        double RelativePhi = 0.0;
        for(int i=0;i<BPE;i++)
        {
            double tPhi0 = _e->NeB[i]->dPhi[j][0] * FlowWeight[i];
            if(tPhi0 < 0.0)
            {
                OutTable[i]  = 1;
                RelativePhi += tPhi0*_e->rVolume;
            }
        }

        double tSumPhi0 = RelativePhi + _e->sPhiL3[j][0];
        if(tSumPhi0 >= 0.0)
        {
            _e->sPhiL3[j][0] = tSumPhi0;
            _e->sPhiL3[j][1] = tSumPhi0 * _e->PhiL3[j][1];
            _e->sPhiL3[j][2] = tSumPhi0 * _e->PhiL3[j][2];
            tMass           += _e->sPhiL3[j][1];
        }
        else
        {
            // restrict the RelativePhi to sPhiL3[j][0]
            DirtyBit = 1;
            double Restrict = -1.0*_e->sPhiL3[j][0]/RelativePhi;
            ZeroEx(_e->sPhiL3[j],NPHI);
            for(int i=0;i<BPE;i++)
            {
                if(OutTable[i]) ScalerEx(_e->NeB[i]->dPhi[j],Restrict,NPHI);
            }
        }
    }

    ScalerEx(_e->Momentum,tMass/_e->Density,DIM);
    ScalerEx(_e->sPsiL3,tMass/_e->Density,NPSI);
    _e->Mass    = tMass * _e->Volume;
    _e->Density = tMass;

    if(0==DirtyBit) return;
    for(int i=0;i<BPE;i++)
    {
        if(0 == OutTable[i]) continue;
        if(fabs(_e->NeB[i]->dMass) < TOLMASS) continue;
        double bMass = 0.0;
        for(int j=0;j<nm;j++)
        {
            bMass += _e->NeB[i]->dPhi[j][1];
        }
        double sfrac = bMass/_e->NeB[i]->dMass;
        ScalerEx(_e->NeB[i]->dMomentum,sfrac,DIM);
        ScalerEx(_e->NeB[i]->dPsi,sfrac,NPSI);
        _e->NeB[i]->dMass = bMass;
    }
}

void CalculateL3bOut(struct Element * _e, int nm)
{

    /*
     * cut off some overflow material
     */
#if ELEMENT_DEBUG
    _e->Debug = 0.0;
#endif
    if(_e->Density < TOLRHO) return;

    int DirtyBit        = 0;
    int OutTable[BPE]   = {0};
    double _TolMass     = _e->Volume*TOLRHO;

    for(int j=0;j<nm;j++)
    {
        if(_e->sPhiL3[j][0] < TOLVOF) continue;
        double RelativePhi = 0.0;
        for(int i=0;i<BPE;i++)
        {
            double tPhi0 = _e->NeB[i]->dPhi[j][0] * FlowWeight[i];
            if(tPhi0 < 0.0)
            {
                OutTable[i]  = 1;
                RelativePhi += tPhi0*_e->rVolume;
            }
        }

        double tSumPhi0 = RelativePhi + _e->sPhiL3[j][0];
        if(tSumPhi0 < 0.0)
        {
            // restrict the RelativePhi to sPhiL3[j][0]
            DirtyBit = 1;
            double Restrict = -1.0*_e->sPhiL3[j][0]/RelativePhi;
            for(int i=0;i<BPE;i++)
            {
                if(OutTable[i]) ScalerEx(_e->NeB[i]->dPhi[j],Restrict,NPHI);
            }
#if ELEMENT_DEBUG
            _e->Debug += 1.0*(1<<j);
#endif
        }
    }

    if(0==DirtyBit) return;
    for(int i=0;i<BPE;i++)
    {
        if(0 == OutTable[i]) continue;
        if(fabs(_e->NeB[i]->dMass) < _TolMass) continue;
        double bMass = 0.0;
        for(int j=0;j<nm;j++)
        {
            bMass += _e->NeB[i]->dPhi[j][1];
        }
        double sfrac = bMass/_e->NeB[i]->dMass;
        ScalerEx(_e->NeB[i]->dMomentum,sfrac,DIM);
        ScalerEx(_e->NeB[i]->dPsi,sfrac,NPSI);
        _e->NeB[i]->dMass = bMass;
    }

}

int MaterialStep(int Start, int nm,const double Remnant[],const int Sht[],int Positive)
{
    int k = Start;
    if(Positive)
    {
        while(k<nm)
        {
            if(Remnant[Sht[k]] > 0.0)
                break;
            else
                k++;
        }
    } else
    {
        while(k<nm)
        {
            if(Remnant[Sht[k]] < 0.0)
                break;
            else
                k++;
        }
    }

    return k;
}

void MaterialSwap(int mata, int matb, double tpart,struct Element * _e,const double tra[])
{
    // replace flow of mata with matb
    for(int k=0;k<BPE;k++)
    {
        if(tra[k] > 0.0)
        {
            _e->NeB[k]->dPhi[mata][0] += tpart * tra[k] * FlowWeight[k];
            _e->NeB[k]->dPhi[matb][0] -= tpart * tra[k] * FlowWeight[k];
        }
    }
}

void CalculateL3bSimple(struct Element * _e, int nm)
{
#if ELEMENT_DEBUG
    _e->Debug = 0.0;
#endif
    if(_e->Density < TOLRHO) return;

    int DirtyBit        = 0;
    int OutTable[BPE]   = {0};
    double Outflow[NMT],Remnant[NMT];

    for(int j=1;j<nm;j++)
    {
        if(_e->sPhiL3[j][0] < TOLVOF) continue;
        double RelativePhi = 0.0;
        for(int i=0;i<BPE;i++)
        {
            double tPhi0 = _e->NeB[i]->dPhi[j][0] * FlowWeight[i];
            if(tPhi0 < 0.0)
            {
                OutTable[i]  = 1;
                RelativePhi += tPhi0;
            }
        }

        double tSumPhi0 = RelativePhi + _e->sPhiL3[j][0]*_e->Volume;

        Outflow[j] = RelativePhi;
        Remnant[j] = tSumPhi0;

        if(tSumPhi0 < 0.0)
        {
            // restrict the RelativePhi to sPhiL3[j][0]
            DirtyBit = 1;
#if ELEMENT_DEBUG
            _e->Debug += 1.0*(1<<j);
#endif
        }
    }

    if(0==DirtyBit) return;
    double OutRatio[NMT][BPE] = {0.0};
    for(int j=0;j<nm;j++)
    {
        if((_e->sPhiL3[j][0] < TOLVOF) || (fabs(Outflow[j]) < TOLVOF*_e->Volume))
            ZeroEx(OutRatio[j],BPE);
        else
        {
            for(int i=0;i<BPE;i++)
            {
                if(OutTable[i])
                {
                    OutRatio[j][i] = FlowWeight[i] * _e->NeB[i]->dPhi[j][0]/Outflow[j];
                } else
                    OutRatio[j][i] = 0.0;
            }
        }
    }

    if(Remnant[1]>0.0 && Remnant[2]<0.0)
    {
        MaterialSwap(2,1,Min(-Remnant[2],Remnant[1]),_e,OutRatio[2]);
    } else if(Remnant[1]<0.0 && Remnant[2]>0.0)
    {
        MaterialSwap(1,2,Min(-Remnant[1],Remnant[2]),_e,OutRatio[1]);
    }

    for(int i=0;i<BPE;i++)
    {
        if(0==OutTable[i]) continue;
        double bMass = 0.0;
        for(int j=0;j<nm;j++)
        {
            if(_e->NeB[i]->dPhi[j][0] * FlowWeight[i] > 0.0)
            {
                fprintf(stdout,"Error in Adjust!\n");
                exit(0);
            }
            _e->NeB[i]->dPhi[j][1] = _e->PhiL3[j][1] * _e->NeB[i]->dPhi[j][0];
            _e->NeB[i]->dPhi[j][2] = _e->PhiL3[j][2] * _e->NeB[i]->dPhi[j][0];
            bMass += _e->NeB[i]->dPhi[j][1];
        }

        _e->NeB[i]->dMass = bMass;
        ScalerMoveEx(_e->NeB[i]->dMomentum,_e->Momentum,bMass/_e->Mass,DIM);
        ScalerMoveEx(_e->NeB[i]->dPsi,_e->PsiL3,bMass,NPSI);
    }
}


void CalculateL3bAdjust(struct Element * _e, int nm)
{
#if ELEMENT_DEBUG
    _e->Debug = 0.0;
#endif
    if(_e->Density < TOLRHO) return;

    int DirtyBit        = 0;
    int OutTable[BPE]   = {0};
    double Outflow[NMT],Remnant[NMT];

    for(int j=0;j<nm;j++)
    {
        if(_e->sPhiL3[j][0] < TOLVOF) continue;
        double RelativePhi = 0.0;
        for(int i=0;i<BPE;i++)
        {
            double tPhi0 = _e->NeB[i]->dPhi[j][0] * FlowWeight[i];
            if(tPhi0 < 0.0)
            {
                OutTable[i]  = 1;
                RelativePhi += tPhi0;
            }
        }

        double tSumPhi0 = RelativePhi + _e->sPhiL3[j][0]*_e->Volume;

        Outflow[j] = RelativePhi;
        Remnant[j] = tSumPhi0;

        if(tSumPhi0 < 0.0 && j!=0)
        {
            // restrict the RelativePhi to sPhiL3[j][0]
            DirtyBit = 1;
#if ELEMENT_DEBUG
            _e->Debug += 1.0*(1<<j);
#endif
        }
    }

    if(0==DirtyBit) return;

    double OutRatio[NMT][BPE] = {0.0};
    for(int j=0;j<nm;j++)
    {
        if((_e->sPhiL3[j][0] < TOLVOF) || (fabs(Outflow[j]) < TOLVOF*_e->Volume))
            ZeroEx(OutRatio[j],BPE);
        else
        {
            for(int i=0;i<BPE;i++)
            {
                if(OutTable[i])
                {
                    OutRatio[j][i] = FlowWeight[i] * _e->NeB[i]->dPhi[j][0]/Outflow[j];
                } else
                    OutRatio[j][i] = 0.0;
            }
        }
    }

    int Sht[NMT];
    double Adjusted[NMT]={0.0};
    Shuffle(Sht,nm);
    int Incr = MaterialStep(0,nm,Remnant,Sht,0); // over-empty material index
    int Decr = MaterialStep(0,nm,Remnant,Sht,1); // candidate to replace over-empty space

    while(Incr<nm && Decr <nm)
    {
        double IncrPart = -(Remnant[Sht[Incr]] + Adjusted[Sht[Incr]]);
        double DecrPart =   Remnant[Sht[Decr]] + Adjusted[Sht[Decr]];
        if(DecrPart >= IncrPart)
        {
            Adjusted[Sht[Incr]] += IncrPart;
            Adjusted[Sht[Decr]] -= IncrPart;
            MaterialSwap(Sht[Incr],Sht[Decr],IncrPart,_e,OutRatio[Sht[Incr]]);
            Incr = MaterialStep(Incr+1,nm,Remnant,Sht,0);
        }
        else
        {
            Adjusted[Sht[Incr]] +=  DecrPart;
            Adjusted[Sht[Decr]] -=  DecrPart;
            MaterialSwap(Sht[Incr],Sht[Decr],DecrPart,_e,OutRatio[Sht[Incr]]);
            Decr = MaterialStep(Decr+1,nm,Remnant,Sht,1);
        }
    }

    // compress the element

    // recalculate flow on boundary
    for(int i=0;i<BPE;i++)
    {
        if(0==OutTable[i]) continue;


        double bMass = 0.0;
        for(int j=0;j<nm;j++)
        {
            if(_e->NeB[i]->dPhi[j][0] * FlowWeight[i] > 0.0)
            {
                fprintf(stdout,"Error in Adjust!\n");
                exit(0);
            }
            _e->NeB[i]->dPhi[j][1] = _e->PhiL3[j][1] * _e->NeB[i]->dPhi[j][0];
            _e->NeB[i]->dPhi[j][2] = _e->PhiL3[j][2] * _e->NeB[i]->dPhi[j][0];
            bMass += _e->NeB[i]->dPhi[j][1];
        }

        _e->NeB[i]->dMass = bMass;
        ScalerMoveEx(_e->NeB[i]->dMomentum,_e->Momentum,bMass/_e->Mass,DIM);
        ScalerMoveEx(_e->NeB[i]->dPsi,_e->PsiL3,bMass,NPSI);

        /*double sfrac = bMass/_e->NeB[i]->dMass;
        ScalerEx(_e->NeB[i]->dMomentum,sfrac,DIM);
        ScalerEx(_e->NeB[i]->dPsi,sfrac,NPSI);
        _e->NeB[i]->dMass = bMass;*/
    }

}

void CalculateL3bIn(struct Element * _e, int nm)
{

    /*
     * The information of momentum or velocity in elements are
     * stored in boundary.
     * here, _e->Momentum is cleared and its definition is changed]
     * to the increment of momentum because of advection.
     */
    int InTable[BPE] = {0};
    for(int i=0;i<BPE;i++)
    {
        if(FlowWeight[i]*_e->NeB[i]->dMass < 0.0) continue;
        InTable[i] = 1;
        _e->Mass += FlowWeight[i]*_e->NeB[i]->dMass;
        ScalerAddition(_e->Momentum,_e->NeB[i]->dMomentum,FlowWeight[i]);
        ScalerAdditionEx(_e->sPsiL3,_e->NeB[i]->dPsi,FlowWeight[i],NPSI);
    }

    double Zeta = 0.0;
    for(int i = 0; i < nm; i++)
    {
        for (int j = 0; j < BPE; j++)
        {
            if(InTable[j]==0) continue;
            ScalerAdditionEx(_e->sPhiL3[i], _e->NeB[j]->dPhi[i], FlowWeight[j] * _e->rVolume, NPHI);
        }

        _e->sPhiL3[i][0] = Wind(_e->sPhiL3[i][0],0.0,1.0);
        if(i != 0)
            Zeta += _e->sPhiL3[i][0];
    }

    /* this is similar with sale2d
    Compress or expand each material proportionally with
    volume fraction of the material.  In essence, this assumes
    the compressibility of each material is the same. . . */

    if((Zeta > 1.0-TOLVOF) || (_e->sPhiL3[0][0] < TOLVOF))
    {
        _e->sPhiL3[0][0] = 0.0;
        for (int i = 1; i < nm; i++)
            _e->sPhiL3[i][0] = _e->sPhiL3[i][0] / Zeta ;
    }
    else
    {
        _e->sPhiL3[0][0] = 1.0 - Zeta;
    }

    /*
     * for some extreme cases, this method is not valid
     * i.e, the volume Phi[k][0] limits to 0 when residual mass Phi[k][1] doesn't limit to 0.
     * the density will be infinity.
     */
    double tMass = 0.0;
    for(int i = 1; i < nm; i++)
    {
        if (_e->sPhiL3[i][0] > TOLVOF)
        {
            _e->PhiL3[i][0] = 1.0;
            _e->PhiL3[i][1] = _e->sPhiL3[i][1] / _e->sPhiL3[i][0];
            _e->PhiL3[i][2] = _e->sPhiL3[i][2] / _e->sPhiL3[i][0];
        }
        else
        {
            _e->PhiL3[i][0] = 1.0;
            _e->PhiL3[i][1] = 0.0;
            _e->PhiL3[i][2] = 0.0;
            ZeroEx(_e->sPhiL3[i],NPHI);
        }
        _e->VOF[i]  = _e->sPhiL3[i][0];
        tMass      += _e->sPhiL3[i][1];
    }
    tMass *= _e->Volume;

    if(tMass < TOLRHO*_e->Volume)
    {
        for(int k=0;k<nm;k++)
        {
            ZeroEx(_e->sPhiL3[k],NPHI);
            ZeroEx(_e->PhiL3[k],NPHI);
            _e->PhiL3[k][0] = 1.0;
        }
        Zero(_e->Momentum);
        ZeroEx(_e->sPsiL3,NPSI);
        ZeroEx(_e->PsiL3,NPSI);
        ZeroEx(_e->VOF,NMT);
        _e->sPhiL3[0][0] = 1.0;
        _e->VOF[0] = 1.0;
        _e->Mass = 0.0;
        // need improve !!
    }
    else
    {
        _e->Mass = tMass;
        _e->VOF[0] = _e->sPhiL3[0][0];
        ScalerMoveEx(_e->PsiL3,_e->sPsiL3,1.0/_e->Mass,NPSI);
    }

    // the average density is sum by weight of materials
    // Sync-Backward
    _e->Density = _e->Mass * _e->rVolume;
    // Copy Psi[0] to damage
    _e->PsiL3[0]  = Wind(_e->PsiL3[0],0.0,1.0);
    _e->sPsiL3[0] = _e->PsiL3[0] * _e->Mass;
    _e->Damage    = _e->PsiL3[0];

    // stress dev, is stored in sPsiL3[1:6]
}


void ANEOSCalL3c(struct Element *_e, const struct ANEOSTable _a[], const struct StateReference _s[], int nm)
{
    /*
     * this subroutine is used to calculate the pressure/module
     * and other state-related parameters that are updated implicitly
     * need to be completed
     */

    // for the vacuum/void ThetaL3 = 0.0
    ZeroEx(_e->ThetaL3[0], NTHETA);

    for (int i = 1; i < nm; i++)
    {
        if (_e->VOF[i] > TOLVOF)
        {
            ANEOSCalL3t(_e, &_a[i], i);
            /*
             * Check the min_pressure and the cs_sound
             */
            _e->ThetaL3[i][1] = Max(_e->ThetaL3[i][1],_s[i].minPint*(1.0-_e->Damage) + _s[i].minPdam*_e->Damage);
            _e->ThetaL3[i][2] = Max(1.0,_e->ThetaL3[i][2]);

        }
        else
        {
//            _e->VOF[i] = 0.0;
            ZeroEx(_e->ThetaL3[i], NTHETA);
        }
    }

    /*
     * pressure and temperature (+ Csound) is weighted by VOF
     * The Shear module is the harmonic average weighted by VOF
     */
    _e->Temperature = 0.0;
    _e->Pressure    = 0.0;
    _e->Csound      = 0.0;
    _e->ShearModule = 0.0;


    if(_e->VOF[0] < TOLVOF)
    {
        for (int i = 1; i < nm; i++)
        {
            _e->Temperature += _e->VOF[i] * _e->ThetaL3[i][0];
            _e->Pressure    += _e->VOF[i] * _e->ThetaL3[i][1];
            double tK = _s[i].GaCoef * _e->PhiL3[i][1] * _e->ThetaL3[i][2] * _e->ThetaL3[i][2];
            if(tK > TOL)
                _e->ShearModule += _e->VOF[i] / tK;
            if(_e->ThetaL3[i][2] > TOL)
                _e->Csound += _e->VOF[i]/_e->ThetaL3[i][2];
        }
        if(_e->ShearModule > 0.0)
            _e->ShearModule = 1.0 / _e->ShearModule;
        if(_e->Csound > 0.0)
            _e->Csound = 1.0/_e->Csound;
    }
    else
    {
        double tV = 0.0;
        for(int i=1;i<nm;i++) tV += _e->VOF[i];
        for (int i = 1; i < nm; i++)
        {
            _e->Temperature += _e->VOF[i] * _e->ThetaL3[i][0];
//            _e->Pressure    += _e->VOF[i] * _e->ThetaL3[i][1];
            double tK = _s[i].GaCoef * _e->PhiL3[i][1] * _e->ThetaL3[i][2] * _e->ThetaL3[i][2];
            if(tK > TOL)
                _e->ShearModule += _e->VOF[i] / tK;
            if(_e->ThetaL3[i][2] > TOL)
                _e->Csound += _e->VOF[i]/_e->ThetaL3[i][2];
        }
        if(_e->ShearModule > 0.0)
            _e->ShearModule = tV / _e->ShearModule;
        if(_e->Csound > 0.0)
            _e->Csound = tV /_e->Csound;

        _e->Pressure = 0.0;
    }
}


void CalculateL3c(struct Element *_e, const struct EosTable _et[], const struct StateReference _s[],const struct MeshInfo * _minfo)
{
    /*
     * this subroutine is used to calculate the pressure/module
     * and other state-related parameters that are updated implicitly
     * Add support for tillotson Eos
     */

    // for the vacuum/void ThetaL3 = 0.0
    ZeroEx(_e->ThetaL3[0], NTHETA);

    int nm = _minfo->MaterialNum;
    for (int i = 1; i < nm; i++)
    {
        if (_e->VOF[i] > TOLVOF)
        {
            _et[i].CalL3t(_e, &_et[i],i);
        }
        else
        {
            ZeroEx(_e->ThetaL3[i], NTHETA);
        }
        /*
        * Check the min_pressure and the cs_sound
        */
        _e->ThetaL3[i][1] = Max(_e->ThetaL3[i][1],_s[i].minPint*(1.0-_e->Damage) + _s[i].minPdam*_e->Damage);
        _e->ThetaL3[i][2] = Max(1.0,_e->ThetaL3[i][2]);
    }

    /*
     * pressure and temperature (+ Csound) is weighted by VOF
     * The Shear module is the harmonic average weighted by VOF
     */


    if((1==_minfo->PartPressure && _e->VOF[0]>= 1.0-TOLVOF)||(0==_minfo->PartPressure && _e->VOF[0]>=TOLVOF))
    {
        _e->Temperature = 0.0;
        _e->Pressure    = 0.0;
        _e->Csound      = 1.0;
        _e->ShearModule = 0.0;
    }
    else if(_e->VOF[0] > TOLVOF)
    {
        _e->Temperature = 0.0;
        _e->Pressure    = 0.0;

        double tV = 1.0 - _e->VOF[0];

        for(int i = 1; i < nm; i++)
        {
            _e->Temperature += _e->VOF[i] * _e->ThetaL3[i][0];
            _e->Pressure    += _e->VOF[i] * _e->ThetaL3[i][1];
        }
        _e->Temperature = _e->Temperature/tV;
        _e->Csound = 1.0;
        _e->ShearModule = 0.45 * _e->Density;
    }
    else
    {
        _e->Temperature = 0.0;
        _e->Pressure    = 0.0;
        _e->Csound      = 0.0;
        _e->ShearModule = 0.0;

        for (int i = 1; i < nm; i++)
        {
            _e->Temperature += _e->VOF[i] * _e->ThetaL3[i][0];
            _e->Pressure    += _e->VOF[i] * _e->ThetaL3[i][1];
            double tK = _s[i].GaCoef * _e->PhiL3[i][1] * _e->ThetaL3[i][2] * _e->ThetaL3[i][2];
            if(tK > TOL)
                _e->ShearModule += _e->VOF[i] / tK;
            if(_e->ThetaL3[i][2] > TOL)
                _e->Csound += _e->VOF[i]/_e->ThetaL3[i][2];
        }
        if(_e->ShearModule > 0.0)
            _e->ShearModule = 1.0 / _e->ShearModule;
        if(_e->Csound > 0.0)
            _e->Csound = 1.0/_e->Csound;
    }
}

void CalculateL3c_gas(struct Element *_e, const struct AirEOSTable _a[], int nm)
{
    /*
     * this subroutine is used to calculate the pressure/module
     * and other state-related parameters that are updated implicitly
     * need to be completed
     */

    // for the vacuum/void ThetaL3 = 0.0
    ZeroEx(_e->ThetaL3[0], NTHETA);

    for (int i = 1; i < nm; i++)
    {
        if (_e->VOF[i] > TOLVOF)
        {
            CalculateL3atm(_e,&_a[i],i);
        }
        else
        {
            _e->VOF[i] = 0.0;
            ZeroEx(_e->ThetaL3[i], NTHETA);
        }
    }

    /*
     * pressure and temperature (+ Csound) is weighted by VOF
     * The Shear module is the harmonic average weighted by VOF
     */
    _e->Temperature = 0.0;
    _e->Pressure    = 0.0;
    _e->Csound      = 0.0;
    _e->ShearModule = 0.0;

    if(_e->VOF[0] > 1.0 - TOLVOF)
    {
        _e->Temperature = 0.0;
        _e->Pressure    = 0.0;
        _e->Csound      = 0.0;
        _e->ShearModule = 0.0;
    }
    else
    {
        for (int i = 1; i < nm; i++)
        {
            _e->Temperature += _e->VOF[i] * _e->ThetaL3[i][0];
            _e->Pressure    += _e->VOF[i] * _e->ThetaL3[i][1];
            _e->Csound      += _e->VOF[i] * _e->ThetaL3[i][2];
        }

        _e->Temperature /= 1.0 - _e->VOF[0];
        _e->Pressure    /= 1.0 - _e->VOF[0];
        _e->Csound      /= 1.0 - _e->VOF[0];
        _e->ShearModule /= 1.0 - _e->VOF[0];

        if(_e->VOF[0] > TOLVOF || _e->Mass < TOLMASS)
        {
            _e->Pressure    = 0.0;
            _e->ShearModule = 0.0;
        }
    }

}

// version of not vectorized operation

void MapE2Vb(struct Vertex * iv)
{
    iv->CsoundL1     = 0.0;
    iv->State = 0;
    double TotalMass = 0.0;
    double deltaMomentum[DIM] = {0.0};
    int VoidVertex = 0;
    for(int j=0;j<EPV;j++)
    {
        if(iv->NeE[j]->VOF[0]>1.0-TOLVOF)
        {
            VoidVertex += 1;
        }

        TotalMass    += iv->NeE[j]->Mass * iv->NeE[j]->CutRatio[PCVE[j]];
        iv->CsoundL1 += iv->NeE[j]->Mass * iv->NeE[j]->CutRatio[PCVE[j]] * iv->NeE[j]->Csound;
        ScalerAddition(deltaMomentum,iv->NeE[j]->Momentum,iv->NeE[j]->CutRatio[PCVE[j]]);
    }

    if((TotalMass < TOLRHO * iv->Volume) || (iv->VOF[0] > 1.0 - TOLVOF))
    {
        Zero(iv->Velocity);
        iv->MassL1   = 0.0;
        iv->CsoundL1 = 0.0;
        iv->State   = 1;
    }
    else
    {
        ScalerMove(iv->Velocity,deltaMomentum,1.0/TotalMass);
        iv->MassL1    = TotalMass;
        iv->CsoundL1 /= TotalMass;
        iv->State    = 0;
    }

    if(VoidVertex >= 1 && 0==iv->State)
    {
        iv->State = 2;
    }

}

void MapV2Eb(struct Element * ie, const struct MeshInfo *_info)
{
    /*
     * a weaken version of MapV2Ea which only process the velocity
     * NOTE: Velocity is updated in CalculateL1
     * VOF and other variables should not be changed !!
     * When all the variables in Vertex is given
     * map the variables into element.
     * This subroutine utilizes PhiL3
     * Within an element, the velocity for different materials is consistent,
     * then the momentum stored in PhiL3 is removed.
     */

    /*
     * New version 0409
     * ie->Velocity is removed
     * ie->Momentum is the momentum in this element.
     */
    Zero(ie->Momentum);
    for(int j=0;j<VPE;j++)
    {
        ScalerAddition(ie->Momentum,ie->NeV[j]->Velocity,ie->NeV[j]->MassL1 * ie->NeV[j]->SubVolumeRatio[PCVE[j]]);
    }

    // At the Last, whether to crop the velocity ??
    // In ths subroutine, (similar to iSALE2d-Dellen) restrict the velocity in a interval
    // to avoid some surprise values.
//    double nVel = Length(ie->Velocity);
//    if(nVel > _info->MaxVelocity)
//        ScalerMove(ie->Velocity,ie->Velocity,_info->MaxVelocity/nVel);

}

void VelocityUpdate(struct Vertex * _vertex)
{
    // __name: CalculateL1a --> VelocityUpdate
    // calculate the acceleration of a vertex and
    // store it in VelocityL1

    if(1 == _vertex->State)
    {
        Zero(_vertex->VelocityL1);
    } else if(2== _vertex->State)
    {
        Copy(_vertex->BodyF,_vertex->VelocityL1);
    }
    else
    {
        double VelocityLt1[DIM] = {0.0};
        double VelocityLt2[DIM] = {0.0};
        for (int k = 0; k < EPV; k++)
        {
            struct Element *kNeE = _vertex->NeE[k];
            double * kSubFace = kNeE->SubFace[PCVE[k]];
            MatrixTimeSym(VelocityLt2,kNeE->PsiL3+1,kSubFace);
            ScalerAddition(VelocityLt2, kSubFace, -1.0 * (kNeE->Pressure + kNeE->ArtificialPressure));
            Addition(VelocityLt1, VelocityLt2);
        }
        for(int k=0;k<DIM;k++)
            _vertex->VelocityL1[k] = VelocityLt1[k]/_vertex->MassL1 + _vertex->BodyF[k];
        if(_vertex->Damping>0.0)
            ScalerAddition(_vertex->VelocityL1,_vertex->Velocity,-1.0*_vertex->Damping);
    }
}

double CalculateL1b(struct Vertex * iv, double dt0, double MaxVel)
{
    double dt = 1.1 * dt0;
    if(iv->State) return dt;
    double Vel = Length(iv->Velocity);
    if(Vel > MaxVel)
    {
        Scaler(iv->Velocity,MaxVel/Vel);
        Vel = MaxVel;
    }

    double Cwave = Vel + iv->CsoundL1;
    if(Cwave > TOLVEL)
        dt = Min(COURANT1*iv->dx_min/Cwave,dt);
    return dt;
}


void CalculateL1c(struct Vertex * _vertex, double dt)
{
    ScalerMove(_vertex->DistanceL1,_vertex->Velocity,0.5*dt);
    ScalerAddition(_vertex->Velocity,_vertex->VelocityL1,dt);
    ScalerAddition(_vertex->DistanceL1,_vertex->Velocity,0.5*dt);
}

void CalculateL1cABF(struct Vertex * _vertex, double dt ,const double ABFAn[])
{
    ScalerMove(_vertex->DistanceL1,_vertex->Velocity,0.5*dt);
    for(int i=0;i<DIM;i++)
    {
        _vertex->Velocity[i] += ABFAn[OABF-1] * _vertex->VelocityL1[i] * dt;
        for(int j=0;j<OABF-1;j++)
            _vertex->Velocity[i] += dt*ABFAn[OABF - 2 - j]*_vertex->MtVelocityL1[j][i];
    }
    ScalerAddition(_vertex->DistanceL1,_vertex->Velocity,0.5*dt);
    for(int k=OABF-2;k>0;k--)
        Copy(_vertex->MtVelocityL1[k-1],_vertex->MtVelocityL1[k]);
    Copy(_vertex->VelocityL1,_vertex->MtVelocityL1[0]);
}

void CalculateABFAn(const double dt[], double ABFAn[])
{
    if(OABF==3)
    {
        double a = dt[2], b = dt[1], c = dt[0];
        ABFAn[2] = (2.0*c*c + 3.0*c*(a+2.0*b) + 6.0*b*(a+b))/(6.0*b*(a + b));
        ABFAn[1] = -1.0*c*(3.0*a + 3.0*b + 2.0*c)/(6.0 * a * b);
        ABFAn[0] = c*(3.0*b + 2.0*c)/(6.0*a*(a+b));
    }
    else if(OABF==1)
    {
        ABFAn[0] = 1.0;
    }
    // other order of Adams-Bashforth method is not implemented.
}


double vofMSTACS(double nVOFD, double Gf, double CourantD)
{
    double nVOFCDSMSTACS,nVOFfSTOIC;
    if(CourantD < 1.0/3.0)
        nVOFCDSMSTACS = Min(nVOFD/CourantD,1.0);
    else
        nVOFCDSMSTACS = Min(3.0*nVOFD,1.0);

    if(nVOFD < 1.0/5.0)
        nVOFfSTOIC = 3.0 * nVOFD;
    else if(nVOFD < 1.0/2.0)
        nVOFfSTOIC = (1.0 + nVOFD)/2.0;
    else
        nVOFfSTOIC = Min((6.0*nVOFD + 3.0)/8.0,1.0);
    return Gf*nVOFCDSMSTACS + (1.0 - Gf)*nVOFfSTOIC;
}

double vofCICSAM(double nVOFD, double Gf, double CourantD)
{
    double nVOFCBC, nVOFUQ;
    nVOFCBC = Min(nVOFD / CourantD, 1.0);
    nVOFUQ  = Min(CourantD * nVOFD + (1.0 - CourantD) * (6.0 * nVOFD + 3.0) / 8.0, nVOFCBC);
    return Gf * nVOFCBC + (1.0 - Gf) * nVOFUQ;
}

double vofSTACS(double nVOFD, double Gf)
{
    double nVOFfSUP,nVOFfSTOIC;
    if(nVOFD < 1.0/3.0)
        nVOFfSUP = 2.0*nVOFD;
    else if(nVOFD < 1.0/2.0)
        nVOFfSUP = (1.0 + nVOFD)/2.0;
    else
        nVOFfSUP = Min(3.0*nVOFD/2.0,1.0);
    if(nVOFD < 1.0/5.0)
        nVOFfSTOIC = 3.0 * nVOFD;
    else if(nVOFD < 1.0/2.0)
        nVOFfSTOIC = (1.0 + nVOFD)/2.0;
    else
        nVOFfSTOIC = Min((6.0*nVOFD + 3.0)/8.0,1.0);

    return Gf*nVOFfSUP + (1.0 - Gf)*nVOFfSTOIC;
}

double vofMHRIC(double nVOFD, double Gf)
{
    // MHRIC
    double nVOFs1,nVOFs2;
    nVOFs1 = Min(2.0*nVOFD,1.0);
    nVOFs2 = Min((6.0*nVOFD + 3.0)/8.0,nVOFs1);
    return Gf*nVOFs1 + (1.0 - Gf)*nVOFs2;
}