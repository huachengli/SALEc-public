//
// Created by huacheng on 6/23/21.
//

#include "Variable.h"
#include "Setup.h"
#include "Iteration.h"
#include "Ghost.h"
#include "Parallelism.h"
#include "Reconstruction.h"
#include "InputParser.h"
#include "Tracer.h"


int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);

    InputFile * ifp = OpenInputFile("../SALEc.inp");

    struct Mesh pMesh;
    struct MeshInfo pMeshInfo;
    struct ProcessInfo procInfo;
    TracerInfo pTracerInfo;

    SetMeshInfo(ifp,&pMeshInfo,&procInfo,&pMesh);

    struct EosTable ETB[NMT];
    struct StateReference SRF[NMT];
    EosInit(&pMeshInfo,ETB,SRF);
    LoadStateRef(ifp,&pMeshInfo,SRF);

    GetMeshGeometry(&pMesh, &pMeshInfo);
    BuildMeshGeometry(&pMesh, &pMeshInfo);

    Calculate_SV(&pMesh);
    ParaSearchDomainV(&pMesh,&procInfo);
    ParaSearchDomainE(&pMesh,&procInfo);
    ParaSortDomainE(&pMesh,&procInfo);
    InitProcessInfo(&pMesh,&procInfo);
    ParaPackRecvSendInit(&procInfo);


    InitMesh(ifp,&pMesh,ETB);
    SetVelStavlity(&pMesh,TestVelStavlity);

#ifdef SALEC_USE_CHECKPOINT
    InitChk(ifp,&pMesh,&procInfo,&pMeshInfo);
#endif

    for(int i=0;i<pMesh.ndb;i++) // bid is used in interface reconstruction
        InsertBoundaryBID(pMesh.DomainB[i]);

    SetBC(ifp,&pMesh);
    SetDamping(&pMeshInfo,&pMesh);

#ifdef SALEC_USE_TRACER
    InitTracerInfo(&pTracerInfo,&pMeshInfo);
    InitTracer(&pTracerInfo,ifp,&pMesh);
#endif

    char _prefix[MaxStrLen];
    GetValueS(ifp,"output.prefix",_prefix,"ParaTest");
    double StepTime  = GetValueD(ifp,"cycle.time0","0.0");
    double StartTime = StepTime;
    double EndTime   = GetValueD(ifp,"cycle.time1","5.0");
    double dt        = GetValueD(ifp,"cycle.dt0","1.0e-3");
    double dt_max    = GetValueD(ifp,"cycle.maxdt","1.0e-2");
    double OutTime   = GetValueD(ifp,"output.time","1.0");

    int  MaxSteps = GetValueI(ifp,"cycle.maxstep","100");
    int  InitStep = GetValueI(ifp,"output.init","0");
    int  kOutStep = 0;

    CloseInputFile(ifp);
    MPI_Barrier(SALEC_CART_COMM);

    double st = MPI_Wtime();

#ifdef MPI_DEBUG_HOOK
    int local_rank;
    MPI_Comm_rank(SALEC_CART_COMM,&local_rank);
    int debugvalue = 1;
    while(debugvalue && (0==local_rank))
        sleep(5);
#endif

    for(int  iStep = 0;iStep<=MaxSteps && StepTime<= EndTime + dt_max;iStep++)
    {

        
        for(int i=pMesh.nde1;i<pMesh.nde;i++)
        {
            CalculateL3bOut(pMesh.DomainE[i],pMeshInfo.MaterialNum);
        }

        
        for(int i=pMesh.nde1;i<pMesh.nde;i++)
        {
            CalculateL3b(pMesh.DomainE[i],&pMeshInfo);
            CalculateL3c(pMesh.DomainE[i], ETB, SRF, &pMeshInfo);
        }
        ParaPackBegin(&procInfo,IdPres);
        ParaPackBegin(&procInfo,IdMass);
        ParaPackBegin(&procInfo,IdVOF);
        ParaPackBegin(&procInfo,IdCsound);
        ParaPackBegin(&procInfo,IdMomentum);

#ifdef SALEC_USE_TRACER
        TracerMove(&pTracerInfo,&pMesh,0.5*dt);
        TracerInfoExchangeBegin(&pTracerInfo,&procInfo,&pMesh);
#endif

        for(int i=0;i<pMesh.nde1;i++)
        {
            CalculateL3bOut(pMesh.DomainE[i],pMeshInfo.MaterialNum);
        }

        
        for(int i=0;i<pMesh.nde1;i++)
        {
            CalculateL3b(pMesh.DomainE[i],&pMeshInfo);
            CalculateL3c(pMesh.DomainE[i], ETB, SRF, &pMeshInfo);
        }

        
        for(int i=0;i<pMesh.ndv1;i++)
        {
            MapE2Vc(pMesh.DomainV[i]);
        }
        ParaPackEnd(&procInfo,IdVOF);
        FlushCondEId(&pMesh,IdFCEL1);

#ifdef SALEC_USE_TRACER
        TracerInfoExchangeEnd(&pTracerInfo,&procInfo);
        TracerExchangeBegin(&pTracerInfo,&procInfo,&pMesh);
#endif

        for(int i=pMesh.ndv1;i<pMesh.ndv;i++)
        {
            MapE2Vc(pMesh.DomainV[i]);
        }

        
        for(int i=pMesh.nde1;i<pMesh.nde;i++)
        {
            CutGradVOF(pMesh.DomainE[i],pMeshInfo.MaterialNum);
            ConstructPlane(pMesh.DomainE[i]);
        }
        ParaPackBegin(&procInfo,IdGradVOF);

        for(int i=0;i<pMesh.nde1;i++)
        {
            CutGradVOF(pMesh.DomainE[i],pMeshInfo.MaterialNum);
            ConstructPlane(pMesh.DomainE[i]);
        }
        ParaPackEnd(&procInfo,IdGradVOF);
        PostPackEnd(&procInfo,IdGradVOF);
        FlushCondEId(&pMesh,IdFCEL2);
        
        for(int i=0;i<pMesh.ndv1;i++)
        {
            MapE2Vb(pMesh.DomainV[i]);
        }
        ParaPackEnd(&procInfo,IdMass);
        ParaPackEnd(&procInfo,IdCsound);
        ParaPackEnd(&procInfo,IdMomentum);
        FlushCondEId(&pMesh,IdFCEL3);

        for(int i=pMesh.ndv1;i<pMesh.ndv;i++)
        {
            MapE2Vb(pMesh.DomainV[i]);
        }

        /* ******************************************
         * the "real" end of a calculation cycle,
         * the output or check point should be added here
         * *******************************************
         */
#ifdef MPI_DEBUG_HOOK
        MPI_Barrier(SALEC_CART_COMM);
#endif

        if((fabs(StepTime-kOutStep*OutTime - StartTime) <= 0.5*dt_max) || (iStep<=InitStep))
        {
            char fname[200];
            sprintf(fname,"%s.proc%d.%d.vts",_prefix,procInfo.rank,kOutStep);
            TestMesh(&pMesh,fname);
            if(procInfo.rank==0)
            {
                fprintf(stdout,"[*%3d]Step:%04d, dt = %1.10e; t = %5.5e\n",kOutStep,iStep,dt,StepTime);
            }

#ifdef SALEC_USE_TRACER
            OutputTracer(&pTracerInfo,&pMesh,kOutStep,procInfo.rank);
#endif

            if(fabs(StepTime-kOutStep*OutTime - StartTime) > 0.5*dt_max)
                StartTime -= OutTime;

            kOutStep++;
        } else if(0==iStep%100)
        {
            if(procInfo.rank==0)
            {
                fprintf(stdout,"[    ]Step:%04d, dt = %1.10e; t = %5.5e\n",iStep,dt,StepTime);
            }
        }

#ifdef SALEC_USE_CHECKPOINT
        if((iStep%pMeshInfo.chkStep==0) && (pMeshInfo.chkStep>0))
        {
            UpdateChk(&pMesh,&pMeshInfo,&procInfo,iStep);
//            CompareChk(&pMesh,&pMeshInfo,&procInfo,iStep);
        }
#endif

#ifdef SALEC_USE_PERFMONITER
 // performance monitor between progress

        char t_message[400],t_pname[100];
        int t_length;
        MPI_Get_processor_name(t_pname,&t_length);

        if(procInfo.rank ==0 && iStep%100 == 0)
        {
            sprintf(t_message,"before-barrier:%d:rdt=%f on %s\n",procInfo.rank,MPI_Wtime()-st,t_pname);
            RunPoint_s(stdout,t_message);
        }
        MPI_Barrier(SALEC_CART_COMM);

        if(procInfo.rank ==0 && iStep%100 == 0)
        {
            RunPoint_s(stdout,"s_mpi:barrier-allreduce");
            fprintf(stdout,"rdt=%f\n",MPI_Wtime()-st);
        }
#endif
        /*
         * ******************************************************
         * start a new cycle
         * ******************************************************
         */
        double tmpdt = dt_max;
        for(int i=0;i<pMesh.ndv;i++)
            tmpdt = Min(CalculateL1b(pMesh.DomainV[i],dt,pMeshInfo.MaxVelocity),tmpdt);
        MPI_Allreduce(&tmpdt, &dt, 1, MPI_DOUBLE, MPI_MIN, SALEC_CART_COMM);

#ifdef SALEC_USE_PERFMONITER
        if(procInfo.rank ==0 && iStep%100 == 0)
        {
            RunPoint_s(stdout,"s_mpi:complete-allreduce");
            fprintf(stdout,"rdt=%f\n",MPI_Wtime()-st);
        }
#endif

#ifdef MPI_DEBUG_HOOK
        MPI_Barrier(SALEC_CART_COMM);
#endif

#ifdef SALEC_USE_TRACER
        TracerExchangeEnd(&pTracerInfo,&procInfo,&pMesh);
        TracerMove(&pTracerInfo,&pMesh, 0.5*dt);
#endif
        for(int i=pMesh.nde1;i<pMesh.nde;i++)
        {
            CalculateGradM(pMesh.DomainE[i],pMeshInfo.MaterialNum);
            RheologyUpdate(pMesh.DomainE[i],dt,SRF,pMeshInfo.MaterialNum);
            BlockVibration(pMesh.DomainE[i],dt,SRF,pMeshInfo.MaterialNum,StepTime);
            FailureUpdate(pMesh.DomainE[i],dt,SRF,&pMeshInfo);
            CalculateABV2(pMesh.DomainE[i],pMeshInfo.ARTVIS,pMeshInfo.ARTVIS2);
        }
        ParaPackBegin(&procInfo,IdPresABV);
        ParaPackBegin(&procInfo,IdPsiL3);

        
        for(int i=0;i<pMesh.nde1;i++)
        {
            CalculateGradM(pMesh.DomainE[i],pMeshInfo.MaterialNum);
            RheologyUpdate(pMesh.DomainE[i],dt,SRF,pMeshInfo.MaterialNum);
            BlockVibration(pMesh.DomainE[i],dt,SRF,pMeshInfo.MaterialNum,StepTime);
            FailureUpdate(pMesh.DomainE[i],dt,SRF,&pMeshInfo);
            CalculateABV2(pMesh.DomainE[i],pMeshInfo.ARTVIS,pMeshInfo.ARTVIS2);
        }

        
        for(int i=0;i<pMesh.ndv1;i++)
        {
            VelocityUpdate(pMesh.DomainV[i]);
            CalculateL1c(pMesh.DomainV[i],dt);
        }
        ParaPackEnd(&procInfo,IdPresABV);
        ParaPackEnd(&procInfo,IdPres);
        ParaPackEnd(&procInfo,IdPsiL3);
        FlushCondEId(&pMesh,IdFCEL4);

        
        for(int i=pMesh.ndv1;i<pMesh.ndv;i++)
        {
            VelocityUpdate(pMesh.DomainV[i]);
            CalculateL1c(pMesh.DomainV[i],dt);
        }

        
        for(int i=0;i<pMesh.ndb;i++)
            CalculateBFV(pMesh.DomainB[i]);

        for(int k=pMesh.nde1;k<pMesh.nde;k++)
        {
            CalculateCourant(pMesh.DomainE[k]);
            EnergyUpdate(pMesh.DomainE[k], &pMeshInfo, dt);
        }
        ParaPackBegin(&procInfo,IdCourant);
        ParaPackBegin(&procInfo,IdMomentum);
        ParaPackBegin(&procInfo,IdPhiL3);

        
        for(int k=0;k<pMesh.nde1;k++)
        {
            CalculateCourant(pMesh.DomainE[k]);
            EnergyUpdate(pMesh.DomainE[k], &pMeshInfo, dt);
        }

        ParaPackEnd(&procInfo,IdPhiL3);
        ParaPackEnd(&procInfo,IdMomentum);
        ParaPackEnd(&procInfo,IdCourant);
        FlushCondEId(&pMesh,IdFCEL5);
        FlushVnode(&pMesh);

        for(int i=0;i<pMesh.ndb;i++)
        {
            CalculateVOFM(pMesh.DomainB[i],pMeshInfo.CISSAMK1,pMeshInfo.MaterialNum);
        }

        StepTime += dt;
    }

    MPI_Barrier(SALEC_CART_COMM);
    double ed = MPI_Wtime();
    if(procInfo.rank==0)
    {
        fprintf(stdout,"Time used in MAIN LOOP: %8.7e\n",ed-st);
        fprintf(stdout,"Writing VTM files \n");
        for(int k=0;k<kOutStep;k++)
            write_vtm(procInfo.npgs,k,_prefix);
    }

    fflush(stdout);
   /*
    DataClean(ATB,&pMesh,pMeshInfo.MaterialNum);
    InfoClean(&pMeshInfo,&procInfo);
   */

#ifdef SALEC_USE_TRACER
    TracerClean(&pTracerInfo);
#endif

    MPI_Barrier(SALEC_CART_COMM);
    MPI_Finalize();
    return 0;
}


