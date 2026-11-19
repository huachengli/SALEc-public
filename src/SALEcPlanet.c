//
// Created by huacheng on 11/8/22.
//

#include "Variable.h"
#include "Setup.h"
#include "Iteration.h"
#include "Ghost.h"
#include "Parallelism.h"
#include "Reconstruction.h"
#include "InputParser.h"
#include "Tracer.h"
#include "Planet.h"

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);

    // Yue
    char inpfile[256];
    if(argc == 1)
    {
        sprintf(inpfile,"../SALEc.inp");
    }
    else if(argc == 2)
    {
        sprintf(inpfile,"%s",argv[1]);
    }else
    {
        fprintf(stderr,"invalid args for %s\n",__FILE__);
        exit(0);
    }

    InputFile * ifp = OpenInputFile(inpfile);
    //  Yue
    // InputFile *ifp = OpenInputFile("../SALEc.inp");  Yue

    struct Mesh pMesh;
    struct MeshInfo pMeshInfo;
    struct ProcessInfo procInfo;
    TracerInfo pTracerInfo;

    SetMeshInfo(ifp, &pMeshInfo, &procInfo, &pMesh);
    SetOutputDir();

    struct EosTable ETB[NMT];
    struct StateReference SRF[NMT];
    EosInit(&pMeshInfo, ETB, SRF);
    LoadStateRef(ifp, &pMeshInfo, SRF);

    GetMeshGeometry(&pMesh, &pMeshInfo);
    BuildMeshGeometry(&pMesh, &pMeshInfo);

    Calculate_SV(&pMesh);
    ParaSearchDomainV(&pMesh,&procInfo);
    ParaSearchDomainE(&pMesh,&procInfo);
    ParaSortDomainE(&pMesh,&procInfo);
    InitProcessInfo(&pMesh,&procInfo);
    ParaPackRecvSendInit(&procInfo);

    // ParaDebug();
    InitPlanetMesh(ifp,&pMesh,ETB,SRF);
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

    struct CycleControl StepControl;
    LoadCycleControl(ifp,&StepControl);
    int real_step = 0;

    while(StepControl.cycle_continue)
    {
        CheckCycleControl(&StepControl);
        if(StepControl.out_vts || StepControl.out_init)
        {
            char fname[200];
            sprintf(fname,"vts/%s.proc%d.%d.vts",StepControl.prefix,procInfo.rank,StepControl.k_out_step);
            TestMesh(&pMesh,fname);
            if(procInfo.rank==0)
            {
                fprintf(stdout,"[*%3d]Step:%04d, dt = %1.10e; t = %5.5e\n",
                        StepControl.k_out_step,StepControl.i_step,StepControl.dt,StepControl.step_time);
                write_vtm(procInfo.npgs,StepControl.k_out_step,StepControl.prefix,pTracerInfo.Prefix);
            }
#ifdef SALEC_USE_TRACER
            OutputTracer(&pTracerInfo,&pMesh,StepControl.k_out_step,procInfo.rank);
            OutputTracerVtp(&pTracerInfo,&pMesh,StepControl.k_out_step,procInfo.rank);
            FlushEjecta(&pTracerInfo,StepControl.k_out_step,procInfo.rank);
#endif

            if(!StepControl.out_vts)
                StepControl.start_time -= StepControl.out_time;

            StepControl.k_out_step ++;
        } else if(StepControl.out_log)
        {
            if(procInfo.rank==0)
            {
                fprintf(stdout,"[    ]Step:%04d, dt = %1.10e; t = %5.5e\n",
                        StepControl.i_step,StepControl.dt,StepControl.step_time);
            }
        }
#ifdef SALEC_USE_TRACER
        OutputEjecta(&pTracerInfo,&pMesh,StepControl.step_time,procInfo.rank);
#endif

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
        TracerMove(&pTracerInfo,&pMesh,0.5*StepControl.dt);
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

#ifdef SALEC_USE_CHECKPOINT
        if(((StepControl.i_step)%pMeshInfo.chkStep==0) && (pMeshInfo.chkStep>0))
        {
            UpdateChk(&pMesh,&pMeshInfo,&procInfo,StepControl.i_step);
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
        double tmpdt = StepControl.max_dt;
        for(int i=0;i<pMesh.ndv;i++)
            tmpdt = Min(CalculateL1b(pMesh.DomainV[i],StepControl.dt,pMeshInfo.MaxVelocity),tmpdt);
        MPI_Allreduce(&tmpdt, &(StepControl.dt), 1, MPI_DOUBLE, MPI_MIN, SALEC_CART_COMM);

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
        TracerMove(&pTracerInfo,&pMesh, 0.5*StepControl.dt);
#endif
        for(int i=pMesh.nde1;i<pMesh.nde;i++)
        {
            CalculateGradM(pMesh.DomainE[i],pMeshInfo.MaterialNum);
            RheologyUpdate(pMesh.DomainE[i],StepControl.dt,SRF,pMeshInfo.MaterialNum);
            BlockVibration(pMesh.DomainE[i],StepControl.dt,SRF,pMeshInfo.MaterialNum,StepControl.step_time);
            FailureUpdate(pMesh.DomainE[i],StepControl.dt,SRF,&pMeshInfo);
            CalculateABV2(pMesh.DomainE[i],pMeshInfo.ARTVIS,pMeshInfo.ARTVIS2);
        }
        ParaPackBegin(&procInfo,IdPresABV);
        ParaPackBegin(&procInfo,IdPsiL3);


        for(int i=0;i<pMesh.nde1;i++)
        {
            CalculateGradM(pMesh.DomainE[i],pMeshInfo.MaterialNum);
            RheologyUpdate(pMesh.DomainE[i],StepControl.dt,SRF,pMeshInfo.MaterialNum);
            BlockVibration(pMesh.DomainE[i],StepControl.dt,SRF,pMeshInfo.MaterialNum,StepControl.step_time);
            FailureUpdate(pMesh.DomainE[i],StepControl.dt,SRF,&pMeshInfo);
            CalculateABV2(pMesh.DomainE[i],pMeshInfo.ARTVIS,pMeshInfo.ARTVIS2);
        }


        for(int i=0;i<pMesh.ndv1;i++)
        {
            VelocityUpdate(pMesh.DomainV[i]);
            CalculateL1c(pMesh.DomainV[i],StepControl.dt);
        }
        ParaPackEnd(&procInfo,IdPresABV);
        ParaPackEnd(&procInfo,IdPres);
        ParaPackEnd(&procInfo,IdPsiL3);
        FlushCondEId(&pMesh,IdFCEL4);


        for(int i=pMesh.ndv1;i<pMesh.ndv;i++)
        {
            VelocityUpdate(pMesh.DomainV[i]);
            CalculateL1c(pMesh.DomainV[i],StepControl.dt);
        }


        for(int i=0;i<pMesh.ndb;i++)
            CalculateBFV(pMesh.DomainB[i]);

        for(int k=pMesh.nde1;k<pMesh.nde;k++)
        {
            CalculateCourant(pMesh.DomainE[k]);
            EnergyUpdate(pMesh.DomainE[k], &pMeshInfo, StepControl.dt);
        }
        ParaPackBegin(&procInfo,IdCourant);
        ParaPackBegin(&procInfo,IdMomentum);
        ParaPackBegin(&procInfo,IdPhiL3);


        for(int k=0;k<pMesh.nde1;k++)
        {
            CalculateCourant(pMesh.DomainE[k]);
            EnergyUpdate(pMesh.DomainE[k], &pMeshInfo, StepControl.dt);
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

        UpdateCycleControl(&StepControl);
    }

#ifdef SALEC_USE_TRACER
    TracerClean(&pTracerInfo);
#endif

    MPI_Barrier(SALEC_CART_COMM);

    MPI_Finalize();
}
