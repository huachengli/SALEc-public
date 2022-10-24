//
// Created by huacheng on 8/16/21.
//

#include "Output.h"
#include "InputParser.h"

int main(int argc,char *argv[])
{
    InputFile * ifp = OpenInputFile("../SALEc.inp");
    char _prefix[MaxStrLen];
    GetValueS(ifp,"output.prefix",_prefix,"ParaTest");

    int npgx = GetValueI(ifp,"processor.npgx","2");
    int npgy = GetValueI(ifp,"processor.npgy","2");
    int npgz = GetValueI(ifp,"processor.npgz","2");
    int npgs = npgx*npgy*npgz;

    int  MaxStep,MinStep;
    if(2==argc)
    {
        MaxStep = strtol(argv[1],NULL,10);
    } else if(3==argc)
    {
        MinStep = strtol(argv[1],NULL,10);
        MaxStep = strtol(argv[2],NULL,10);
    } else
    {
        MaxStep = GetValueI(ifp,"cycle.maxstep","100");
    }
    CloseInputFile(ifp);

    for(int k=MinStep;k<MaxStep;k++)
    {
        fprintf(stdout,"Writing VTM files ...\n");
        write_vtm(npgs,k,_prefix);
        fprintf(stdout,"Write %s.%04d.vtm\n",_prefix,k);
    }
    fprintf(stdout,"End of Writing VTM files\n");
    return 0;
}