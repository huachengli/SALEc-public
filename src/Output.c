//
// Created by huacheng on 2021/2/3.
//

#include "Output.h"
#ifdef USE_GZDIR
#include <zlib.h>
#include <assert.h>
#define CHUNK 16384
#endif
void vts_file_header(struct Mesh * E, FILE *fp)
{
    const char format[] =
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"StructuredGrid\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
            "  <StructuredGrid WholeExtent=\"%s\">\n"
            "    <Piece Extent=\"%s\">\n";

    char Wextent[64], Pextent[64], header[1024];

    snprintf(Pextent, 64, "%d %d %d %d %d %d",
             0, E->nxp - 1,
             0, E->nyp - 1,
             0, E->nzp - 1);
    snprintf(Wextent, 64, "%d %d %d %d %d %d",
             E->noffset, E->nxp - 1 - E->noffset,
             E->noffset, E->nyp - 1 - E->noffset,
             E->noffset, E->nzp - 1 - E->noffset);

    snprintf(header, 1024, format, Wextent, Pextent);

    fputs(header, fp);
}

void vts_file_trailer(struct Mesh * E, FILE *fp)
{
    const char trailer[] =
            "    </Piece>\n"
            "  </StructuredGrid>\n"
            "</VTKFile>\n";

    fputs(trailer, fp);
}

void vtk_point_data_header(struct Mesh *E, FILE *fp)
{
    fputs("      <PointData>\n", fp);
}

void vtk_point_data_trailer(struct Mesh *E, FILE *fp)
{
    fputs("      </PointData>\n", fp);
}

void vtk_cell_data_header(struct Mesh *E, FILE *fp)
{
    fputs("      <CellData>\n", fp);
}

void vtk_cell_data_trailer(struct Mesh *E, FILE *fp)
{
    fputs("      </CellData>\n", fp);
}

void vtk_output_temp(struct Mesh *E, FILE *fp)
{
    float* floattemp = (float *)malloc(E->ne*sizeof(float));
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"temperature\" format=\"%s\">\n",E->vtk_format);
    for(int i=0;i < E->ne;i++)
        floattemp[i] =  (float) E->Elements[i].Temperature;
    if(0==strcmp(E->vtk_format,"binary"))
        write_binary_array(E->ne,floattemp,fp);
    else
        write_ascii_array(E->ne,1,floattemp,fp);
    fputs("        </DataArray>\n", fp);
    free(floattemp);
}

void vtk_output_density(struct Mesh *E, FILE *fp)
{
    float* floattemp = (float *)malloc(E->ne*sizeof(float));
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"density\" format=\"%s\">\n",E->vtk_format);
    for(int i=0;i < E->ne;i++)
        floattemp[i] =  (float) E->Elements[i].Density;
    if(0==strcmp(E->vtk_format,"binary"))
        write_binary_array(E->ne,floattemp,fp);
    else
        write_ascii_array(E->ne,1,floattemp,fp);

    fputs("        </DataArray>\n", fp);
    free(floattemp);
}


void vtk_output_sound(struct Mesh *E, FILE *fp)
{
    float* floattemp = (float *)malloc(E->ne*sizeof(float));
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"sound\" format=\"%s\">\n",E->vtk_format);
    for(int i=0;i < E->ne;i++)
        floattemp[i] =  (float) E->Elements[i].Csound;
    if(0==strcmp(E->vtk_format,"binary"))
        write_binary_array(E->ne,floattemp,fp);
    else
        write_ascii_array(E->ne,1,floattemp,fp);

    fputs("        </DataArray>\n", fp);
    free(floattemp);
}

void vtk_output_module(struct Mesh *E, FILE *fp)
{
    float* floattemp = (float *)malloc(E->ne*sizeof(float));
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Shear\" format=\"%s\">\n",E->vtk_format);
    for(int i=0;i < E->ne;i++)
        floattemp[i] =  (float) E->Elements[i].ShearModule;
    if(0==strcmp(E->vtk_format,"binary"))
        write_binary_array(E->ne,floattemp,fp);
    else
        write_ascii_array(E->ne,1,floattemp,fp);

    fputs("        </DataArray>\n", fp);
    free(floattemp);
}

void vtk_output_debug(struct Mesh *E, FILE *fp)
{
#if ELEMENT_DEBUG
    float* floattemp = (float *)malloc(E->ne*sizeof(float));
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"DEBUG\" format=\"%s\">\n",E->vtk_format);
    for(int i=0;i < E->ne;i++)
        floattemp[i] =  (float) E->Elements[i].Debug;
    if(0==strcmp(E->vtk_format,"binary"))
        write_binary_array(E->ne,floattemp,fp);
    else
        write_ascii_array(E->ne,1,floattemp,fp);

    fputs("        </DataArray>\n", fp);
    free(floattemp);
#endif
}

void vtk_output_sInvSte(struct Mesh *E, FILE *fp)
{
    float* floattemp = (float *)malloc(E->ne*sizeof(float));
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"sInvSte2\" format=\"%s\">\n",E->vtk_format);
    for(int i=0;i < E->ne;i++)
        floattemp[i] =  (float) E->Elements[i].sInvSte2;
    if(0==strcmp(E->vtk_format,"binary"))
        write_binary_array(E->ne,floattemp,fp);
    else
        write_ascii_array(E->ne,1,floattemp,fp);

    fputs("        </DataArray>\n", fp);
    free(floattemp);
}


void vtk_output_vof(struct Mesh *E, FILE *fp)
{
    float* floattemp = (float *)malloc(E->ne*sizeof(float));
    for(int k=0;k<E->nvof;k++)
    {
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"VOF-%d\" format=\"%s\">\n",k,E->vtk_format);
        for(int i=0;i < E->ne;i++)
            floattemp[i] =  (float) E->Elements[i].VOF[k];
        if(0==strcmp(E->vtk_format,"binary"))
            write_binary_array(E->ne,floattemp,fp);
        else
            write_ascii_array(E->ne,1,floattemp,fp);
        fputs("        </DataArray>\n", fp);
    }
    free(floattemp);
}

void vtk_output_stress(struct Mesh *E, FILE *fp)
{
    float* floattemp = (float *)malloc(E->ne*sizeof(float));
    for(int k=0;k<DIM;k++)
        for(int j=k;j<DIM;j++)
        {
            if(k>j) continue;
            fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Stress[%d][%d]\" format=\"%s\">\n",k,j,E->vtk_format);
            for(int i=0;i < E->ne;i++)
                floattemp[i] =  (float) E->Elements[i].PsiL3[1+ SymIndex(k,j)];
            if(0==strcmp(E->vtk_format,"binary"))
                write_binary_array(E->ne,floattemp,fp);
            else
                write_ascii_array(E->ne,1,floattemp,fp);
            fputs("        </DataArray>\n", fp);
        }
    free(floattemp);
}

void vtk_output_strain(struct Mesh *E, FILE *fp)
{
    float* floattemp = (float *)malloc(E->ne*sizeof(float));
    for(int k=0;k<DIM;k++)
        for(int j=k;j<DIM;j++)
        {
            if(k>j) continue;
            fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Strain[%d][%d]\" format=\"%s\">\n",k,j,E->vtk_format);
            for(int i=0;i < E->ne;i++)
                floattemp[i] =  (float) E->Elements[i].Strain[k][j];
            if(0==strcmp(E->vtk_format,"binary"))
                write_binary_array(E->ne,floattemp,fp);
            else
                write_ascii_array(E->ne,1,floattemp,fp);
            fputs("        </DataArray>\n", fp);
        }
    free(floattemp);
}


void vtk_output_gradvof(struct Mesh *E, FILE *fp)
{
    float* floattemp = (float *)malloc(E->ne*sizeof(float));
    for(int k=0;k<DIM;k++)
    {
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"GradVOF[0][%d]\" format=\"%s\">\n",k,E->vtk_format);
        for(int i=0;i < E->ne;i++)
            floattemp[i] =  (float) E->Elements[i].GradVOF[0][k];
        if(0==strcmp(E->vtk_format,"binary"))
            write_binary_array(E->ne,floattemp,fp);
        else
            write_ascii_array(E->ne,1,floattemp,fp);
        fputs("        </DataArray>\n", fp);
    }
    free(floattemp);
}


void vtk_output_pressure(struct Mesh *E, FILE *fp)
{
    float* floattemp = (float *)malloc(E->ne*sizeof(float));
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"pressure\" format=\"%s\">\n",E->vtk_format);
    for(int i=0;i < E->ne;i++)
        floattemp[i] =  (float) E->Elements[i].Pressure;

    if(0==strcmp(E->vtk_format,"binary"))
        write_binary_array(E->ne,floattemp,fp);
    else
        write_ascii_array(E->ne,1,floattemp,fp);

    fputs("        </DataArray>\n", fp);
    free(floattemp);
}

void vtk_output_damage(struct Mesh *E, FILE *fp)
{
    float* floattemp = (float *)malloc(E->ne*sizeof(float));
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Damage\" format=\"%s\">\n",E->vtk_format);
    for(int i=0;i < E->ne;i++)
        floattemp[i] =  (float) E->Elements[i].PsiL3[0];

    if(0==strcmp(E->vtk_format,"binary"))
        write_binary_array(E->ne,floattemp,fp);
    else
        write_ascii_array(E->ne,1,floattemp,fp);

    fputs("        </DataArray>\n", fp);
    free(floattemp);
}



void write_ascii_array(int nn, int perLine, float *array, FILE *fp)
{
    int i;
    switch(perLine) {
        case 1:
            for(i=0; i<nn; i++)
                fprintf(fp, "%.4e\n", array[i]);
            break;
        case 3:
            for(i=0; i < nn/3; i++)
                fprintf(fp,"%.4e %.4e %.4e\n",array[3*i],array[3*i+1],array[3*i+2]);
            break;
        case 6:
            for(i=0; i < nn/6; i++)
                fprintf(fp,"%.4e %.4e %.4e %.4e %.4e %.4e\n",
                        array[6*i],array[6*i+1],array[6*i+2],
                        array[6*i+3],array[6*i+4],array[6*i+5]);
            break;
        default:
            break;
    }
}

void zlibcompress(unsigned char* in, int nn, unsigned char** out, int *nn2)
/* function to compress "in" to "out" reducing size from nn to nn2 */
{
#ifdef USE_GZDIR
    int ntemp=0;

    /* in and out of z-stream */
    unsigned char inz[CHUNK];
    unsigned char outz[CHUNK];

    /* compression level */
    int level = Z_DEFAULT_COMPRESSION;
    int ret,flush;
    int i,j,k;

    /* zlib compression stream */
    z_stream strm;

    /* hope compressed data will be <= uncompressed */
    *out = malloc(sizeof(unsigned char)*nn);

    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    /* zlib init */
    ret = deflateInit(&strm, level);
    if (ret == Z_OK){
        i=0;     // position in "in" array
        do{
            j=0; // position in "inz"
            do{
                inz[j++]=in[i++];
            } while((j<CHUNK) && (i<nn)); // stopps if "inz"-buffer is full or "in" array empty
            strm.avail_in=j;              // set number of input chars

            flush = (i==nn) ? Z_FINISH : Z_NO_FLUSH; // done?
            strm.next_in = inz;           // set input buffer

            do{
                strm.avail_out = CHUNK;   // set number of max output chars
                strm.next_out = outz;     // set output buffer

                /* zlib compress */
                ret = deflate(&strm, flush);
                assert(ret != Z_STREAM_ERROR);

                /* zlib changed strm.avail_out=CHUNK
                 to the number of chars we can NOT use
                 in outz */

                for (k=0;k<CHUNK-strm.avail_out;k++){
                    (*out)[ntemp+k]=outz[k];
                }

                /* increase position in "out" */
                ntemp+=(CHUNK-strm.avail_out);
            }while(strm.avail_out==0);
            assert(strm.avail_in == 0);

        }while (flush != Z_FINISH);
    }
    else{fprintf(stderr,"Error during compression init\n");}

    // now we know how short "out" should be!
    *nn2=ntemp;
    *out = realloc(*out,sizeof(unsigned char)*ntemp);

    (void)deflateEnd(&strm);
#endif
    return;
}

void FloatToUnsignedChar(float * floatarray, int nn, unsigned char * chararray)
{
    /* simple float to unsigned chararray routine via union
    nn=length(intarray) chararray has to be BIG ENOUGH! */
    int i;
    union FloatToUnsignedChars
    {
        float input;
        unsigned char output[4];
    } floattransform;

    for (i=0; i<nn; i++){
        floattransform.input=floatarray[i];
        chararray[4*i]=floattransform.output[0];
        chararray[4*i+1]=floattransform.output[1];
        chararray[4*i+2]=floattransform.output[2];
        chararray[4*i+3]=floattransform.output[3];
    }
}

void base64plushead(unsigned char * in, int nn, int orinn, unsigned char* out)
{
    /* writing vtk compatible zlib compressed base64 encoded data to "out" */
    int i;
    unsigned char * b64head;
    int b64bodylength;
    unsigned char * b64body;
    /* header of data */
    unsigned char * charhead = malloc(sizeof(unsigned char)*16);
    /* - consists of "1" (number of pieces) */
    /* - original datalength in byte */
    /* - original datalength in byte */
    /* - new datalength after z-lib compression */
    int * headInts= malloc(sizeof(int)*4);
    headInts[0]=1;
    headInts[1]=orinn;
    headInts[2]=orinn;
    headInts[3]=nn;
    // transform to unsigned char
    IntToUnsignedChar(headInts,4,charhead);

    // base64: 16byte -> 24byte
    b64head =  malloc(sizeof(unsigned char)*24);
    // fills b64head
    base64(charhead, 16, b64head);

    // base64 data
    b64bodylength = 4*ceil((double) nn/3.0);
    b64body = malloc(sizeof(unsigned char)*b64bodylength);
    // writes base64 data to b64body
    base64(in,nn,b64body);

    // combines header and body
    for (i=0; i<24 ; i++){
        out[i]=b64head[i];
    }

    for (i=0; i<b64bodylength ; i++){
        out[24+i]=b64body[i];
    }

    if(b64body){free(b64body);}
    if(b64head){free(b64head);}
    if(headInts){free(headInts);}
    if(charhead){free(charhead);}
}

void IntToUnsignedChar(int * intarray, int nn, unsigned char * chararray)
{
    /* simple int - to unsigned chararray routine via union
    nn=length(intarray) chararray has to be BIG ENOUGH! */
    int i;
    union IntToUnsignedChars
    {
        int input;
        unsigned char output[4];
    } inttransform;

    for (i=0; i<nn; i++){
        inttransform.input=intarray[i];
        chararray[4*i]=inttransform.output[0];
        chararray[4*i+1]=inttransform.output[1];
        chararray[4*i+2]=inttransform.output[2];
        chararray[4*i+3]=inttransform.output[3];
    }
}

void base64(unsigned char * in, int nn, unsigned char* out)
{
    /*takes *in*-array and "in"-length-"nn" and fills "out"-array
    with base64(in) "out" needs to be big enough!!!
    length(out) >= 4* |^ nn/3.0 ^| */
    char cb64[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    int len;
    int i;

    for (i=0; i < nn; i+=3){

        len = (3 < nn-i ? 3 : nn-i);
        if (len >= 3){
            /* normal base64 encoding */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) | ((in[i+1] & 0xf0) >> 4) ];
            out[i/3*4+2] = cb64[ ((in[i+1] & 0x0f) << 2) | ((in[i+2] & 0xc0) >> 6)];
            out[i/3*4+3] = cb64[ in[i+2] & 0x3f ];
        } else if (len == 2){
            /* at the end of array fill up with '=' */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) | ((in[i+1] & 0xf0) >> 4) ];
            out[i/3*4+2] = cb64[((in[i+1] & 0x0f) << 2)];
            out[i/3*4+3] = (unsigned char) '=';
        } else if (len == 1){
            /* at the end of array fill up with '=' */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) ];
            out[i/3*4+2] = (unsigned char) '=';
            out[i/3*4+3] = (unsigned char) '=';
        }
    }
}


void write_binary_array(int nn, float* array, FILE * f)
{
    /* writes vtk-data array of floats and performs zip and base64 encoding */
    int chararraylength=4*nn;	/* nn floats -> 4*nn unsigned chars */
    unsigned char * chararray = malloc (chararraylength * sizeof(unsigned char));
    int compressedarraylength = 0;
    unsigned char * compressedarray;
    unsigned char ** pointertocompressedarray= &compressedarray;
    int base64plusheadlength;
    unsigned char * base64plusheadarray;

    FloatToUnsignedChar(array,nn,chararray);

    /* compression routine */
    zlibcompress(chararray,chararraylength,pointertocompressedarray,&compressedarraylength);

    /* special header for zip compressed and bas64 encoded data
    header needs 4 int32 = 16 byte -> 24 byte due to base64 (4*16/3) */
    base64plusheadlength = 24 + 4*ceil((double) compressedarraylength/3.0);
    base64plusheadarray  = malloc(sizeof(unsigned char)* base64plusheadlength);

    /* fills base64plusheadarray with everything ready for simple writing */
    base64plushead(compressedarray,compressedarraylength, chararraylength, base64plusheadarray);

    fwrite(base64plusheadarray,sizeof(unsigned char),base64plusheadlength,f);
    fprintf(f,"\n");
    free(chararray);
    free(base64plusheadarray);
    free(compressedarray);
}

void vtk_output_velo(struct Mesh *E, FILE *fp)
{
    int nodes=E->nv;
    float* floatvel = (float *) malloc(nodes*3*sizeof(float));
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"%s\">\n", E->vtk_format);
    for(int i=0;i<E->nv;i++)
    {
        floatvel[3*i + 0] = (float) E->Vertexes[i].Velocity[0];
        floatvel[3*i + 1] = (float) E->Vertexes[i].Velocity[1];
        floatvel[3*i + 2] = (float) E->Vertexes[i].Velocity[2];
    }
    if(0==strcmp(E->vtk_format,"binary"))
        write_binary_array(nodes*3,floatvel,fp);
    else
        write_ascii_array(nodes*3,3,floatvel,fp);
    fputs("        </DataArray>\n", fp);
    free(floatvel);
}

void vtk_output_acce(struct Mesh *E, FILE *fp)
{
    int nodes=E->nv;
    float* floatvel = (float *) malloc(nodes*3*sizeof(float));
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Acceleration\" NumberOfComponents=\"3\" format=\"%s\">\n", E->vtk_format);
    for(int i=0;i<E->nv;i++)
    {
        floatvel[3*i + 0] = (float) E->Vertexes[i].VelocityL1[0];
        floatvel[3*i + 1] = (float) E->Vertexes[i].VelocityL1[1];
        floatvel[3*i + 2] = (float) E->Vertexes[i].VelocityL1[2];
    }
    if(0==strcmp(E->vtk_format,"binary"))
        write_binary_array(nodes*3,floatvel,fp);
    else
        write_ascii_array(nodes*3,3,floatvel,fp);
    fputs("        </DataArray>\n", fp);
    free(floatvel);
}


void vtk_output_coord(struct Mesh *E, FILE *fp)
{
    /*
     * Output Cartesian coordinates as most VTK visualization softwares
     * assume it.
     */

    int nodes = E->nv;
    float* floatpos = (float*) malloc(nodes*3*sizeof(float));

    fputs("      <Points>\n", fp);
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"%s\">\n", E->vtk_format);

    for(int i=0;i<nodes;i++)
    {
        for(int j=0;j<DIM;j++)
            floatpos[3*i + j] = (float) E->Vertexes[i].Position[j];
    }

    if(0==strcmp(E->vtk_format,"binary"))
        write_binary_array(nodes*3,floatpos,fp);
    else
        write_ascii_array(nodes*3,3,floatpos,fp);
    fputs("        </DataArray>\n", fp);
    fputs("      </Points>\n", fp);
    free(floatpos);
}

void vtk_output_all(struct Mesh *E,FILE * fp)
{
    /* first, write volume data to vts file */
    vts_file_header(E, fp);
    /* write node-based field */
    vtk_point_data_header(E, fp);
    vtk_output_velo(E, fp);
    vtk_output_acce(E, fp);
    vtk_point_data_trailer(E, fp);

    /* write element-based field */
    vtk_cell_data_header(E, fp);

    /* TODO: comp_el, heating */
    vtk_output_temp(E, fp);
    vtk_output_pressure(E,fp);
    vtk_output_damage(E,fp);
    vtk_output_density(E,fp);
    vtk_output_sound(E,fp);
    vtk_output_module(E,fp);
    vtk_output_debug(E,fp);
    vtk_output_sInvSte(E,fp);
    vtk_output_vof(E,fp);
    vtk_output_gradvof(E,fp);
    vtk_output_stress(E,fp);
    vtk_output_strain(E,fp);
    vtk_cell_data_trailer(E, fp);

    /* write coordinate */
    vtk_output_coord(E, fp);
    vts_file_trailer(E, fp);
}


void vtk_output(struct Mesh *E,FILE * fp)
{
    /* first, write volume data to vts file */
    vts_file_header(E, fp);
    /* write node-based field */
    vtk_point_data_header(E, fp);
    vtk_output_velo(E, fp);
    vtk_output_acce(E, fp);
    vtk_point_data_trailer(E, fp);

    /* write element-based field */
    vtk_cell_data_header(E, fp);

    vtk_output_temp(E, fp);
    vtk_output_pressure(E,fp);
    vtk_output_damage(E,fp);
    vtk_output_density(E,fp);
    vtk_output_sound(E,fp);
    vtk_output_module(E,fp);
    vtk_output_sInvSte(E,fp);
    vtk_output_vof(E,fp);
    vtk_cell_data_trailer(E, fp);

    /* write coordinate */
    vtk_output_coord(E, fp);
    vts_file_trailer(E, fp);
}


void write_vtm(int nproc, int cycles, const char _prefix[])
{
    const char header[] =
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
            "  <vtkMultiBlockDataSet>\n";

    char fname[50];
    sprintf(fname,"%s.%04d.vtm",_prefix,cycles);
    FILE * fp = fopen(fname,"w");
    fputs(header, fp);

    for(int k=0; k<nproc; k++) {
        fprintf(fp, "    <DataSet index=\"%d\" file=\"%s.proc%d.%d.vts\"/>\n",
                k,_prefix, k, cycles);
    }
    fputs("  </vtkMultiBlockDataSet>\n",fp);
    fputs("</VTKFile>",fp);
    fclose(fp);
}

void ChkWrite(struct Mesh * _mesh,int _xrank, const char _chkPrefix[])
{
    char fname[50];
    sprintf(fname,"%s.proc%04d.chk",_chkPrefix,_xrank);
    FILE *fp = fopen(fname,"wb");
    fwrite(&(_mesh->ne),sizeof(int),1,fp);
    size_t dsize = sizeof(double);
    char ChkNum = 'K';
    for(int k=0;k<_mesh->ne;k++)
    {
        struct Element * ie = _mesh->Elements + k;
        double * VibVelocity = ie->PsiL3 + 8;
        fwrite(&k,sizeof(int),1,fp);
        fwrite(&(ie->PsiL3[0]),dsize,NPSI,fp);
        fwrite(&(ie->Momentum[0]),dsize,DIM,fp);
        fwrite(&(ie->VOF[0]),dsize,NMT,fp);
        fwrite(&(ie->PhiL3[0][0]),dsize,NMT*NPHI,fp);
        fwrite(&(ie->VibPressure),dsize,1,fp);
        fwrite(&(*VibVelocity),dsize,1,fp);
        fwrite(&(ie->Center[0]),dsize,DIM,fp);
        fwrite(&(ChkNum),sizeof(char),1,fp);
    }
    fclose(fp);
}

void ChkLoad(struct Mesh * _mesh,int _xrank, const char _chkPrefix[])
{
    char fname[50];
    sprintf(fname,"%s.proc%04d.chk",_chkPrefix,_xrank);
    FILE *fp = fopen(fname,"rb");
    if(NULL == fp)
    {
        fprintf(stdout,"can not open chkfile: %s\n",fname);
        exit(0);
    }

    size_t dsize = sizeof(double);
    char ChkNum = 'K';
    int ChkElements = 0;

    if(1!=fread(&(ChkElements),sizeof(int),1,fp) || ChkElements!=_mesh->ne)
    {
        fprintf(stdout,"%d Elements in Chk file, %d Elements in mesh\n",ChkElements,_mesh->ne);
        exit(0);
    }

    int ChkIndex = 0;
    for(int k=0;k<_mesh->ne;k++)
    {
        struct Element * ie = _mesh->Elements + k;
        double * VibVelocity = ie->PsiL3 + 8;
        size_t rChunksSize = 0;
        rChunksSize += fread(&ChkIndex,sizeof(int),1,fp);
        rChunksSize += fread(&(ie->PsiL3[0]),dsize,NPSI,fp);
        rChunksSize += fread(&(ie->Momentum[0]),dsize,DIM,fp);
        rChunksSize += fread(&(ie->VOF[0]),dsize,NMT,fp);
        rChunksSize += fread(&(ie->PhiL3[0][0]),dsize,NMT*NPHI,fp);
        rChunksSize += fread(&(ie->VibPressure),dsize,1,fp);
        rChunksSize += fread(&(*VibVelocity),dsize,1,fp);
        rChunksSize += fread(&(ie->Center[0]),dsize,DIM,fp);
        rChunksSize += fread(&ChkNum,sizeof(char),1,fp);

        const size_t rChunkSizeExpected = (4+NPSI+DIM+NMT*(NPHI+1)) + DIM;
        if(ChkNum!='K' || ChkIndex!=k || rChunksSize != rChunkSizeExpected)
        {
            fprintf(stdout,"CheNum or ChkIndex unexpected in Chkfile\n "
                           "ChunkSize=%ld(%ld expected),ChkNum=%c(K expected)\n",
                           rChunksSize,rChunkSizeExpected,ChkNum);
            exit(0);
        }

        ie->Density = 0.0;
        for(int i=0;i<NMT;i++)
        {
            ScalerMove(ie->sPhiL3[i],ie->PhiL3[i],ie->VOF[i]);
            ie->Density += ie->sPhiL3[i][1];
        }
        ie->Mass = ie->Density*ie->Volume;

        ScalerMove(ie->sPsiL3,ie->PsiL3,ie->Mass);
        ie->Damage = ie->PsiL3[0];

    }
    fclose(fp);
}




