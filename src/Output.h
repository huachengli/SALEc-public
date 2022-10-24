//
// Created by huacheng on 2021/2/3.
//

#ifndef SALEC_OUTPUT_H
#define SALEC_OUTPUT_H

#include "Variable.h"
#include <stdio.h>
#include <string.h>
//#define USE_GZDIR 1
void write_ascii_array(int nn, int perLine, float *array, FILE *fp);
void zlibcompress(unsigned char* in, int nn, unsigned char** out, int *nn2);
void FloatToUnsignedChar(float * floatarray, int nn, unsigned char * chararray);
void base64plushead(unsigned char * in, int nn, int orinn, unsigned char* out);
void IntToUnsignedChar(int * intarray, int nn, unsigned char * chararray);
void base64(unsigned char * in, int nn, unsigned char* out);
void write_binary_array(int nn, float* array, FILE * f);

void vts_file_header(struct Mesh * E, FILE *fp);
void vts_file_trailer(struct Mesh * E, FILE *fp);
void vtk_point_data_header(struct Mesh *E, FILE *fp);
void vtk_point_data_trailer(struct Mesh *E, FILE *fp);
void vtk_cell_data_header(struct Mesh *E, FILE *fp);
void vtk_output_temp(struct Mesh *E, FILE *fp);
void vtk_cell_data_trailer(struct Mesh *E, FILE *fp);

void vtk_output_temp(struct Mesh *E, FILE *fp);
void vtk_output_coord(struct Mesh *E, FILE *fp);
void vtk_output(struct Mesh *E,FILE * fp);
void vtk_output_all(struct Mesh *E,FILE * fp);
void write_vtm(int nproc, int cycles, const char _prefix[]);

void ChkWrite(struct Mesh * _mesh,int _xrank, const char _chkPrefix[]);
void ChkLoad(struct Mesh * _mesh,int _xrank, const char _chkPrefix[]);
#endif //SALEC_OUTPUT_H
