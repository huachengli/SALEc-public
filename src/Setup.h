//
// Created by huacheng on 2020/12/22.
//

#ifndef SALEC_SETUP_H
#define SALEC_SETUP_H

#include <stdlib.h>
#include "Variable.h"
#include "Error.h"
#include "Iteration.h"
#include <string.h>

int GetMeshGeometry(struct Mesh *_mesh, const struct MeshInfo *_info);
int BuildMeshGeometry(struct Mesh *_mesh, const struct MeshInfo *_info);
int Calculate_SV(struct Mesh *_mesh);
void Calculate_SVboundary(struct Boundary *_boundary);

typedef void (*InitConditionV)(struct Vertex *);
typedef void (*InitConditionE)(struct Element *);
//void SetInitVertex(struct Mesh *_mesh, InitConditionV _icf, struct ANEOSTable _t[], int _nm);
void SetInitVertex_gas(struct Mesh *_mesh, InitConditionV _icf, int _nm);
void SetInitElement(struct Mesh * _mesh, InitConditionE _icf);
void SetVelStavlity(struct Mesh *_mesh, InitConditionV _icf);
void MapV2Ea(struct Element * ie,struct ANEOSTable _t[], int _nm);
//void SearchDomain(struct Mesh *_mesh);
void SearchDomain(struct Mesh *_mesh, void (*_cef)(int, struct ConditionE *));
#endif //SALEC_SETUP_H
