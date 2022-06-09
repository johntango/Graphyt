#ifndef EXTRAPOLATE_H
#define EXTRAPOLATE_H

#include "Parameters.hpp"
#include "MaterialPoints.hpp"
#include "Nodes.hpp"
// Task Functions //
CPUGPU void extrapToNodes(Parameters *mpmParameters, MaterialPoints* mpts, Nodes* nds,std::vector<MaterialModels>* materials,int n);
CPUGPU void extrapToMatPoints(Parameters *mpmParameters, MaterialPoints* mpts, Nodes* nds,std::vector<MaterialModels>* materials,int mp);
// Helper Functions //
CPUGPU void   extrapolate_Matpnt_to_Node(Parameters *mpmParameters, MaterialPoints* mpts, Nodes* nds,std::vector<MaterialModels> *materials,int mp, int n);
CPUGPU void   extrapolate_Node_to_Matpnt(Parameters *mpmParameters, MaterialPoints* mpts, Nodes* nds,std::vector<MaterialModels> *materials,int mp, int n);
CPUGPU void   sortParticlesByCellID_CPU(Parameters* mpmParameters,MaterialPoints* mpts, Nodes* nds);
CPUGPU void   sortParticlesByCellID_GPU(Parameters* mpmParameters,MaterialPoints* mpts, Nodes* nds);
CPUGPU double shapeFunctionClassic(MaterialPoints *mpts,Nodes* nds ,int mpIdx, int ndIdx, int N_or_G, int xyz);
CPUGPU double shapeFunctionGIMP(MaterialPoints *mpts,Nodes* nds ,int mpIdx, int ndIdx, int N_or_G, int xyz);
CPUGPU void   getGIMPNodes(Parameters* mpmParameters, MaterialPoints *mpts, Nodes *nds, int mpIdx,int *N);
#endif