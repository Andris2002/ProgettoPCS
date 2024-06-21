#include "DFNstructs.hpp"

#pragma once



using namespace std;
namespace DFN{

struct DFN;
struct Fracture;


bool run_Alg(DFN & dfn);
bool find_trace(Fracture & f1, Fracture & f2, DFN & dfn);
bool Exclude2(Fracture & f1, Fracture & f2);



}
