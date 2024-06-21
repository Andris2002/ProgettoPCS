#pragma once

#include "DFNstructs.hpp"

namespace DFN{
bool ImportDFN(const string& filename,DFN& dfn);
void PrintTraces(const DFN& dfn);
void PrintFractures(const DFN& dfn);
}
