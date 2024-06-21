#include "DFNstructs.hpp"
#include "ImportExport.hpp"
#include "AlgFractures.hpp"
#include "PolygonalMesh.hpp"
#include <fstream>
#include <iostream>
#include <chrono>
#include <thread>

using namespace std;
using namespace DFN;

int main()
{
    DFN::DFN dfn;
    string filename = "./DFN/FR3_data.txt";

    DFN::ImportDFN(filename,dfn);
    DFN::run_Alg(dfn);

    //DFN::Alg_process_cut(dfn) invoca la seconda parte del progetto
    //ESSA NON FUNZIONA PER FILE DIVERSI DA fr3_data.
    //auto valore = DFN::Alg_process_cut(dfn); //disasteriscare solo per provare con FR3_data
    DFN::PrintFractures(dfn);
    DFN::PrintTraces(dfn);

  return 0;
}
