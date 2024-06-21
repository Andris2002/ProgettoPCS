#pragma once

#include "Eigen\Eigen"

using namespace Eigen;
using namespace std;

namespace DFN{
struct DFN
{
    unsigned int FracturesNumber; //Number of fractures defining the DFN.
    vector<unsigned int> FracturesId = {}; //Dynamical vector of size FracturesNumber containing the fractures Ids.
    vector<unsigned int> NumberVerties = {}; //Dynamical vector of size FracturesNumber containing the fractures number of vertices.
    vector<MatrixXd> MatricesVertices = {}; /**Dynamical vector of size FracturesNumber containing the fractures matrix of vertices coordinates.
                                               In particular each matrix had dimension 3*NumberVertices[i] for each i <= FracturesNumber*/
};

struct Traces
{
    unsigned int NumTraces; //Number of traces obtained.
    vector<unsigned int> TracesId = {};//Dynamical vector of size NumTraces containing the traces ids.
    vector<Vector2i> FracturesId = {};//Dynamical vector containig the ids of the points defininig the trace.
    vector<Vector2<Vector3d>> TracesPointsCoordinates = {};/**Dynamical vector whose elemets are bidimensional vectors
                                                              containing the coordinates of the two points defining each trace*/
    vector<bool> TracesTips = {}; //Dynamical vector of size NumTraces containing the information whether a trace is through or not.
};

struct Fractures
{
    unsigned int FractureId;//Fracture Id.
    unsigned int NumTraces;//Number of traces within a fracture.
    unsigned int NumPassingTraces;
    unsigned int NumNonPassingTraces;
    vector<unsigned int> PassingTraces = {};//Dynamical vector of size NumPassingTraces containg the ids of the passing traces.
    vector<unsigned int> NonPassingTraces = {};//Dynamical vector of size NumNonPassingTraces containg the ids of the non passing traces.
    vector<double> LengthsPassingTraces = {};//Dynamical vector of size NumPassingTraces containg the lengths of the passing traces.
    vector<double> LengthsNoNPassingTraces = {};//Dynamical vector of size NumNonPassingTraces containg the lengths of the non passing traces.
};
}
