#pragma once
#include <iostream>
#include "Eigen/Eigen"
#include "AlgFractures.hpp"
#include <cmath>


#define TOLL  1.00e-10

using namespace Eigen;
using namespace std;

namespace DFN{

inline double distance_squared(array<double,3>& P1, array<double,3>& P2){
    double d = (P1[0]-P2[0])*(P1[0]-P2[0]) + (P1[1]-P2[1])*(P1[1]-P2[1]) + (P1[2]-P2[2])*(P1[2]-P2[2]);

    return d;
}

inline array<double,3> vector_product(array<double,3>& P1, array<double,3>& P2){
    array<double,3> a;
    a[0] = P1[1]*P2[2] - P1[2]*P2[1];
    a[1] = P1[2]*P2[0] - P1[0]*P2[2];
    a[2] = P1[0]*P2[1] - P1[1]*P2[0];

    return a;
}

inline void rescale( double r,array<double,3>& P){
    P[0] = P[0]*r;
    P[1] = P[1]*r;
    P[2] = P[2]*r;
}

inline array<double,3> substract(array<double,3>& P1, array<double,3>& P2){
    array<double,3> a;
    a[0] = P1[0] - P2[0];
    a[1] = P1[1] - P2[1];
    a[2] = P1[2] - P2[2];

    return a;

}

inline array<double,3> add(array<double,3>& P1, array<double,3>& P2){
    array<double,3> a;
    a[0] = P1[0] + P2[0];
    a[1] = P1[1] + P2[1];
    a[2] = P1[2] + P2[2];

    return a;

}

inline array<double,3> scalar_vector_product(double k, array<double,3>& P){

    array<double,3> a;
    a[0] = k*P[0];
    a[1] = k*P[1];
    a[2] = k*P[2];

    return a;
}

inline double scalar_product(array<double,3>& P1, array<double,3>& P2){
    double a = P1[0]*P2[0] + P1[1]*P2[1] + P1[2]*P2[2];

    return a;

}

inline bool is_equal(double d1, double d2){

    if (abs(d1 - d2)< TOLL)
        return true;
    return false;
}


struct PolygonalMesh {

    unsigned int numCell0D = 0;
    vector<unsigned int> Cell0DId = {};
    vector<array<double, 3>> Cell0DCoordinates = {};

    unsigned int numCell1D = 0;
    vector<unsigned int> Cell1DId = {};
    vector<array<unsigned int, 2>> boundaryCell0DId = {};

    unsigned int numCell2D = 0;
    vector<vector<unsigned int>> Cell2DVerticesId = {};
    vector<vector<unsigned int>> Cell2DEdgesId = {};
};

struct Fracture
{
    int FractureId;
    /*/unsigned int FractureId;                      //Fracture Id.
    unsigned int NumTraces;                       //Number of traces within a fracture.
    unsigned int NumPassingTraces;
    unsigned int NumNonPassingTraces;
    vector<unsigned int> PassingTraces = {};      //Dynamical vector of size NumPassingTraces containg the ids of the passing traces.
    vector<unsigned int> NonPassingTraces = {};   //Dynamical vector of size NumNonPassingTraces containg the ids of the non passing traces.
    vector<double> LengthsPassingTraces = {};     //Dynamical vector of size NumPassingTraces containg the lengths of the passing traces.
    vector<double> LengthsNoNPassingTraces = {};  //Dynamical vector of size NumNonPassingTraces containg the lengths of the non passing traces.
*/
    unsigned int num_vertices;  //number of Vertices of the given fracture, please fill in when exporting the fracture!!
    vector<array<double,3>> Vertices;
    vector<unsigned int> VerticesId;
    array<double,3> Baricenter;
    array<double,3> NormalVector;
    double Radius;
    bool get_baricenter();
    vector<unsigned int> Traces; //vector with the position of the traces within DFN::Traces that belong to this fracture



};


struct Trace
{
    unsigned int TraceId;
    array<array<double,3>,2> Coordinates;
    array<pair<unsigned int,bool>, 2> AssociatedFractures;
    double TraceLength;

    bool passanteTemp = true; //true se e passante

    unsigned int NumTraces;                                  //Number of traces obtained.
    vector<unsigned int> TracesId = {};                      //Dynamical vector of size NumTraces containing the traces ids.
    vector<Vector2i> FracturesId = {};                       //Dynamical vector containig the ids of the fractures defininig the trace.
    vector<Vector2<Vector3d>> TracesPointsCoordinates = {};  /**Dynamical vector whose elemets are bidimensional vectors
                                                              containing the coordinates of the two points defining each trace*/
    vector<bool> TracesTips = {};                            //Dynamical vector of size NumTraces containing the information whether a trace is through or not.
};

struct DFN
{
    unsigned int FracturesNumber;             //Number of fractures defining the DFN.
    vector<unsigned int> FracturesId = {};    //Dynamical vector of size FracturesNumber containing the fractures Ids.
    vector<unsigned int> NumberVerties = {};  //Dynamical vector of size FracturesNumber containing the fractures number of vertices.
    vector<vector<vector<double>>> MatricesVertices = {};

    unsigned int number_fractures;
    unsigned int number_traces;
    vector<Fracture> Fractures = {};
    vector<Trace> Traces = {};
};
}
