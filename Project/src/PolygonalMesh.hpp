#ifndef POLYGONALMESH_HPP
#define POLYGONALMESH_HPP

#include "DFNstructs.hpp"
#include <vector>
#include <array>
#include <cmath>
#define pi 3.14159265358979323846

using namespace std;

namespace DFN {

/*inline array<double,3> scalarVector_product(double k, array<double,3>& P){

    P[0] = k*P[0];
    P[1] = k*P[1];
    P[2] = k*P[2];

    return P;
}*/
/*inline double calculate_angle(const array<double,3>& baricenter,const array<double,3>& P){

    double distancex = P[0] - baricenter[0];
    double distancey = P[1] - baricenter[1];
    double hypotenuse = sqrt(distancex*distancex + distancey*distancey);
    double angle;
    if( distancex > 0 && distancey >=0){ // sta nel primo quadrante
        angle = acos(distancex / hypotenuse );
    } else if (distancex <= 0 && distancey > 0){ // secondo quadrante
        angle = acos(-distancex / hypotenuse ) + pi/2;
    } else if (distancex < 0 && distancey <= 0){ // terzo quadrante
        angle = acos(-distancex / hypotenuse ) + pi;
    } else { // quarto quadrante
        angle = acos( distancex / hypotenuse ) + 3*pi/2;
    }

    return angle;
}*/

vector<Trace> SortAssociatedTraces(vector<Trace> & v);
pair<vector<unsigned int>, vector<array<double, 3>>> sortingVertices(const vector<array<double, 3>>& vectorPoints);
PolygonalMesh importMesh(Fracture& fracture);
vector<pair<array<double, 3>, pair<unsigned int, unsigned int>>> calculate_intersections(const PolygonalMesh& mesh, const vector<unsigned int>& polygon, const Trace& trace);
vector<pair<Fracture, vector<Trace> > > cutTrace(pair<Fracture, vector<Trace> > pair);
Trace extendTrace(const Trace& trace, const PolygonalMesh& mesh);
PolygonalMesh processFracture(DFN& dfn, Fracture& fracture);
vector<PolygonalMesh> Alg_process_cut(DFN& dfn);

}

#endif
