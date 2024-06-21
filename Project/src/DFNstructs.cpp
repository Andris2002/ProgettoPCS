#include "DFNstructs.hpp"
#include <cmath>

using namespace std;
using namespace DFN;

namespace DFN{


bool Fracture::get_baricenter(){
    if (Vertices.empty())
        return false;

    cout<<num_vertices<<endl;
    //compute baricenter
    for (unsigned int i = 0; i<3; i++){
        Baricenter[i] = (Vertices.at(0)).at(i);
        for (int j = 2; j<num_vertices +1;j++){
            Baricenter[i] = (Baricenter[i]*(j-1) +Vertices.at(j-1).at(i))/j;
        }
    }


    //find the most distant vertice from the baricenter
    double max = 0;
    unsigned int max_vertice = 0;

    for (unsigned int i =0; i < num_vertices ; i++){
        double dist_squared = distance_squared(Baricenter,Vertices[i]);
        if (max < dist_squared){
            max = dist_squared;
            max_vertice = i;
        }
    }


    Radius = sqrt(max);

    //compute normal vector (normal to the plane formed by the poligon
    array<double,3> V1 = substract(Vertices[max_vertice], Baricenter);

    array<double,3> V2 = substract(Vertices[(max_vertice +1)%num_vertices], Baricenter);
    array<double,3> origin ={0,0,0};



    NormalVector = vector_product(V1,V2);
    double length_NormalVector = sqrt(distance_squared(NormalVector,origin));
    rescale(1/length_NormalVector,NormalVector);






    return true;
}
}
