#include "PolygonalMesh.hpp"
#include "DFNstructs.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#define TOLL  1.00e-10




using namespace std;

namespace DFN{

DFN * dfnGlobal = nullptr;
double totVertices = 0;



vector<Trace> SortAssociatedTraces(vector<Trace>& v)
{
    vector<Trace> associatedTraces;
    associatedTraces.reserve(v.size());
    /*
    unordered_map<unsigned int, unsigned int> traceIndexMap; // Mappa gli ID delle tracce agli indici nel vettore dfn.Traces
    // coordinate tracce: array<array<double,3>,2> Coordinates;
    // Popola la mappa degli ID delle tracce
    for (unsigned int i = 0; i < dfn.Traces.size(); i++) {
        traceIndexMap[dfn.Traces[i].TraceId] = i;
    }*/
    vector<Trace> passingTraces;
    vector<Trace> notPassingTraces;
    for (Trace & t : v){

        if (t.passanteTemp)
            passingTraces.push_back(t);
        else
            notPassingTraces.push_back(t);
        //todo: rimuovere v (erase?)
    }
    if (!passingTraces.empty()){
        /*double max =0;
        int maxId = 0;
        for (unsigned int i =0; i< passingTraces.size(); i++){
            if (passingTraces[i].TraceLength >max){
                max = passingTraces[i].TraceLength;
                maxId = i;
            }
        }
        pair<Trace,bool> ritorno;
        ritorno.first = passingTraces[maxId];
        ritorno.second = true;*/

        sort(passingTraces.begin(), passingTraces.end(),
             [](const Trace& a, const Trace& b) {
            return a.TraceLength < b.TraceLength;
             });




    }
    if(!notPassingTraces.empty()){
        /*double max =0;
        int maxId = 0;
        for (unsigned int i =0; i< notPassingTraces.size(); i++){
            if (notPassingTraces[i].TraceLength >max){
                max = notPassingTraces[i].TraceLength;
                maxId = i;
            }
        }
        pair<Trace,bool> ritorno;
        ritorno.first = notPassingTraces[maxId];
        ritorno.second = false;
        return ritorno;*/

        sort(notPassingTraces.begin(), notPassingTraces.end(),
             [](const Trace& a, const Trace& b) {
            return a.TraceLength < b.TraceLength;
             });

    }
    associatedTraces.insert(associatedTraces.end(), notPassingTraces.begin(), notPassingTraces.end());
    associatedTraces.insert(associatedTraces.end(), passingTraces.begin(), passingTraces.end());

    return associatedTraces;

    /*

    // Ordina le tracce passanti per lunghezza decrescente
    sort(passingTraces.begin(), passingTraces.end(),
         [](const Trace& a, const Trace& b) {
             return a.TraceLength > b.TraceLength;
         });

    // Ordina le tracce non passanti per lunghezza decrescente
    sort(notPassingTraces.begin(), notPassingTraces.end(),
         [](const Trace& a, const Trace& b) {
             return a.TraceLength > b.TraceLength;
         });

    // Aggiungi le tracce passanti ordinate all'elenco associato
    associatedTraces.insert(associatedTraces.end(), passingTraces.begin(), passingTraces.end());

    // Aggiungi le tracce non passanti ordinate all'elenco associato
    associatedTraces.insert(associatedTraces.end(), notPassingTraces.begin(), notPassingTraces.end());

    return associatedTraces;*/
}

vector<pair<Fracture,vector<Trace>>> cutTrace(pair<Fracture,vector<Trace>> coppia) {
    //vector<PolygonalMesh> newMeshes;
    //cout << "4" << endl;
    vector<pair<Fracture,vector<Trace>>> v = {};



    auto tracceOrdinate = SortAssociatedTraces(coppia.second);

    auto primaTraccia = tracceOrdinate.back();
    tracceOrdinate.pop_back();

    //classifichiamo i vertici di fracture in base al lato su cui stanno di prima traccia
    auto vettoreTraccia = substract(primaTraccia.Coordinates[1], primaTraccia.Coordinates[0]);
    auto normaleTraccia = vector_product(vettoreTraccia, coppia.first.NormalVector);
    array<double,3> origin = {0,0,0};

    //torviamo i lati sui quali avviene il taglio
    unsigned int taglio[2] = {0,0};
    unsigned int ordineCoordTraccia[2] = {0,1};

    auto diff = substract(coppia.first.Vertices[0], primaTraccia.Coordinates[0]);
    double p = scalar_product(diff, normaleTraccia);

    //ricerca del primo vertici che passi sull-altro lato rispetto a
    for (unsigned int k = 1; k<coppia.first.num_vertices;k++){
        diff = substract(coppia.first.Vertices[k],primaTraccia.Coordinates[0]);
        double lato = p*(scalar_product(diff,normaleTraccia));

        if (lato <= TOLL){

                taglio[0] = k;
                k = coppia.first.num_vertices;
                continue;

        }
    }

    for (unsigned int k = taglio[0] +1; k<coppia.first.num_vertices; k++){
        diff = substract(coppia.first.Vertices[k],primaTraccia.Coordinates[0]);
        double lato = p*(scalar_product(diff,normaleTraccia));

        if (lato >= -TOLL){

            taglio[1] = k;
            k = coppia.first.num_vertices;
            continue;

        }
    }
    if (!taglio[1]) // in questo caso il primo vertice a ritornare sul lato di P* e P* stesso.
        taglio[1] = coppia.first.num_vertices;


    //se e non passante, trovare l'intersezione

    if (!primaTraccia.passanteTemp){

        //intersezione con il primo lato
        auto diff = substract(coppia.first.Vertices[taglio[0]], coppia.first.Vertices[taglio[0] -1]);
        auto normaleAlLato = vector_product(diff, coppia.first.NormalVector);
        double len_normale = distance_squared(origin, normaleAlLato);
        len_normale = sqrt(len_normale);
        rescale(1/len_normale,normaleAlLato);

        diff = substract(primaTraccia.Coordinates[0], coppia.first.Vertices[taglio[0]]);
        double b1 = abs(scalar_product(normaleAlLato,diff));

        diff = substract(primaTraccia.Coordinates[1], coppia.first.Vertices[taglio[0]]);
        double b2 = abs(scalar_product(normaleAlLato,diff));

        if (!is_equal(b1,0) && !is_equal(b2,0)){
            double maxb = fmax(b1,b2);
            double minb = fmin(b1,b2);
            //vertice da prolungare
            unsigned int daProlungare = 0;
            unsigned int nondaProlungare = 1;
            if (is_equal(minb,b2)){
                daProlungare = 1;
                nondaProlungare = 0;
                ordineCoordTraccia[0] = 1;
                ordineCoordTraccia[1] = 0;
            }
            auto diff = substract(primaTraccia.Coordinates[daProlungare], primaTraccia.Coordinates[nondaProlungare]);
            rescale(minb/(maxb-minb), diff);
            primaTraccia.Coordinates[daProlungare] = add(primaTraccia.Coordinates[daProlungare], diff);

        } else if(is_equal(b2,0)){ //in questo caso e il secondo vertice ad essere sul primo lato
            ordineCoordTraccia[0] = 1;
            ordineCoordTraccia[1] = 0;
        }

        //intersezione con il secondo lato
        diff = substract(coppia.first.Vertices[taglio[1]%(coppia.first.num_vertices)], coppia.first.Vertices[taglio[1] -1]);
        normaleAlLato = vector_product(diff, coppia.first.NormalVector);
        len_normale = distance_squared(origin, normaleAlLato);
        len_normale = sqrt(len_normale);
        rescale(1/len_normale,normaleAlLato);

        diff = substract(primaTraccia.Coordinates[0], coppia.first.Vertices[taglio[1]%(coppia.first.num_vertices)]);
        b1 = abs(scalar_product(normaleAlLato,diff));

        diff = substract(primaTraccia.Coordinates[1], coppia.first.Vertices[taglio[1]%(coppia.first.num_vertices)]);
        b2 = abs(scalar_product(normaleAlLato,diff));

        if (!is_equal(b1,0) && !is_equal(b2,0)){
            double maxb = fmax(b1,b2);
            double minb = fmin(b1,b2);
            //vertice da prolungare
            unsigned int daProlungare = 1;
            unsigned int nondaProlungare = 0;
            if (is_equal(minb,b1)){
                daProlungare = 0;
                nondaProlungare = 1;
                ordineCoordTraccia[0] = 1;
                ordineCoordTraccia[1] = 0;
            }
            auto diff = substract(primaTraccia.Coordinates[daProlungare], primaTraccia.Coordinates[nondaProlungare]);
            rescale(minb/(maxb-minb), diff);
            primaTraccia.Coordinates[daProlungare] = add(primaTraccia.Coordinates[daProlungare], diff);

        } else if(is_equal(b1,0)){//in questo caso e il primo vertice ad essere sul secondo lato
            ordineCoordTraccia[0] = 1;
            ordineCoordTraccia[1] = 0;
        }


    }


    vector<array<double,3>> Vertices1;
    vector<unsigned int> VerticesId1;

    vector<array<double,3>> Vertices2;
    vector<unsigned int> VerticesId2;

    //vertici e loro id della prima nuova frattura
    Vertices1.push_back(primaTraccia.Coordinates[ordineCoordTraccia[0]]);
    VerticesId1.push_back(totVertices +1);
    for (unsigned int i= taglio[0];i < taglio[1]; i++){
        Vertices1.push_back(coppia.first.Vertices[i]);
        VerticesId1.push_back(coppia.first.VerticesId[i]);
    }

    Vertices1.push_back(primaTraccia.Coordinates[ordineCoordTraccia[1]]);
    VerticesId1.push_back(totVertices +2);

    //vertici e loro id della seconda nuova frattura
    Vertices2.push_back(primaTraccia.Coordinates[ordineCoordTraccia[0]]);
    VerticesId2.push_back(totVertices +1);

    Vertices2.push_back(primaTraccia.Coordinates[ordineCoordTraccia[1]]);
    VerticesId2.push_back(totVertices +2);

    for (unsigned int i= taglio[1];i < taglio[0] + coppia.first.num_vertices; i++){
        Vertices2.push_back(coppia.first.Vertices[i%(coppia.first.num_vertices)]);
        VerticesId2.push_back(coppia.first.VerticesId[i%(coppia.first.num_vertices)]);
    }

    totVertices = totVertices +2;

    //creazione delle due nuove fratture
    Fracture f1;
    f1.num_vertices = Vertices1.size();
    f1.Vertices = Vertices1;
    f1.VerticesId = VerticesId1;
    f1.get_baricenter();

    Fracture f2;
    f2.num_vertices = Vertices2.size();
    f2.Vertices = Vertices2;
    f2.VerticesId = VerticesId2;
    f2.get_baricenter();

    //assegnazione delle tracce
    vector<Trace> tracce1;
    vector<Trace> tracce2;


    diff = substract(primaTraccia.Coordinates[ordineCoordTraccia[1]],primaTraccia.Coordinates[ordineCoordTraccia[0]]);
    auto indicatore = vector_product(diff,coppia.first.NormalVector);
    double len_indicatore = distance_squared(indicatore,origin);
    len_indicatore = sqrt(len_indicatore);
    rescale(1/len_indicatore,indicatore);


    for (Trace & t: tracceOrdinate){
        diff = substract(t.Coordinates[0], primaTraccia.Coordinates[ordineCoordTraccia[0]]);
        double p1 = scalar_product(diff, indicatore);

        diff = substract(t.Coordinates[1], primaTraccia.Coordinates[ordineCoordTraccia[0]]);
        double p2 = scalar_product(diff, indicatore);

        //la traccia interseca prima traccia
        if (p1*p2 <0 - TOLL){

            auto Intersezione = substract(t.Coordinates[0], t.Coordinates[1]);
            rescale(p2/(p2+p1), Intersezione);
            Intersezione = add(t.Coordinates[1], Intersezione);

            Trace t1;
            t1.TraceId = 0;
            t1.Coordinates[0] = Intersezione;
            t1.Coordinates[1] = t.Coordinates[0];
            t1.TraceLength = sqrt(distance_squared(Intersezione,t.Coordinates[0]));
            t1.passanteTemp = false;

            Trace t2;
            t2.TraceId = 0;
            t2.Coordinates[0] = Intersezione;
            t2.Coordinates[1] = t.Coordinates[1];
            t2.TraceLength = sqrt(distance_squared(Intersezione,t.Coordinates[1]));

            //t1 appartiene alla prima frattura
            if (p1 >0){
                //controlliamo se t1 e passante o non passante
                for (unsigned int i = 0; i < Vertices1.size() -2; i++){
                    auto v1 = substract(t.Coordinates[0],Vertices1[i]);
                    auto v2 = substract(t.Coordinates[0],Vertices1[i +1]);
                    v1 = vector_product(v1,v2);
                    double evaluate = distance_squared(v1,origin);
                    if (is_equal(evaluate, 0)){
                        t1.passanteTemp = true;
                    }
                }

                //controlliamo se t2 e passante o non passante
                for (unsigned int i = 1; i < Vertices2.size() -1; i++){
                    auto v1 = substract(t.Coordinates[1],Vertices2[i]);
                    auto v2 = substract(t.Coordinates[1],Vertices2[i +1]);
                    v1 = vector_product(v1,v2);
                    double evaluate = distance_squared(v1,origin);
                    if (is_equal(evaluate, 0)){
                        t2.passanteTemp = true;
                    }
                }

                tracce1.push_back(t1);
                tracce2.push_back(t2);
            }
            //t1 appartiene alla seconda frattura
            else
            {
                //controlliamo se t2 e passante o non passante
                for (unsigned int i = 0; i < Vertices1.size() -2; i++){
                    auto v1 = substract(t.Coordinates[1],Vertices1[i]);
                    auto v2 = substract(t.Coordinates[1],Vertices1[i +1]);
                    v1 = vector_product(v1,v2);
                    double evaluate = distance_squared(v1,origin);
                    if (is_equal(evaluate, 0)){
                        t2.passanteTemp = true;
                    }
                }

                //controlliamo se t1 e passante o non passante
                for (unsigned int i = 1; i < Vertices2.size() -1; i++){
                    auto v1 = substract(t.Coordinates[0],Vertices2[i]);
                    auto v2 = substract(t.Coordinates[0],Vertices2[i +1]);
                    v1 = vector_product(v1,v2);
                    double evaluate = distance_squared(v1,origin);
                    if (is_equal(evaluate, 0)){
                        t1.passanteTemp = true;
                    }
                }

                tracce1.push_back(t2);
                tracce2.push_back(t1);




            }


        }
        else if (is_equal(p1*p2,0)){ //una delle coordinate della traccia giace proprio su prima traccia
            if (is_equal(p1,0)){// a giacere su prima traccia e la prima coordinata

                if (p2> 0){// la traccia appartiene alla prima frattura
                    //controlliamo se passante o non passante
                    t.passanteTemp = false;
                    for (unsigned int i = 0; i < Vertices1.size() -2; i++){
                        auto v1 = substract(t.Coordinates[1],Vertices1[i]);
                        auto v2 = substract(t.Coordinates[1],Vertices1[i +1]);
                        v1 = vector_product(v1,v2);
                        double evaluate = distance_squared(v1,origin);
                        if (is_equal(evaluate, 0)){
                            t.passanteTemp = true;
                        }
                    }

                    tracce1.push_back(t);

                } else {// la traccia appartiene alla seconda frattura

                    //controlliamo sepassante o non passante
                    t.passanteTemp = false;
                    for (unsigned int i = 1; i < Vertices2.size() -1; i++){
                        auto v1 = substract(t.Coordinates[1],Vertices2[i]);
                        auto v2 = substract(t.Coordinates[1],Vertices2[i +1]);
                        v1 = vector_product(v1,v2);
                        double evaluate = distance_squared(v1,origin);
                        if (is_equal(evaluate, 0)){
                            t.passanteTemp = true;
                        }
                    }

                    tracce2.push_back(t);

                }

            } else {// a giacere su prima traccia e la seconda coordinata
                if (p1> 0){// la traccia appartiene alla prima frattura
                    //controlliamo se passante o non passante
                    t.passanteTemp = false;
                    for (unsigned int i = 0; i < Vertices1.size() -2; i++){
                        auto v1 = substract(t.Coordinates[0],Vertices1[i]);
                        auto v2 = substract(t.Coordinates[0],Vertices1[i +1]);
                        v1 = vector_product(v1,v2);
                        double evaluate = distance_squared(v1,origin);
                        if (is_equal(evaluate, 0)){
                            t.passanteTemp = true;
                        }
                    }

                    tracce1.push_back(t);

                } else {// la traccia appartiene alla seconda frattura

                    //controlliamo sepassante o non passante
                    t.passanteTemp = false;
                    for (unsigned int i = 1; i < Vertices2.size() -1; i++){
                        auto v1 = substract(t.Coordinates[0],Vertices2[i]);
                        auto v2 = substract(t.Coordinates[0],Vertices2[i +1]);
                        v1 = vector_product(v1,v2);
                        double evaluate = distance_squared(v1,origin);
                        if (is_equal(evaluate, 0)){
                            t.passanteTemp = true;
                        }
                    }

                    tracce2.push_back(t);

                }
            }
        }

        else { //la traccia e interamente interna a una delle due fratture
            t.passanteTemp = false;
            if (p1 > 0)
                tracce1.push_back(t);
            else
                tracce2.push_back(t);

        }

    }

    pair<Fracture,vector<Trace>> coppia1;
    pair<Fracture,vector<Trace>> coppia2;

    coppia1.first = f1;
    coppia1.second = tracce1;

    coppia2.first = f2;
    coppia2.second = tracce2;

    vector<pair<Fracture,vector<Trace>>> ritorno;
    ritorno.push_back(coppia1);
    ritorno.push_back(coppia2);


    /*
    if (mesh.Cell2DVerticesId.empty()) {
        cerr << "La mesh non contiene poligoni" << endl;
        newMeshes.push_back(mesh);
        return newMeshes;
    }



    auto& polygon = mesh.Cell2DVerticesId[0];

    // Check for intersection between the segment and the trace
    vector<pair<array<double, 3>, pair<unsigned int,unsigned int>>> intersections = calculate_intersections(mesh,polygon,trace);

    if (intersections.size() == 2) {
        // ho trovato le due intersezioni, ora bisogna capire quali dei vecchi vertici
        // la delimitano, ma in realtà questa info è dentro intersections
        array<double, 3> intersection1 = intersections[0].first;
        array<double, 3> intersection2 = intersections[1].first;
        // Create new intersection point vertices and assign IDs
        unsigned int intersection1Id = mesh.numCell0D++;
        mesh.Cell0DCoordinates.push_back(intersection1);
        unsigned int intersection2Id = mesh.numCell0D++;
        mesh.Cell0DCoordinates.push_back(intersection2);
        mesh.Cell0DId.push_back(intersection1Id);
        mesh.Cell0DId.push_back(intersection2Id);
        // update also the edges
        unsigned int edgeId = mesh.numCell1D++;
        mesh.Cell1DId.push_back(edgeId);
        mesh.boundaryCell0DId.push_back({intersection1Id,intersection2Id});

        // Create new polygons
        vector<unsigned int> newPolygon1, newPolygon2;
        vector<unsigned int> newEdges1, newEdges2;

        pair<unsigned int,unsigned int> vertexIntersection1 = intersections[0].second;// restituisce gli id di v1, v2
        pair<unsigned int,unsigned int> vertexIntersection2 = intersections[1].second;// restituisce gli id di v3, v4

        newPolygon1.push_back(vertexIntersection1.first); // aggiungo gli id di v1
        newPolygon1.push_back(intersection1Id); // aggiungo id di intersezione1
        unsigned int edgeId1 = mesh.numCell1D++;
        newEdges1.push_back(edgeId1);
        mesh.Cell1DId.push_back(edgeId1);
        mesh.boundaryCell0DId.push_back({vertexIntersection1.first,intersection1Id});
        newPolygon1.push_back(intersection2Id);
        newEdges1.push_back(edgeId); // sarebbe il lato traccia


        // Controllo se dopo intersezione 2 c'è v1, altrimenti aggiungo v3
        if (vertexIntersection2.first != newPolygon1[0]) {
            newPolygon1.push_back(vertexIntersection2.first);
            unsigned int edgeId2 = mesh.numCell1D++;
            mesh.Cell1DId.push_back(edgeId2);
            newEdges2.push_back(edgeId2);
            mesh.boundaryCell0DId.push_back({intersection2Id,vertexIntersection2.first});
            unsigned int next = vertexIntersection2.first - 1;
            bool found = false;
            while (!found) {
                // PARTE DA RIVEDERE E BISOGNA AGGIUNGERE I LATI RICORDATI
                // Recupera l'ID del vertice effettivo dal modulo
                unsigned int nextVertexId = (next % mesh.Cell2DVerticesId[0].size()) - 1;
                if (mesh.Cell2DVerticesId[0][nextVertexId] == newPolygon1[0]) {
                    found = true;
                    unsigned int edgeId3 = mesh.numCell1D++;
                    mesh.Cell1DId.push_back(edgeId3);
                    newEdges1.push_back(edgeId3);
                    mesh.boundaryCell0DId.push_back({intersection2Id,vertexIntersection2.first});

                    break;
                } else {
                    newPolygon1.push_back(mesh.Cell2DVerticesId[0][nextVertexId]);
                    for(unsigned int l=0;l<mesh.boundaryCell0DId.size();l++)
                    {
                        array<unsigned int, 2> tempArray = {mesh.Cell2DVerticesId[0][nextVertexId - 1], mesh.Cell2DVerticesId[0][nextVertexId]};

                        // Confronta gli array
                        if (mesh.boundaryCell0DId[l] == tempArray) {
                            unsigned int indice = l;
                            newEdges1.push_back(indice);
                        }
                    }

                    next++;
                    break;  // No need to continue iterating after finding intersections
                }

            }
        }

        newPolygon2.push_back(vertexIntersection1.second);// v2
        unsigned int edgeId4 = mesh.numCell1D++;
        newEdges2.push_back(edgeId4);
        mesh.Cell1DId.push_back(edgeId4);
        mesh.boundaryCell0DId.push_back({vertexIntersection1.second,intersection1Id});
        newPolygon2.push_back(intersection1Id);
        newPolygon2.push_back(intersection2Id);
        newEdges2.push_back(edgeId); // lato traccia

        if (vertexIntersection2.second != newPolygon2[0]) {
            newPolygon2.push_back(vertexIntersection2.second); // aggiungo v4
            unsigned int edgeId22 = mesh.numCell1D++;
            mesh.Cell1DId.push_back(edgeId22);
            newEdges2.push_back(edgeId22);
            mesh.boundaryCell0DId.push_back({intersection2Id,vertexIntersection2.second});
            unsigned int next2 = vertexIntersection2.second + 1;
            bool found2 = false;
            while (!found2) {
                // PARTE DA RIVEDERE E BISOGNA AGGIUNGERE I LATI RICORDATI
                // Recupera l'ID del vertice effettivo dal modulo
                unsigned int nextVertexId2 = (next2 % mesh.Cell2DVerticesId[0].size()) + 1;
                if (mesh.Cell2DVerticesId[0][nextVertexId2] == newPolygon2[0]) {
                    found2 = true;
                    for(unsigned int l=0;l<mesh.boundaryCell0DId.size();l++)
                    {
                        array<unsigned int, 2> tempArray = {mesh.Cell2DVerticesId[0][nextVertexId2 - 1], mesh.Cell2DVerticesId[0][nextVertexId2]};
                        // Confronta gli array
                        if (mesh.boundaryCell0DId[l] == tempArray) {
                            unsigned int indice = l;
                            newEdges1.push_back(indice);
                        }
                    }

                    break;
                } else {

                    newPolygon2.push_back(mesh.Cell2DVerticesId[0][next2]);
                    next2++;


                }

            }
        }
        PolygonalMesh mesh1, mesh2;
        mesh1.numCell0D = newPolygon1.size();
        for (auto id : newPolygon1) {
            mesh1.Cell0DCoordinates.push_back(mesh.Cell0DCoordinates[id]);
            mesh1.Cell0DId.push_back(id);
        }
        mesh1.numCell1D = newEdges1.size();
        mesh1.Cell1DId = newEdges1;
        for (auto edge : newEdges1) {
            mesh1.boundaryCell0DId.push_back(mesh.boundaryCell0DId[edge]);
        }
        mesh1.numCell2D = 1;
        mesh1.Cell2DVerticesId.push_back(newPolygon1);
        mesh1.Cell2DEdgesId.push_back(newEdges1);

        mesh2.numCell0D = newPolygon2.size();
        for (auto id : newPolygon2) {
            mesh2.Cell0DCoordinates.push_back(mesh.Cell0DCoordinates[id]);
            mesh2.Cell0DId.push_back(id);
        }
        mesh2.numCell1D = newEdges2.size();
        mesh2.Cell1DId = newEdges2;
        for (auto edge : newEdges2) {
            mesh2.boundaryCell0DId.push_back(mesh.boundaryCell0DId[edge]);
        }
        mesh2.numCell2D = 1;
        mesh2.Cell2DVerticesId.push_back(newPolygon2);
        mesh2.Cell2DEdgesId.push_back(newEdges2);

        newMeshes.push_back(mesh1);
        newMeshes.push_back(mesh2);
    } else {
        cerr << "Non sono state trovate due intersezioni" << endl;

    }

    */





    return ritorno;

}



PolygonalMesh processFracture(DFN& dfn, Fracture& fracture) {

    totVertices = fracture.num_vertices;
    //cout << " START " << endl;
    // Import the mesh associated with the fracture

    vector<pair<Fracture,vector<Trace>>> underFractures;
    vector<pair<Fracture,vector<Trace>>> underFracturesNew;


    //TROviamo le tracce proprie
    vector<Trace> v;
    for (unsigned int i = 0; i< fracture.Traces.size(); i++){
        v.push_back(dfn.Traces[fracture.Traces[i]]);
    }
    pair<Fracture,vector<Trace>> coppia ;
    coppia.first = fracture;
    coppia.second = v;


    //determiniamo per ogni traccia se è passante o non passante per questa frattura
    for (Trace & t: coppia.second){
        if (t.AssociatedFractures[0].first == coppia.first.FractureId)
            t.passanteTemp = t.AssociatedFractures[0].second;
        else
            t.passanteTemp = t.AssociatedFractures[1].second;
    }

    //inizializziamo il vettore dei vertici
    for (unsigned int i = 0; i< coppia.first.num_vertices; i++){
        coppia.first.VerticesId.push_back(i);
    }

    underFractures.push_back(coppia);

    bool repeat = true;

    while(repeat){

        for (unsigned int i = 0; i< underFractures.size(); i++){
            repeat = false;
            vector<pair<Fracture,vector<Trace>>> coppiaTrovata;
            if (underFractures[i].second.size()){
                coppiaTrovata = cutTrace(underFractures[i]);
                for (pair<Fracture,vector<Trace>> & i: coppiaTrovata){
                    underFracturesNew.push_back(i);
                }
                repeat = true;
            }
            if (!repeat){
                underFracturesNew.push_back(underFractures[i]);
            }



        }
        underFractures = underFracturesNew;
        underFracturesNew.clear();
    }

   // cout<< "finit";
    //creiamo l'oggetto Polygonal Mesh
    PolygonalMesh mesh;

    //Instanziamo le celle 0d
    mesh.numCell0D = totVertices;
    mesh.Cell0DId.reserve(totVertices);
    mesh.Cell0DCoordinates.reserve(totVertices);
    unsigned int maxVert = 0;
    while(maxVert<totVertices){

        for (unsigned int j = 0; j<underFractures.size();j ++){
            for (unsigned int i = 0; i < underFractures[j].first.num_vertices; i++){
                if (underFractures[j].first.VerticesId[i] ==  maxVert){
                    mesh.Cell0DId.push_back(underFractures[j].first.VerticesId[i]);
                    mesh.Cell0DCoordinates.push_back(underFractures[j].first.Vertices[i]);

                }
            }
        }
        maxVert++;
    }

    //Instanziamo le celle 1d
    mesh.numCell1D = mesh.numCell0D*(mesh.numCell0D -1)/2;
    unsigned int count = 0;
    for (unsigned int i=0; i<mesh.numCell0D; i++){
        for (unsigned int j = i+1; j<mesh.numCell0D; j++){
            mesh.Cell1DId.push_back(count);
            array<unsigned int, 2> a;
            a[0] = i;
            a[1] = j;
            mesh.boundaryCell0DId.push_back(a);
            count++;
        }
    }

    //instanziamo le celle 2d
    mesh.numCell2D = underFractures.size();
    for (pair<Fracture,vector<Trace>> & u : underFractures){
        vector<unsigned int> Cell2DVerticeId;
        vector<unsigned int> Cell2DEdgeId;
        for (unsigned int k =0; k< u.first.num_vertices; k++){
            int i = u.first.VerticesId[k];
            int j = u.first.VerticesId[(k+1)%u.first.num_vertices];
            Cell2DVerticeId.push_back(i);
            unsigned int idCell2d;
            if (j>i)
                idCell2d = i*(u.first.num_vertices -1) - i*(i-1)/2 + j-i -1;
            else
                idCell2d = j*(u.first.num_vertices -1) - j*(j-1)/2 + i -j -1;
            Cell2DEdgeId.push_back(idCell2d);


        }
        mesh.Cell2DVerticesId.push_back(Cell2DVerticeId);
        mesh.Cell2DEdgesId.push_back(Cell2DEdgeId);
    }




    /*
    // Find and sort the associated traces
    vector<Trace> associatedTraces = findAndSortAssociatedTraces(dfn, fracture);


    // Process passing traces
    vector<PolygonalMesh> currentMeshes = {mesh};
    for (const Trace& trace : associatedTraces) {
        if (trace.AssociatedFractures[0].second) {
            vector<PolygonalMesh> newMeshes;
            for (PolygonalMesh& currentMesh : currentMeshes) {
                cout << "DENTRO CICLO SOTTO MESH" << endl;
                vector<PolygonalMesh> cutMeshes = cutTrace(currentMesh, trace);
                newMeshes.insert(newMeshes.end(), cutMeshes.begin(), cutMeshes.end());
            }
            currentMeshes = newMeshes;
        }
        else
        {
            // The trace is non-passing, so it needs to be extended
            Trace extendedTrace = extendTrace(trace, mesh);
            vector<PolygonalMesh> newMeshes;
            for (PolygonalMesh& currentMesh : currentMeshes) {
                vector<PolygonalMesh> cutMeshes = cutTrace(currentMesh, extendedTrace);
                newMeshes.insert(newMeshes.end(), cutMeshes.begin(), cutMeshes.end());
            }
            currentMeshes = newMeshes;
        }
    }


    // Process the final set of cut meshes as needed
    for (const PolygonalMesh& cutMesh : currentMeshes) {
        cout << "Processed mesh with " << cutMesh.numCell0D << " vertices, "
             << cutMesh.numCell1D << " edges, and "
             << cutMesh.numCell2D << " faces." << endl;
        // Perform further processing or update the mesh here
    }*/
    return mesh;
}






vector<PolygonalMesh> Alg_process_cut(DFN& dfn)
{
    vector<PolygonalMesh> finalMesh;
    dfnGlobal = & dfn;

    // Itera su ogni frattura
    for (Fracture& fracture : dfn.Fractures)
    {
        PolygonalMesh currentMesh = processFracture(dfn, fracture);

        finalMesh.push_back(currentMesh);

    }
    return finalMesh;
}


}
