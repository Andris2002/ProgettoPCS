#pragma once

#include "DFNstructs.hpp"
#include "ImportExport.hpp"
#include "AlgFractures.hpp"
#include "PolygonalMesh.hpp"

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <iostream>

using ::testing::ElementsAre;
using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::DoubleEq;
using ::testing::UnorderedElementsAre;

#define TOLL  1.00e-10

using namespace DFN;
using namespace std;

namespace DFN{

struct DFN;
struct Fracture;
bool ImportDFN(const string& filename,DFN& dfn);
bool find_trace(Fracture & f1, Fracture & f2, DFN & dfn);
bool run_Alg(DFN & dfn);






TEST(FINDTRACE, FR3_data){

    //controlliamo che le tracce trovate per le fratture contenute nel file siano giuste. Esse dovrebbero essere:
    array<double,3> trace1_1 = {0.8,0,0};
    array<double,3> trace1_2 = {0.8,1,0};

    array<double,3> trace2_1 = {0.316184,0.5,0};
    array<double,3> trace2_2 = {0,0.5,0};

    DFN dfn;

    string filename = "./DFN/FR3_data.txt";
    ImportDFN(filename,dfn);

    run_Alg(dfn);

    ASSERT_EQ(dfn.Traces.size(),2);
    cout<< dfn.Traces[0].Coordinates[0][0];
    EXPECT_THAT(dfn.Traces[0].Coordinates,UnorderedElementsAre(ElementsAre(
                                                                    AllOf(Gt(0.8 - TOLL), Lt(0.8+TOLL)),
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL)),
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL))),
                                                                ElementsAre(
                                                                    AllOf(Gt(0.8 - TOLL), Lt(0.8+TOLL)),
                                                                    AllOf(Gt(1 - TOLL), Lt(1+TOLL)),
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL)))));

    EXPECT_THAT(dfn.Traces[1].Coordinates,UnorderedElementsAre(ElementsAre(
                                                                    AllOf(Gt(3.1618370000000001e-01 - TOLL), Lt(3.1618370000000001e-01+TOLL)),
                                                                    AllOf(Gt(0.5 - TOLL), Lt(0.5+TOLL)),
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL))),
                                                                ElementsAre(
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL)),
                                                                    AllOf(Gt(0.5 - TOLL), Lt(0.5+TOLL)),
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL)))));


}

TEST(FINDTRACE, VerticesOnTrace){
    /*
     * Testiamo la funzione che trova le tracce di due fratture passandole due fratture, non coplanari, con una delle due avente un lato sulla
     * superficie dell'altra
     * */

    Fracture t1;
    Fracture t2;

    array<double,3> t1_v1 = {0,0,0};
    array<double,3> t1_v2 = {0,1,0};
    array<double,3> t1_v3 = {1,1,0};
    array<double,3> t1_v4 = {1,0,0};

    array<double,3> t2_v1 = {0.5,0.25,0};
    array<double,3> t2_v2 = {0.5,0.25,1};
    array<double,3> t2_v3 = {0.5,0.75,1};
    array<double,3> t2_v4 = {0.5,0.75,0};

    t1.Vertices.push_back(t1_v1);
    t1.Vertices.push_back(t1_v2);
    t1.Vertices.push_back(t1_v3);
    t1.Vertices.push_back(t1_v4);

    t2.Vertices.push_back(t2_v1);
    t2.Vertices.push_back(t2_v2);
    t2.Vertices.push_back(t2_v3);
    t2.Vertices.push_back(t2_v4);

    t1.num_vertices = 4;
    t2.num_vertices = 4;

    t1.get_baricenter();
    t2.get_baricenter();

    DFN dfn;
    find_trace(t1,t2, dfn);

    ASSERT_EQ(dfn.Traces.size(),1);

    cout<<"first coordinate: "<< dfn.Traces[0].Coordinates[0][0]<<" "<<dfn.Traces[0].Coordinates[0][1]<<" "<<dfn.Traces[0].Coordinates[0][2]<<
        " second coordinate: "<<dfn.Traces[0].Coordinates[1][0]<<" "<<dfn.Traces[0].Coordinates[1][1]<<" "<<dfn.Traces[0].Coordinates[1][2]<<endl;


    EXPECT_THAT(dfn.Traces[0].Coordinates,UnorderedElementsAre(ElementsAre(
                                                                    AllOf(Gt(0.5 - TOLL), Lt(0.5+TOLL)),
                                                                    AllOf(Gt(0.25 - TOLL), Lt(0.25+TOLL)),
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL))),
                                                                ElementsAre(
                                                                    AllOf(Gt(0.5 - TOLL), Lt(0.5+TOLL)),
                                                                    AllOf(Gt(0.75 - TOLL), Lt(0.75+TOLL)),
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL)))));




}

TEST(FINDTRACE, VerticesOnTraceCoplanar){
    /*
     * Testiamo la funzione che trova le tracce di due fratture passandole due fratture,che sono coplanari e con due lati che coincidono
     * */

    Fracture t1;
    Fracture t2;

    array<double,3> t1_v1 = {0,0,0};
    array<double,3> t1_v2 = {0,1,0};
    array<double,3> t1_v3 = {1,1,0};
    array<double,3> t1_v4 = {1,0,0};

    array<double,3> t2_v1 = {1,0.25,0};
    array<double,3> t2_v2 = {2,0.25,0};
    array<double,3> t2_v3 = {2,0.75,0};
    array<double,3> t2_v4 = {1,0.75,0};

    t1.Vertices.push_back(t1_v1);
    t1.Vertices.push_back(t1_v2);
    t1.Vertices.push_back(t1_v3);
    t1.Vertices.push_back(t1_v4);

    t2.Vertices.push_back(t2_v1);
    t2.Vertices.push_back(t2_v2);
    t2.Vertices.push_back(t2_v3);
    t2.Vertices.push_back(t2_v4);

    t1.num_vertices = 4;
    t2.num_vertices = 4;

    t1.get_baricenter();
    t2.get_baricenter();

    DFN dfn;
    find_trace(t1,t2, dfn);

    ASSERT_EQ(dfn.Traces.size(),1);

    cout<<"first coordinate: "<< dfn.Traces[0].Coordinates[0][0]<<" "<<dfn.Traces[0].Coordinates[0][1]<<" "<<dfn.Traces[0].Coordinates[0][2]<<
        " second coordinate: "<<dfn.Traces[0].Coordinates[1][0]<<" "<<dfn.Traces[0].Coordinates[1][1]<<" "<<dfn.Traces[0].Coordinates[1][2]<<endl;


    EXPECT_THAT(dfn.Traces[0].Coordinates,UnorderedElementsAre(ElementsAre(
                                                                    AllOf(Gt(1 - TOLL), Lt(1+TOLL)),
                                                                    AllOf(Gt(0.25 - TOLL), Lt(0.25+TOLL)),
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL))),
                                                                ElementsAre(
                                                                    AllOf(Gt(1 - TOLL), Lt(1+TOLL)),
                                                                    AllOf(Gt(0.75 - TOLL), Lt(0.75+TOLL)),
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL)))));




}


TEST(FINDTRACE, hexagons){
    /*
     * il test viene eseguito su due esagoni con baricentro nell'origine e lato unitario, uno giacente sul piano z=0, l'altro sul piano x=0
     * */
    double h = sqrt(3)/2;
    Fracture t1;
    Fracture t2;

    array<double,3> t1_v1 = {1,0,0};
    array<double,3> t1_v2 = {0.5,h,0};
    array<double,3> t1_v3 = {-0.5,h,0};
    array<double,3> t1_v4 = {-1,0,0};
    array<double,3> t1_v5 = {-0.5,-h,0};
    array<double,3> t1_v6 = {0.5,-h,0};

    array<double,3> t2_v1 = {0,1,0};
    array<double,3> t2_v2 = {0,0.5,h};
    array<double,3> t2_v3 = {0,-0.5,h};
    array<double,3> t2_v4 = {0,-1,0};
    array<double,3> t2_v5 = {0,-0.5,-h};
    array<double,3> t2_v6 = {0,0.5,-h};

    t1.Vertices.push_back(t1_v1);
    t1.Vertices.push_back(t1_v2);
    t1.Vertices.push_back(t1_v3);
    t1.Vertices.push_back(t1_v4);
    t1.Vertices.push_back(t1_v5);
    t1.Vertices.push_back(t1_v6);

    t2.Vertices.push_back(t2_v1);
    t2.Vertices.push_back(t2_v2);
    t2.Vertices.push_back(t2_v3);
    t2.Vertices.push_back(t2_v4);
    t2.Vertices.push_back(t2_v5);
    t2.Vertices.push_back(t2_v6);

    t1.num_vertices = 6;
    t2.num_vertices = 6;

    t1.get_baricenter();
    t2.get_baricenter();

    DFN dfn;
    find_trace(t1,t2, dfn);

    ASSERT_EQ(dfn.Traces.size(),1);

    cout<<"first coordinate: "<< dfn.Traces[0].Coordinates[0][0]<<" "<<dfn.Traces[0].Coordinates[0][1]<<" "<<dfn.Traces[0].Coordinates[0][2]<<
        " second coordinate: "<<dfn.Traces[0].Coordinates[1][0]<<" "<<dfn.Traces[0].Coordinates[1][1]<<" "<<dfn.Traces[0].Coordinates[1][2]<<endl;


    EXPECT_THAT(dfn.Traces[0].Coordinates,UnorderedElementsAre(ElementsAre(
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL)),
                                                                    AllOf(Gt(h - TOLL), Lt(h+TOLL)),
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL))),
                                                                ElementsAre(
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL)),
                                                                    AllOf(Gt(-h - TOLL), Lt(-h+TOLL)),
                                                                    AllOf(Gt(0 - TOLL), Lt(0+TOLL)))));




}


TEST(EXCLUDE2, FR10_data){

    DFN dfn1;
    string filename = "./DFN/FR10_data.txt";
    ImportDFN(filename,dfn1);

    DFN dfn2;
    ImportDFN(filename,dfn2);

    unsigned int n = 10;
    unsigned int i = 0;
    unsigned int j = 0;
    dfn1.Traces.reserve(n*(n-1));
    for ( i = 0; i<n; i++){

        j = i+1;

        for( j ; j<n;j++ ){

                find_trace(dfn1.Fractures[i],dfn1.Fractures[j], dfn1);
            }

    }

    run_Alg(dfn2);

    ASSERT_EQ(dfn1.Traces.size(),dfn2.Traces.size());

    for (unsigned int k = 0; k< dfn2.Traces.size(); k++){

        EXPECT_THAT(dfn1.Traces[k].Coordinates,UnorderedElementsAre(ElementsAre(
                            AllOf(Gt(dfn2.Traces[k].Coordinates[0][0] - TOLL), Lt(dfn2.Traces[k].Coordinates[0][0]+TOLL)),
                            AllOf(Gt(dfn2.Traces[k].Coordinates[0][1] - TOLL), Lt(dfn2.Traces[k].Coordinates[0][1]+TOLL)),
                            AllOf(Gt(dfn2.Traces[k].Coordinates[0][2] - TOLL), Lt(dfn2.Traces[k].Coordinates[0][2]+TOLL))),
                        ElementsAre(
                            AllOf(Gt(dfn2.Traces[k].Coordinates[1][0] - TOLL), Lt(dfn2.Traces[k].Coordinates[1][0]+TOLL)),
                            AllOf(Gt(dfn2.Traces[k].Coordinates[1][1] - TOLL), Lt(dfn2.Traces[k].Coordinates[1][1]+TOLL)),
                            AllOf(Gt(dfn2.Traces[k].Coordinates[1][2] - TOLL), Lt(dfn2.Traces[k].Coordinates[1][2]+TOLL)))));

    }

}


}
