#include "DFNstructs.hpp"
#include <iostream>
#include <fstream>

using namespace std;

namespace DFN{

bool ImportDFN(const string& filename,DFN& dfn)
{
    ifstream file;
    file.open(filename);

    if(file.fail())
    {
        cerr << "File not opened correctly" << endl;
        return false;
    }

    string line;
    unsigned int numberfractures;
    string delimitator;

    getline(file,line);
    getline(file,line);

    istringstream converter(line);
    converter >> numberfractures; //I extract the number of fractures in the DFN.
    dfn.FracturesNumber = numberfractures;
    cout << "The DFN has: " << numberfractures << " fractures." << endl;

    dfn.FracturesId.reserve(numberfractures);
    dfn.NumberVerties.reserve(numberfractures);

    unsigned int k = 0;
    while(k < numberfractures){

        unsigned int fractureid = 0;
        unsigned int numvertices = 0;
        Fracture fracture;

        getline(file,line);
        getline(file,line);
        istringstream converter(line);

        if(converter >> fractureid >> delimitator >> numvertices) /**I extract the id of each fracture and its
                                                                     number of vertices and I add them in the
                                                                     respective vectors.*/
        {
            dfn.FracturesId.push_back(fractureid);
            fracture.FractureId = fractureid;
            dfn.NumberVerties.push_back((numvertices));
            fracture.num_vertices = numvertices;
            fracture.Vertices.reserve(numvertices);
        }

        getline(file,line);

        vector<vector<double>> Matrix;
        for(unsigned int i = 0; i < 3; i++ ) /**This for loop extracts the coordinates of each vertice of a fracture
                                                and it adds them firstly in a matrix and then in a vector of matrices. */
        {
            getline(file,line);
            istringstream vertexStream(line);
            vector<double> vertices;

            for(unsigned int j = 0; j < numvertices; j++)
            {
                double value;
                if(vertexStream >> value)
                {
                    vertices.push_back(value);
                    cout << value << " ";
                    vertexStream >> delimitator;
                }

            }
            Matrix.push_back(vertices);
            dfn.MatricesVertices.push_back(Matrix);

            cout << endl;
        }
        cout << endl;

        for(unsigned int t = 0; t < numvertices; t++)
        {
            array<double,3> columnVertices;
            for(unsigned int w = 0; w <3; w++)
            {
                columnVertices[w] = Matrix[w][t];
            }
            //fracture.Vertices[t].reserve(3);
            fracture.Vertices.push_back(columnVertices);
        }

        cout << "The coordinates of the vertices of the fracture with id " << dfn.FracturesId[k] << " are:" << endl;
        for(unsigned int i = 0; i < fracture.Vertices.size(); i++) /**This loop adds the coordinates of the vertices of a fracture
                                                                       in a vector in the struct fracture*/
        {
            cout << "[ ";
            for(unsigned int j = 0; j < 3; j++)
            {
                cout << fracture.Vertices[i][j] << " ";
            }
            cout << "]";
        }
        cout << endl;
        cout << endl;

        cout <<"computing baricenter"<<endl;
        fracture.get_baricenter();
        dfn.Fractures.push_back(fracture);


        k++;
    }

    file.close();
    return true;
}

void PrintTraces(const DFN& dfn)
{
    ofstream outputfile("ResultsTraces.txt");

    string header;
    string second_header;
    outputfile << header << "# Traces" << endl;
    outputfile << dfn.number_traces << endl;

    for(unsigned int i = 0; i < dfn.number_traces; i++)
    {
        outputfile << second_header << "# TraceId;    FractureId1;    FractureId2;    X1;    Y1;    Z1;    X2;    Y2;    Z2" << endl;
        outputfile << "\t" <<dfn.Traces[i].TraceId << "\t" << dfn.Traces[i].AssociatedFractures[0].first << "\t" << dfn.Traces[i].AssociatedFractures[1].first;
        outputfile << "\t" << dfn.Traces[i].Coordinates[0][0] << "\t" << dfn.Traces[i].Coordinates[0][1] << "\t" << dfn.Traces[i].Coordinates[0][2] << "\t" <<dfn.Traces[i].Coordinates[1][0] << "\t" << dfn.Traces[i].Coordinates[1][1] << "\t" << dfn.Traces[i].Coordinates[1][2] << endl;
        outputfile << endl;
    }

    outputfile.close();
}

void PrintFractures(const DFN& dfn)
{
    ofstream outputfile("ResultsFractures.txt");

    outputfile << "# Number of Fractures" << endl;
    outputfile << dfn.FracturesNumber << endl;

    for(unsigned int i = 0; i < dfn.FracturesNumber; i++)
    {
        outputfile << "# FractureId;" <<"\t" << "NumTraces" << endl;
        outputfile << " \t" << dfn.Fractures[i].FractureId << ";" << "\t " << dfn.Fractures[i].Traces.size() << endl;
        outputfile << "# TraceId;" << "\t" << "Tips;" << "\t" << "Length" << endl;
        for(unsigned int j = 0; j < dfn.Fractures[i].Traces.size(); j++)
        {
            unsigned int id = dfn.Fractures[i].Traces[j];
            unsigned int tip;
            if(i == dfn.Traces[j].AssociatedFractures[0].first)
            {
                tip = dfn.Traces[j].AssociatedFractures[0].second;
                if(tip == true)
                {
                    tip = false;
                }
                else if(tip == false)
                {
                    tip = true;
                }
            }
            else
            {
                tip = dfn.Traces[j].AssociatedFractures[1].second;
                if(tip == true)
                {
                    tip = false;
                }
                else if(tip == false)
                {
                    tip = true;
                }
            }
            outputfile << "\t" << id << ";\t "<<tip << ";\t " << dfn.Traces[j].TraceLength << endl;
        }
        outputfile << endl;
    }

    outputfile.close();
}
}
