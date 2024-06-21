#include "AlgFractures.hpp"
#include <cmath>

#define TOLL  1.00e-10


using namespace DFN;
using namespace std;

namespace DFN{


bool find_trace_coplanar(Fracture & f1, Fracture & f2, DFN & dfn){

    /*Cicliamo sugli n vertici di f1. Per ogni vertice P[i],  si calcola il prodotto scalare tra il vettore formato da P[i] e P[(i+1)%n] e i vettori formati da
     * P[i] e il j-esimo vertice Q[j] di f2. Se calcolando i prodotti scalari, qualcuno di questi due vettori risulta parallelo, significa che
     * Q[j] si trova sulla retta formata da P[i] e P[(i+1)%n]. Allora calcoliamo anche P[i] - Q[(j+1)%m]. Se anche esso risulta parallelo a
     * P[i] - P[(i+1)%n], abbiamo trovato due lati dei due poligoni che giacciono sulla stessa retta. A questo punto applichiamo il passo 5
     * di find_trace(...) per capire se oltre a trovarsi sulla stessa retta, c'è anche intersezione tra i lati e se si, in quali punti.
     * */

    auto n1 = f1.num_vertices;
    auto n2 = f2.num_vertices;

    for (unsigned int i =0; i < n1  ; i++){

        array<double,3> edge_f1 = substract(f1.Vertices[(i+1)%n1], f1.Vertices[i]);
        double len1 = distance_squared(f1.Vertices[(i+1)%n1], f1.Vertices[i]);

        for (unsigned int j=0; j<n2 ; j++){

            array<double,3> edge_f2 = substract(f2.Vertices[j], f1.Vertices[i]);
            double p = scalar_product(edge_f1,edge_f2);
            double len2 = distance_squared(f2.Vertices[j], f1.Vertices[i]);

            if (is_equal(p*p, len1*len2))// in tale caso i due vettori sono paralleli
            {
                array<double,3> edge_next = substract(f2.Vertices[(j+1)%n2], f1.Vertices[i]);
                p = scalar_product(edge_f1,edge_next);
                len2 = distance_squared(f2.Vertices[(j+1)%n2], f1.Vertices[i]);


                if (is_equal(p*p, len1*len2))//abbiamo trovato due lati che giacciono sulla stessa retta.
                {
                    len1 = sqrt(len1);
                    rescale(1/len1, edge_f1);
                    auto diff = substract(f2.Vertices[j],f1.Vertices[i]);
                    double I1s = scalar_product(diff,edge_f1);

                    diff = substract(f2.Vertices[(j+1)%n2],f1.Vertices[i]);
                    double I2s = scalar_product(diff,edge_f1);

                    if (I1s<0){
                        if (I2s<0)
                            return false;
                        else if (I2s<len1){
                            Trace t;
                            t.TraceId = dfn.number_traces;
                            t.Coordinates[0] = f1.Vertices[i];
                            t.Coordinates[1] = f2.Vertices[(j+1)%n2];
                            t.AssociatedFractures[0].first = f1.FractureId;
                            t.AssociatedFractures[0].second = false;
                            t.AssociatedFractures[1].first = f2.FractureId;
                            t.AssociatedFractures[1].second = false;
                            t.TraceLength = sqrt(distance_squared(f1.Vertices[i],f2.Vertices[(j+1)%n2]));
                            dfn.Traces.push_back(t);
                            dfn.number_traces = dfn.number_traces + 1;
                            f1.Traces.push_back(t.TraceId);
                            f2.Traces.push_back(t.TraceId);

                            //cout<< "found!"<<endl;
                            return true;


                        }
                        else{
                            Trace t;
                            t.TraceId = dfn.number_traces;
                            t.Coordinates[0] = f1.Vertices[i];
                            t.Coordinates[1] = f1.Vertices[(i+1)%n1];
                            t.AssociatedFractures[0].first = f1.FractureId;
                            t.AssociatedFractures[0].second = true;
                            t.AssociatedFractures[1].first = f2.FractureId;
                            t.AssociatedFractures[1].second = false;
                            t.TraceLength = sqrt(distance_squared(f1.Vertices[i],f1.Vertices[(i+1)%n1]));
                            dfn.Traces.push_back(t);
                            dfn.number_traces = dfn.number_traces + 1;
                            f1.Traces.push_back(t.TraceId);
                            f2.Traces.push_back(t.TraceId);

                            //cout<< "found!"<<endl;
                            return true;


                        }

                    }
                    else if (I1s< len1){

                        if (I2s<0){
                            Trace t;
                            t.TraceId = dfn.number_traces;
                            t.Coordinates[0] = f1.Vertices[i];
                            t.Coordinates[1] = f2.Vertices[j];
                            t.AssociatedFractures[0].first = f1.FractureId;
                            t.AssociatedFractures[0].second = false;
                            t.AssociatedFractures[1].first = f2.FractureId;
                            t.AssociatedFractures[1].second = false;
                            t.TraceLength = sqrt(distance_squared(f1.Vertices[i],f2.Vertices[j]));
                            dfn.Traces.push_back(t);
                            dfn.number_traces = dfn.number_traces + 1;
                            f1.Traces.push_back(t.TraceId);
                            f2.Traces.push_back(t.TraceId);

                            //cout<< "found!"<<endl;
                            return true;

                        }
                        else if (I2s<len1){
                            Trace t;
                            t.TraceId = dfn.number_traces;
                            t.Coordinates[0] = f2.Vertices[j];
                            t.Coordinates[1] = f2.Vertices[(j+1)%n2];
                            t.AssociatedFractures[0].first = f1.FractureId;
                            t.AssociatedFractures[0].second = false;
                            t.AssociatedFractures[1].first = f2.FractureId;
                            t.AssociatedFractures[1].second = true;
                            t.TraceLength = sqrt(distance_squared(f2.Vertices[j],f2.Vertices[(j+1)%n2]));
                            dfn.Traces.push_back(t);
                            dfn.number_traces = dfn.number_traces + 1;
                            f1.Traces.push_back(t.TraceId);
                            f2.Traces.push_back(t.TraceId);

                            //cout<< "found!"<<endl;
                            return true;

                        }
                        else{
                            Trace t;
                            t.TraceId = dfn.number_traces;
                            t.Coordinates[0] = f2.Vertices[j];
                            t.Coordinates[1] = f1.Vertices[(i+1)%n1];
                            t.AssociatedFractures[0].first = f1.FractureId;
                            t.AssociatedFractures[0].second = false;
                            t.AssociatedFractures[1].first = f2.FractureId;
                            t.AssociatedFractures[1].second = false;
                            t.TraceLength = sqrt(distance_squared(f2.Vertices[j],f1.Vertices[(i+1)%n1]));
                            dfn.Traces.push_back(t);
                            dfn.number_traces = dfn.number_traces + 1;
                            f1.Traces.push_back(t.TraceId);
                            f2.Traces.push_back(t.TraceId);

                            //cout<< "found!"<<endl;
                            return true;


                        }

                    }
                    else if (I1s >= len1){

                        if (I2s<0){
                            Trace t;
                            t.TraceId = dfn.number_traces;
                            t.Coordinates[0] = f1.Vertices[i];
                            t.Coordinates[1] = f1.Vertices[(i+1)%n1];
                            t.AssociatedFractures[0].first = f1.FractureId;
                            t.AssociatedFractures[0].second = true;
                            t.AssociatedFractures[1].first = f2.FractureId;
                            t.AssociatedFractures[1].second = false;
                            t.TraceLength = sqrt(distance_squared(f1.Vertices[i],f1.Vertices[(i+1)%n1]));
                            dfn.Traces.push_back(t);
                            dfn.number_traces = dfn.number_traces + 1;
                            f1.Traces.push_back(t.TraceId);
                            f2.Traces.push_back(t.TraceId);

                            //cout<< "found!"<<endl;
                            return true;
                        }
                        else if (I2s<len1){//e se I1s==I2s??
                            Trace t;
                            t.TraceId = dfn.number_traces;
                            t.Coordinates[0] = f2.Vertices[(j+1)%n2];
                            t.Coordinates[1] = f1.Vertices[(i+1)%n1];
                            t.AssociatedFractures[0].first = f1.FractureId;
                            t.AssociatedFractures[0].second = false;
                            t.AssociatedFractures[1].first = f2.FractureId;
                            t.AssociatedFractures[1].second = false;
                            t.TraceLength = sqrt(distance_squared(f2.Vertices[(j+1)%n2],f1.Vertices[(i+1)%n1]));
                            dfn.Traces.push_back(t);
                            dfn.number_traces = dfn.number_traces + 1;
                            f1.Traces.push_back(t.TraceId);
                            f2.Traces.push_back(t.TraceId);

                            //cout<< "found!"<<endl;
                            return true;



                        }
                        else
                            return false;


                    }

                }
            }


        }


    }
    return false;
}





bool find_trace(Fracture & f1, Fracture & f2, DFN & dfn){
    /*
     * Il primo controllo che viene fatto, è vedere se i due poligoni sono coplanari. In tale caso la possibile traccia comune, siccome i poligoni
     * non si possono sovrapporre, sarà il lato di uno dei due. La ricerca di questo lato è svolta dalla funzione find_trace_coplanar()
     * altrimenti l'algoritmo è svolto nei seguenti cinque passi:
    Passo 1: si prende il primo vertice di f1 che non sia nel piano di f2 (cioe se il primo vertice è nel piano di f2, si prende il
    secondo, esso non può essere adiacente, altrimenti le due fratture giacciono sullo stesso piano). Sia questo vertice P*.
    Passo 2: Si classificano i vertici di f1 in base a se stanno sullo stesso lato di P* rispetto al piano di f2 o meno. In particolare si trovano il primo vertice
    che passi sull'altro lato e il primo che ritorni sullo stesso lato di P* (in caso P*stesso).
    Passo 3: si trovano le coordinate dell'intersezione dei lati del poligono f1 con il piano del poligono f2. Esso sarà un segmento individuato dai
    vertici I1 I2.
    Passo4: trovato tale segmento, si classificano i  vertici di f2, in base a quale lato stanno della retta formata dal segmento rispetto
    al primo vertice. In particolare si trova il primo vertice che passi sull'altro lato e il primo che ritorni sullo stesso lato del segmento
    Passo5: si trovano le intersezioni dei lati del poligono f2 con i vertici del segmento. In questo modo si può capire se il segmento è
    passante o non passante e ottenere la traccia.
*/

    // controlliamo coplanarita (se i vettori normali sono paralleli, le fratture sono coplanari per il fatto che in caso contrario la possibilità
    //di intersezione è già stata esclusa in Exclude(). Se Exclude non viene utilizzato, bisogna aggiungere la possibilità di fratture parallele
    //ma non coplanari
    double h = scalar_product(f1.NormalVector,f2.NormalVector);
    if (is_equal(abs(h),1))
        return find_trace_coplanar(f1,f2,dfn);


    //cout<<"step one"<< endl;

    array<array<double,3>,2> Intersections; //le due intersizioni dei lati di f1 con il piano di f2
    array<double,3> * Intersection_pointers[2];//Per tener traccia se si sono trovati I1 e I2 o ancora no
    Intersection_pointers[0] = nullptr;
    Intersection_pointers[1] = nullptr;

    unsigned int taglio[2] = {0,0};

    auto diff = substract(f1.Vertices[0], f2.Baricenter);
    double p = scalar_product(diff, f2.NormalVector);

    //ricerca del primo vertice del poligono che passi sull'altro lato di P*. Tale vertice deve essere il primo a soddisfare
    //  p*(scalar_product(substract(f1.Vertices[1], f2.Baricenter),f2.NormalVector)) < 0
    if (is_equal(p,0)){// in questo caso il primo vertice del primo poligono è sul piano del secondo poligono
        taglio[0] = 0;
        Intersections[0] = f1.Vertices[0];
        Intersection_pointers[0] = &Intersections[0];
        diff = substract(f1.Vertices[1], f2.Baricenter);
        p = scalar_product(diff, f2.NormalVector);
    }
    else{
        for (unsigned int k = 1; k<f1.num_vertices;k++){
            diff = substract(f1.Vertices[k],f2.Baricenter);
            double lato = p*(scalar_product(diff,f2.NormalVector));

            if (lato <= TOLL){
                if (is_equal(lato,0)){//abbiamo trovato un vertice che sia esattamente sul piano di f2
                    Intersections[0] = f1.Vertices[k];
                    Intersection_pointers[0] = &Intersections[0];
                    k = f1.num_vertices;//!
                    continue;
                } else {
                    taglio[0] = k;
                    k = f1.num_vertices;
                    continue;
                }
            }
        }
        if (!taglio[0]){// i due poligoni non si intersecano
            return false;
        }
    }


    //cout<<"step two"<<endl;
    //ricerva del primo vertice di f1 che ritorni sullo stesso lato di P*. Esso deve avere indice maggiore di taglio1 e deve essere il primo a soddisfare:
    //  p*(scalar_product(substract(f1.Vertices[1], f2.Baricenter),f2.NormalVector)) > 0
    for (unsigned int k = taglio[0] +1; k<f1.num_vertices; k++){
        diff = substract(f1.Vertices[k],f2.Baricenter);
        double lato = p*(scalar_product(diff,f2.NormalVector));
        if (lato >= -TOLL){
            if (is_equal(lato,0)){//abbiamo trovato un vertice che giace sul piano di f2
                Intersections[1] = f1.Vertices[k];
                Intersection_pointers[1] = &Intersections[1];
                k = f1.num_vertices;//!
                continue;
            } else {
                taglio[1] = k;
                k = f1.num_vertices;
                continue;
            }
        }
    }
    if (!taglio[1]) // in questo caso il primo vertice a ritornare sul lato di P* e P* stesso.
        taglio[1] = f1.num_vertices;

    //se si hanno due punti, uno su un lato di un piano, l'altro sull'altro lato, si può trovare il punto di intersezione del segmento che
    //congiunge i due vertici con il piano stesso. Noi avendo classificato i vertici del poligono in base a quale lato stanno del piano,
    //è quello che andiamo a fare ora: trovare l'intersezione dei segmenti unenti i vertici di f2 (cioè i lati di f2) con il piano di f2

    //cout<<"step three"<<endl;
    for (unsigned int i = 0; i<2;i++){
        if (Intersection_pointers[i] == nullptr){

            diff = substract(f1.Vertices[taglio[i] - 1], f2.Baricenter);
            double b1 = abs(scalar_product(diff,f2.NormalVector));

             //prendiamo il modulo perchè a taglio[1] potrebbe essere assegnato il valore f1.num_vertices
            diff = substract(f1.Vertices[taglio[i]%f1.num_vertices], f2.Baricenter);
            double b2 = abs(scalar_product(diff,f2.NormalVector));

            //double d = sqrt(distance_squared(f1.Vertices[taglio[i] - 1],f1.Vertices[taglio[i]%f1.num_vertices]));

            double alpha = 1/(b1/b2 + 1);

            Intersections[i] = substract(f1.Vertices[taglio[i] - 1], f1.Vertices[taglio[i]%f1.num_vertices]);
            rescale(alpha,Intersections[i]);
            Intersections[i] = add(Intersections[i], f1.Vertices[taglio[i]%f1.num_vertices]);


        }
    }

    //Passo4
    //cout<<"step four"<<endl;

    auto versor_intersection = substract(Intersections[1],Intersections[0]);// n e il vettore normale al segmento delle intersezioni
    double modulus_intersection = sqrt(distance_squared(Intersections[1],Intersections[0]));
    rescale(1/modulus_intersection,versor_intersection);
    auto n = vector_product(f2.NormalVector,versor_intersection);

    array<array<double,3>,2> IntersectionS; //S per segmento
    array<double,3> * IntersectionS_pointers[2];
    IntersectionS_pointers[0] = nullptr;
    IntersectionS_pointers[1] = nullptr;

    unsigned int taglioS[2] = {0,0};

    diff = substract(f2.Vertices[0], Intersections[0]);
    p = scalar_product(diff, n);

    if (is_equal(p,0)){// in questo caso il primo vertice del primo poligono è sul piano del secondo poligono
        taglioS[0] = 0;
        IntersectionS[0] = f2.Vertices[0];
        IntersectionS_pointers[0] = &IntersectionS[0];
        diff = substract(f2.Vertices[1], Intersections[0]);
        p = - scalar_product(diff, n);
    }  else{
        for (unsigned int k = 1; k<f2.num_vertices;k++){
            diff = substract(f2.Vertices[k],Intersections[0]);
            double lato = p*(scalar_product(diff,n));

            if (lato <= TOLL){
                if (is_equal(lato,0)){//abbiamo trovato un vertice che sia esattamente sul piano di f2
                    IntersectionS[0] = f2.Vertices[k];
                    IntersectionS_pointers[0] = &IntersectionS[0];
                    k = f2.num_vertices;//!
                    continue;
                } else {
                    taglioS[0] = k;
                    k = f2.num_vertices;
                    continue;
                }
            }
        }
        if (!taglioS[0]){// i due poligoni non si intersecano
            return false;
        }
    }

    for (unsigned int k = taglioS[0] +1; k<f2.num_vertices; k++){
        diff = substract(f2.Vertices[k],Intersections[0]);
        double lato = p*(scalar_product(diff,n));
        if (lato >= -TOLL){
            if (is_equal(lato,0)){//abbiamo trovato un vertice che giace sul piano di f2
                IntersectionS[1] = f2.Vertices[k];
                IntersectionS_pointers[1] = &IntersectionS[1];
                k = f2.num_vertices;//!
                continue;
            } else {
                taglioS[1] = k;
                k = f2.num_vertices;
                continue;
            }
        }
    }

    if (!taglioS[1]) // in questo caso il primo vertice a ritornare sul lato di P* e P* stesso.
        taglioS[1] = f2.num_vertices;


    for (unsigned int i = 0; i<2;i++){
        if (IntersectionS_pointers[i] == nullptr){

            diff = substract(f2.Vertices[taglioS[i] - 1], Intersections[0]);
            double b1 = abs(scalar_product(diff,n));

            //prendiamo il modulo perchè a taglio[1] potrebbe essere assegnato il valore f1.num_vertices
            diff = substract(f2.Vertices[taglioS[i]%f2.num_vertices], Intersections[0]);
            double b2 = abs(scalar_product(diff,n));

            //double d = sqrt(distance_squared(f2.Vertices[taglioS[i] - 1],f2.Vertices[taglioS[i]%f2.num_vertices]));

            double alpha = 1/(b1/b2 + 1);

            IntersectionS[i] = substract(f2.Vertices[taglioS[i] - 1], f2.Vertices[taglioS[i]%f2.num_vertices]);
            rescale(alpha,IntersectionS[i]);
            IntersectionS[i] = add(IntersectionS[i], f2.Vertices[taglioS[i]%f2.num_vertices]);


        }
    }

    /*
    ora bisogna individuare i vertici finali delle intersezione tra i poligoni, che saranno due tra i quattro vertici Intersections[0,1] e
    intersectionS[0,1]
    */

    //cout<< "step five"<<endl;

    diff = substract(IntersectionS[0],Intersections[0]);
    double I1s = scalar_product(diff,versor_intersection);

    diff = substract(IntersectionS[1],Intersections[0]);
    double I2s = scalar_product(diff,versor_intersection);

    if (is_equal(I1s,0) && is_equal(I2s,modulus_intersection)){
        Trace t;
        t.TraceId = dfn.number_traces;
        t.Coordinates[0] = Intersections[0];
        t.Coordinates[1] = Intersections[1];
        t.AssociatedFractures[0].first = f1.FractureId;
        t.AssociatedFractures[0].second = true;
        t.AssociatedFractures[1].first = f2.FractureId;
        t.AssociatedFractures[1].second = true;
        t.TraceLength = sqrt(distance_squared(IntersectionS[1],Intersections[0]));
        dfn.Traces.push_back(t);
        dfn.number_traces = dfn.number_traces + 1;
        f1.Traces.push_back(t.TraceId);
        f2.Traces.push_back(t.TraceId);

        return true;
    }


    if (I1s<0){
        if (I2s<0)
            return false;
        else if (I2s<modulus_intersection){
            Trace t;
            t.TraceId = dfn.number_traces;
            t.Coordinates[0] = Intersections[0];
            t.Coordinates[1] = IntersectionS[1];
            t.AssociatedFractures[0].first = f1.FractureId;
            t.AssociatedFractures[0].second = false;
            t.AssociatedFractures[1].first = f2.FractureId;
            t.AssociatedFractures[1].second = false;
            t.TraceLength = sqrt(distance_squared(IntersectionS[1],Intersections[0]));
            dfn.Traces.push_back(t);
            dfn.number_traces = dfn.number_traces + 1;
            f1.Traces.push_back(t.TraceId);
            f2.Traces.push_back(t.TraceId);

            //cout<< "found!"<<endl;


        }
        else{
            Trace t;
            t.TraceId = dfn.number_traces;
            t.Coordinates[0] = Intersections[0];
            t.Coordinates[1] = Intersections[1];
            t.AssociatedFractures[0].first = f1.FractureId;
            t.AssociatedFractures[0].second = true;
            t.AssociatedFractures[1].first = f2.FractureId;
            t.AssociatedFractures[1].second = false;
            t.TraceLength = sqrt(distance_squared(Intersections[1],Intersections[0]));
            dfn.Traces.push_back(t);
            dfn.number_traces = dfn.number_traces + 1;
            f1.Traces.push_back(t.TraceId);
            f2.Traces.push_back(t.TraceId);

            //cout<< "found!"<<endl;


        }

    }
    else if (I1s< modulus_intersection){

        if (I2s<0){
            Trace t;
            t.TraceId = dfn.number_traces;
            t.Coordinates[0] = Intersections[0];
            t.Coordinates[1] = IntersectionS[0];
            t.AssociatedFractures[0].first = f1.FractureId;
            t.AssociatedFractures[0].second = false;
            t.AssociatedFractures[1].first = f2.FractureId;
            t.AssociatedFractures[1].second = false;
            t.TraceLength = sqrt(distance_squared(IntersectionS[0],Intersections[0]));
            dfn.Traces.push_back(t);
            dfn.number_traces = dfn.number_traces + 1;
            f1.Traces.push_back(t.TraceId);
            f2.Traces.push_back(t.TraceId);

           //cout<< "found!"<<endl;

        }
        else if (I2s<modulus_intersection){//e se I1s==I2s??
            Trace t;
            t.TraceId = dfn.number_traces;
            t.Coordinates[0] = IntersectionS[0];
            t.Coordinates[1] = IntersectionS[1];
            t.AssociatedFractures[0].first = f1.FractureId;
            t.AssociatedFractures[0].second = false;
            t.AssociatedFractures[1].first = f2.FractureId;
            t.AssociatedFractures[1].second = true;
            t.TraceLength = sqrt(distance_squared(IntersectionS[1],IntersectionS[0]));
            dfn.Traces.push_back(t);
            dfn.number_traces = dfn.number_traces + 1;
            f1.Traces.push_back(t.TraceId);
            f2.Traces.push_back(t.TraceId);

            //cout<< "found!"<<endl;

        }
        else{
            Trace t;
            t.TraceId = dfn.number_traces;
            t.Coordinates[0] = IntersectionS[0];
            t.Coordinates[1] = Intersections[1];
            t.AssociatedFractures[0].first = f1.FractureId;
            t.AssociatedFractures[0].second = false;
            t.AssociatedFractures[1].first = f2.FractureId;
            t.AssociatedFractures[1].second = false;
            t.TraceLength = sqrt(distance_squared(Intersections[1],IntersectionS[0]));
            dfn.Traces.push_back(t);
            dfn.number_traces = dfn.number_traces + 1;
            f1.Traces.push_back(t.TraceId);
            f2.Traces.push_back(t.TraceId);

            //cout<< "found!"<<endl;


        }

    }
    else if (I1s >= modulus_intersection){

        if (I2s<0){
            Trace t;
            t.TraceId = dfn.number_traces;
            t.Coordinates[0] = Intersections[0];
            t.Coordinates[1] = Intersections[1];
            t.AssociatedFractures[0].first = f1.FractureId;
            t.AssociatedFractures[0].second = true;
            t.AssociatedFractures[1].first = f2.FractureId;
            t.AssociatedFractures[1].second = false;
            t.TraceLength = sqrt(distance_squared(Intersections[1],Intersections[0]));
            dfn.Traces.push_back(t);
            dfn.number_traces = dfn.number_traces + 1;
            f1.Traces.push_back(t.TraceId);
            f2.Traces.push_back(t.TraceId);

            //cout<< "found!"<<endl;
        }
        else if (I2s<modulus_intersection){//e se I1s==I2s??
            Trace t;
            t.TraceId = dfn.number_traces;
            t.Coordinates[0] = IntersectionS[1];
            t.Coordinates[1] = Intersections[1];
            t.AssociatedFractures[0].first = f1.FractureId;
            t.AssociatedFractures[0].second = false;
            t.AssociatedFractures[1].first = f2.FractureId;
            t.AssociatedFractures[1].second = false;
            t.TraceLength = sqrt(distance_squared(IntersectionS[1],Intersections[1]));
            dfn.Traces.push_back(t);
            dfn.number_traces = dfn.number_traces + 1;
            f1.Traces.push_back(t.TraceId);
            f2.Traces.push_back(t.TraceId);

            //cout<< "found!"<<endl;



        }
        else
            return false;


    }



    return true;

}




//prima ottimizzazione: due poligoni non possono intersecarsi se le palle centro il loro baricentro, raggio la distanza tra il baricentro
//e il vertice piu distante da quest'utlimo non si intersecano.
//ritrorna vero se la possibilita di intersezione e da escludersi
bool Exclude1(Fracture & f1, Fracture & f2){

    double d = distance_squared(f1.Baricenter,f2.Baricenter);
    if (d < f1.Radius + f2.Radius)
        return false;
    return true;

}


bool Exclude2(Fracture & f1, Fracture & f2){

    /*
     * Prima esclusione: aprossimiamo f1 ed f2 come due sfere centro il baricentro, raggio la distanza del vertice più lontano dal baricentro.
     * Se la somma dei due raggi è minore della distanza tra i baricentri, la possibilità di intersezione è esclusa.
     * Seconda esclusione: approssimiamo la due fratture con dei cerchi centro il baricentro, raggio la distanza del vertice più lontano dal
     * baricentro. Troviamo il vettore t che identifica la direzione della retta di intersezione dei due piani formati da f1 e f2.
     * Notiamo che gli unici punti di intersezione possono avvenire su questa retta. (Se le fratture sono coplanari, la possibilità di intersezione
     * non viene esclusa.)
     * Consideriamo le proiezioni sul piano ortogonale a t dei due segmenti identificati da (Baricenter1 +-  Radius1*n1xt) e
     * (Baricenter2 +- Radius2*( n2xt), dove x sta per prodotto vettoriale, n1 ed n2 sono i versori normali ai piani dei due poligoni.
     * Se su questo piano il primo segmento non interseca la retta formata dal secondo, sicuramente i due poligoni non si intersecano.
     * Terza esclusione: se invece c'è intersezione, significa che il cerchio del primo poligono si interseca con la retta t.
     * Allora troviamo l'intersezione della retta del secondo segmento con la retta t. Se questo punto giace al di fuori del secondo segmento, non
     * c'è intersezione.
     * Quarta esclusione: se infine i segmenti ottenuti dall'intersezione della retta t con i due cerchi di f1 e f2 non si intersecano, non c'è
     * intersezione.
    */

    //prima esclusione
    double d = distance_squared(f1.Baricenter,f2.Baricenter);
    d = sqrt(d);
    if (d > f1.Radius + f2.Radius + TOLL)
        return true;

    auto t = vector_product(f1.NormalVector,f2.NormalVector);

    //caso dei segmenti coplanari
    array<double,3> origin ={0,0,0};
    double normt = distance_squared(t, origin);

    if ( normt < TOLL ){ //i piani di f1 e f2 risultano paralleli
        auto diff = substract(f1.Baricenter, f2.Baricenter);
        double h = scalar_product(diff, f1.NormalVector);

        if (is_equal(h,0)) // f1 e f2 sono coplanari, la possibilità di intersezione non viene esclusa
            return false;
        else
            return true;
    }

    //normalizziamo t
    normt = sqrt(normt);
    rescale(1/normt, t);

    //seconda esclusione
    auto h1 = vector_product(t, f1.NormalVector);
    auto h2 = vector_product(t, f2.NormalVector);

    rescale(f1.Radius, h1);
    rescale(f2.Radius, h2);

    //i due estremi del primo segmento
    auto P1_1 = add(f1.Baricenter, h1);
    auto P1_2 = substract(f1.Baricenter, h1);

    //i due estremi del secondo segmento
    auto P2_1 = add(f2.Baricenter, h2);
    auto P2_2 = substract(f2.Baricenter, h2);

    //i due estremi del primo segmento relativi al secondo
    auto P1_1_Bari2 = substract(P1_1,f2.Baricenter);
    auto P1_2_Bari2 = substract(P1_2,f2.Baricenter);

    double p1 = scalar_product(f2.NormalVector, P1_1_Bari2);
    double p2 = scalar_product(f2.NormalVector, P1_2_Bari2);

    if (p1*p2> -TOLL)//i due segmenti non si intersecano, o si toccano solo gli estremi
        return true;

    //terza esclusione
    //troviamo l'intersezione del primo segmento con la retta t
    double alpha = abs(p2)/(abs(p1)+abs(p2));
    auto I1 = substract(P1_1, P1_2);
    rescale(alpha, I1);
    I1 = add(P1_2, I1);

    //troviamo l'intersezione del secondo segmento con la retta t
    auto f2bari_I1 = substract(f2.Baricenter, I1);
    double dist_on_t = scalar_product(t,f2bari_I1); //non prendere il valore assoluto è giusto!
    auto I2 = scalar_vector_product(dist_on_t, t);
    I2 = add(I1,I2);
    double dist_I2q = distance_squared(f2.Baricenter,I2);
    if (dist_I2q > (f2.Radius)*(f2.Radius) + TOLL)
        return true;
/*
    //quarta esclusione
    double dist_I1q = distance_squared(f1.Baricenter,I1);
    double portion_of_f1_on_t = sqrt(f1.Radius*f1.Radius - dist_I1q);
    double portion_of_f2_on_t = sqrt(f2.Radius*f2.Radius - dist_I2q);

    if (abs(dist_on_t) > portion_of_f1_on_t + portion_of_f2_on_t + TOLL)
        return true;

*/
    return false;


}




bool run_Alg(DFN & dfn){

    unsigned int count = 0;

    unsigned int n = dfn.Fractures.size();
    unsigned int i = 0;
    unsigned int j = 0;
    dfn.Traces.reserve(n*(n-1));
    for (  i = 0; i<n; i++){

        j = i+1;

        for( j ; j<n;j++ ){
            if (!Exclude2(dfn.Fractures[i],dfn.Fractures[j])){
                //cout<< "computing trace of "<< i<<j<<endl;
                count++;
                find_trace(dfn.Fractures[i],dfn.Fractures[j], dfn);
            }


        }

    }
    for(int i = 0; i<dfn.Traces.size();i++){
        cout<< "Trace"<<i<<":";
        cout<<"first coordinate: "<< dfn.Traces[i].Coordinates[0][0]<<" "<<dfn.Traces[i].Coordinates[0][1]<<" "<<dfn.Traces[i].Coordinates[0][2]<<
            " second coordinate: "<<dfn.Traces[i].Coordinates[1][0]<<" "<<dfn.Traces[i].Coordinates[1][1]<<" "<<dfn.Traces[i].Coordinates[1][2]<<endl;
    }

    cout<< count<<endl;
    return true;


}



}

