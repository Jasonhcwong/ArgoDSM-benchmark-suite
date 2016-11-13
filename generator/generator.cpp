// simple graph generator

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <unordered_set>

//#define DEBUG


/*
 *  Return a random element from an unordered set of vertices
 */
int get_random_element(std::unordered_set<int>* u_set, int num_vertices) {

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> rand_vertex(0, num_vertices-1);

    int elem;
    int stub = rand_vertex(rng);
    if (u_set->count(stub) == 0) {
        u_set->insert(stub);
        std::unordered_set<int>::iterator it = u_set->find(stub);
        if (it == u_set->end()) elem = *(u_set->begin());
        else elem = *(++it);
        u_set->erase(stub);
    }
    else elem = stub;
    return elem;
}


/*
 *  Decide whether an edge is safe to be added to the graph
 */
bool is_legit(int current, int neighbor, int max_degree, std::vector<std::vector<int>*>* graph) {
    bool legit = true;
    if (current == neighbor)
        legit = false;

    if ((int)graph->at(current)->size() == max_degree) legit = false;

    for (int k = 0; k < (int)graph->at(current)->size(); k++)
        if (graph->at(current)->at(k) == neighbor)
            legit = false;

    for (int k = 0; k < (int)graph->at(neighbor)->size(); k++)
        if (graph->at(neighbor)->at(k) == current)
            legit = false;

    return legit;
}


/*
 *  Output the graph to the desired text file
 */
void write_graph_to_file(std::string filename, int weighted, std::vector<std::vector<int>*>* graph,  std::vector<std::vector<int>*>* weights) {
    std::ofstream output(filename);
    for(int j = 0; j < (int)graph->size(); j++) {
        int degree = (int)graph->at(j)->size();
        if (degree > 0)
            for(int i = 0; i < degree; i++) {
                output << j << " " << graph->at(j)->at(i);
                if (weighted != 0)
                    output << " " << weights->at(j)->at(i);
                output << "\n";
            }
    }
    output.close();
}


/*
 *  Generate a graph file with N vertices, M edges, K max degree with weight 0-20, partitioned in L number of components
 */
int main(int argc, char** argv) {

    if (argc != 7) {
        std::cout << "./generator <num_vertices> <num_edges> <max_degree> <max_weight> <num_components> <output_file>" << std::endl;
        return 1;
    }

    int num_vertices     = atoi(argv[1]);
    int num_edges        = atoi(argv[2]);
    int max_degree       = atoi(argv[3]);
    int weighted         = atoi(argv[4]);
    int num_components   = atoi(argv[5]);
    std::string filename = argv[6];

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> rand_vertex(0, num_vertices-1);
    std::uniform_int_distribution<int> rand_weight(1, 20);


    std::vector<std::vector<int>*>* graph = new std::vector<std::vector<int>*>(num_vertices);
    std::vector<std::vector<int>*>* weights = new std::vector<std::vector<int>*>(num_vertices);
    std::unordered_set<int>* vertices = new std::unordered_set<int>();
    std::vector<std::unordered_set<int>*>* components = new std::vector<std::unordered_set<int>*>(num_components);

    for (int i = 0; i < (int)graph->size(); i++) {
        graph->at(i) = new std::vector<int>();
        weights->at(i) = new std::vector<int>();
    }

    for (int i = 0; i < num_vertices; i++)
        vertices->insert(i);

    for (int i = 0; i < num_components; i++) {
        components->at(i) = new std::unordered_set<int>();
        int elem = get_random_element(vertices, num_vertices);
        components->at(i)->insert(elem);
        vertices->erase(elem);
    }

    int current, neighbor;
    int count = 0;
    bool legit;

    while (count < num_edges) {

        int i = count % num_components;

        if (vertices->size() > 0) {
            // Connect all vertices in each component
            current = get_random_element(vertices, num_vertices);
            neighbor = get_random_element(components->at(i), num_vertices);
            legit = is_legit(current, neighbor, max_degree, graph);
            if(legit) vertices->erase(current);
        } else {
            // Insert remaining edges, not necessary to connect a component
            current = get_random_element(components->at(i), num_vertices);
            neighbor = get_random_element(components->at(i), num_vertices);
            legit = is_legit(current, neighbor, max_degree, graph);
        }

        if (legit) {
#ifdef DEBUG
            std::cout << "Edge " << count << ": components[" << i << "]: " << current << " - " << neighbor << "\n";
#endif
            components->at(i)->insert(current);
            graph->at(current)->push_back(neighbor);
            weights->at(current)->push_back(rand_weight(rng));
            count++;
        }

    }

#ifdef DEBUG
    for (int i = 0; i < num_components; i++) {
        std::cout << "\n";
        std::cout << i << ": { ";
        for (auto elem: *(components->at(i)))
            std::cout << elem << " ";
        std::cout << "}\n";
    }
#endif

    write_graph_to_file(filename, weighted, graph, weights);

    std::cout << "\nGenerated graph: ";
    std::cout << "\nVertices: " << num_vertices;
    std::cout << "\nEdges: " << num_edges;
    std::cout << "\nMax degree: " << max_degree;
    std::cout << "\nComponents: " << num_components;
    std::cout << "\nFile: " << filename;

    for (int i = 0; i < (int)graph->size(); i++) {
        graph->at(i) = new std::vector<int>();
        weights->at(i) = new std::vector<int>();
    }

    for (int i = 0; i < (int)components->size(); i++)
        delete components->at(i);

    delete vertices;
    delete components;
    delete graph;
    delete weights;

    return 0;
}
