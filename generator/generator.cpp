// simple graph generator

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <unordered_set>


int main(int argc, char** argv) {

    if (argc != 6) {
        std::cout << "./generator <vertices> <edges> <max_degree> <max_weight> <output_file>" << std::endl;
        return 1;
    }

    int num_vertices         = atoi(argv[1]);
    int edges            = atoi(argv[2]);
    int max_degree       = atoi(argv[3]);
    int weighted         = atoi(argv[4]);
    std::string filename = argv[5];

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> rand_degree(0, max_degree);
    std::uniform_int_distribution<int> rand_vertex(0, num_vertices-1);
    std::uniform_int_distribution<int> rand_weight(1, 20);

    std::vector<std::vector<int>*>* graph   = new std::vector<std::vector<int>*>(num_vertices);
    std::vector<std::vector<int>*>* weights = new std::vector<std::vector<int>*>(num_vertices);

    for (int i = 0; i < (int)graph->size(); i++) {
        graph->at(i) = new std::vector<int>;
        weights->at(i) = new std::vector<int>;
    }

    std::unordered_set<int> vertices;

    for (int i = 0; i < num_vertices; i++)
        vertices.insert(i);

    int current, neighbor, elem;
    int count = 0;

    while (count < edges) {

        elem = rand_vertex(rng);
        if (vertices.count(elem) == 0) {
            vertices.insert(elem);
            std::unordered_set<int>::iterator it = vertices.find(elem);
            if (it == vertices.end()) current = (*vertices.begin());
            else current = (*it++);
            vertices.erase(elem);
        }
        else current = elem;

        elem = rand_vertex(rng);
        if (vertices.count(elem) == 0) {
            vertices.insert(elem);
            std::unordered_set<int>::iterator it = vertices.find(elem);
            if (it == vertices.end()) neighbor = (*vertices.begin());
            else neighbor = (*it++);
            vertices.erase(elem);
        }
        else neighbor = elem;

        bool legit = true;
        if (current == neighbor) legit = false;
        if ((int)graph->at(current)->size() == max_degree) legit = false;

        for(int k = 0; k < (int)graph->at(current)->size(); k++)
            if (graph->at(current)->at(k) == neighbor)
                legit = false;

        for(int k = 0; k < (int)graph->at(neighbor)->size(); k++)
            if (graph->at(neighbor)->at(k) == current)
                legit = false;

        if (legit) {
            graph->at(current)->push_back(neighbor);
            weights->at(current)->push_back(rand_weight(rng));
            count++;
        }
    }

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


    std::cout << "\nA graph with " << num_vertices << " vertices, " << edges << " edges and maximum degree of " << max_degree << "\nhas been written to " << filename << std::endl;
    return 0;
}
