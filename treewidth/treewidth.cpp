#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/random/linear_congruential.hpp>
#include <getopt.h>

using namespace std;
using namespace boost;

typedef boost::adjacency_list<> Graph;
typedef boost::erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;

void print_usage() {
    std::cout << "Usage: program [options]\n"
              << "Options:\n"
              << "  --graph <input graph>                               Input file of graph\n"
              << "  --eliminationOrdering=min-deg|min-fill|max-card     Set elimination ordering\n"
              << "  --output <targetFolder>                             Folder where results are stored\n"
              << "  --evaluation                                        Perform evaluation\n"
              << "  --help                                              Display this help message\n";
}

void print_graph(Graph& g) {
      // Output the graph
        auto vertex_idMap = get(vertex_index, g);
        graph_traits <Graph>::vertex_iterator i, end;
        graph_traits <Graph>::adjacency_iterator ai, a_end;

        for (tie(i, end) = vertices(g); i != end; ++i) {
            cout << vertex_idMap[*i] << ": ";

            for (tie(ai, a_end) = adjacent_vertices(*i, g); ai != a_end; ++ai) {
                cout << vertex_idMap[*ai];
                if (boost::next(ai) != a_end)
                    cout << ", ";
            }
            cout << endl;
        }
}

vector<string> get_tokens_from_line(string s) {
    stringstream ss(s);
    string token;
    vector<string> tokens;
    char delimiter = ',';

    while (getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }

    return tokens;
}

void add_edge(Graph& g, const string& v1, const string& v2, const string& label,
              map<string, Vertex>& vertex_map) {
    Vertex vertex1, vertex2;

    if (vertex_map.find(v1) == vertex_map.end()) {
        vertex1 = boost::add_vertex(g);
        vertex_map[v1] = vertex1;
    } else {
        vertex1 = vertex_map[v1];
    }

    if (vertex_map.find(v2) == vertex_map.end()) {
        vertex2 = add_vertex(g);
        vertex_map[v2] = vertex2;
    } else {
        vertex2 = vertex_map[v2];
    }

    Edge e;
    bool inserted;
    tie(e, inserted) = add_edge(vertex1, vertex2, g);
    if (inserted && !label.empty()) {
        //put(g, e, label);
    }
}

// Function to add a vertex with a label
void add_vertex(Graph& g, const string& v, const string& label,
                map<string, Vertex>& vertex_map) {
    if (vertex_map.find(v) == vertex_map.end()) {
        Vertex vertex = add_vertex(g);
        vertex_map[v] = vertex;
        //put(vertex_name, g, vertex, label);
    } else {
        Vertex vertex = vertex_map[v];
        //put(vertex_name, g, vertex, label);
    }
}


void read_input(string file, Graph& g) {
    // Open the CSV file
    ifstream infile(file);
    if (!infile.is_open()) {
        cerr << "Error opening file." << endl;
        exit(1);
    }

    map<string, Vertex> vertex_map;
    string line;
    bool reading_edges = true;

    // Read the file line by line
    while (getline(infile, line)) {
        vector<string> tokens = get_tokens_from_line(line);

        if (tokens[1].empty()) {
            reading_edges = false;
        }

        string v1, v2, label;

        if (reading_edges) {
            // Read edge
            v1 = tokens[0];
            v2 = tokens[1];
            if (tokens.size() == 3) {
                label = tokens[2];
            } else {
                label = "";
            }
            add_edge(g, v1, v2, label, vertex_map);
        } else {
            // Read vertex
            v1 = tokens[0];
            if (tokens.size() == 2) {
                label = tokens[1];
            } else {
                label = "";
            }
            add_vertex(g, v1, label, vertex_map);
        }
    }

    infile.close();
}

void run_experiment() {
    cout << "Experiments for n and k" << endl;
}

int main(int argc, char* argv[]) {

    int opt;
    int option_index = 0;
    string graph_input;
    string eliminationOrdering = "min-deg"; // default value
    string output;
    bool treeDecomposition = false;
    bool evaluation = false;

    struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"graph", required_argument, 0, 'g'},
        {"eliminationOrdering", required_argument, 0, 'e'},
        {"treeDecomposition", no_argument, 0, 't'},
        {"output", required_argument, 0, 'o'},
        {"evaluation", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "hg:e:v", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'h':
                print_usage();
                return 0;
            case 'g':
                graph_input = optarg;
                break;
            case 'e':
                eliminationOrdering = optarg;
                if (eliminationOrdering != "min-deg" && eliminationOrdering != "min-fill" && eliminationOrdering != "max-card") {
                    cerr << "Invalid value for --eliminationOrdering: " << eliminationOrdering << "\n";
                    return 1;
                }
                break;
            case 't':
                treeDecomposition = true;
                break;
            case 'o':
                output = optarg;
                break;
            case 'v':
                evaluation = true;
                break;
            default:
                print_usage();
                return 1;
        }
    }

    // Verify required options
    if (!evaluation && graph_input.empty()) {
        cerr << "Error: --graph <input graph> is required\n";
        print_usage();
        return 1;
    }

    // Process the options
    cout << "Input graph: " << graph_input << endl;
    cout << "Elimination ordering: " << eliminationOrdering << endl;
    cout << "Tree decomposition generated: " << (treeDecomposition? "yes" : "no") << endl;
    cout << "Evaluation: " << (evaluation ? "enabled" : "disabled") << endl;
    cout << "Results stored: " << (output.empty()? "---" : output) << endl << endl;

    if (evaluation) {
        cout << "Start evaluation:" << endl;
        run_experiment();
        return 0;
    } else {
        cout << "Reading input graph:" << endl;
        Graph g;
        map<string, Vertex> vertex_map;
        read_input(graph_input, g);
        print_graph(g);
    }

    return 0;

    /*

    minstd_rand gen;
    // Create graph with 100 nodes and edges with probability 0.05
    int n = 2;
    float p = 0.5;
    Graph g(ERGen(gen, n, p), ERGen(), n);

    auto vertex_idMap = get(vertex_index, g);
    graph_traits <Graph>::vertex_iterator i, end;
    graph_traits <Graph>::adjacency_iterator ai, a_end;

    for (tie(i, end) = vertices(g); i != end; ++i) {
        cout << vertex_idMap[*i] << ": ";

        for (tie(ai, a_end) = adjacent_vertices(*i, g); ai != a_end; ++ai) {
            cout << vertex_idMap[*ai];
            if (boost::next(ai) != a_end)
                cout << ", ";
        }
        cout << endl;
    }

    */


    return 0;
}
