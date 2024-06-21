#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/random/linear_congruential.hpp>
#include <getopt.h>

using namespace std;
using namespace boost;

struct VertexProperties {
    vector<string> nodes;
};

typedef boost::property<boost::vertex_name_t, std::string, VertexProperties> VertexProperty;
typedef boost::adjacency_list<mapS, listS, undirectedS, property<vertex_name_t,string>> Graph;
typedef boost::erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<Graph>::adjacency_iterator AdjacencyIterator;


map<Vertex, string> lookup;


string join(Graph g, const vector<Vertex>& v, const string& delim) {
    ostringstream s;
    for (const auto& i : v) {
        if (&i != &v[0]) {
            s << delim;
        }
        s << get(vertex_name, g, i);
    }
    return s.str();
}

void print_usage() {
    std::cout << "Usage: program [options]\n"
              << "Options:\n"
              << "  --graph <input graph>                               Input file of graph\n"
              << "  --heuristic=min-deg|min-fill|max-card               Set elimination ordering\n"
              << "  --output <targetFolder>                             Folder where results are stored\n"
              << "  --evaluation                                        Perform evaluation\n"
              << "  --help                                              Display this help message\n";
}

void print_graph(Graph& g) {
      // Output the graph
        //auto vertex_idMap = get(vertex_index, g);
        graph_traits <Graph>::vertex_iterator i, end;
        graph_traits <Graph>::adjacency_iterator ai, a_end;

        if (num_vertices(g) == 0) {
            cout << "Empty graph" << endl;
        } else {
            for (tie(i, end) = vertices(g); i != end; ++i) {
                cout << get(vertex_name, g, *i) << ": ";

                for (tie(ai, a_end) = adjacent_vertices(*i, g); ai != a_end; ++ai) {
                    cout << get(vertex_name,g, *ai);
                    //cout << lookup[*ai];
                    if (boost::next(ai) != a_end)
                        cout << ", ";
                }
                cout << endl;
            }
        }
}

void update_lookup_table(Graph& g) {
    graph_traits <Graph>::vertex_iterator i, end;

    for (tie(i, end) = vertices(g); i != end; i++) {

    }
    vertices(g);
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
    tie(e, inserted) = add_edge(vertex2, vertex1, g);
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
        put(vertex_name, g, vertex, label);
    } else {
        Vertex vertex = vertex_map[v];
        put(vertex_name, g, vertex, label);
    }
}


/*Graph::vertex_descriptor add_vertex(const VertexProperties& p, Graph& g) {
    Graph::vertex_descriptor v = boost::add_vertex(g);
    g[v] = p;
    return v;
}*/

void read_input(string file, Graph& g) {
    // Open the CSV file
    ifstream infile(file);
    if (!infile.is_open()) {
        cerr << "Error opening file." << endl;
        exit(1);
    }

    string line;
    bool reading_edges = true;
    map<string,Vertex> vertex_map;

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
            if (tokens.size() == 3) {
                label = tokens[2];
            } else {
                label = "";
            }
            add_vertex(g, v1, label, vertex_map);
        }
    }

    map<string,Vertex>::iterator it;
    for (it = vertex_map.begin(); it != vertex_map.end(); it++){
            //How do I access each element?
            lookup[it->second] = it->first;
    }


    infile.close();
}


void remove_vertex_and_reconnect_graph(Graph& g, Vertex v) {
    pair<AdjacencyIterator, AdjacencyIterator> neighbours_of_v = adjacent_vertices(v, g);

    for (AdjacencyIterator i = neighbours_of_v.first; i != neighbours_of_v.second; i++) {
        pair<AdjacencyIterator, AdjacencyIterator> neighbours_of_v_tmp = adjacent_vertices(v, g);
        for (AdjacencyIterator ai = neighbours_of_v_tmp.first; ai != neighbours_of_v_tmp.second; ai++) {
            if (*i != *ai) {
                add_edge(*i, *ai, g);
            }
        }
    }

    clear_vertex(v, g);
    remove_vertex(v, g);
}


Vertex find_min_deg_vertex(Graph& g) {
    graph_traits <Graph>::vertex_iterator i, end;
    Vertex min_deg_vertex;
    int min_degree = INT_MAX;
    for (tie(i,end) = vertices(g); i != end; i++) {
        int out_deg = out_degree(*i, g);
        if (out_deg < min_degree) {
            min_degree = out_deg;
            min_deg_vertex = *i;
        }
    }
    return min_deg_vertex;
}

Vertex find_min_fill_in_vertex(Graph& g) {
    graph_traits <Graph>::vertex_iterator v, end;
    Vertex min_fill_in_vertex;
    int min_fill_in_edges = INT_MAX;

    for (tie(v,end) = vertices(g); v != end; v++) {
        int fill_in_edges = 0;
        // assume we remove v, how many fill in edges do we need to add
        //
        pair<AdjacencyIterator, AdjacencyIterator> neighbours_of_v = adjacent_vertices(*v, g);
        for (AdjacencyIterator adjacent_v = neighbours_of_v.first; adjacent_v != neighbours_of_v.second; adjacent_v++) {
            pair<AdjacencyIterator, AdjacencyIterator> neighbours_of_current_iters = adjacent_vertices(*adjacent_v, g);
            vector<Vertex> neighbours_of_current;
            copy(neighbours_of_current_iters.first, neighbours_of_current_iters.second, std::back_inserter(neighbours_of_current));
            pair<AdjacencyIterator, AdjacencyIterator> neighbours_of_v_tmp = adjacent_vertices(*v, g);
            for (AdjacencyIterator adjacent_to_orig_v = neighbours_of_v_tmp.first; adjacent_to_orig_v != neighbours_of_v_tmp.second; adjacent_to_orig_v++) {
                if (adjacent_to_orig_v != adjacent_v) {
                    if (find(neighbours_of_current.begin(), neighbours_of_current.end(), *adjacent_to_orig_v) == neighbours_of_current.end()) {
                        fill_in_edges += 1;
                    }
                }
            }
        }

        fill_in_edges = fill_in_edges / 2;

        if (fill_in_edges < min_fill_in_edges) {
            min_fill_in_edges = fill_in_edges;
            min_fill_in_vertex = *v;
        }
    }

    return min_fill_in_vertex;
}

Vertex find_max_neighbours_in_ordering_vertex(Graph& g, vector<Vertex>& ordering) {
    graph_traits <Graph>::vertex_iterator v, end;
    Vertex max_neighbours_in_ordering_vertex;
    int max_num_of_neighbours_in_ordering = 0;

    for (tie(v,end) = vertices(g); v != end; v++) {
        if (find(ordering.begin(), ordering.end(), *v) != ordering.end()) {
            // vertex v is already in ordering
            continue;
        }

        int neighbours_in_ordering_vertex = 0;
        pair<AdjacencyIterator, AdjacencyIterator> neighbours_of_current_iters = adjacent_vertices(*v, g);
        for (AdjacencyIterator adjacent_to_current = neighbours_of_current_iters.first; adjacent_to_current != neighbours_of_current_iters.second; adjacent_to_current++) {
            if (find(ordering.begin(), ordering.end(), *adjacent_to_current) != ordering.end()) {
                neighbours_in_ordering_vertex++;
            }
        }
        // get num of neighbours of v in ordering

        if (neighbours_in_ordering_vertex > max_num_of_neighbours_in_ordering) {
            max_num_of_neighbours_in_ordering = neighbours_in_ordering_vertex;
            max_neighbours_in_ordering_vertex = *v;
        }
    }

    return max_neighbours_in_ordering_vertex;
}

vector<Vertex> min_degree_elimination_heuristic(Graph g) {
    vector<Vertex> ordering;

    int vertices_size = num_vertices(g);

    for (int i = 0; i < vertices_size; i++) {
        Vertex min_deg_vertex = find_min_deg_vertex(g);
        ordering.push_back(min_deg_vertex);
        remove_vertex_and_reconnect_graph(g, min_deg_vertex);
        cout << endl << endl;
        print_graph(g);
        cout << endl << endl;
    }

    return ordering;
}

vector<Vertex> min_fill_in_eliminiation_heuristic(Graph g) {
    vector<Vertex> ordering;

    int vertices_size = num_vertices(g);

    for (int i = 0; i < vertices_size; i++) {
        Vertex min_fill_in_vertex = find_min_fill_in_vertex(g);
        ordering.push_back(min_fill_in_vertex);
        remove_vertex_and_reconnect_graph(g, min_fill_in_vertex);
        cout << endl << endl;
        print_graph(g);
        cout << endl << endl;
    }

    return ordering;
}


vector<Vertex> max_card_heuristic(Graph g) {
    vector<Vertex> ordering;

    int vertices_size = num_vertices(g);

    // select some vertex as initial first element
    Vertex init_vertex = find_min_deg_vertex(g);
    ordering.push_back(init_vertex);

    for (int i = 1; i < vertices_size; i++) {
        // select vertex with has highest num of neighbours in ordering
        Vertex vertex = find_max_neighbours_in_ordering_vertex(g, ordering);
        ordering.push_back(vertex);
    }

    return ordering;
}

/*Vertex find_decomposition_node(Graph& g, Graph& decomposition, Vertex vertex) {
    graph_traits <Graph>::vertex_iterator v, end;
    for (tie(v,end) = vertices(decomposition); v != end; v++) {
        get(vertex_name, decomposition)[source(*v,decomposition)];
        VertexProperties vp = decomposition[*v];
        vector<string> labels = vp.nodes;
        pair<AdjacencyIterator, AdjacencyIterator> adjacent_to_vertex = adjacent_vertices(vertex, g);
        bool all_contained = true;
        for (AdjacencyIterator ai = adjacent_to_vertex.first; ai != adjacent_to_vertex.second; ai++) {
            if (find(labels.begin(), labels.end(), lookup[*ai]) != labels.end()) {
                all_contained = false;
            }
        }
        if (all_contained) {
            return *v;
        }
    }
    return vertex;
}


Graph create_tree_decomposition(Graph& g, Graph& decomposition, vector<Vertex> ordering) {
    if (num_vertices(g) == 1) {
        // return basic node new tree node,label with current node
        auto vertex_pair = boost::vertices(g);
        Graph::vertex_descriptor v = *vertex_pair.first;
        VertexProperties vp;
        vp.nodes.push_back(lookup[v]);
        add_vertex(vp, decomposition);
    }
    Vertex vertex = ordering.front();
    ordering.erase(ordering.begin());
    remove_vertex_and_reconnect_graph(g, vertex);
    decomposition = create_tree_decomposition(g, decomposition, ordering);
    Vertex decomposition_node = find_decomposition_node(g, decomposition, vertex);
    // new decomp_node, make labels with vertex and all its neighbours
    // add edge between decomposition_node and new decomp_node
    VertexProperties vp;
    vp.nodes.push_back(lookup[vertex]);
    pair<AdjacencyIterator, AdjacencyIterator> neighbours_of_current = adjacent_vertices(vertex, g);
    for (AdjacencyIterator ai = neighbours_of_current.first; ai != neighbours_of_current.second; ai++) {
        vp.nodes.push_back(lookup[*ai]);
    }
    Vertex new_node = add_vertex(vp, decomposition);
    add_edge(decomposition_node, new_node, decomposition);
    return decomposition;
}

Graph create_tree_decomposition(Graph g, vector<Vertex> ordering) {
    Graph tree_decomposition;
    create_tree_decomposition(g, tree_decomposition, ordering);
    return tree_decomposition;
}
*/

void run_experiment() {
    minstd_rand gen;
    // Create graph with 100 nodes and edges with probability 0.05
    vector<int> n_values = {10,100,1000};
    vector<double> p_values = {1/4,1/2,3/4};

    for (const int& n: n_values) {
        for (const double& p: p_values) {
            // run tests on given n,p setting
            for (int i = 0; i < 100; i++) {
                Graph g(ERGen(gen, n, p), ERGen(), n);

                /*vector<Vertex> min_deg_ordering = min_degree_elimination_heuristic(g);
                vector<Vertex> min_fill_in_ordering = min_fill_in_eliminiation_heuristic(g);
                vector<Vertex> max_card_ordering = max_card_heuristic(g);

                Graph min_deg_tree_decomposition = create_tree_decomposition(g, min_deg_ordering);
                Graph min_fill_in_tree_decomposition = create_tree_decomposition(g, min_fill_in_ordering);
                Graph max_card_tree_decomposition = create_tree_decomposition(g, max_card_ordering);
                */

            //
            }
        }
    }

}

int main(int argc, char* argv[]) {

    int opt;
    int option_index = 0;
    string graph_input;
    string elimination_ordering_type = "min-deg"; // default value
    string output;
    bool treeDecomposition = false;
    bool evaluation = false;

    struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"graph", required_argument, 0, 'g'},
        {"heuristic", required_argument, 0, 'e'},
        {"treeDecomposition", no_argument, 0, 't'},
        {"output", required_argument, 0, 'o'},
        {"evaluation", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "hg:e:vto:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'h':
                print_usage();
                return 0;
            case 'g':
                graph_input = optarg;
                break;
            case 'e':
                elimination_ordering_type = optarg;
                if (elimination_ordering_type != "min-deg" && elimination_ordering_type != "min-fill" && elimination_ordering_type != "max-card") {
                    cerr << "Invalid value for --eliminationOrdering: " << elimination_ordering_type << "\n";
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
    cout << "Elimination ordering: " << elimination_ordering_type << endl;
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
        read_input(graph_input, g);
        print_graph(g);

        vector<Vertex> elimination_ordering;

        // min-deg|min-fill|max-card
        if(elimination_ordering_type == "min-deg") {
            elimination_ordering = min_degree_elimination_heuristic(g);
        } else if (elimination_ordering_type == "min-fill") {
            elimination_ordering = min_fill_in_eliminiation_heuristic(g);
        } else if (elimination_ordering_type == "max-card") {
            elimination_ordering = max_card_heuristic(g);
        } else {
            cerr << "Error: Unknown elimination ordering heuristic";
            return 1;
        }

        cout << "Elimination ordering:" << endl;
        cout << join(g, elimination_ordering, ",") << endl;


        // TODO create a tree decomposition from ordering#
        Graph decomposition;

        if (treeDecomposition) {
            //decomposition = create_tree_decomposition(g, elimination_ordering);
        }


        if (!output.empty()) {
            // store resulting tree decomposition in
        }

    }

    return 0;
}
