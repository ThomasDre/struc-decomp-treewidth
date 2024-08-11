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
#include <boost/graph/copy.hpp>
#include <getopt.h>
#include <string>
#include <chrono>
#include <iomanip>
#include <ctime>


using namespace std;
using namespace boost;

struct VertexProperties {
    string label;
    vector<string> nodes;
};


typedef boost::property<boost::vertex_name_t, std::string, VertexProperties> VertexProperty;
typedef boost::adjacency_list<mapS, listS, undirectedS, VertexProperties> Graph;
typedef boost::erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<Graph>::adjacency_iterator AdjacencyIterator;

map<Vertex, Vertex> lookup;


/**
--------------------------------------------------
-               Utils and Usage                  -
--------------------------------------------------
**/

// helper function
string join(Graph& g, const vector<Vertex>& v, const string& delim) {
    ostringstream s;
    int index = 0;
    for (const Vertex i : v) {
        auto prop_map = get(vertex_bundle, g);
        VertexProperties vp = prop_map[i];
        s << vp.label;
        if (index < (v.size() - 1)) {
            s << delim;
        }
        index++;
    }
    return s.str();
}

// helper function
std::string vector_to_string(const vector<string>& vec) {
    stringstream result;
    for (size_t i = 0; i < vec.size(); ++i) {
        result << vec[i];
        if (i < vec.size() - 1) {
            result << ";";
        }
    }
    return result.str();
}

void print_usage() {
    std::cout << "Usage: program [options]\n"
              << "Options:\n"
              << "  --graph <input graph>                               Input file of graph\n"
              << "  --heuristic=min-deg|min-fill|max-card               Set elimination ordering\n"
              << "  --evaluation                                        Perform evaluation\n"
              << "  --help                                              Display this help message\n";
}

/**
Method used to create a look up table that stores which vertex of original graph corresponds to new vertex in copied/updated graph

This is necessary to keep original vertices of the input graph, vertices of the elimination ordering and vertices of temporary graphs
(e.g. during the eliminate and reconnect phase) in synch.
**/
void create_lookup_table(Graph& orig, Graph& g) {
    graph_traits <Graph>::vertex_iterator i,end;
    vector<Vertex> orig_vertices;
    for (tie(i,end) = vertices(orig); i != end; i++) {
        orig_vertices.push_back(*i);
    }
    vector<Vertex> g_vertices;
    for (tie(i,end) = vertices(g); i != end; i++) {
        g_vertices.push_back(*i);
    }

    for (int i = 0; i <  orig_vertices.size(); i++) {
        Vertex orig_v = orig_vertices[i];
        Vertex g_v = g_vertices[i];
        lookup[g_v] = orig_v;
    }
}

/**
Method that updates lookup table of the vertices of a given graph when the addresses of the nodes change (upon copying, recreating the graph, etc.)
**/
void update_prop_map(Graph& orig_g, Graph& g) {
    auto prop_map = get(vertex_bundle, orig_g);
    graph_traits <Graph>::vertex_iterator i, end;

    vector<VertexProperties> vps;
    for (tie(i,end) = vertices(orig_g); i!= end; i++) {
        VertexProperties vp_i = prop_map[*i];
        vps.push_back(vp_i);
    }

    int index = 0;
    for (tie(i,end) = vertices(g); i != end; i++) {
        VertexProperties vp_i = vps[index];
        prop_map[*i] = vp_i;
        index++;
    }
}







/**
--------------------------------------------------
-     Parsing Input And Storing Results          -
--------------------------------------------------
**/

/**
Stores the given elimination ordering in a txt file
**/
string store_ordering(vector<Vertex>& elimination_ordering, Graph& g, string type) {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << put_time(&tm, "%d-%m-%Y %H-%M-%S");
    string filename = "results/eliminationOrderings/" + type + "_" + oss.str() + ".txt";
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return "";
    }

    auto prop_map = get(vertex_bundle, g);

    for (const auto& vertex : elimination_ordering) {
        VertexProperties vp_vertex = prop_map[vertex];
        outFile << vp_vertex.label << " ";
    }

    outFile.close();

    return filename;
}

/**
Stores the given tree decomposition in a csv file
**/
string store_decomposition(Graph& decomposition, Graph& g) {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << put_time(&tm, "%d-%m-%Y %H-%M-%S");
    string filename = "results/decompositions/tree_" + oss.str() + ".csv";
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return "";
    }

    graph_traits <Graph>::vertex_iterator i,end;

    auto prop_map = get(vertex_bundle, decomposition);
    map<Vertex,string> node_map;


    int index = 1;
    for (tie(i,end) = vertices(decomposition); i != end; i++) {
        string nodeLabel = "N" + to_string(index);
        node_map[*i] = nodeLabel;
        index++;
    }

    graph_traits <Graph>::edge_iterator edge_i, edge_end;

    for (tie(edge_i, edge_end) = edges(decomposition); edge_i != edge_end; edge_i++) {
        Vertex source_v = source(*edge_i, decomposition);
        Vertex target_v = target(*edge_i, decomposition);
        string source_label = node_map[source_v];
        string target_label = node_map[target_v];
        outFile << source_label << "," << target_label << "," << endl;
    }

    index = 1;
    for (tie(i,end) = vertices(decomposition); i != end; i++) {
        VertexProperties vp_i = prop_map[*i];
        vector<string> labels = vp_i.nodes;

        string labels_string = vector_to_string(labels);
        string nodeLabel = "N" + index;

        outFile << ("N" + to_string(index)) << ",," << labels_string << endl;
        index++;
    }

    return filename;
}


void save_data(const std::string &filename, const std::vector<int> &x) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (size_t i = 0; i < x.size(); ++i) {
            file << x[i] << endl;
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << "\n";
    }
}



// print the given graph
void print_graph(Graph& g) {
        graph_traits <Graph>::vertex_iterator i, end;
        graph_traits <Graph>::adjacency_iterator ai, a_end;

        if (num_vertices(g) == 0) {
            cout << "Empty graph" << endl;
        } else {
            for (tie(i, end) = vertices(g); i != end; ++i) {
                string line = "";
                auto prop_map = get(vertex_bundle, g);
                VertexProperties vp_i = prop_map[*i];
                line.append(vp_i.label);
                line.append(": ");

                for (tie(ai, a_end) = adjacent_vertices(*i, g); ai != a_end; ++ai) {
                    VertexProperties vp_ai = prop_map[*ai];
                    line.append(vp_ai.label);
                    if (boost::next(ai) != a_end) {
                        line.append(",");
                    }
                }
                cout << line << endl;
            }
        }
}

// print content of given csv file
void print_csv(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open the file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::cout << line << std::endl;
    }

    file.close();
}

/**
Return the individual parts (tokens) of a line in a csv fomat
**/
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

// Function to add a edge with a label
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
    VertexProperties vp;
    vp.label = label;
    if (vertex_map.find(v) == vertex_map.end()) {
        Vertex vertex = add_vertex(vp, g);
        vertex_map[v] = vertex;
        // put(vertex_name, g, vertex, label);
    } else {
        Vertex vertex = vertex_map[vp, v];
        auto prop_map = get(vertex_bundle, g);
        prop_map[vertex] = vp;
        // put(vertex_name, g, vertex, label);
    }
}


Graph::vertex_descriptor add_vertex(const VertexProperties& p, Graph& g) {
    Graph::vertex_descriptor v = boost::add_vertex(p, g);
    return v;
}

/**
Method to parse given csv file
**/
void read_input(string file, Graph& g) {
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
                if (!label.empty() && label[label.size() - 1] == '\r') {
                    label.erase(label.size()-1);
                }
            } else {
                label = "";
            }
            add_edge(g, v1, v2, label, vertex_map);
        } else {
            // Read vertex
            v1 = tokens[0];
            if (tokens.size() == 3) {
                label = tokens[2];
                if (!label.empty() && label[label.size() - 1] == '\r') {
                    label.erase(label.size()-1);
                }
            } else {
                label = "";
            }
            add_vertex(g, v1, label, vertex_map);
        }
    }

    infile.close();
}






/**
--------------------------------------------------
-        Tree decomposition specific function    -
--------------------------------------------------
**/



/**
Returns the tree width of a given tree decomposition.
Note: width = |size of biggest bag| - 1
**/
int get_tree_width(Graph& decomposition) {
    graph_traits <Graph>::vertex_iterator i, end;
    auto prop_map = get(vertex_bundle, decomposition);
    int max_bag_size = 0;

    for (tie(i,end) = vertices(decomposition); i != end; i++) {
        VertexProperties vp_i = prop_map[*i];
        vector<string> labels = vp_i.nodes;
        if (labels.size() > max_bag_size) {
            max_bag_size = labels.size();
        }
    }

    return max_bag_size - 1;
}


/**
Helper method used during tree decomposition conversion, to find node of given decomposition tree whose bag contains all neighbours of given vertex
**/
Vertex find_decomposition_node(Graph& g, Graph& decomposition, Vertex vertex) {
    graph_traits <Graph>::vertex_iterator v, end;
    auto prop_map = get(vertex_bundle, g);
    for (tie(v,end) = vertices(decomposition); v != end; v++) {
        VertexProperties vp = decomposition[*v];
        vector<string> labels = vp.nodes;
        pair<AdjacencyIterator, AdjacencyIterator> adjacent_to_vertex = adjacent_vertices(vertex, g);
        bool all_contained = true;
        for (AdjacencyIterator ai = adjacent_to_vertex.first; ai != adjacent_to_vertex.second; ai++) {
            VertexProperties vp_ai = prop_map[*ai];
            if (find(labels.begin(), labels.end(), vp_ai.label) == labels.end()) {
                all_contained = false;
            }
        }
        if (all_contained) {
            return *v;
        }
    }
    return vertex;
}

/**
The given vertex v is removed from the graph g and the graph g is reconnected by adding fill-in-edges
**/
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

/**
Returns the vertex with the minimum degree
**/
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

/**
Returns the vertex which produces the minimum number of fill-in-edges upon its removal
**/
Vertex find_min_fill_in_vertex(Graph& g) {
    graph_traits <Graph>::vertex_iterator v, end;
    Vertex min_fill_in_vertex;
    int min_fill_in_edges = INT_MAX;

    for (tie(v,end) = vertices(g); v != end; v++) {
        // for all vertices in G: |V|
        int fill_in_edges = 0;
        // assume we remove v, how many fill in edges do we need to add
        //
        pair<AdjacencyIterator, AdjacencyIterator> neighbours_of_v = adjacent_vertices(*v, g);
        for (AdjacencyIterator adjacent_v = neighbours_of_v.first; adjacent_v != neighbours_of_v.second; adjacent_v++) {
            // for all neighbours of current v: |V|
            pair<AdjacencyIterator, AdjacencyIterator> neighbours_of_current_iters = adjacent_vertices(*adjacent_v, g);
            vector<Vertex> neighbours_of_current;
            copy(neighbours_of_current_iters.first, neighbours_of_current_iters.second, std::back_inserter(neighbours_of_current));
            pair<AdjacencyIterator, AdjacencyIterator> neighbours_of_v_tmp = adjacent_vertices(*v, g);
            for (AdjacencyIterator adjacent_to_orig_v = neighbours_of_v_tmp.first; adjacent_to_orig_v != neighbours_of_v_tmp.second; adjacent_to_orig_v++) {
                // for all neighbours of current v: |V|
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

/**
Returns the vertex which has the maximum number of its neighbors already contained in the given elimination ordering list.
**/
Vertex find_max_neighbours_in_ordering_vertex(Graph& g, vector<Vertex>& ordering) {
    graph_traits <Graph>::vertex_iterator v, end;
    Vertex max_neighbours_in_ordering_vertex;
    int max_num_of_neighbours_in_ordering = -1;

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

        if (neighbours_in_ordering_vertex > max_num_of_neighbours_in_ordering) {
            max_num_of_neighbours_in_ordering = neighbours_in_ordering_vertex;
            max_neighbours_in_ordering_vertex = *v;
        }
    }

    return max_neighbours_in_ordering_vertex;
}








/**
--------------------------------------------------
-                 Heuristics                     -
--------------------------------------------------
**/


/**
Creates an elimination ordering based on the minimum-degree heuristic
**/
vector<Vertex> min_degree_elimination_heuristic(Graph g, Graph& orig) {
    update_prop_map(orig, g);
    create_lookup_table(orig,g);
    vector<Vertex> ordering;

    int vertices_size = num_vertices(g);

    for (int i = 0; i < vertices_size; i++) {
        Vertex min_deg_vertex = find_min_deg_vertex(g);
        ordering.push_back(lookup[min_deg_vertex]);
        remove_vertex_and_reconnect_graph(g, min_deg_vertex);
    }

    return ordering;
}

/**
Creates an elimination ordering based on the minimum-fill-in-edges heuristic
**/
vector<Vertex> min_fill_in_eliminiation_heuristic(Graph g, Graph& orig) {
    update_prop_map(orig, g);
    create_lookup_table(orig,g);
    vector<Vertex> ordering;

    int vertices_size = num_vertices(g);

    for (int i = 0; i < vertices_size; i++) {
        Vertex min_fill_in_vertex = find_min_fill_in_vertex(g);
        ordering.push_back(lookup[min_fill_in_vertex]);
        remove_vertex_and_reconnect_graph(g, min_fill_in_vertex);
    }

    return ordering;
}

/**
Creates an elimination ordering based on the maximum cardinality heuristic
**/
vector<Vertex> max_card_heuristic(Graph g, Graph& orig) {
    update_prop_map(orig, g);
    create_lookup_table(orig,g);
    vector<Vertex> ordering;

    int vertices_size = num_vertices(g);

    // select some vertex as initial first element
    Vertex init_vertex = find_min_deg_vertex(g);
    ordering.push_back(lookup[init_vertex]);

    for (int i = 1; i < vertices_size; i++) {
        // select vertex which has highest num of neighbours in ordering
        Vertex vertex = find_max_neighbours_in_ordering_vertex(orig, ordering);
        ordering.push_back(vertex);
    }

    return ordering;
}


/**
--------------------------------------------------
-             Tree Decomposition                 -
--------------------------------------------------
**/

Graph create_tree_decomposition(Graph& g, Graph& help, Graph& decomposition, vector<Vertex> ordering) {
    auto prop_map = get(vertex_bundle, g);
    if (num_vertices(g) == 1) {
        // return basic node new tree node,label with current node
        auto vertex_pair = boost::vertices(g);
        Graph::vertex_descriptor v = *vertex_pair.first;
        VertexProperties vp;
        VertexProperties vp_v = prop_map[v];
        vp.nodes.push_back(vp_v.label);
        add_vertex(vp, decomposition);
        return decomposition;
    }
    Vertex vertex = ordering.front();
    ordering.erase(ordering.begin());


    Graph g_copy = Graph(g);
    create_lookup_table(g_copy,g);

    remove_vertex_and_reconnect_graph(g, vertex);
    decomposition = create_tree_decomposition(g, help, decomposition, ordering);
    Vertex decomposition_node = find_decomposition_node(g_copy, decomposition, lookup[vertex]);
    VertexProperties vp_decomposition_node = prop_map[decomposition_node];
    // new decomp_node, make labels with vertex and all its neighbours
    // add edge between decomposition_node and new decomp_node
    VertexProperties vp;
    //vp.nodes.push_back(vp_decomposition_node.label);
    VertexProperties vp_vertex = prop_map[lookup[vertex]];
    vp.nodes.push_back(vp_vertex.label);
    pair<AdjacencyIterator, AdjacencyIterator> neighbours_of_current = adjacent_vertices(lookup[vertex], g);
    for (AdjacencyIterator ai = neighbours_of_current.first; ai != neighbours_of_current.second; ai++) {
        VertexProperties vp_ai = prop_map[*ai];
        vp.nodes.push_back(vp_ai.label);
    }
    Vertex new_node = add_vertex(vp, decomposition);
    add_edge(decomposition_node, new_node, decomposition);
    return decomposition;
}

/**
Generates a tree decomposition from given ordering for the given input graph
**/
Graph create_tree_decomposition(Graph& g, Graph help, vector<Vertex> ordering) {
    Graph tree_decomposition;
    create_tree_decomposition(g, help, tree_decomposition, ordering);
    return tree_decomposition;
}

/**
Runs the experiments for given settings of number of nodes (10,100,1000) and connectivity probability (0.25,0.5,0.75) by generating instances of Erdos-Renyi based tree instances
**/
void run_experiment() {
    minstd_rand gen;
    vector<int> n_values = {10,100,1000};
    vector<double> p_values = {0.25,0.5,0.75};

    for (const int& n: n_values) {
        for (const double& p: p_values) {
            // run tests on given n,p setting
            vector<int> min_deg_dps;
            vector<int> min_fill_in_dps;
            vector<int> max_card_dps;

            cout << "Run experiment for n=" << n << ", p=" << p << endl;

            for (int i = 0; i < 100; i++) {
                Graph g(ERGen(gen, n, p), ERGen(), n);
                Graph g1 = Graph(g);
                Graph g2 = Graph(g);

                //std::chrono::steady_clock::time_point begin_t1 = std::chrono::steady_clock::now();
                vector<Vertex> min_deg_ordering = min_degree_elimination_heuristic(g, g);
                //std::chrono::steady_clock::time_point end_t1 = std::chrono::steady_clock::now();
                //std::cout << "Time difference min_deg_ordering= " << (std::chrono::duration_cast<std::chrono::microseconds>(end_t1 - begin_t1).count())/1000000.0 << "[s]" << std::endl;

                //std::chrono::steady_clock::time_point begin_t2 = std::chrono::steady_clock::now();
                vector<Vertex> min_fill_in_ordering = min_fill_in_eliminiation_heuristic(g1, g1);
                //std::chrono::steady_clock::time_point end_t2 = std::chrono::steady_clock::now();
                //std::cout << "Time difference min_fill_in_ordering= " << (std::chrono::duration_cast<std::chrono::microseconds>(end_t2 - begin_t2).count())/1000000.0  << "[s]" << std::endl;

                //std::chrono::steady_clock::time_point begin_t3 = std::chrono::steady_clock::now();
                vector<Vertex> max_card_ordering = max_card_heuristic(g2, g2);
                //std::chrono::steady_clock::time_point end_t3 = std::chrono::steady_clock::now();
                //std::cout << "Time difference max_card_ordering= " << (std::chrono::duration_cast<std::chrono::microseconds>(end_t3 - begin_t3).count())/1000000.0  << "[s]" << std::endl;

                //std::chrono::steady_clock::time_point begin_dec1 = std::chrono::steady_clock::now();
                Graph min_deg_tree_decomposition = create_tree_decomposition(g, g, min_deg_ordering);
                //std::chrono::steady_clock::time_point end_dec1 = std::chrono::steady_clock::now();
                //std::cout << "Time difference min_deg_composition= " << (std::chrono::duration_cast<std::chrono::microseconds>(end_dec1 - begin_dec1).count())/1000000.0  << "[s]" << std::endl;

                //std::chrono::steady_clock::time_point begin_dec2 = std::chrono::steady_clock::now();
                Graph min_fill_in_tree_decomposition = create_tree_decomposition(g1, g1, min_fill_in_ordering);
                //std::chrono::steady_clock::time_point end_dec2 = std::chrono::steady_clock::now();
                //std::cout << "Time difference min_fill_in_composition= " << (std::chrono::duration_cast<std::chrono::microseconds>(end_dec2 - begin_dec2).count())/1000000.0  << "[s]" << std::endl;

                //std::chrono::steady_clock::time_point begin_dec3 = std::chrono::steady_clock::now();
                Graph max_card_tree_decomposition = create_tree_decomposition(g2, g2, max_card_ordering);
                //std::chrono::steady_clock::time_point end_dec3 = std::chrono::steady_clock::now();
                //std::cout << "Time difference max_card_composition= " << (std::chrono::duration_cast<std::chrono::microseconds>(end_dec3 - begin_dec3).count())/1000000.0  << "[ss]" << std::endl;

                int tw_min_deg = get_tree_width(min_deg_tree_decomposition);
                int tw_min_fill_in = get_tree_width(min_fill_in_tree_decomposition);
                int tw_max_card = get_tree_width(max_card_tree_decomposition);

                min_deg_dps.push_back(tw_min_deg);
                min_fill_in_dps.push_back(tw_min_fill_in);
                max_card_dps.push_back(tw_max_card);
            }

            string min_deg_results = "min_deg_" + to_string(n) + "_" + to_string(p) + "_results.dat";
            string min_fill_in_results = "min_fill_in_" + to_string(n) + "_" + to_string(p) + "_results.dat";
            string max_card_results = "max_card_" + to_string(n) + "_" + to_string(p) + "_results.dat";

            save_data("results/eval/" + min_deg_results, min_deg_dps);
            save_data("results/eval/" + min_fill_in_results, min_fill_in_dps);
            save_data("results/eval/" + max_card_results, max_card_dps);
        }
    }
}


int main(int argc, char* argv[]) {

    int opt;
    int option_index = 0;
    string graph_input;
    string elimination_ordering_type = "min-deg"; // default value
    bool evaluation = false;

    struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"graph", required_argument, 0, 'g'},
        {"heuristic", required_argument, 0, 'e'},
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
    cout << "Evaluation mode: " << (evaluation ? "yes" : "no") << endl;


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
            elimination_ordering = min_degree_elimination_heuristic(g, g);
        } else if (elimination_ordering_type == "min-fill") {
            elimination_ordering = min_fill_in_eliminiation_heuristic(g, g);
        } else if (elimination_ordering_type == "max-card") {
            elimination_ordering = max_card_heuristic(g, g);
        } else {
            cerr << "Error: Unknown elimination ordering heuristic";
            return 1;
        }

        // store ordering
        string ordering_file_name = store_ordering(elimination_ordering, g, elimination_ordering_type);

        cout << endl << "Elimination ordering:" << endl;
        cout << join(g, elimination_ordering, ",") << endl;
        cout << "**********************" << endl;


        Graph decomposition;
        decomposition = create_tree_decomposition(g, g, elimination_ordering);
        string tree_decomposition_file_name = store_decomposition(decomposition, g);

        cout << "Tree Width:" << endl;
        cout << get_tree_width(decomposition) << endl;
        cout <<  "**********************" << endl;

        cout << "Tree decomposition:" << endl;
        print_csv(tree_decomposition_file_name);
        cout << "**********************" << endl;

        if (!ordering_file_name.empty()) {
            cout << "Ordering saved: " << ordering_file_name << endl;
        }
        if (!tree_decomposition_file_name.empty()) {
            cout << "Tree decomposition saved: " << tree_decomposition_file_name << endl;
        }
    }

    return 0;
}
