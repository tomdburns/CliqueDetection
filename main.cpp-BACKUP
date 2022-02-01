// Calculates the correspondence graph between
// two imported graphs and outputs it to a file

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <stdio.h>  /* printf, fopen */
#include <stdlib.h> /* Exit Commands */
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
using namespace std;


class Atom {
    public:
        string ident;
        double x_coord;
        double y_coord;
        double z_coord;
        Atom(){
            ident = "NaN";
            x_coord = 0.0;
            y_coord = 0.0;
            z_coord = 0.0;
        }
        Atom(string aLabel, double aX,
            double aY, double aZ){
            ident = aLabel;
            x_coord = aX;
            y_coord = aY;
            z_coord = aZ;
            }
};


class Node {
    public:
        string node_1;
        string node_2;
        string ident;
        Node(string aNode, string bNode, string nIdent) {
            node_1 = aNode;
            node_2 = bNode;
            ident = nIdent;
        }
};


class Edge {
    public:
        string vert1_a;
        string vert1_b;
        string vert2_a;
        string vert2_b;
        double dist_a;
        double dist_b;
        Edge (string v1a, string v1b, string v2a, string v2b, double adist, double bdist) {
            vert1_a = v1a;
            vert1_b = v1b;
            vert2_a = v2a;
            vert2_b = v2b;
            dist_a = adist;
            dist_b = bdist;
        }
};


class Conn {
    // This class contains connectivity info for each node
    public:
        int number;
        vector<string> connects;
        Conn() {
            number = 0;
        }
};


string import_file(string filename)
{
    // Imports raw data from the input file
    ifstream InputFile;
    char data;
    string raw;
    InputFile.open(filename);
    if (InputFile.is_open()) {
        while (!InputFile.eof()) {
            InputFile >> data;
            raw += data;
            //cout << data << " ";
        }
    }
    else {
        exit (EXIT_FAILURE);
    }

    InputFile.close();
    //cout << data << endl;
    //char raw_data[raw.size() + 1];
    //strcpy(raw_data, raw.c_str());
    return raw;
}


vector<string> split(string str, string sep)
{   // Finds all of the locations in a string delimited by
    // the variable delimit
    char* cstr=const_cast<char*>(str.c_str());
    char* current;
    std::vector<std::string> arr;
    current=strtok(cstr,sep.c_str());
    while(current!=NULL){
        arr.push_back(current);
        current=strtok(NULL,sep.c_str());
    }
    return arr;

    //bool run = true;
    //vector<string> splits;
    // the first set is counting the number of instances
    // the delimiter is found
    //int found = 0;
    //int last = 0;
    //string current = data;
    //while (run) {
    //    current.find(delimit);
    //    run = false;
    //}

    //return 0;
}


vector<string> format_atoms(string data)
{   // Formats the raw data and returns it

    // Gather all of the Raw data and parse it into
    // individual graph sections
    string del1 = "AtomList";
    string del2 = "Edges";
    string del3 = "Ref";
    string del4 = "Vertices";
    string del5 = "Original";

    int end1 = data.find(del1) + del1.length();
    int start2 = data.find(del2); // - end1;
    int end2 = start2 + del2.length();
    int start3 = data.find(del3); // - end2;
    int end3 = start3 + del3.length();
    int start4 = data.find(del4); // - end3;
    int end4 = start4 + del4.length();
    int start5 = data.find(del5); // - end4;
    int end5 = start5 + del5.length();
    int datasize = data.length();

    string raw_atoms = data.substr(end1, start2 - end1);
    string raw_edges = data.substr(end2, start3 - end2);
    string raw_refs  = data.substr(end3, start4 - end3);
    string raw_verts = data.substr(end4, start5 - end4);
    string raw_orig  = data.substr(end5, datasize + 1);

    // Formats Atom List
    vector<string> new_atoms = split(raw_atoms, ",");
    string atoms[new_atoms.size()];
    vector<string> atom_list;
    for (size_t i = 0; i < new_atoms.size(); i++) {
        vector<string> small_atoms = split(new_atoms[i].c_str(), "\'");
        if (i == 0) {
            atoms[i] = small_atoms[1].c_str();
        }
        else {
            atoms[i] = small_atoms[0].c_str();
        }
        atom_list.push_back(atoms[i]);
    }

    return atom_list;
}


map<string, Atom> format_origin(string data)
{   // Formats the raw data and returns it

    // Gather all of the Raw data and parse it into
    // individual graph sections
    string del1 = "AtomList";
    string del2 = "Edges";
    string del3 = "Ref";
    string del4 = "Vertices";
    string del5 = "Original";

    int end1 = data.find(del1) + del1.length();
    int start2 = data.find(del2); // - end1;
    int end2 = start2 + del2.length();
    int start3 = data.find(del3); // - end2;
    int end3 = start3 + del3.length();
    int start4 = data.find(del4); // - end3;
    int end4 = start4 + del4.length();
    int start5 = data.find(del5); // - end4;
    int end5 = start5 + del5.length();
    int datasize = data.length();

    string raw_atoms = data.substr(end1, start2 - end1);
    string raw_edges = data.substr(end2, start3 - end2);
    string raw_refs  = data.substr(end3, start4 - end3);
    string raw_verts = data.substr(end4, start5 - end4);
    string raw_orig  = data.substr(end5, datasize + 1);

    // Format the Original Data
    vector<string> new_orig = split(raw_orig, "}");
    map<string, Atom> original;
    for (size_t i = 0; i < new_orig.size(); i++) {
        vector<string> new_1 = split(new_orig[i].c_str(), ":");
        vector<string> raw_coord = split(new_1[3].c_str(), ",");
        double x, y, z;
        for (size_t j = 0; j < raw_coord.size(); j++) {
            size_t offset = 0;
            //cout << j << " >> " << raw_coord[j] << endl;
            if (j == 0) {
                vector<string> string_coord = split(raw_coord[j].c_str(), "[");
                string coord = string_coord[0];
                x = stod(coord, &offset);
            }
            else if (j == 2) {
                vector<string> string_coord = split(raw_coord[j].c_str(), "]");
                string coord = string_coord[0];
                z = stod(coord, &offset);
            }
            else {
                y = stod(raw_coord[j], 0);
            }
        }
        vector<string> new_2 = split(new_orig[i].c_str(), "\'");
        string atom_label = new_2[1].c_str();
        string atom_type = new_2[5].c_str();
        Atom atomic_prop(atom_type, x, y, z);
        original.insert(make_pair(atom_label, atomic_prop));
    }

    return original;
}


map<string, string> format_refs(string data)
{   // Formats the raw data and returns it

    // Gather all of the Raw data and parse it into
    // individual graph sections
    string del1 = "AtomList";
    string del2 = "Edges";
    string del3 = "Ref";
    string del4 = "Vertices";
    string del5 = "Original";

    int end1 = data.find(del1) + del1.length();
    int start2 = data.find(del2); // - end1;
    int end2 = start2 + del2.length();
    int start3 = data.find(del3); // - end2;
    int end3 = start3 + del3.length();
    int start4 = data.find(del4); // - end3;
    int end4 = start4 + del4.length();
    int start5 = data.find(del5); // - end4;
    int end5 = start5 + del5.length();
    int datasize = data.length();

    string raw_atoms = data.substr(end1, start2 - end1);
    string raw_edges = data.substr(end2, start3 - end2);
    string raw_refs  = data.substr(end3, start4 - end3);
    string raw_verts = data.substr(end4, start5 - end4);
    string raw_orig  = data.substr(end5, datasize + 1);

    // Formats the Refs Dict
    vector<string> new_refs = split(raw_refs, ",");
    map<string, string> refs;
    for (size_t i = 0; i < new_refs.size(); i++) {
        vector<string> small = split(new_refs[i].c_str(), "\'");
        if (i == 0) {
            refs.insert(make_pair(small[1], small[3]));
        }
        else {
            refs.insert(make_pair(small[0], small[2]));
        }
    }

    return refs;
}


map<string, double> format_edges(string data)
{   // Formats the raw data and returns it

    // Gather all of the Raw data and parse it into
    // individual graph sections
    string del1 = "AtomList";
    string del2 = "Edges";
    string del3 = "Ref";
    string del4 = "Vertices";
    string del5 = "Original";

    int end1 = data.find(del1) + del1.length();
    int start2 = data.find(del2); // - end1;
    int end2 = start2 + del2.length();
    int start3 = data.find(del3); // - end2;
    int end3 = start3 + del3.length();
    int start4 = data.find(del4); // - end3;
    int end4 = start4 + del4.length();
    int start5 = data.find(del5); // - end4;
    int end5 = start5 + del5.length();
    int datasize = data.length();

    string raw_atoms = data.substr(end1, start2 - end1);
    string raw_edges = data.substr(end2, start3 - end2);
    string raw_refs  = data.substr(end3, start4 - end3);
    string raw_verts = data.substr(end4, start5 - end4);
    string raw_orig  = data.substr(end5, datasize + 1);

    // Formats the Edges Dict - This one is trickier
    // because it needs to convert the string to a double
    vector<string> new_edges = split(raw_edges,  ",");
    map<string, double> edges;
    for (size_t i = 0; i < new_edges.size(); i++) {
        vector<string> small = split(new_edges[i].c_str(), "\'");
        if (i == 0) {
            size_t offset = 0;
            vector<string> raw_valvect = split(small[2].c_str(), ":");
            string raw_val = raw_valvect[0].c_str();
            double val = stod(raw_val, &offset);
            edges.insert(make_pair(small[1], val));
        }
        else {
            size_t offset = 0;
            vector<string> raw_valvect = split(small[1].c_str(), ":");
            string raw_val = raw_valvect[0].c_str();
            double val = stod(raw_val, &offset);
            edges.insert(make_pair(small[0], val));
        }
    }

    return edges;
}


map<string, string> format_vertices(string data)
{   // Formats the raw data and returns it

    // Gather all of the Raw data and parse it into
    // individual graph sections
    string del1 = "AtomList";
    string del2 = "Edges";
    string del3 = "Ref";
    string del4 = "Vertices";
    string del5 = "Original";

    int end1 = data.find(del1) + del1.length();
    int start2 = data.find(del2); // - end1;
    int end2 = start2 + del2.length();
    int start3 = data.find(del3); // - end2;
    int end3 = start3 + del3.length();
    int start4 = data.find(del4); // - end3;
    int end4 = start4 + del4.length();
    int start5 = data.find(del5); // - end4;
    int end5 = start5 + del5.length();
    int datasize = data.length();

    string raw_atoms = data.substr(end1, start2 - end1);
    string raw_edges = data.substr(end2, start3 - end2);
    string raw_refs  = data.substr(end3, start4 - end3);
    string raw_verts = data.substr(end4, start5 - end4);
    string raw_orig  = data.substr(end5, datasize + 1);

    // Get vertex identities and store in in a map (dict)
    vector<string> new_verts = split(raw_verts, ",");
    map<string, string> verts; //[new_verts.size()];
    for (size_t i = 0; i < new_verts.size(); i++){
        vector<string> small = split(new_verts[i].c_str(), "\'");
        if (i == 0) {
            verts.insert(make_pair(small[1], small[3]));
        }
        else {
            verts.insert(make_pair(small[0], small[2]));
        }
    }

    return verts;
}


vector<Node> cartesian_product(map<string, string> set1, map<string, string> set2)
{   // Calculates the cartesian product of two sets - only accepts vertices
    // which have the same atom type
    vector<Node> node_list;
    for(map<string, string>::iterator it=set1.begin(); it!=set1.end(); ++it){
        string outer_tag = it->first;
        string outer_idn = it->second;
        for(map<string, string>::iterator in=set2.begin(); in!=set2.end(); ++in){
            string inner_tag = in->first;
            string inner_idn = in->second;
            if (inner_idn == outer_idn) {
                Node node(outer_tag, inner_tag, outer_idn);
                node_list.push_back(node);
            }
        }
    }
    return node_list;
}


int string_to_int(string raw_val)
{   // Function that returns the integer id of a node
    // by removing the V and converting the rest to an int
    size_t offset = 0;
    vector<string> vect_val = split(raw_val, "V");
    string raw_str = vect_val[0].c_str();
    int result = stoi(raw_str, &offset);
    return result;
}


string make_label(int val1, int val2)
{   // Makes the labels for the edges, makes sure
    // smaller number comes first
    string one, two;
    if (val1 == val2) {
        cout << "ERROR VAL1 = VAL2" << endl;
        exit(EXIT_FAILURE);
    }
    long double val1_d = 1.0 * val1;
    long double val2_d = 1.0 * val2;
    string val1_sraw = to_string(val1_d);
    string val2_sraw = to_string(val2_d);
    vector<string> val1_v = split(val1_sraw, ".");
    vector<string> val2_v = split(val2_sraw, ".");
    string val1_s = val1_v[0].c_str();
    string val2_s = val2_v[0].c_str();
    if (val1 > val2) {
        one = val2_s;
        two = val1_s;
    }
    else {
        one = val1_s;
        two = val2_s;
    }
    string tag = (one + "-" + two);
    return tag;
}


vector<Edge> correspondence_edges(vector<Node> node_list, map<string, double> edge_a, map<string, double> edge_b)
{   // Finds the edges in the correspondence plot
    bool dynamic = false;
    double tol;
    vector<Edge> edges;
    if (dynamic) {
        tol = 0.2; // Dynamic
    }
    else {
        tol = 0.4; // Same cutoff pete used
    }
    for (size_t i = 0; i < node_list.size(); i++) {
        string node1_a = node_list[i].node_1;
        int val_1a = string_to_int(node1_a);
        string node1_b = node_list[i].node_2;
        int val_1b = string_to_int(node1_b);
        for (size_t j = 0; j < node_list.size(); j++){
            if (i == j){
                continue;
            }
            bool keep = false;
            string node2_a = node_list[j].node_1;
            string node2_b = node_list[j].node_2;
            int val_2a = string_to_int(node2_a);
            int val_2b = string_to_int(node2_b);
            if (val_1a == val_2a || val_1b == val_2b){
                continue;
            }
            string label_a = make_label(val_1a, val_2a);
            string label_b = make_label(val_1b, val_2b);
            long double dist_a = edge_a[label_a];
            long double dist_b = edge_b[label_b];
            if (dynamic) {
                double max_d = dist_a + (dist_a * tol);
                double min_d = dist_a - (dist_a * tol);
                if (dist_b > min_d && dist_b < max_d) {
                    keep = true;
                    //cout << dist_a << " " << dist_b << endl;
                }
            }
            else {
                long double diff = fabs(dist_a - dist_b);
                if (diff < tol) {
                    keep = true;
                    //cout << dist_a << " " << dist_b << endl;
                }
            }
            if (!keep) {
                continue;
            }
            Edge this_edge(node1_a, node1_b, node2_a, node2_b, dist_a, dist_b);
            edges.push_back(this_edge);
        }
    }
    return edges;
}


map<string, Conn> connectivity(vector<Node> nodes, vector<Edge> edges)
{   // Build a connectivity table for the graph
    map<string, Conn> connects;
    vector<string> done;
    for (size_t i = 0; i < nodes.size(); i++) {
        string node_id = nodes[i].node_1 + "-" + nodes[i].node_2;
        if (find(done.begin(), done.end(), node_id) != done.end()) {
            continue;
        }
        done.push_back(node_id);
        Conn node_info;
        node_info.number = 0;
        for (size_t j = 0; j < edges.size(); j++) {
            string node_ida = (edges[j].vert1_a + "-" + edges[j].vert1_b);
            string node_idb = (edges[j].vert2_a + "-" + edges[j].vert2_b);
            vector<string> curr = node_info.connects;
            if (node_id == node_ida) {
                if (find(curr.begin(), curr.end(), node_idb) != curr.end()) {
                    continue;
                }
                else {
                    node_info.number++;
                    node_info.connects.push_back(node_idb);
                }
                //cout << "NODE 1" << endl;
            }
            else if (node_id == node_idb) {
                if (find(curr.begin(), curr.end(), node_ida) != curr.end()) {
                    continue;
                }
                else {
                    node_info.number++;
                    node_info.connects.push_back(node_ida);
                }
            }
        }
        connects.insert(make_pair(node_id, node_info));
    }
    return connects;
}


vector<string> maximum_clique(vector<Node> nodes, vector<Edge> edges)
{   // Performs the maximum clique detection
    int min_clique = 4; // This is the smallest clique that will
                        // be considered
    // We must first build the connectivity table
    map<string, Conn> conn = connectivity(nodes, edges);

    // Then we rank their nodes based off connectivity
    vector<int> values;
    int max_conn = 0;
    map<int, vector<string>> classes;
    for(map<string, Conn>::iterator it=conn.begin(); it!=conn.end(); ++it) {
        string node = it->first;
        Conn data = it->second;
        //cout << node << " " << data.number << endl;
        if (find(values.begin(), values.end(), data.number) != values.end()) {
            continue;
        }
        values.push_back(data.number);
        if (data.number > max_conn) {
            max_conn = data.number;
        }
    }
    map<int, vector<string>> ranks;
    for (size_t i = max_conn; i > (min_clique - 1); i--){
        bool found = false;
        bool point = false;
        for(map<int, vector<string>>::iterator it=ranks.begin(); it!=ranks.end(); ++it) {
            int j = it->first;
            if (i == j) {
                found = true;
            }
        }
        vector<string> f_nodes;
        for(map<string, Conn>::iterator it=conn.begin(); it!=conn.end(); ++it) {
            string node = it->first;
            Conn data = it->second;
            if (data.number == i) {
                f_nodes.push_back(node);
                point = true;
                if (found) {
                    ranks[i].push_back(node);
                }
            }
        if (!found && point) {
            ranks.insert(make_pair(i, f_nodes));
        }
        }
    }

    int largest_clique = 0;
    vector<string> cliques;
    // Now that I have ranked all of the nodes we can
    // iterate over everything and perform the maximum clique location
    for (size_t i = max_conn; i > (min_clique - 1); i--) {
        bool found = false;
        if ((i + 1) < largest_clique && i <= (min_clique + 1)) {
            break;
        }
        for(map<int, vector<string>>::iterator it=ranks.begin(); it!=ranks.end(); ++it) {
            int j = it->first;
            if (i == j) {
                found = true;
            }
        }
        if (found) {
            for (size_t k = 0; k < ranks[i].size(); k++) {
                string head_a = ranks[i][k];
                for (size_t h = max_conn; h > (min_clique - 1); h--) {
                    for (size_t p = 0; p < ranks[h].size(); p++) {
                        string head_b = ranks[h][p];
                        if (head_a == head_b || conn[head_b].number < min_clique || conn[head_b].number < largest_clique) {
                            continue;
                        }
                        bool vetted = false;
                        int v_count = 0;
                        int v_size = conn[head_b].number;
                        while (!vetted) {
                            if (v_count >= v_size) {
                                vetted = true;
                                break;
                            }
                            vector<string> branch = {head_a, head_b};
                            vector<bool> success = {false, false};
                            int branch_count = 2;
                            for (size_t stick = v_count; stick < conn[head_b].connects.size(); stick++) {
                                string head_c = conn[head_b].connects[stick];
                                vector<string> sub_cons = conn[head_c].connects;
                                bool sub_success = true;
                                for (size_t brnch = 0; brnch < branch.size(); brnch++){
                                    if(std::find(sub_cons.begin(), sub_cons.end(), branch[brnch]) != sub_cons.end()) {
                                        success[brnch] = true;
                                    }
                                }
                                for (size_t brnch = 0; brnch < success.size(); brnch++) {
                                    if (!success[brnch]) {
                                        sub_success = false;
                                    }
                                    success[brnch] = false;
                                }
                                if (sub_success) {
                                    branch.push_back(head_c);
                                    success.push_back(false);
                                    branch_count++;
                                }
                            }
                            v_count++;
                            if (branch_count > largest_clique) {
                                string branch_id = "";
                                for (size_t brnch = 0; brnch < branch.size(); brnch++) {
                                    branch_id = branch_id + "," + branch[brnch];
                                }
                                cout << branch_count << " " << branch_id << endl;
                                cliques = {branch_id};
                                largest_clique = branch_count;
                            } else if (branch_count == largest_clique) {
                                string branch_id = "";
                                for (size_t brnch = 0; brnch < branch.size(); brnch++) {
                                    branch_id = branch_id + "," + branch[brnch];
                                }
                                bool save = true;
                                if(std::find(cliques.begin(), cliques.end(), branch_id) != cliques.end()) {
                                    save = false;
                                }
                                if (save) {
                                    cliques.push_back(branch_id);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return cliques;
}


void print_cliques(vector<string> cliques, map<string, Atom> origin_a, map<string, string> refs, string loc_a, string loc_b)
{   // prints the cliques to xyz files
    //string source = "C:\\Users\\tomda\\Documents\\2019\\March\\March 19th 2019\\";
    string name;
    vector<string> loc_a_vect = split(loc_a, "-");
    vector<string> loc_b_vect = split(loc_b, "-");
    string name_a = loc_a_vect[0].c_str();
    string name_b = loc_b_vect[0].c_str();
    for (size_t i = 0; i < cliques.size(); i++) {
        ofstream myfile;
        //name = source + "Clique_" + to_string(i) + ".xyz";
        long double i_dbl = 1.0 * i;
        string i_str = to_string(i_dbl);
        vector<string> i_vec = split(i_str, ".");
        string i_s = i_vec[0].c_str();
        name = name_a + "_and_" + name_b + "-";
        name = name + "Clique_" + i_s + ".xyz";
        vector<string> data = split(cliques[i], ",");
        int atom_count = 0;
        for (size_t j = 0; j < data.size(); j++) {
            atom_count++;
        }
        myfile.open (name);
        myfile << atom_count << "\n";
        myfile << "Maximum Clique " << i << " \n";
        for (size_t j = 0; j < data.size(); j++) {
            vector<string> labels = split(data[j], "-");
            string label = refs[labels[0]];
            Atom atom = origin_a[label];
//            cout << label << " " << origin_a[label].ident << endl;;
            myfile << " " << origin_a[label].ident << "\t";
            myfile << origin_a[label].x_coord << "\t";
            myfile << origin_a[label].y_coord << "\t";
            myfile << origin_a[label].z_coord << "\n";
        }
//        cout << atom_count << endl;
        myfile.close();
    }
}


// Main execution of the script
int main(int argc, char*argv[])
{
    // Declare variables
    string data_a, data_b;
//    string loc_a = "C:\\Users\\tomda\\Documents\\2019\\March\\March 19th 2019\\Test.cgraph";
//    string loc_b = "C:\\Users\\tomda\\Documents\\2019\\March\\March 19th 2019\\Test2.cgraph";
//    string atom_types[] = {"C", "H", "N", "O"};

//    char loc_a = argv[1];
//    char loc_b = argv[2];
    string loc_a = argv[1];
    string loc_b = argv[2];
    if (argc < 3) {
        cout << "Need 2 additional command line arguments" << endl;
        exit (EXIT_FAILURE);
    }

    // Load data from plot 1
    data_a = import_file(loc_a);
    map<string, string> verts_a = format_vertices(data_a);
    map<string, double> edges_a = format_edges(data_a);
    map<string, string> refs_a = format_refs(data_a);
    map<string, Atom> origin_a = format_origin(data_a);
    vector<string> atoms_a = format_atoms(data_a);

    // Load data from plot2
    data_b = import_file(loc_b);
    map<string, string> verts_b = format_vertices(data_b);
    map<string, double> edges_b = format_edges(data_b);
    map<string, string> refs_b = format_refs(data_b);
    map<string, Atom> origin_b = format_origin(data_b);
    vector<string> atoms_b = format_atoms(data_b);

    // Find the nodes of the correspondence graph
    vector<Node> corr_nodes = cartesian_product(verts_a, verts_b);

    // Find the edges of the correspondence graph
    vector<Edge> corr_edges = correspondence_edges(corr_nodes, edges_a, edges_b);

    // Now that we have our correspondence graph we can start the
    // maximum clique search
    vector<string> cliques = maximum_clique(corr_nodes, corr_edges);

    cout << "\nCliques Located:" << endl;
    for (size_t i = 0; i < cliques.size(); i++) {
        cout << cliques[i] << endl;
    }

    print_cliques(cliques, origin_a, refs_a, loc_a, loc_b);

    return 0;
}
