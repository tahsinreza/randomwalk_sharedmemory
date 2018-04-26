//#include <atomic>
//#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "omp.h"

#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/discrete_distribution.hpp>

//#include "graph.hpp"
#include "random_walk.hpp"
#include "util.hpp"

enum DISTRIBUTION {UNIFORM, RANDOM, BIASED, ALL};

template <typename Vertex, typename Edge, typename EdgeListTuple, typename Metadata>
Edge read_edge_list_file(const std::string input_filename, EdgeListTuple& edge_list,
  Vertex& max_vertex, Edge skip_lines_count = 0) {
  std::ifstream input_file(input_filename, std::ifstream::in);
  Edge edge_count(0);
  Vertex s(0), t(0);
  Metadata a(0), b(0), c(0), d(0);
  std::string line;

  while(std::getline(input_file, line)) {
    std::istringstream iss(line);
    edge_count++;
    //if (edge_count >= skip_lines_count || skip_lines_count == 0) {
      iss >> s >> t >> a >> b >> c >> d;
      edge_list.push_back(std::forward_as_tuple(s, t, a, b, c, d));

      auto temp_max_vertex = s >= t ? s : t;

      if (max_vertex < temp_max_vertex) {
        max_vertex = temp_max_vertex;
      }

    //}
  }
  input_file.close();
  edge_count-= skip_lines_count;
  assert(edge_count > 0);
  return edge_count;
}

template <typename EdgeListTuple, typename EdgeList>
void generate_edge_list(EdgeListTuple& edge_list, EdgeList& edges) {
  for (size_t e = 0; e < edge_list.size(); e++) {
    edges.push_back(std::get<1>(edge_list[e]));
  }
}

template <typename EdgeListTuple, typename EdgeList>
void generate_edge_metadata(EdgeListTuple& edge_list, EdgeList& edges) {
  for (size_t e = 0; e < edge_list.size(); e++) {
    //edges.push_back(std::get<1>(edge_list[e]));
    std::cout 
      << std::get<0>(edge_list[e]) << " " 
      << std::get<1>(edge_list[e]) << " "
      << std::get<2>(edge_list[e]) << " "
      << std::get<3>(edge_list[e]) << " "
      << std::get<4>(edge_list[e]) << " "
      << std::get<5>(edge_list[e]) << " "      
      << std::endl;
  }
}

template <typename EdgeListTuple, typename VertexList, typename Vertex>          
Vertex generate_vertex_list(EdgeListTuple& edge_list, VertexList& vertices,
  VertexList& vertex_degree, const Vertex max_vertex) {

  //std::ofstream vertex_file("vertex_file_tmp", std::ofstream::out);

  Vertex vertex_count = 0;
  //Vertex max_vertex = std::get<0>(edge_list[edge_list.size() - 1]);
  Vertex l = 0; // edge list index
  Vertex degree = 0;
  Vertex current_vertex = vertex_count;
  Vertex source;
  //Vertex target;

  do {
    auto edge = edge_list[l];
    source = std::get<0>(edge);
    //target = std::get<1>(edge);
    if (source == current_vertex) {
      degree++;
      l++;
    } else {
      //std::cout << current_vertex << std::endl;
      vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
                                                     vertices[vertices.size() - 1]));
      vertex_degree.push_back(degree);

      //VertexData v_data = get_random_uint(rnd_a, rnd_b, rnd_eng);
      //vertex_data.push_back(v_data);

      // write vertex info to file
      //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " " << "\n";

      //vertex_data_file << current_vertex << " " << v_data << "\n";

      // update vertices array
      degree = 0;
      vertex_count++;
      current_vertex = vertex_count;
    }
  } while(current_vertex <= max_vertex);

  // add the last dummy vertex to the vertex list
  vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
                                                 vertices[vertices.size() - 1]));
  //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " "
  //  << "\n";
  //vertex_file.close();
  //vertex_data_file.close();

  return vertex_count;
}

void usage(uint64_t walker_count, uint64_t walker_hop_count, 
  size_t thread_count, std::string dist_opt) {  
  std::cerr << "Usage: -d <string> -w <uint> -s <uint> -t <uint> [edgelist_input_filename walks_output_filepath]\n"
    //[vertex_input_file edge_input_file vertex_data_input_file \
    //pattern_input_file vertex_rank_output_file]\n"
    << " -d <string>   - distribution type (uniform, random, or biased, Default is " << dist_opt << ")\n"
    << " -h            - print help and exit\n"
    << " -w <uint>     - number of walkers (nonzero unsigned integer, Default is " << walker_count << ")\n"
    << " -s <uint>     - walker max step count (nonzero unsigned integer, Default is " << walker_hop_count << ")\n"
    << " -t <uint>     - number of CPU threads (nonzero unsigned integer, Default is " << thread_count << ")\n"
    << "[file ...] - list of input and output files (required)\n\n";
}

void parse_cmd_line(int argc, char** argv, 
  std::string& vertex_input_filename, std::string& edge_input_filename, 
  std::string& vertex_data_input_filename, std::string& pattern_input_filename, 
  std::string& vertex_rank_output_filename, std::string& walks_output_filename, 
  uint64_t& walker_count, uint64_t& walker_hop_count, size_t& thread_count, 
  size_t default_thread_count, DISTRIBUTION& wlkr_dist) {

  std::cout << "CMD line:";
  for (int i=0; i<argc; ++i) {
    std::cout << " " << argv[i];
  }
  std::cout << std::endl;
  std::cout << std::endl;

  std::vector<std::string> dist_opts = {"uniform", "random", "biased", "all"};
  std::string dist_opt = dist_opts[wlkr_dist];
  //uint64_t walker_count_d = walker_count;
  //uint64_t walker_hop_count_d = walker_hop_count;

  if (argc < 3) { 
    std::cerr << "Too few arguments."<<std::endl;
    usage(walker_count, walker_hop_count, default_thread_count, dist_opt);
    exit(-1);
  }

  bool prn_help = false;
  char c;

  while ((c = getopt(argc, argv, "d:w:s:t:h ")) != -1) {
     switch (c) {
       case 'h':
         prn_help = true;
         break;
       case 'd':
         prn_help = true;
         for (size_t i = 0; i < dist_opts.size(); i++) {
           if(boost::iequals(dist_opts[i], optarg)) {
             wlkr_dist = static_cast<DISTRIBUTION>(i);
             prn_help = false;
             break;
           }
         }
         break;
      case 'w':
         walker_count = std::stoull(optarg);
         break;
      case 's':
         walker_hop_count = std::stoull(optarg);
         break;
      case 't':
         thread_count = std::stoull(optarg);
         break;
      default:
         std::cerr << "Unrecognized option: " << c << "." <<std::endl;
         prn_help = true;
         break;
     }
   }

   if (prn_help) {
     usage(walker_count, walker_hop_count, default_thread_count, dist_opt);
     exit(-1);
   }

   // optind is initialized to 1
   /*vertex_input_filename = argv[optind];
   edge_input_filename = argv[optind + 1];
   vertex_data_input_filename = argv[optind + 2];
   //pattern_input_filename = argv[optind + 3];
   //vertex_rank_output_filename = argv[optind + 4];*/

   edge_input_filename = argv[optind];
   walks_output_filename = argv[optind + 1]; 

   /*std::cout << "\n" //<< vertex_input_filename << " "
     << edge_input_filename << " "
     //<< vertex_data_input_filename << " "
     //<< pattern_input_filename << " "
     //<< vertex_rank_output_filename << " "
     << walks_output_filename << " "
     //<< walker_count << " "
     //<< wlkr_dist 
     << "\n" << std::endl;*/
}

int main(int argc, char** argv) {
  
  // parse commandline input
  std::string vertex_input_filename = "dummy";
  std::string edge_input_filename;
  std::string vertex_data_input_filename = "dummy";
  std::string edge_data_input_filename = "dummy"; 
  std::string pattern_input_filename = "dummy";
  std::string vertex_rank_output_filename = "dummy";
  std::string walks_output_filename; // filepath
  uint64_t walker_count = 10;
  uint64_t walker_hop_count = 5;
  size_t thread_count = 0;
  DISTRIBUTION wlkr_dist = RANDOM;

  size_t host_thread_count = 0;

  {
  #pragma omp parallel
  { 
    host_thread_count = omp_get_num_threads();
  }
  //std::cout << "Number of available threads: " << host_thread_count 
  //  << std::endl;
  }

  parse_cmd_line(argc, argv, vertex_input_filename, edge_input_filename,
    vertex_data_input_filename, pattern_input_filename,
    vertex_rank_output_filename, walks_output_filename,
    walker_count, walker_hop_count, thread_count, host_thread_count, 
    wlkr_dist);   
 
  //std::cout << "Application ... " << std::endl;

  std::chrono::time_point<std::chrono::steady_clock> start_time;
  std::chrono::time_point<std::chrono::steady_clock> end_time;
  double elapsed_time;

  // host info
  char hostname[HOST_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  std::cout << "Hostname: " << hostname << std::endl;
   
  //size_t host_thread_count = 0;
   
  //{
  //#pragma omp parallel
  //{ 
  //  host_thread_count = omp_get_num_threads();
  //}
  std::cout << "Number of available threads on " << hostname << ": " 
    << host_thread_count << std::endl;
  //}
   
  if (thread_count > 0) {
    omp_set_dynamic(0);
    omp_set_num_threads(thread_count);
    {
    #pragma omp parallel
    {
      host_thread_count = omp_get_num_threads();
    }
    std::cout << "Number of threads to be used: " << host_thread_count << std::endl;
    } 
  }

  std::cout << std::endl;
 
  //////////////////////////////////////////////////////////////////////////////

  // graph construction
  typedef uint64_t Vertex;
  typedef uint64_t Edge;
  typedef uint64_t VertexData;
  typedef uint64_t EdgeData;
  typedef uint64_t Metadata;

  //typedef std::vector<std::tuple<Vertex, Vertex>> EdgeListTuple;
  typedef std::vector<std::tuple<Vertex, Vertex, Vertex, Vertex, Vertex, Vertex>> EdgeListTuple;
  EdgeListTuple edge_list(0);
  typedef std::vector<Vertex> EdgeList;
  EdgeList edges(0);
  typedef std::vector<Edge> VertexList;
  VertexList vertices(0);
  VertexList vertex_degree(0);

  Vertex vertex_count = 0;
  Edge edge_count = 0;
  Vertex max_vertex = 0;
  Edge max_degree = 0;
 
  // build in-memory CSR graph  

  Edge skip_lines_count = 0; // TODO: read from the commandline 

  std::cout << "Building in-memory CSR graph ... " << std::endl;

  std::chrono::time_point<std::chrono::steady_clock> global_start_time =
  std::chrono::steady_clock::now();

  // read an edgelist file
  std::cout << "Reading edgelist file ..." << std::endl;
  edge_count = read_edge_list_file<Vertex, Edge, EdgeListTuple, Metadata>
    (edge_input_filename, edge_list, max_vertex, skip_lines_count);
  std::cout << "Size of edge list: " << edge_list.size() << std::endl;
  std::cout << "Max vertex: " << max_vertex << std::endl;

  // sort edges by source
  std::cout << "Sorting edges by source vertex ..." << std::endl;
  std::stable_sort(edge_list.begin(), edge_list.end(),
    [](const std::tuple<Vertex, Vertex, Vertex, Vertex, Vertex, Vertex>& a,
       const std::tuple<Vertex, Vertex, Vertex, Vertex, Vertex, Vertex>& b) -> bool {
         return std::get<0>(a) < std::get<0>(b);
       });

  // generate vetex list
  //std::cout << "Generating vertex list ..." << std::endl;
  std::cout << "Creating CSR row pointers ..." << std::endl;
  vertex_count = generate_vertex_list(edge_list, vertices, vertex_degree, max_vertex);
  std::cout << "Size of vertex list: " << vertices.size() << std::endl;
  std::cout << "Size of vertex degree list: " << vertex_degree.size() << std::endl;

  // sort targets in increasing order
  std::cout << "Sorting neighbours ..." << std::endl;
  #pragma omp parallel for
  for (size_t v = 0; v < vertices.size() - 1; v++) {
    size_t start = vertices[v];
    size_t end = vertices[v + 1];
    std::stable_sort(edge_list.begin() + start,
                     edge_list.begin() + end,
       [](const std::tuple<Vertex, Vertex, Vertex, Vertex, Vertex, Vertex>& a,
          const std::tuple<Vertex, Vertex, Vertex, Vertex, Vertex, Vertex>& b) -> bool {
            return std::get<1>(a) < std::get<1>(b);
          });
  }

  // generate edge list
  std::cout << "Generating edge list ..." << std::endl;    
  generate_edge_list(edge_list, edges);
  std::cout << "Size of edge list: " << edges.size() << std::endl;

  std::cout << "CSR Graph generation completed." << std::endl;
  std::cout << "Number of vertices: " << vertex_count << std::endl;
  std::cout << "Number of edges: " << edge_count << std::endl;  
  std::cout << "Max vertex: " << max_vertex << std::endl;

  // Test
  //for (Vertex v = 0; v < vertex_count; v++) {
    //std::cerr << v << " " << vertex_degree[v] << " " << vertices[v] << std::endl; 
  //  for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
  //    uint64_t v_nbr = edges[e];
  //    assert(v_nbr <= max_vertex);
  //
  //  }  
  //} 
 
  //std::cout << "Generating edge metadata ..." << std::endl;
  //generate_edge_metadata(edge_list, edges); 
   
  // Test

  std::chrono::time_point<std::chrono::steady_clock> global_end_time =
    std::chrono::steady_clock::now();
  double global_elapsed_time = getElapsedTimeSecond(global_start_time, 
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds." 
    << std::endl;

  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////

  // random walk

  global_start_time = std::chrono::steady_clock::now();

  std::cout << "Random walk ..." << std::endl; 

  typedef random_walk<Vertex, Edge, VertexList, EdgeList, uint64_t> RandomWalk;
  RandomWalk rw(walker_count, walker_hop_count, host_thread_count); 
  rw.do_random_walk(vertices, vertex_degree, edges, walks_output_filename);

  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  std::cout << std::endl;
   
  //////////////////////////////////////////////////////////////////////////////
 
  return 0;  
} // end of main
