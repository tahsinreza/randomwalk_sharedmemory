#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "omp.h"

#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>

template <typename RandomNumberEngine, typename Uint>
Uint uinform_random_uint(RandomNumberEngine& rnd_eng, Uint a, Uint b) {
  std::uniform_int_distribution<Uint> uni_int_dist(a, b);
  return uni_int_dist(rnd_eng);
}

template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList, typename Uint>
class random_walk {
  public: 
    random_walk(Uint _walker_count, Uint _walker_hop_count, 
      size_t _host_thread_count) :  
      walker_count(_walker_count), walker_hop_count(_walker_hop_count), 
      host_thread_count(_host_thread_count), rnd_eng_per_thread(0), 
      walker_memory_per_thread(0), walks_output_file_per_thread(_host_thread_count) {
      std::cout << "Random walk initialized with " 
        << walker_count << " walkers." << std::endl;

      // setup parallel random number engine
      std::cout << "Setting up parallel random number engine ... " << std::endl;
      for (size_t tc = 0; tc < host_thread_count; tc++) {
        boost::random::random_device rnd_dev;
        rnd_eng_per_thread.push_back(boost::random::mt19937_64(rnd_dev));
      }

      // setup memory used by walkers
      std::cout << "Setting up parallel walkers ... " << std::endl;
      for (size_t tc = 0; tc < host_thread_count; tc++) {
        std::vector<Vertex> walker_memory(0);
        walker_memory_per_thread.push_back(walker_memory);
      }
    }
  
    ~random_walk() {}

    void do_random_walk(VertexList& vertices, VertexList& vertex_degree,
      EdgeList& edges, std::string walks_output_filepath) {

      // setup parallel output files
      for (size_t tc = 0; tc < host_thread_count; tc++) {
        std::string walks_output_filename = walks_output_filepath +
          "/thread_" + std::to_string(tc);
        walks_output_file_per_thread[tc].open(walks_output_filename, 
          std::ofstream::out);
      }
  
      std::cout << "Random walk in progress ..." << std::endl; 

      #pragma omp parallel for schedule(static, 1)
      for(size_t walker_ID = 0; walker_ID < walker_count; walker_ID++) {
        size_t thread_ID = omp_get_thread_num();
        //std::cout << "Thread# " << thread_ID << " " 
        //  << " Walker# " << walker_ID << std::endl;
       
        bool walker_alive = false;
        Uint current_walker_hop_count = 0;
    
        walker_memory_per_thread[thread_ID].clear();
     
        //boost::random::uniform_int_distribution<Vertex>
        //  uni_int_dist(0, vertices.size() - 1);
        //Vertex rnd_num = uni_int_dist(rnd_eng_per_thread[thread_ID]); 
        //Vertex rnd_num = get_uinform_random_uint
        //  (rnd_eng_per_thread[thread_ID], 0, vertices.size() - 1);
        //assert(rnd_num <= vertices.size() - 1);
        //Vertex current_vertex = vertices[rnd_num];
        //assert(current_vertex <= max_vertex);        
    
        Vertex current_vertex = source_vertex_random(vertices, thread_ID);
        walker_alive = true; 
        walker_memory_per_thread[thread_ID].push_back(current_vertex); 
        current_walker_hop_count++;

        auto max_trials = 5;
        auto trials = 0;  

        do {
          if (vertex_degree[current_vertex] > 0) {
            current_vertex = next_vertex_random(current_vertex, vertices, 
              edges, thread_ID);
            walker_memory_per_thread[thread_ID].push_back(current_vertex);
            current_walker_hop_count++;
            trials = 0; 
          } else {
            trials++;   
          }
        } while (current_walker_hop_count < walker_hop_count && 
          trials < max_trials); 

        //output_walks(walks_output_filepath, walker_ID, thread_ID);
        output_walks(walks_output_file_per_thread[thread_ID], walker_ID, thread_ID);
      } // for 
 
      // close files
      for (auto& f : walks_output_file_per_thread) {
        f.close();
      }
 
      std::cout << "Random walk completed." << std::endl;
    }

  private:
    Vertex source_vertex_random(VertexList& vertices, size_t thread_ID) {
      Vertex rnd_num = uinform_random_uint
        (rnd_eng_per_thread[thread_ID], static_cast<Vertex>(0), 
        static_cast<Vertex>(vertices.size() - 2)); // important // TODO: use vertex_count
      assert(rnd_num <= vertices.size() - 2);
      Vertex current_vertex = rnd_num;    
      //assert(current_vertex <= max_vertex);  
      return current_vertex;  
    }

    Vertex next_vertex_random(Vertex& current_vertex, VertexList& vertices, 
      EdgeList& edges, size_t thread_ID) {
      Edge rnd_num = uinform_random_uint
        (rnd_eng_per_thread[thread_ID], vertices[current_vertex], 
        vertices[current_vertex + 1] - 1);
      //std::cout << rnd_num << " " << edges.size() << std::endl;
      assert(rnd_num <= edges.size() - 1); 
      Vertex next_vertex = edges[rnd_num]; 
      //assert(next_vertex <= max_vertex);  
      return next_vertex;   
    } 

    void output_walks(std::ofstream& walks_output_file, size_t walker_ID, 
      size_t thread_ID) {
      //std::string walks_output_filename = walks_output_filepath + 
      //  "/walker_" + std::to_string(walker_ID);
      //std::ofstream walks_output_file(walks_output_filename, 
      //  std::ofstream::out);
      walks_output_file << walker_ID << ", ";
      for (size_t v = 0; v < walker_memory_per_thread[thread_ID].size(); v++) {
        walks_output_file << walker_memory_per_thread[thread_ID][v] << ", ";  
      }  
      walks_output_file <<"\n";       
      //walks_output_file.close(); 
    } 

    void output_walks(std::string walks_output_filepath, size_t walker_ID, 
      size_t thread_ID) {
      std::string walks_output_filename = walks_output_filepath + 
        "/walker_" + std::to_string(walker_ID);
      std::ofstream walks_output_file(walks_output_filename, 
        std::ofstream::out);
      for (size_t v = 0; v < walker_memory_per_thread[thread_ID].size(); v++) {
        walks_output_file << walker_memory_per_thread[thread_ID][v] << ", ";  
      }  
      walks_output_file <<"\n";       
      walks_output_file.close(); 
    } 

    size_t host_thread_count;
    Uint walker_count;
    Uint walker_hop_count;
    std::vector<boost::random::mt19937_64> rnd_eng_per_thread;
    std::vector<std::vector<Vertex>> walker_memory_per_thread;
    std::vector<std::ofstream> walks_output_file_per_thread;
};
