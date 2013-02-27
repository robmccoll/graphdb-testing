/*  _________________________________________________________________________
 *
 *  MTGL: The MultiThreaded Graph Library
 *  Copyright (c) 2008 Sandia Corporation.
 *  This software is distributed under the BSD License.
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *  For more information, see the README file in the top MTGL directory.
 *  _________________________________________________________________________
 */

/****************************************************************************/
/*! \file persist_graph.cpp

    \brief Creates or reads in a graph or vector and then mmaps the object's
           data to make the object persist in memory.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 10/7/2010
*/
/****************************************************************************/

#include <string>
#include <iomanip>
#include <cstdio>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/stinger_graph_adapter.hpp>
#include <mtgl/mtgl_io.hpp>
#include <mtgl/generate_rmat_graph.hpp>
#include <mtgl/dynamic_array.hpp>

#define MAX_NUM_WORDS 10
#define MAX_NUM_OBJECTS 100
#define WORD_SIZE 256

using namespace mtgl;

typedef stinger_graph<int64_t, int64_t, int64_t> StingGraph;

struct ObjectInfo {
  char name[WORD_SIZE];
  unsigned long order;
  unsigned long size;
  unsigned long type;
  size_t mem_size;
  void* mem_ptr;
};

int parse_words(const std::string& input, char words[][WORD_SIZE])
{
  int num_words = 0;

  // Find the beginning of the first word.
  size_t beg_pos = 0;
  for ( ; beg_pos < input.size() && isspace(input[beg_pos]); ++beg_pos);

  while (beg_pos < input.size() && num_words < MAX_NUM_WORDS)
  {
    // Find the end of the word.
    size_t end_pos = beg_pos;
    for ( ; end_pos < input.size() && !isspace(input[end_pos]); ++end_pos);

    std::strcpy(words[num_words],
                input.substr(beg_pos, end_pos - beg_pos).c_str());

    ++num_words;

    // Find the beginning of the next word.
    beg_pos = end_pos;
    for ( ; beg_pos < input.size() && isspace(input[beg_pos]); ++beg_pos);
  }

  return num_words;
}

// Count until 0 and chop of every tenth
int num_digits(unsigned long num)
{
  // Handle 0.
  if (num == 0) return 1;

  int n = 0;

  while (num > 0)
  {
    ++n;
    num /= 10;
  }

  return n;
}

template <typename T>
void add_vector(T* data, unsigned long size, char* vector_name,
                ObjectInfo object_list[], int& num_objects)
{
  // Make sure the name is unique.
  int i = 0;
  for ( ; i < num_objects && strcmp(vector_name, object_list[i].name) != 0;
       ++i);

  if (i < num_objects)
  {
    std::cout << "Duplicate name: " << vector_name << std::endl;
    return;
  }

  unsigned long mapped_size = 4 * sizeof(unsigned long) + size * sizeof(T);

  int mmap_prots = PROT_READ | PROT_WRITE;
  int mmap_flags = MAP_SHARED;

#ifdef __MTA__
  mmap_flags |= MAP_ANON;
  int fd = open(vector_name, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
#else
  int fd = shm_open(vector_name, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
#endif

  if (fd < 0)
  {
    std::cout << "Error opening file: "<< vector_name << std::endl;
    return;
  }

#ifndef __MTA__
  if (ftruncate(fd, mapped_size) != 0)
  {
    std::cout << "Error resizing memory region: "<< vector_name << "    size: "
              << mapped_size << std::endl;
    return;
  }
#endif

  void* mapped_mem = mmap(0, mapped_size, mmap_prots, mmap_flags, fd, 0);

  if (mapped_mem == MAP_FAILED)
  {
    close(fd);

#ifdef __MTA__
    remove(vector_name);
#else
    shm_unlink(vector_name);
#endif

    std::cout << "Error mapping memory: "<< vector_name << std::endl;
    return;
  }

  unsigned long* ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);

  // Write the size and object type to the mapped memory header.
  ul_mapped_mem[0] = mmap_traits<T>::type;
  ul_mapped_mem[1] = mapped_size;
  ul_mapped_mem[2] = size;

  // Get a pointers to the location in the mapped memory for the vector.
  T* mmap_data = reinterpret_cast<T*>(ul_mapped_mem + 4);

  // Copy the vector to the mapped memory.
  memcpy(mmap_data, data, size * sizeof(T));

  close(fd);

  strcpy(object_list[num_objects].name, vector_name);
  object_list[num_objects].order = size;
  object_list[num_objects].size = 0;
  object_list[num_objects].type = mmap_traits<T>::type;
  object_list[num_objects].mem_size = mapped_size;
  object_list[num_objects].mem_ptr = mapped_mem;

  ++num_objects;

  std::cout << "Successfully mmapped " << vector_name << ": "
            << size << std::endl;
}

template <typename T>
void add_vector(char* file_name, char* vector_name, ObjectInfo object_list[],
                int& num_objects)
{
  // Read in the array.
  unsigned long size;
  T* v = read_array<T>(file_name, size);

  add_vector(v, size, vector_name, object_list, num_objects);
}

template <typename Graph>
void add_graph(Graph& g, char* file_name, char* graph_name,
               ObjectInfo object_list[], int& num_objects)
{
  // Make sure the name is unique.
  int i = 0;
  for ( ; i < num_objects && strcmp(graph_name, object_list[i].name) != 0; ++i);

  if (i < num_objects)
  {
    std::cout << "Duplicate name: " << graph_name << std::endl;
    return;
  }

  // Determine which input format is used: automatically generated rmat or
  // file-based input.
  int use_rmat = 1;
  int llen = strlen(file_name);

  for (i = 0; i < llen; ++i)
  {
    if (!(file_name[i] >= '0' && file_name[i] <= '9'))
    {
      // String contains non-numeric characters; it must be a filename.
      use_rmat = 0;

      break;
    }
  }

  if (use_rmat)
  {
    generate_rmat_graph(g, atoi(file_name), 8);
  }
  else if (llen > 3 && strcmp(&file_name[llen-3], "mtx") == 0)
  {
    // Matrix-market input.
    dynamic_array<double> vals;

    read_matrix_market(g, file_name, vals);

    if (vals.size() > 0)
    {
      int glen = strlen(graph_name);
      char vals_name[WORD_SIZE];

      strcpy(vals_name, graph_name);
      strcpy(&vals_name[glen], ".vals");

      add_vector(vals.get_data(), vals.size(), vals_name,
                 object_list, num_objects);
    }
  }
  else if (llen > 4 && strcmp(&file_name[llen-4], "srcs") == 0)
  {
    // Snapshot input: <file>.srcs, <file>.dests
    char* srcs_fname = file_name;
    char dests_fname[WORD_SIZE];

    strcpy(dests_fname, srcs_fname);
    strcpy(&dests_fname[llen - 4], "dests");

    read_binary(g, srcs_fname, dests_fname);
  }
  else if (llen > 6 && strcmp(&file_name[llen-6], "dimacs") == 0)
  {
    // DIMACS input.
    dynamic_array<long> weights;

    read_dimacs(g, file_name, weights);

    int glen = strlen(graph_name);
    char weights_name[WORD_SIZE];

    strcpy(weights_name, graph_name);
    strcpy(&weights_name[glen], ".weights");

    add_vector(weights.get_data(), weights.size(), weights_name,
               object_list, num_objects);
  }
  else
  {
    std::cout << "Invalid filename: '" << file_name << "'" << std::endl;
    return;
  }

  if (num_vertices(g) == 0)
  {
    std::cout << "Error reading graph." << std::endl;
    return;
  }

  int mmap_prots = PROT_READ | PROT_WRITE;
  int mmap_flags = MAP_SHARED;

#ifdef __MTA__
  mmap_flags |= MAP_ANON;
  int fd = open(graph_name, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
#else
  int fd = shm_open(graph_name, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
#endif

  if (fd < 0)
  {
    std::cout << "Error opening file: "<< graph_name << std::endl;
    return;
  }

  unsigned long mapped_size = g.get_mmap_size();

#ifndef __MTA__
  if (ftruncate(fd, mapped_size) != 0)
  {
    std::cout << "Error resizing memory region: "<< graph_name << "    size: "
              << mapped_size << std::endl;
    return;
  }
#endif

  void* mapped_mem = mmap(0, mapped_size, mmap_prots, mmap_flags, fd, 0);

  if (mapped_mem == MAP_FAILED)
  {
    close(fd);

#ifdef __MTA__
    remove(graph_name);
#else
    shm_unlink(graph_name);
#endif

    std::cout << "Error mapping memory: "<< graph_name << std::endl;
    return;
  }

  g.write_mmap(mapped_mem);

  close(fd);

  strcpy(object_list[num_objects].name, graph_name);
  object_list[num_objects].order = num_vertices(g);
  object_list[num_objects].size = num_edges(g);
  object_list[num_objects].type = mmap_traits<Graph>::type;
  object_list[num_objects].mem_size = mapped_size;
  object_list[num_objects].mem_ptr = mapped_mem;

  ++num_objects;

  std::cout << "Successfully mmapped " << graph_name << ": "
            << num_vertices(g) << ", " << num_edges(g) << std::endl;
}

template <typename Graph>
void add_graph(char* file_name, char* graph_name, ObjectInfo object_list[],
               int& num_objects)
{
  Graph g;

  add_graph(g, file_name, graph_name, object_list, num_objects);
}

template <>
void add_graph<stinger_graph_adapter<StingGraph> >(
  char* file_name, char* graph_name, ObjectInfo object_list[], int& num_objects)
{
  StingGraph sg(1);
  stinger_graph_adapter<StingGraph> g(sg);

  add_graph(g, file_name, graph_name, object_list, num_objects);
}

template <typename Graph>
void print_graph(Graph& g, ObjectInfo& object_info)
{
  bool print_it = true;

  // Prompt the user if the graph is large.
  if (object_info.order > 50 || object_info.size > 200)
  {
    std::cout << "This graph is large.  Do you really want to print? "
              << "(y/n) ";

    std::string response;
    std::getline(std::cin, response);

    if (response != "y" && response != "Y") print_it = false;
  }

  if (print_it)
  {
    // Initialize the graph from memory.
    g.read_mmap(object_info.mem_ptr);

    // Print the graph.
    std::cout << "Graph: (" << num_vertices(g) << ", " << num_edges(g)
              << ")" << std::endl;
    print(g);
  }
}

template <typename Graph>
void print_graph(ObjectInfo& object_info)
{
  Graph g;

  print_graph(g, object_info);
}

template <>
void print_graph<stinger_graph_adapter<StingGraph> >(ObjectInfo& object_info)
{
  StingGraph sg(1);
  stinger_graph_adapter<StingGraph> g(sg);

  print_graph(g, object_info);
}

template <typename T>
void print_vector(ObjectInfo& object_info)
{
  bool print_it = true;

  // Prompt the user if the vector is large.
  if (object_info.order > 50)
  {
    std::cout << "This vector is large.  Do you really want to print? "
              << "(y/n) ";

    std::string response;
    std::getline(std::cin, response);

    if (response != "y" && response != "Y") print_it = false;
  }

  if (print_it)
  {
    // Initialize the vector from memory.
    unsigned long* temp = reinterpret_cast<unsigned long*>(object_info.mem_ptr);
    T* vec = reinterpret_cast<T*>(temp + 4);

    // Print the vector.
    std::cout << "Size: " << object_info.order << std::endl;
    for (unsigned long i = 0; i < object_info.order; ++i)
    {
      std::cout << vec[i] << std::endl;
    }
  }
}

int main(int argc, char* argv[])
{
  typedef compressed_sparse_row_graph<directedS> DGraph;
  typedef compressed_sparse_row_graph<undirectedS> UGraph;
  typedef compressed_sparse_row_graph<bidirectionalS> BGraph;
  typedef stinger_graph_adapter<StingGraph> SGraph;

  char words[MAX_NUM_WORDS][WORD_SIZE];
  ObjectInfo object_list[MAX_NUM_OBJECTS];
  int num_objects = 0;

  std::cout << "> ";
  std::string input;
  std::getline(std::cin, input);

  while (input != "quit")
  {
    int num_words = parse_words(input, words);

    std::cout << std::endl;

    if (num_words == 0)
    {
      // Do nothing.  Just print the prompt again.
    }
    if (strcmp(words[0], "add") == 0 && (num_words == 3 || num_words == 4))
    {
      if (num_words == 3 || strcmp(words[3], "csr_directed") == 0 ||
          strcmp(words[3], "csrd") == 0)
      {
        add_graph<DGraph>(words[1], words[2], object_list, num_objects);
      }
      else if (strcmp(words[3], "csr_undirected") == 0 ||
               strcmp(words[3], "csru") == 0)
      {
        add_graph<UGraph>(words[1], words[2], object_list, num_objects);
      }
      else if (strcmp(words[3], "csr_bidirectional") == 0 ||
               strcmp(words[3], "csrb") == 0)
      {
        add_graph<BGraph>(words[1], words[2], object_list, num_objects);
      }
      else if (strcmp(words[3], "stinger") == 0 ||
               strcmp(words[3], "sting") == 0)
      {
        add_graph<SGraph>(words[1], words[2], object_list, num_objects);
      }
      else if (strcmp(words[3], "char_vector") == 0 ||
               strcmp(words[3], "cvec") == 0)
      {
        add_vector<char>(words[1], words[2], object_list, num_objects);
      }
      else if (strcmp(words[3], "int_vector") == 0 ||
               strcmp(words[3], "ivec") == 0)
      {
        add_vector<int>(words[1], words[2], object_list, num_objects);
      }
      else if (strcmp(words[3], "unsigned_int_vector") == 0 ||
               strcmp(words[3], "uivec") == 0)
      {
        add_vector<unsigned int>(words[1], words[2], object_list, num_objects);
      }
      else if (strcmp(words[3], "long_vector") == 0 ||
               strcmp(words[3], "lvec") == 0)
      {
        add_vector<long>(words[1], words[2], object_list, num_objects);
      }
      else if (strcmp(words[3], "unsigned_long_vector") == 0 ||
               strcmp(words[3], "ulvec") == 0)
      {
        add_vector<unsigned long>(words[1], words[2], object_list, num_objects);
      }
      else if (strcmp(words[3], "double_vector") == 0 ||
               strcmp(words[3], "dvec") == 0)
      {
        add_vector<double>(words[1], words[2], object_list, num_objects);
      }
      else
      {
        std::cout << "Invalid command." << std::endl;
      }
    }
    else if ((strcmp(words[0], "delete") == 0 || strcmp(words[0], "del") == 0 ||
              strcmp(words[0], "remove") == 0 || strcmp(words[0], "rm") == 0) &&
             num_words == 2)
    {
      // Find the object.
      int i = 0;
      for ( ; i < num_objects && strcmp(words[1], object_list[i].name) != 0;
           ++i);

      if (i >= num_objects)
      {
        std::cout << "Invalid object name." << std::endl << std::endl
                  << "> ";
        std::getline(std::cin, input);
        continue;
      }

      // Delete the memory associated with the object.
      if (munmap((caddr_t) object_list[i].mem_ptr,
                 object_list[i].mem_size) != 0)
      {
        std::cout << "Error unmapping memory: " << object_list[i].name
                  << std::endl;
      }


      // Remove the filename associated with the object.
#ifdef __MTA__
      if (remove(object_list[i].name) != 0)
#else
      if (shm_unlink(object_list[i].name) != 0)
#endif
      {
        std::cout << "Error deleting file: " << object_list[i].name
                  << std::endl;
      }

      // Move all objects after the deleted object forward one position in
      // the object list.
      for ( ; i < num_objects - 1; ++i)
      {
        strcpy(object_list[i].name, object_list[i+1].name);
        object_list[i].order = object_list[i+1].order;
        object_list[i].size = object_list[i+1].size;
        object_list[i].type = object_list[i+1].type;
        object_list[i].mem_size = object_list[i+1].mem_size;
        object_list[i].mem_ptr = object_list[i+1].mem_ptr;
      }

      --num_objects;
    }
    else if ((strcmp(words[0], "list") == 0 || strcmp(words[0], "ls") == 0) &&
             num_words == 1)
    {
      if (num_objects > 0)
      {
        // Find the maximum width of each of the first three columns.
        int col_1_size = num_digits(object_list[0].order);
        int col_2_size = num_digits(object_list[0].size);
        int col_3_size = strlen(object_list[0].name);

        for (int i = 1; i < num_objects; ++i)
        {
          int my_col_1_size = num_digits(object_list[i].order);
          int my_col_2_size = num_digits(object_list[i].size);
          int my_col_3_size = strlen(object_list[i].name);

          if (my_col_1_size > col_1_size) col_1_size = my_col_1_size;
          if (my_col_2_size > col_2_size) col_2_size = my_col_2_size;
          if (my_col_3_size > col_3_size) col_3_size = my_col_3_size;
        }

        for (int i = 0; i < num_objects; ++i)
        {
          std::string object_type;

          if (object_list[i].type == mmap_traits<DGraph>::type)
          {
            object_type = "CSR Directed";
          }
          else if (object_list[i].type == mmap_traits<UGraph>::type)
          {
            object_type = "CSR Undirected";
          }
          else if (object_list[i].type == mmap_traits<BGraph>::type)
          {
            object_type = "CSR Bidirectional";
          }
          else if (object_list[i].type == mmap_traits<SGraph>::type)
          {
            object_type = "STINGER";
          }
          else if (object_list[i].type == mmap_traits<char>::type)
          {
            object_type = "Char Vector";
          }
          else if (object_list[i].type == mmap_traits<int>::type)
          {
            object_type = "Int Vector";
          }
          else if (object_list[i].type == mmap_traits<unsigned int>::type)
          {
            object_type = "Unsigned Int Vector";
          }
          else if (object_list[i].type == mmap_traits<long>::type)
          {
            object_type = "Long Vector";
          }
          else if (object_list[i].type == mmap_traits<unsigned long>::type)
          {
            object_type = "Unsigned Long Vector";
          }
          else if (object_list[i].type == mmap_traits<double>::type)
          {
            object_type = "Double Vector";
          }
          else
          {
            object_type = "Unknown";
          }

          std::cout << std::setw(col_1_size) << std::right
                    << object_list[i].order << "   "
                    << std::setw(col_2_size) << std::right
                    << object_list[i].size << "   "
                    << std::setw(col_3_size) << std::left
                    << object_list[i].name << "   "
                    << std::setw(0) << object_type << std::endl;
        }
      }
    }
    else if (strcmp(words[0], "print") == 0 && num_words == 2)
    {
      // Find the object.
      int i = 0;
      for ( ; i < num_objects && strcmp(words[1], object_list[i].name) != 0;
           ++i);

      if (i >= num_objects)
      {
        std::cout << "Invalid object name." << std::endl << std::endl;
      }
      else if (object_list[i].type == mmap_traits<DGraph>::type)
      {
        print_graph<DGraph>(object_list[i]);
      }
      else if (object_list[i].type == mmap_traits<UGraph>::type)
      {
        print_graph<UGraph>(object_list[i]);
      }
      else if (object_list[i].type == mmap_traits<BGraph>::type)
      {
        print_graph<BGraph>(object_list[i]);
      }
      else if (object_list[i].type == mmap_traits<SGraph>::type)
      {
        print_graph<SGraph>(object_list[i]);
      }
      else if (object_list[i].type == mmap_traits<char>::type)
      {
        print_vector<char>(object_list[i]);
      }
      else if (object_list[i].type == mmap_traits<int>::type)
      {
        print_vector<int>(object_list[i]);
      }
      else if (object_list[i].type == mmap_traits<unsigned int>::type)
      {
        print_vector<unsigned int>(object_list[i]);
      }
      else if (object_list[i].type == mmap_traits<long>::type)
      {
        print_vector<long>(object_list[i]);
      }
      else if (object_list[i].type == mmap_traits<unsigned long>::type)
      {
        print_vector<unsigned long>(object_list[i]);
      }
      else if (object_list[i].type == mmap_traits<double>::type)
      {
        print_vector<double>(object_list[i]);
      }
      else
      {
        std::cout << "Unknown object type." << std::endl << std::endl;
      }
    }
    else if (strcmp(words[0], "help") == 0)
    {
      if (num_words == 1)
      {
        std::cout << "Type 'help <command>' for help on a specific command."
                  << std::endl << std::endl
                  << "Available commands:" << std::endl
                  << "    add" << std::endl
                  << "    delete (del, remove, rm)" << std::endl
                  << "    help" << std::endl
                  << "    list (ls)" << std::endl
                  << "    print" << std::endl
                  << "    quit" << std::endl;
      }
      else if (num_words == 2)
      {
        if (strcmp(words[1], "add") == 0)
        {
          std::cout << "add: Reads or creates a graph or vector and stores "
                    << "it in memory." << std::endl
                    << "Usage: add {<p> | <filename>} <name> <object_type"
                    << " = csr_directed>" << std::endl << std::endl
                    << "  Generates the rmat graph of size <p> or reads the "
                    << "graph or vector from" << std::endl
                    << "  disk given by <filename> into a object of type "
                    << "<object_type>.  Stores" << std::endl
                    << "  the resulting object into memory with the label "
                    << "<name>." << std::endl << std::endl
                    << "  Specifying <p> generates an r-mat graph with 2^p "
                    << "vertices, and specifying" << std::endl
                    << "  filename requests that a DIMACS, MatrixMarket, "
                    << "sources/dests, " << std::endl
                    << "  or vector file be read to build the object."
                    << std::endl << std::endl
                    << "  DIMACS files must end with suffix '.dimacs'."
                    << std::endl
                    << "  MatrixMarket files must end with suffix '.mtx'."
                    << std::endl
                    << "  Sources/dests files must end with suffixes '.srcs' "
                    << " and '.dests'." << std::endl
                    << "  Vector files can have any extension."
                    << std::endl << std::endl
                    << "  Some file formats include extra data such as edge "
                    << "or vertex properties." << std::endl
                    << "  The data is stored in a vector using <name> as the "
                    << "base with an added" << std::endl
                    << "  extension.  DIMACS files have an edge weight that "
                    << "is stored in the vector" << std::endl
                    << "  '<name>.weights'.  Matrix market files have an "
                    << "optional edge value that" << std::endl
                    << "  is stored in the vector '<name>.vals'.  If the file "
                    << "has no edge value," << std::endl
                    << "  then the vector will not be created."
                    << std::endl << std::endl
                    << "  The supported graph types are:" << std::endl
                    << "    csr_directed (csrd)       "
                    << "compressed_sparse_row_graph<directedS>" << std::endl
                    << "    csr_undirected (csru)     "
                    << "compressed_sparse_row_graph<undirectedS>" << std::endl
                    << "    csr_bidirectional (csrb)  "
                    << "compressed_sparse_row_graph<bidirectionalS>"
                    << std::endl << std::endl
                    << "  The supported vector types are:" << std::endl
                    << "    char_vector (cvec)              char*" << std::endl
                    << "    int_vector (ivec)               int*" << std::endl
                    << "    unsigned_int_vector (uivec)     unsigned int*"
                    << std::endl
                    << "    long_vector (lvec)              long*" << std::endl
                    << "    unsigned_long_vector (ulvec)    unsigned long*"
                    << std::endl
                    << "    double_vector (dvec)            double*"
                    << std::endl << std::endl
                    << "  The default object type is csr_directed if no type "
                    << "is given." << std::endl;
        }
        else if (strcmp(words[1], "delete") == 0 ||
                 strcmp(words[1], "del") == 0 ||
                 strcmp(words[1], "remove") == 0 ||
                 strcmp(words[1], "rm") == 0)
        {
          std::cout << "delete (del, remove, rm): Deletes an object from "
                    << "memory." << std::endl
                    << "Usage: delete <name>" << std::endl << std::endl
                    << "  Deletes the object given by <name> from memory."
                    << std::endl;
        }
        else if (strcmp(words[1], "help") == 0)
        {
          std::cout << "help: Descirbes the usage of this program or its "
                    << "commands." << std::endl
                    << "Usage: help <command>" << std::endl << std::endl
                    << "  Gives detailed help information for <command>."
                    << std::endl;
        }
        else if (strcmp(words[1], "list") == 0 || strcmp(words[1], "ls") == 0)
        {
          std::cout << "list (ls): Lists all the objects currently in memory."
                    << std::endl
                    << "Usage: list" << std::endl << std::endl
                    << "  The first column gives the number of vertices for "
                    << "graphs and the size for" << std::endl
                    << "  vectors.  The second column gives the number of "
                    << "edges for graphs and is" << std::endl
                    << "  meaningless for vectors.  The third column gives "
                    << "the name of the object." << std::endl
                    << "  The fourth column gives the object type."
                    << std::endl;
        }
        else if (strcmp(words[1], "print") == 0)
        {
          std::cout << "print: Prints an object to the screen."
                    << std::endl
                    << "Usage: print <name>" << std::endl << std::endl
                    << "  Prints the object given by <name> to the screen.  "
                    << "If the object is a graph" << std::endl
                    << "  and has more than 50 vertices or 200 edges or if "
                    << "the object is a vector" << std::endl
                    << "  with more than 50 entries, the user is prompted "
                    << "to confirm he really wants" << std::endl
                    << "  to print the object." << std::endl;
        }
        else if (strcmp(words[1], "quit") == 0)
        {
          std::cout << "quit: Removes all objects from memory and quits."
                    << std::endl
                    << "Usage: quit" << std::endl;
        }
        else
        {
          std::cout << "Invalid command." << std::endl;
        }
      }
      else
      {
        std::cout << "Invalid command." << std::endl;
      }
    }
    else
    {
      std::cout << "Invalid command." << std::endl;
    }


    std::cout << std::endl << "> ";
    std::getline(std::cin, input);
  }

  for (int i = 0; i < num_objects; ++i)
  {
    // Delete the memory associated with the object.
    if (munmap((caddr_t) object_list[i].mem_ptr, object_list[i].mem_size) != 0)
    {
      std::cout << "Error unmapping memory: " << object_list[i].name
                << std::endl;
    }

#ifdef __MTA__
    if (remove(object_list[i].name) != 0)
#else
    if (shm_unlink(object_list[i].name) != 0)
#endif
    {
      std::cout << "Error deleting file: " << object_list[i].name << std::endl;
    }
  }

  return 0;
}
