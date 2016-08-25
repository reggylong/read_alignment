
#include <array>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <sys/time.h>
#include <sys/resource.h>
#include <thread>
#include <mutex>
#include <unordered_set>
#include <cassert>
#include <cmath>
#include "semaphore.h"

using namespace std;

static const size_t LENGTH = 15;
static const size_t SPACE = pow(4, LENGTH);
static const long THREADS = 6;

static inline uint32_t chtoi(char c) {
  switch (toupper(c)) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    default:
      throw invalid_argument("Expected [AGCT], got " + c);
  }
}

static inline uint32_t ktoi(const string &str) {
  uint32_t base = 1;
  assert(str.length() == LENGTH);
  uint32_t total = 0;
  for (int i = str.length() - 1; i >= 0; i--) {
    total += base * chtoi(str[i]);
    base *= 4;
  }
  return total; 
}

static inline void add(vector<vector<size_t>> &m, const string &kmer, long pos) {
  uint32_t index = ktoi(kmer);
  m[index].emplace_back(pos);
}

static inline void add(unordered_set<size_t> &s, vector<size_t> &v) {
  for (size_t i = 0; i < v.size(); i++) {
    s.insert(v[i]);
  }
}

void process_chromosome(string ref, vector<vector<size_t>> &m, long pos) {
  string kmer = "";

  for (size_t i = 0; i < ref.length() - LENGTH; i += LENGTH) {
    add(m, ref.substr(i, LENGTH), pos);
    pos += LENGTH;
  }
}

long get_job(ifstream &f, string &ref) {
  string line;
  long size = 0;
  while (getline(f, line)) {
    if (line[0] == '>') {
      cout << line << endl;
      break;
    }
    if (line.find('N') != string::npos) {
      size += line.length();
    } else ref.append(line);
  }
  return ref.size();

}



void map_reads(char *fastqname, vector<vector<size_t>> &m) {

  ifstream f;
  f.open(fastqname);
  string line;
  long counter = 0;
  double average = 0.0;
  long unmapped = 0;
  while (getline(f, line)) {
    if (counter % 100000 == 0) {
      cout << "Examined " << counter << " lines" << endl;
    }
    if (counter % 4 != 1) {
      counter++;
      continue;
    }
    
    unordered_set<size_t> locations;
    for (size_t i = 0; i < line.length() - LENGTH; i++) {
      string kmer = line.substr(i, LENGTH); 
      uint32_t index = ktoi(kmer);
      add(locations, m[index]);
    }
    average += locations.size();
    if (locations.size() == 0) {
      unmapped++;
    }
    counter++;
  }

  average /= (counter / 4);
  cout << "Average # of mapped locations: " << average << endl;
  cout << "Unmapped reads: " << unmapped << endl;
}

void distribute_work(char *fname, vector<vector<size_t>> &m) {

  ifstream f;
  f.open(fname);

  string line;
  // 1 index
  long pos = 1;
  getline(f, line); // remove first header
  cout << line << endl;

  while (!f.eof()) {
    string ref;
    ref.reserve(10000000);
    long offset = get_job(f, ref);
    process_chromosome(ref, m, pos);
    pos += offset;
  }
  f.close();

  cout << "Size of map: " << m.size() << endl;
  double average = 0.0;
  long count = 0;
  for (auto iter = m.begin(); iter != m.end(); ++iter) {
    average += iter->size();
    if (!iter->empty()) count++;
  }
  average /= count;
  cout << "Average # of elements: " << average << endl;
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  cout << "Memory usage " << r_usage.ru_maxrss * 1024 << endl;
}


int main(int argc, char *argv[]) {
  ios::sync_with_stdio(false);

  if (argc < 3) {
    cout << "Usage: align [fastq file] [reference genome path]" << endl;
    exit(1);
  }

  vector<vector<size_t>> m;
  m.resize(SPACE, vector<size_t>(0));
  
  cout << "Preallocated vector with capacity: " << m.capacity() << endl;
  distribute_work(argv[2], m);
  //map_reads(argv[1], m);

  return 0;
}
