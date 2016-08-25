
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
#include <set>
#include <cassert>
#include <cmath>
#include <list>
#include <algorithm>
#include "semaphore.h"

#define DEBUG false 

using namespace std;

static const size_t LENGTH = 14;
static const size_t SPACE = pow(4, LENGTH);
static const long THREADS = 6;

static __attribute__ ((noinline)) uint32_t chtoi(char c) {
  // A & 3 = 01
  // C & 3 = 11
  // G & 3 = 11
  // T & 3 = 10
  if (c == 'G') return 2;
  return (uint32_t) c & 3;
}

static __attribute__ ((noinline)) bool ktoi(const string &read, const long pos, uint32_t &index) {
  uint32_t base = 1;
  uint32_t total = 0;
  if (DEBUG) assert(pos + LENGTH - 1 < read.length());
  for (int i = pos + LENGTH - 1; i >= pos; i--) {
    if (toupper(read[i]) == 'N') return false;
    total += base * chtoi(toupper(read[i]));
    base *= 4;
  }
  index = total;
  return true;
}

static __attribute__ ((noinline)) bool add(vector<vector<size_t>> &m, const string &read, long pos) {
  uint32_t index = 0;
  if (ktoi(read, pos, index)) {
    m[index].emplace_back(pos);
    return true;
  }
  return false;
}

// Importantly, we do not make distinctions between chromosomes
void process_chromosome(string &ref, vector<vector<size_t>> &m, long pos) {
  string kmer = "";

  for (size_t i = pos; i < ref.length() - LENGTH; i += LENGTH) {
    add(m, ref, pos);
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

  static vector<size_t> locations;
  locations.reserve(10000);

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
    
    for (size_t i = 0; i < line.length() - LENGTH; i++) {
      uint32_t index = 0;
      if (ktoi(line, i, index)) locations.insert(locations.end(), m[index].begin(), m[index].end());
    }
    sort(locations.begin(), locations.end());
    locations.erase(unique(locations.begin(), locations.end()), locations.end());
    average += locations.size();
    if (locations.size() == 0) {
      unmapped++;
    }
    counter++;
    locations.clear();
  }

  average /= (counter / 4);
  cout << "Average # of mapped locations: " << average << endl;
  cout << "Unmapped reads: " << unmapped << endl;
}

void distribute_work(char *fname, vector<vector<size_t>> &m) {

  ifstream f;
  f.open(fname);

  string line;
  // 0 index
  long pos = 0;
  getline(f, line); // remove first header
  cout << line << endl;

  string ref;
  ref.reserve(500000000);
  while (!f.eof()) {
    long offset = get_job(f, ref);
    process_chromosome(ref, m, pos);
    pos += offset;
    ref.clear();
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
  map_reads(argv[1], m);

  return 0;
}
