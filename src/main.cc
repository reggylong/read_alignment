
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
#define VERBOSE false

using namespace std;

static const size_t LENGTH = 14;
static const size_t LOCATIONS = 5;
static const size_t SPACE = pow(4, LENGTH);
static const long THREADS = 6;
string complete_ref;

char reverse_bp(char c) {
  switch (toupper(c)) {
    case 'A':
      return 'T';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'T':
      return 'A';
    default:
      return 'N';
  }
}

void reverse_read(const string &read, string &reversed) {
  for (int_fast32_t i = read.length() - 1; i >= 0; i--) {
    reversed += reverse_bp(read[i]);
  }
}

static uint_fast32_t chtoi(char c) {
  // A & 3 = 01
  // C & 3 = 11
  // G & 3 = 11
  // T & 3 = 10
  if (c == 'G') return 2;
  return (uint_fast32_t) c & 3;
}

static bool ktoi(const string &read, const long pos, uint_fast32_t &index) {
  uint_fast32_t base = 1;
  uint_fast32_t total = 0;
  for (int i = pos + LENGTH - 1; i >= pos; i--) {
    if (toupper(read[i]) == 'N') return false;
    total += base * chtoi(toupper(read[i]));
    base *= 4;
  }
  index = total;
  return true;
}

static bool add(vector<vector<uint_fast32_t>> &m, const string &read, size_t pos) {
  if (pos == 33491329) {
    cout << "HELLO WORLD" << endl;
  }
  uint_fast32_t index = 0;
  if (ktoi(read, pos, index)) {
    m[index].emplace_back(pos);
    return true;
  }
  return false;
}

// Importantly, we do not make distinctions between chromosomes
void process_reference(vector<vector<uint_fast32_t>> &m) {
  string kmer = "";

  for (size_t i = 0; i < complete_ref.length() - LENGTH; i += LENGTH) {
    add(m, complete_ref, i);
  }
}



uint_fast32_t hamming(const string &read, uint_fast32_t pos) {
  uint_fast32_t dist = 0;
  for (size_t i = 0; i < read.length(); i++) {
    if (read[i] != complete_ref[pos + i]) dist++;
  }
  return dist;
}

void compute_locations(const string &read, vector<vector<uint_fast32_t>> &m, vector<uint_fast32_t> &locations) {
  size_t i = 0;
  uint_fast32_t index = 0;
  while (i < read.size() && (i < LENGTH || locations.size() < LOCATIONS)) {
    if (ktoi(read, i, index)) {
      for (uint_fast32_t j = 0; j < m[index].size(); j++) {
        locations.push_back(m[index][j] - i);
      }
    }
    i++;
  }
}

uint_fast32_t argmin(vector<uint_fast32_t> &v) {
  uint_fast32_t i = 0;
  uint_fast32_t minimum = UINT_FAST32_MAX;
  for (uint_fast32_t j = 0; j < v.size(); j++) {
    if (v[j] < minimum) {
      i = j;
      minimum = v[j];
    }
  }
  return i;
}

void map_reads(char *fastqname, vector<vector<uint_fast32_t>> &m) {

  vector<uint_fast32_t> locations;
  vector<uint_fast32_t> reversed_locations;

  ifstream f;
  f.open(fastqname);
  string line;
  long counter = 0;
  double average = 0.0;
  double dist = 0.0;
  long unmapped = 0;
  while (getline(f, line)) {
    if (counter % 100000 == 0) {
      cout << "Examined " << counter << " lines" << endl;
    }
    if (counter % 4 != 1) {
      counter++;
      continue;
    }

    if (VERBOSE) {
      cout << "[map_reads] Examining line: " << line << endl;
    }
    compute_locations(line, m, locations);
    string reversed_read;
    reverse_read(line, reversed_read);

    compute_locations(reversed_read, m, reversed_locations);
    locations.insert(locations.end(), reversed_locations.begin(), reversed_locations.end());

    sort(locations.begin(), locations.end());
    locations.erase(unique(locations.begin(), locations.end()), locations.end());
    average += locations.size();
    if (locations.size() == 0) {
      unmapped++;
    } else {
      vector<uint_fast32_t> dists;
      for (uint_fast32_t k = 0; k < locations.size(); k++) {
        dists.emplace_back(hamming(line, locations[k]));
      }
      uint_fast32_t index = argmin(dists);
      if (DEBUG) {
        for (uint_fast32_t i = 0; i < locations.size(); i++) {
          cout << complete_ref.substr(locations[i], 100) << " : " << locations[i] << endl;
        }
      }
      dist += dists[index];
      cout << "Index of min: " << locations[index] << endl;
    }
    counter++;


    locations.clear();
  }

  average /= (counter / 4);
  dist /= (counter / 4);
  cout << "Average hamming dist: " << dist << endl;
  cout << "Average # of mapped locations: " << average << endl;
  cout << "Unmapped reads: " << unmapped << endl;
}

void read_reference(char *fname, vector<vector<uint_fast32_t>> &m) {

  ifstream f;
  f.open(fname);

  string line;
  // 0 index
  getline(f, line); // remove first header
  cout << line << endl;

  while (getline(f, line)) {
    if (line[0] == '>') {
      cout << line << endl;
    } else {
      complete_ref.append(line);
    }
  }

  f.close();

}

void map_stats(vector<vector<uint_fast32_t>> &m) {
  cout << "Size of map: " << m.size() << endl;
  double average = 0.0;
  long count = 0;
  for (auto iter = m.begin(); iter != m.end(); ++iter) {
    average += iter->size();
    if (!iter->empty()) count++;
  }
  average /= count;
  cout << "Average # of elements per kmer: " << average << endl;
}

int main(int argc, char *argv[]) {
  ios::sync_with_stdio(false);

  if (argc < 3) {
    cout << "Usage: align [fastq file] [reference genome path]" << endl;
    exit(1);
  }

  complete_ref.reserve(4000000000);
  vector<vector<uint_fast32_t>> m;
  m.resize(SPACE, vector<uint_fast32_t>(0));
  
  cout << "Preallocated vector with capacity: " << m.capacity() << endl;
  read_reference(argv[2], m);
  cout << "Reference genome has length: " << complete_ref.size() << endl;

  process_reference(m);
  map_stats(m);

  map_reads(argv[1], m);

  cout << complete_ref.substr(33491329, 100) << endl;
  return 0;
}
