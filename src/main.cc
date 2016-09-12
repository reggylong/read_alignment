
#include <array>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
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
#define STATS true

using namespace std;

static const size_t LENGTH = 14;
static const size_t LOCATIONS = 100;
static const size_t MIN_COST = 5;
static const size_t CUTOFF = 100000;
static const size_t SPACE = pow(4, LENGTH);
static const long THREADS = 6;
vector<uint_fast32_t> counts(101, 0);
string complete_ref;

static uint_fast32_t edit_dist(string str1, string str2) {
  static vector<vector<uint_fast32_t>> matrix_edits(str1.length() + 1, vector<uint_fast32_t>(str2.length() + 1));

  for (uint32_t i = 0; i <= str1.length(); i++) {
    matrix_edits[i][0] = i;
  }
  for (uint32_t j = 0; j <= str2.length(); j++) {
    matrix_edits[0][j] = j;
  }
  matrix_edits[0][0] = 0;
  for (uint32_t i = 1; i <= str1.length(); i++) {
    for (uint32_t j = 1; j <= str2.length(); j++) {
      if (str1[i-1] == str2[j-1]) {
        matrix_edits[i][j] = matrix_edits[i-1][j-1];
      } else {
        matrix_edits[i][j] = min(matrix_edits[i-1][j-1] + 1, min(matrix_edits[i-1][j] + 1, matrix_edits[i][j-1] + 1));
      }
    }
  }
  return matrix_edits[str1.length()][str2.length()];
}

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
    if (pos + i >= complete_ref.length()) return dist + (read.length() - i);
    if (read[i] != complete_ref[pos + i]) dist++;
  }
  return dist;
}

void compute_locations(const string &read, vector<vector<uint_fast32_t>> &m, vector<uint_fast32_t> &locations) {
  size_t i = 0;
  uint_fast32_t index = 0;
  while (i < read.size()) {// && (i < LENGTH || locations.size() < LOCATIONS)) {
    if (ktoi(read, i, index) && m[index].size() < CUTOFF) {
      for (uint_fast32_t j = 0; j < m[index].size(); j++) {
        if (m[index][j] - i + read.length() >= complete_ref.size()) continue; 
        locations.push_back(m[index][j] - i);
      }
    }
    i++;
  }
}

int_fast32_t argmin(const vector<uint_fast32_t> &v) {
  int_fast32_t i = 0;
  uint_fast32_t minimum = UINT_FAST32_MAX;
  for (uint_fast32_t j = 0; j < v.size(); j++) {
    if (v[j] < minimum) {
      i = j;
      minimum = v[j];
    }
  }
  return i;
}

int_fast32_t compute_hamming(string &read, vector<uint_fast32_t> &locations, vector<uint_fast32_t> &dists) {
  if (locations.empty()) {
    return -1;
  }
  for (uint_fast32_t i = 0; i < locations.size(); i++) {
    dists.emplace_back(hamming(read, locations[i]));
  }
  int_fast32_t index = argmin(dists);
  if (VERBOSE) cout << read << endl;
  if (VERBOSE) cout << complete_ref.substr(locations[index], 100) << endl;
  return index;
}

void map_reads(char *fastqname, vector<vector<uint_fast32_t>> &m) {

  static vector<uint_fast32_t> locations;
  static vector<uint_fast32_t> reversed_locations;
  locations.reserve(10000);
  reversed_locations.reserve(10000);

  ifstream f;
  f.open(fastqname);
  string line;
  long counter = 0;
  double average = 0.0;
  double dist = 0.0;
  double edit = 0.0;

  long unmapped = 0;
  string name = "";
  while (getline(f, line)) {
    if (counter % 1000000 == 0) {
      cout << "Examined " << counter << " lines" << endl;
    }
    if (counter % 4 != 1) {
      if (counter % 4 == 0) name = line;
      counter++;
      continue;
    }

    if (VERBOSE) {
      cout << "[map_reads] Examining line: " << line << endl;
    }
    vector<uint_fast32_t> dists;
    vector<uint_fast32_t> reversed_dists;

    compute_locations(line, m, locations);
    string reversed_read;
    reverse_read(line, reversed_read);

    compute_locations(reversed_read, m, reversed_locations);

    sort(locations.begin(), locations.end());
    sort(reversed_locations.begin(), reversed_locations.end());

    locations.erase(unique(locations.begin(), locations.end()), locations.end());
    reversed_locations.erase(unique(reversed_locations.begin(), reversed_locations.end()), reversed_locations.end());

    average += locations.size();
    average += reversed_locations.size();

    if (locations.empty() && reversed_locations.empty()) {
      unmapped++;
      counts[line.length()]++;
    } else {
      int_fast32_t locations_index = compute_hamming(line, locations, dists);
      uint_fast32_t value = (locations_index >= 0) ? dists[locations_index] : UINT_FAST32_MAX;
      int_fast32_t reversed_index = compute_hamming(reversed_read, reversed_locations, reversed_dists);
      uint_fast32_t reversed_value = (reversed_index >= 0) ? reversed_dists[reversed_index] : UINT_FAST32_MAX; 

      uint_fast32_t h_dist = 0;
      uint_fast32_t e_dist = 0;
      // compute edit distance
      if (value < reversed_value) {
        if (STATS) e_dist = edit_dist(line, complete_ref.substr(locations[locations_index], line.length()));
        h_dist = value;
      } else {
        if (STATS) e_dist = edit_dist(reversed_read, complete_ref.substr(reversed_locations[reversed_index], reversed_read.length()));
        h_dist = reversed_value;
      }

      if (VERBOSE) cout << "[map_reads] " << name << " " << h_dist << endl;
      if (STATS && VERBOSE) cout << "[map_reads] " << name << " " << e_dist << endl;
      dist += h_dist;
      edit += e_dist;
      if (STATS) {
        counts[e_dist]++;
      } else {
        counts[h_dist]++;
      }
    }
    counter++;
    locations.clear();
    reversed_locations.clear();
  }

  average /= (counter / 4);
  dist /= (counter / 4);
  edit /= (counter / 4);
  cout << "Average hamming dist: " << dist << endl;
  if (STATS) cout << "Average edit dist: " << edit << endl;
  cout << "Average # of mapped locations: " << average << endl;
  cout << "Unmapped reads: " << unmapped << endl;
}

void read_reference(char *fname, vector<vector<uint_fast32_t>> &m) {

  ifstream f;
  f.open(fname);

  string line;
  // 0 index
  getline(f, line); // remove first header

  cout << "[read_reference] " << line << endl;
  while (getline(f, line)) {
    if (line[0] == '>') {
      cout << "[read_reference] " << line << endl;
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
    cout << "Usage: align [fastq file] [reference genome path] [optional: stats]" << endl;
    exit(1);
  }

  complete_ref.reserve(4000000000);
  vector<vector<uint_fast32_t>> m;
  m.resize(SPACE, vector<uint_fast32_t>(0));

  cout << "Preallocated vector with capacity: " << m.capacity() << endl;
  read_reference(argv[2], m);
  cout << "Reference genome has length: " << complete_ref.size() << endl;
  process_reference(m);
  cout << "Finished processing reference" << endl;
  map_stats(m);

  map_reads(argv[1], m);

  if (argc > 3) {
    ofstream f;
    f.open(argv[3]);
    double total = 0;
    for (auto n : counts) {
      total += n;
    }
    for (uint_fast32_t i = 0; i < counts.size(); i++) {
      if (counts[i] > 0) f << i << " " << counts[i] / total << " " << counts[i] << endl;
    }
    f.close();
  }

  return 0;
}
