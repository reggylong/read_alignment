
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>
#include <sys/time.h>
#include <sys/resource.h>
#include <thread>
#include <mutex>
#include "semaphore.h"

using namespace std;

static const size_t LENGTH = 11;
static const long THREADS = 6;

static inline void add(map<string, vector<size_t> > &m, const string &kmer, long pos) {
  if (m.find(kmer) == m.end()) {
    vector<size_t> v;
    v.emplace_back(pos);
    m[kmer] = v;
  } else {
    m[kmer].emplace_back(pos);
  }
}


void process_chromosome(string ref, map<string, vector<size_t> > &m, long pos) {
  s.signal(on_thread_exit);
  string kmer = "";

  for (size_t i = 0; i < ref.length() - LENGTH; i += LENGTH) {
    add(m, ref.substr(i, LENGTH), pos);
    pos += LENGTH;
  }
}

long get_job(ifstream &f, string &ref) {
  string line;
  while (getline(f, line)) {
    if (line[0] == '>') {
      cout << line << endl;
      break;
    }
    if (ref.find('N') != string::npos) continue;
    ref.append(line);
  }
  return ref.size();

}

void distribute_work(char *fname, map<string, vector<size_t> > &m) {

  ifstream f;
  f.open(fname);

  string line;
  // 1 index
  long pos = 1;
  long prev_pos = 1;
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
  for (auto& kv : m) {
    cout << kv.first << endl;
  }
}


int main(int argc, char *argv[]) {
  ios::sync_with_stdio(false);

  if (argc < 2) {
    cout << "Usage: align [reference genome path]" << endl;
    exit(1);
  }

  map<string, vector<size_t> > m;
  distribute_work(argv[1], m);

  return 0;
}
