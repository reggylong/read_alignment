
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>
#include <sys/time.h>
#include <sys/resource.h>
#include <condition_variable>
#include <mutex>
#include <thread>

#include "semaphore.h"

using namespace std;

static const size_t LENGTH = 11;
static const long THREADS = 4;
semaphore s(THREADS);
mutex write_m;

static inline void add(map<string, vector<size_t> > &m, const string &kmer, long pos) {
  if (m.find(kmer) == m.end()) {
    vector<size_t> v;
    v.emplace_back(pos);
    m[kmer] = v;
  } else {
    m[kmer].emplace_back(pos);
  }
}


void read(stringstream &ss, map<string, vector<size_t> > &m, long pos) {
  char c;
  string kmer = "";
  while (ss.get(c)) {
    if (kmer.length() < LENGTH) {
      kmer += c;
      if (kmer.length() == LENGTH) {
        add(m, kmer, pos - LENGTH);
      } else {
        kmer.erase(0, 1);
        kmer += c;
        add(m, kmer, pos - LENGTH);
      }
    }
    pos++;
  }
  s.signal();
}

void get_single_job(ifstream &f, stringstream &ss, long &pos) {
  string line;
  while (getline(f, line)) {
    if (line[0] == '<') {
      break;
    }
    ss << line;
    pos += line.length();
  }
}

void distribute_work(char *fname, map<string, vector<size_t> > &m) {
  ifstream f;
  f.open(fname);
  string kmer = "";

  vector<thread> threads;
  string line;
  getline(f, line);
  // 1 index
  long pos = 1;
  long prev_pos = 1;
  while (!f.eof()) {
    stringstream ss;
    get_single_job(f, ss, pos);
    s.wait();
    threads.emplace_back(thread(read, m, ss, prev_pos));
    prev_pos = pos;
  }
  f.close();
  for (size_t i = 0; i < threads.size(); i++) {
    threads[i].join();
  }
  cout << "Size of map: " << m.size() << endl;
  //size_t usage = (11 + sizeof(size_t) + sizeof(_Rb_tree_node_base)) * m.size();
  //cout << "(Under) Estimated Memory Usage: " << usage << endl;
  int who = RUSAGE_SELF;
  struct rusage usage;
  int ret;
  ret = getrusage(who, &usage);
  cout << "Max memory usage: " << usage.ru_maxrss * 1024 << " bytes" << endl;;
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
