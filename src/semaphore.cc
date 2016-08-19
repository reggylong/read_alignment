

#include "semaphore.h"
#include <mutex>
#include <condition_variable>


using namespace std;

semaphore::semaphore(long value) : value(value) {}

void semaphore::wait() {
  unique_lock<mutex> lock(m);
  while (value == 0) {
    cv.wait(lock);
  }
  value--;
}

void semaphore::signal() {
  unique_lock<mutex> lock(m);
  value++;
  if (value == 1) cv.notify_all();
}
