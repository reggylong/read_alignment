
#ifndef _SEMAPHORE_
#define _SEMAPHORE_

#include <mutex>
#include <condition_variable>

class semaphore {
  public:
    semaphore(long value = 0);
    void wait();
    void signal();
  
  private:
    long value;
    std::mutex m;
    std::condition_variable cv;

    semaphore(const semaphore& orig) = delete;
    const semaphore& operator=(const semaphore& rhs) const = delete;

};
#endif
