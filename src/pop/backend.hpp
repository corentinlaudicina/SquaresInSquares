#pragma once

// Batch relaxation behind an interface: packings in, relaxed packings out. The CPU backend runs one
// packing per task on a persistent thread pool with one Workspace per thread; every packing is
// relaxed independently, so the result does not depend on the thread count. A GPU backend, if one
// is ever approved, implements the same interface (float32 relax; legalization stays on the CPU).

#include "pop/relax_capped.hpp"
#include "system.hpp"

#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace pack {

class RelaxBackend {
public:
    virtual ~RelaxBackend() = default;
    virtual void relax_batch(std::vector<Packing>& batch, const CappedParams& prm,
                             std::vector<CappedStats>& stats) = 0;
    // Cost only (no gradient, no relaxation): the GA's pre-filter.
    virtual void cost_batch(const std::vector<Packing>& batch, Penalty pen, std::vector<CostInfo>& out) = 0;
};

// A fixed pool; parallel_for(n, f) calls f(index, worker) for index in [0, n) and returns when all ran.
class ThreadPool {
public:
    explicit ThreadPool(int threads);
    ~ThreadPool();
    int size() const { return static_cast<int>(workers_.size()) + 1; }  // the caller works too
    void parallel_for(int n, const std::function<void(int, int)>& f);

private:
    void loop(int worker);
    std::vector<std::thread> workers_;
    std::mutex mu_;
    std::condition_variable wake_, done_;
    const std::function<void(int, int)>* job_ = nullptr;
    int n_ = 0, next_ = 0, running_ = 0;
    long generation_ = 0;
    bool stop_ = false;
};

class CpuBackend : public RelaxBackend {
public:
    explicit CpuBackend(int threads = 0);  // <= 0: hardware concurrency
    void relax_batch(std::vector<Packing>& batch, const CappedParams& prm, std::vector<CappedStats>& stats) override;
    void cost_batch(const std::vector<Packing>& batch, Penalty pen, std::vector<CostInfo>& out) override;
    int threads() const { return pool_.size(); }
    void parallel_for(int n, const std::function<void(int, int)>& f) { pool_.parallel_for(n, f); }

private:
    ThreadPool pool_;
    std::vector<Workspace> ws_;
};

// Unit-square Systems <-> genomes.
Packing to_packing(const System& sys);
System to_system(const Packing& p);

}  // namespace pack
