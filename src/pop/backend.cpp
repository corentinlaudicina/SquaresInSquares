#include "pop/backend.hpp"

#include <algorithm>

namespace pack {

ThreadPool::ThreadPool(int threads) {
    for (int w = 1; w < std::max(1, threads); ++w) workers_.emplace_back([this, w] { loop(w); });
}

ThreadPool::~ThreadPool() {
    {
        std::lock_guard<std::mutex> lk(mu_);
        stop_ = true;
    }
    wake_.notify_all();
    for (std::thread& t : workers_) t.join();
}

void ThreadPool::loop(int worker) {
    long seen = 0;
    for (;;) {
        std::unique_lock<std::mutex> lk(mu_);
        wake_.wait(lk, [&] { return stop_ || generation_ != seen; });
        if (stop_) return;
        seen = generation_;
        ++running_;
        while (next_ < n_) {
            const int i = next_++;
            lk.unlock();
            (*job_)(i, worker);
            lk.lock();
        }
        if (--running_ == 0) done_.notify_all();
    }
}

void ThreadPool::parallel_for(int n, const std::function<void(int, int)>& f) {
    std::unique_lock<std::mutex> lk(mu_);
    job_ = &f;
    n_ = n;
    next_ = 0;
    ++generation_;
    wake_.notify_all();
    ++running_;
    while (next_ < n_) {  // the caller is worker 0
        const int i = next_++;
        lk.unlock();
        f(i, 0);
        lk.lock();
    }
    --running_;
    done_.wait(lk, [&] { return running_ == 0; });
    job_ = nullptr;
}

CpuBackend::CpuBackend(int threads)
    : pool_(threads > 0 ? threads : static_cast<int>(std::max(1u, std::thread::hardware_concurrency()))),
      ws_(static_cast<size_t>(pool_.size())) {}

void CpuBackend::relax_batch(std::vector<Packing>& batch, const CappedParams& prm, std::vector<CappedStats>& stats) {
    stats.assign(batch.size(), CappedStats{});
    pool_.parallel_for(static_cast<int>(batch.size()), [&](int i, int w) {
        stats[static_cast<size_t>(i)] = relax_capped(batch[static_cast<size_t>(i)], prm, ws_[static_cast<size_t>(w)]);
    });
}

void CpuBackend::cost_batch(const std::vector<Packing>& batch, Penalty pen, std::vector<CostInfo>& out) {
    out.assign(batch.size(), CostInfo{});
    pool_.parallel_for(static_cast<int>(batch.size()), [&](int i, int w) {
        out[static_cast<size_t>(i)] = cost_grad(batch[static_cast<size_t>(i)], pen, ws_[static_cast<size_t>(w)], nullptr);
    });
}

Packing to_packing(const System& sys) {
    Packing p;
    p.h = sys.L;
    p.xyt.reserve(3 * sys.particles.size());
    for (const Particle& q : sys.particles) {
        p.xyt.push_back(q.x());
        p.xyt.push_back(q.y());
        p.xyt.push_back(q.theta());
    }
    return p;
}

System to_system(const Packing& p) {
    System s;
    s.species.push_back(make_square(1.0));
    for (int i = 0; i < p.n(); ++i) {
        const size_t k = 3 * static_cast<size_t>(i);
        s.add_particle(p.xyt[k], p.xyt[k + 1], p.xyt[k + 2], 0);
    }
    s.L = p.h;
    s.recount();
    return s;
}

}  // namespace pack
