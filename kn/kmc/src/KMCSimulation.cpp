#include "KMCSimulation.h"

#include <chrono>
#include <utility>
#include <mpi.h>

namespace kmc {

constexpr size_t kEventListSize = Al_const::kNumFirstNearestNeighbors;
KMCSimulation::KMCSimulation(cfg::Config config,
                             unsigned long long int log_dump_steps,
                             unsigned long long int config_dump_steps,
                             unsigned long long int maximum_number,
                             const std::set<std::string> &type_set,
                             unsigned long long int steps,
                             double energy,
                             double time,
                             const std::string &json_parameters_filename,
                             size_t lru_size)
    : config_(std::move(config)),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_number_(maximum_number),
      steps_(steps),
      energy_(energy),
      time_(time),
      vacancy_index_(cfg::GetVacancyIndex(config_)),
      lru_cache_barrier_predictor_(json_parameters_filename,
                                   config_, type_set),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())) {
  MPI_Init(nullptr, nullptr);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  mpi_rank_ = static_cast<size_t>(mpi_rank);
  mpi_size_ = static_cast<size_t>(mpi_size);

  event_list_.reserve(kEventListSize);
  if (mpi_rank_ == 0) {
    std::cout << "Using " << mpi_size_ << " processes." << std::endl;
  }
  if (mpi_size_ > 12) {
    std::cout << "Cannot use more than 12 precesses. Terminating.\n";
    MPI_Finalize();
    exit(0);
  }
}
KMCSimulation::~KMCSimulation() {
  MPI_Finalize();
}
void KMCSimulation::BuildEventListSerial() {
  event_list_.clear();
  total_rate_ = 0;
  for (const auto neighbor_index :
      config_.GetAtomList()[vacancy_index_].GetFirstNearestNeighborsList()) {
    std::pair<size_t, size_t> atom_id_jump_pair(vacancy_index_, neighbor_index);
    KMCEvent event(atom_id_jump_pair, lru_cache_barrier_predictor_.GetBarrierAndDiff(config_, atom_id_jump_pair));

    total_rate_ += event.GetForwardRate();
    event_list_.push_back(std::move(event));
  }
#ifndef NDEBUG
  for (const auto &event : event_list_) {
    std::cerr << event.GetForwardBarrier() << '\t';
  }
  std::cerr << '\n';
  for (const auto &event : event_list_) {
    std::cerr << event.GetForwardRate() << '\t';
  }
  std::cerr << '\n';
#endif
}
void KMCSimulation::BuildEventListParallel() {
  event_list_.clear();
  total_rate_ = 0;

  auto quotient = kEventListSize / mpi_size_;
  auto remainder = kEventListSize % mpi_size_;
  // quotient part
  for (size_t j = 0; j < quotient; ++j) {
    const auto i = j * mpi_size_ + mpi_rank_;
    const auto
        neighbor_index = config_.GetAtomList()[vacancy_index_].GetFirstNearestNeighborsList()[i];
    std::pair<size_t, size_t> atom_id_jump_pair(vacancy_index_, neighbor_index);
    KMCEvent event(atom_id_jump_pair, lru_cache_barrier_predictor_.GetBarrierAndDiff(config_, atom_id_jump_pair));

    const double this_rate = event.GetForwardRate();
    double sum_rates = 0;
    MPI_Reduce(&this_rate, &sum_rates, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mpi_rank_ == 0) {
      std::vector<KMCEvent> collected_event_ctor_pair_list(mpi_size_);
      MPI_Gather(&event, sizeof(KMCEvent), MPI_BYTE,
                 collected_event_ctor_pair_list.data(),
                 sizeof(KMCEvent), MPI_BYTE, 0, MPI_COMM_WORLD);

      total_rate_ += sum_rates;
      std::copy(collected_event_ctor_pair_list.begin(),
                collected_event_ctor_pair_list.end(),
                std::back_inserter(event_list_));
    } else {
      MPI_Gather(&event, sizeof(KMCEvent), MPI_BYTE,
                 nullptr, 0, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    // if (mpi_rank_ == 0) {
    //   std::vector<KMCEvent> collected_list;
    //   double sum_rates_;
    //
    //   boost::mpi::reduce(world_, event.GetRate(), sum_rates_, std::plus<>(), 0);
    //   boost::mpi::gather(world_, event, collected_list, 0);
    //
    //   total_rate_ += sum_rates_;
    //   std::copy(collected_list.begin(), collected_list.end(), std::back_inserter(event_list_));
    // } else {
    //   boost::mpi::reduce(world_, event.GetRate(), std::plus<>(), 0);
    //   boost::mpi::gather(world_, event, 0);
    // }
  }
  // remainder part
  if (remainder) {
    auto i = quotient * mpi_size_ + mpi_rank_;
    auto neighbor_index = config_.GetAtomList()[vacancy_index_].GetFirstNearestNeighborsList()[i];
    const std::pair<size_t, size_t> atom_id_jump_pair(vacancy_index_, neighbor_index);
    KMCEvent event(atom_id_jump_pair, lru_cache_barrier_predictor_.GetBarrierAndDiff(config_, atom_id_jump_pair));

    if (mpi_rank_ == 0) {
      total_rate_ += event.GetForwardRate();
      event_list_.push_back(event);
      for (size_t j = 1; j < remainder; ++j) {
        MPI_Recv(&event, sizeof(KMCEvent), MPI_BYTE,
                 static_cast<int>(j), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        total_rate_ += event.GetForwardRate();
        event_list_.push_back(event);
      }
    } else {
      if (i < kEventListSize) {
        MPI_Send(&event, sizeof(KMCEvent), MPI_BYTE,
                 0, 0, MPI_COMM_WORLD);
      }
    }
  }
  MPI_Bcast(event_list_.data(), sizeof(KMCEvent) * kEventListSize, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&total_rate_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // boost::mpi::broadcast(world_, event_list_, 0);
  // boost::mpi::broadcast(world_, total_rate_, 0);

  /* calculate relative and cumulative probability */
  double cumulative_probability = 0.0;
  for (auto &event : event_list_) {
    event.CalculateProbability(total_rate_);
    cumulative_probability += event.GetProbability();
    event.SetCumulativeProvability(cumulative_probability);
  }
// #ifndef NDEBUG
//   if (mpi_rank_ == 0) {
//     for (const auto &event : event_list_) {
//       std::cerr << event.GetForwardBarrier() << '\t';
//     }
//     std::cerr << '\n';
//     for (const auto &event : event_list_) {
//       std::cerr << event.GetForwardRate() << '\t';
//     }
//     std::cerr << '\n';
//   }
// #endif
}

size_t KMCSimulation::SelectEvent() const {
  static std::uniform_real_distribution<double> distribution(0.0, 1.0);

  const double random_number = distribution(generator_);
  auto it = std::lower_bound(event_list_.begin(),
                             event_list_.end(),
                             random_number,
                             [](const auto &lhs, double value) {
                               return lhs.GetCumulativeProvability() < value;
                             });
  // If not find (maybe generated 1), which rarely happens, returns the last event
  if (it == event_list_.cend()) {
    it--;
  }
  return static_cast<size_t>(std::distance(event_list_.begin(), it));
}
void KMCSimulation::Simulate() {
  std::ofstream ofs("kmc_log.txt", std::ofstream::out | std::ofstream::app);
  if (mpi_rank_ == 0) {
    ofs << "steps\ttime\tenergy\tEa\tdE\n";
    ofs.precision(8);
  }

  while (steps_ < maximum_number_) {
    // log and config file
    if (mpi_rank_ == 0) {
      if (steps_ % log_dump_steps_ == 0) {
        ofs << steps_ << '\t' << time_ << '\t' << energy_
            << '\t' << one_step_barrier_ << '\t' << one_step_change_ << std::endl;
      }
      if (steps_ % config_dump_steps_ == 0) {
        cfg::Config::WriteConfig(config_, std::to_string(steps_) + ".cfg", 2);
      }
    }
    // world_.barrier();
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_size_ == 1) {
      BuildEventListSerial();
    } else {
      BuildEventListParallel();
    }
    size_t event_index;
    if (mpi_rank_ == 0) {
      event_index = SelectEvent();
    }
    // world_.barrier();
    MPI_Bcast(&event_index, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    // boost::mpi::broadcast(world_, event_index, 0);
    // world_.barrier();
    const auto &executed_invent = event_list_[event_index];
    cfg::AtomsJump(config_, executed_invent.GetAtomIDJumpPair());

    // std::pair<size_t,size_t> CheckingPair{0,0};

    if (mpi_rank_ == 0) {
      // update time and energy
      static std::uniform_real_distribution<double> distribution(0.0, 1.0);
      constexpr double kPrefactor = 1e14;
      auto one_step_time = log(distribution(generator_)) / (total_rate_ * kPrefactor);
      time_ -= one_step_time;
      one_step_change_ = executed_invent.GetEnergyChange();
      energy_ += one_step_change_;
      one_step_barrier_ = executed_invent.GetForwardBarrier();
      // CheckTimeAndFix(-one_step_time);
    }
    ++steps_;
    // world_.barrier();
  }

}

} // namespace kmc