#include "KMCSimulation.h"

#include <chrono>
#include <utility>
#include <boost/serialization/vector.hpp>
#include "BarrierPredictor.h"

namespace kmc {

constexpr size_t kEventListSize = Al_const::kNumFirstNearestNeighbors;
KMCSimulation::KMCSimulation(
    cfg::Config config,
    unsigned long long int log_dump_steps,
    unsigned long long int config_dump_steps,
    unsigned long long int maximum_number,
    std::unordered_map<std::string, double> type_category_hashmap,
    unsigned long long int steps,
    double energy,
    double time)
    : config_(std::move(config)),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_number_(maximum_number),
      steps_(steps),
      energy_(energy),
      time_(time),
      vacancy_index_(cfg::GetVacancyIndex(config_)),
      mpi_rank_(static_cast<size_t>(world_.rank())),
      mpi_size_(static_cast<size_t>(world_.size())),
      barrier_predictor_(config_, std::move(type_category_hashmap)),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())) {
  event_list_.reserve(kEventListSize);
  if (world_.rank() == 0) {
    std::cout << "Using " << world_.size() << " processes." << std::endl;
  }
}
KMCSimulation::~KMCSimulation() = default;
void KMCSimulation::BuildEventListSerial() {
  event_list_.clear();
  total_rate_ = 0;
  for (const auto neighbor_index :
      config_.GetAtomList()[vacancy_index_].GetFirstNearestNeighborsList()) {
    std::pair<size_t, size_t> jump_pair(vacancy_index_, neighbor_index);
    KMCEvent event(jump_pair, barrier_predictor_.GetBarrierAndDiff(config_, jump_pair));

    total_rate_ += event.GetRate();
    event_list_.push_back(std::move(event));
  }
}
void KMCSimulation::BuildEventListParallel() {
  event_list_.clear();
  total_rate_ = 0;

  auto quotient = kEventListSize / mpi_size_;
  auto remainder = kEventListSize % mpi_size_;
  // quotient part
  for (size_t j = 0; j < quotient; ++j) {
    auto i = j * mpi_size_ + mpi_rank_;
    auto neighbor_index = config_.GetAtomList()[vacancy_index_].GetFirstNearestNeighborsList()[i];
    std::pair<size_t, size_t> jump_pair(vacancy_index_, neighbor_index);
    KMCEvent event(jump_pair, barrier_predictor_.GetBarrierAndDiff(config_, jump_pair));

    if (mpi_rank_ == 0) {
      std::vector<KMCEvent> collected_list;
      double sum_rates_;

      boost::mpi::reduce(world_, event.GetRate(), sum_rates_, std::plus<>(), 0);
      boost::mpi::gather(world_, event, collected_list, 0);

      total_rate_ += sum_rates_;
      std::copy(collected_list.begin(), collected_list.end(), std::back_inserter(event_list_));
    } else {
      boost::mpi::reduce(world_, event.GetRate(), std::plus<>(), 0);
      boost::mpi::gather(world_, event, 0);
    }
  }
  // remainder part
  if (remainder) {
    auto i = quotient * mpi_size_ + mpi_rank_;
    auto neighbor_index = config_.GetAtomList()[vacancy_index_].GetFirstNearestNeighborsList()[i];
    std::pair<size_t, size_t> jump_pair(vacancy_index_, neighbor_index);
    KMCEvent event(jump_pair, barrier_predictor_.GetBarrierAndDiff(config_, jump_pair));

    if (mpi_rank_ == 0) {
      total_rate_ += event.GetRate();
      event_list_.push_back(event);
      for (int j = 1; j < static_cast<int>(remainder); ++j) {
        world_.recv(j, 0, event);
        total_rate_ += event.GetRate();
        event_list_.push_back(event);
      }
    } else {
      if (i < kEventListSize)
        world_.send(0, 0, event);
    }
  }

  boost::mpi::broadcast(world_, event_list_, 0);
  boost::mpi::broadcast(world_, total_rate_, 0);

  /* calculate relative and cumulative probability */
  double cumulative_provability = 0.0;
  for (auto &event : event_list_) {
    event.CalculateProbability(total_rate_);
    cumulative_provability += event.GetProbability();
    event.SetCumulativeProvability(cumulative_provability);
  }
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
void KMCSimulation::CheckAndFix([[maybe_unused]] double one_step_time) {}
void KMCSimulation::Simulate() {
  std::ofstream ofs("kmc_log.txt", std::ofstream::out | std::ofstream::app);
  if (mpi_rank_ == 0) {
    ofs << "steps    time    energy    barrier\n";
    ofs.precision(8);
  }

  while (steps_ < maximum_number_) {
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
    boost::mpi::broadcast(world_, event_index, 0);
    // world_.barrier();
    const auto &executed_invent_pair = event_list_[event_index].GetJumpPair();
    cfg::AtomsJump(config_, executed_invent_pair);

    // std::pair<size_t,size_t> CheckingPair{0,0};

    if (mpi_rank_ == 0) {
      // update time and energy
      static std::uniform_real_distribution<double> distribution(0.0, 1.0);
      auto one_step_time = log(distribution(generator_)) / (total_rate_ * 1e14);
      time_ -= one_step_time;
      energy_ += event_list_[event_index].GetEnergyChange();

      CheckAndFix(-one_step_time);

      // log and config file
      if (steps_ % log_dump_steps_ == 0) {
        ofs << steps_ << " " << time_ << " " << energy_ << " "
            << event_list_[event_index].GetBarrier() << std::endl;
      }
      if (steps_ % config_dump_steps_ == 0)
        cfg::Config::WriteConfig(config_, std::to_string(steps_) + ".cfg", true);
    }
    steps_++;
    // world_.barrier();
  }

}

} // namespace kmc