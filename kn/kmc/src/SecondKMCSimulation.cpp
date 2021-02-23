#include "SecondKMCSimulation.h"

#include <utility>
#include <chrono>
#include <utility>
namespace kmc {

constexpr size_t kFirstEventListSize = Al_const::kNumFirstNearestNeighbors;
constexpr size_t kSecondEventListSize = 11;

SecondKMCSimulation::SecondKMCSimulation(cfg::Config config,
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
                                   config_, type_set, lru_size),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())) {
  event_list_.resize(kFirstEventListSize * kSecondEventListSize);

  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank_);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size != kFirstEventListSize * kSecondEventListSize) {
    std::cout << "Must use " << kFirstEventListSize * kSecondEventListSize
              << " precesses. Terminating.\n";
    MPI_Finalize();
    exit(0);
  }
  MPI_Comm_group(MPI_COMM_WORLD, &world_group_);

  const int first_rank[12] =
      {0 * kSecondEventListSize, 1 * kSecondEventListSize, 2 * kSecondEventListSize,
       3 * kSecondEventListSize, 4 * kSecondEventListSize, 5 * kSecondEventListSize,
       6 * kSecondEventListSize, 7 * kSecondEventListSize, 8 * kSecondEventListSize,
       9 * kSecondEventListSize, 10 * kSecondEventListSize, 11 * kSecondEventListSize};
  MPI_Group_incl(world_group_, 12, first_rank, &first_group_);
  // MPI_Group_excl(world_group, 12, first_rank, &nonfirst_group);
  // if they are in first_rank
  MPI_Comm_create(MPI_COMM_WORLD, first_group_, &first_comm_);

  if (MPI_COMM_NULL != first_comm_) {
    MPI_Comm_rank(first_comm_, &first_group_rank_);
  }

  int color = world_rank_ / static_cast<int>(kSecondEventListSize);
  MPI_Comm_split(MPI_COMM_WORLD, color, world_rank_, &second_comm_);
  MPI_Comm_group(second_comm_, &second_group_);
  // MPI_Group_rank(second_group, &second_group_rank_);
  MPI_Comm_rank(second_comm_, &second_group_rank_);
  if (world_rank_ == 0) {
    std::cout << "Using " << world_size << " processes." << std::endl;
  }
}
SecondKMCSimulation::~SecondKMCSimulation() {
  if (MPI_GROUP_NULL != first_group_) MPI_Group_free(&first_group_);
  if (MPI_GROUP_NULL != second_group_) MPI_Group_free(&second_group_);
  if (MPI_COMM_NULL != first_comm_) MPI_Comm_free(&first_comm_);
  if (MPI_COMM_NULL != second_comm_) MPI_Comm_free(&second_comm_);
  MPI_Finalize();
}

// get the energy change and the probability from j to k, pjk by the reference.
// And return the index of j and k. Only applied to 12 sub-primary processes.
// And then pass to others, now we will have 12 different numbers for 12 different second groups.
std::pair<size_t, size_t> SecondKMCSimulation::BuildProbabilityListParallel(
    double &first_probability, double &first_energy_change) {
  std::pair<size_t, size_t> first_jump_pair{0, 0};
  if (first_comm_ != MPI_COMM_NULL) {
    first_total_rate_ = 0;

    const auto first_neighbor_index =
        config_.GetAtomList()[vacancy_index_].GetFirstNearestNeighborsList()[static_cast<size_t>(first_group_rank_)];
    first_jump_pair = {vacancy_index_, first_neighbor_index};
    KMCEvent event(first_jump_pair, lru_cache_barrier_predictor_.GetBarrierAndDiff(config_, first_jump_pair));

    const double this_rate = event.GetRate();
    MPI_Allreduce(&this_rate, &first_total_rate_, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    first_probability = this_rate / first_total_rate_;
    first_energy_change = event.GetEnergyChange();
    // MPI_Allgather(&first_probability_, 1, MPI_DOUBLE, probability_list_.data(), 1, MPI_DOUBLE, first_comm_);
  }

  MPI_Bcast(&first_jump_pair, sizeof(std::pair<size_t, size_t>), MPI_BYTE, 0, second_comm_);
  MPI_Bcast(&first_probability, 1, MPI_DOUBLE, 0, second_comm_);
  MPI_Bcast(&first_energy_change, 1, MPI_DOUBLE, 0, second_comm_);

  return first_jump_pair;
}

// Return the indexed of the corresponding second neighbors
std::vector<size_t> SecondKMCSimulation::GetSecondNeighborsIndexes() {
  std::vector<size_t> second_neighbors_indexes;
  second_neighbors_indexes.reserve(kSecondEventListSize);
  if (first_comm_ != MPI_COMM_NULL) {
    const auto first_neighbor_index =
        config_.GetAtomList()[vacancy_index_].
            GetFirstNearestNeighborsList()[static_cast<size_t>(first_group_rank_)];
    for (const auto second_neighbor_index
        : config_.GetAtomList()[first_neighbor_index].GetFirstNearestNeighborsList()) {
      if (second_neighbor_index != vacancy_index_) {
        second_neighbors_indexes.push_back(second_neighbor_index);
      }
    }
  } else{
    second_neighbors_indexes.resize(kSecondEventListSize);
  }
  MPI_Bcast(static_cast<void *>(second_neighbors_indexes.data()),
            kSecondEventListSize,
            MPI_UNSIGNED_LONG,
            0,
            second_comm_);
// #ifdef NDEBUG
//   std::cout << "size :" << second_neighbors_indexes.size() << '\n';
// #endif
  return second_neighbors_indexes;
}
void SecondKMCSimulation::BuildEventListParallel() {
  double first_probability, first_energy_change;

  auto first_jump_pair = BuildProbabilityListParallel(first_probability, first_energy_change);
  auto second_neighbors_indexes = GetSecondNeighborsIndexes();


  cfg::AtomsJump(config_, first_jump_pair);
  second_total_rate_ = 0.0;
  const auto
      second_neighbor_index = second_neighbors_indexes[static_cast<size_t>(second_group_rank_)];
  std::pair<size_t, size_t> second_jump_pair(vacancy_index_, second_neighbor_index);
  KMCEvent second_event
      (second_jump_pair,
       lru_cache_barrier_predictor_.GetBarrierAndDiff(config_, second_jump_pair));

  second_event.SetRate(second_event.GetRate() * first_probability);
  second_event.SetEnergyChange(second_event.GetEnergyChange() + first_energy_change);
  const double this_rate = second_event.GetRate();
  MPI_Allreduce(&this_rate, &second_total_rate_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allgather(&second_event,
                sizeof(KMCEvent),
                MPI_BYTE,
                event_list_.data(),
                sizeof(KMCEvent),
                MPI_BYTE,
                MPI_COMM_WORLD);
  /* calculate relative and cumulative probability */
  double cumulative_provability = 0.0;
  for (auto &event : event_list_) {
    event.CalculateProbability(second_total_rate_);
    cumulative_provability += event.GetProbability();
    event.SetCumulativeProvability(cumulative_provability);
  }
  cfg::AtomsJump(config_, first_jump_pair);
}

size_t SecondKMCSimulation::SelectEvent() const {
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
void SecondKMCSimulation::Simulate() {
  std::ofstream ofs("kmc_log.txt", std::ofstream::out | std::ofstream::app);
  if (world_rank_ == 0) {
    ofs << "steps\ttime\tenergy\tEa\tdE\n";
    ofs.precision(8);
  }

  while (steps_ < maximum_number_) {
    // log and config file
    if (world_rank_ == 0) {
      if (steps_ % log_dump_steps_ == 0) {
        ofs << steps_ << '\t' << time_ << '\t' << energy_
            << '\t' << one_step_barrier_ << '\t' << one_step_change_ << std::endl;
      }
      if (steps_ % config_dump_steps_ == 0) {
        cfg::Config::WriteConfig(config_, std::to_string(steps_) + ".cfg", true);
      }
    }
    // world_.barrier();
    MPI_Barrier(MPI_COMM_WORLD);

    BuildEventListParallel();

    size_t event_index;
    if (world_rank_ == 0) {
      event_index = SelectEvent();
// #ifdef NDEBUG
//       std::cout << "event choose " << event_index << '\n';
// #endif
    }
    // world_.barrier();
    MPI_Bcast(&event_index, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    // boost::mpi::broadcast(world_, event_index, 0);
    // world_.barrier();
    std::pair<size_t, size_t> first_jump_pair{vacancy_index_,
                                              config_.GetAtomList()[vacancy_index_].GetFirstNearestNeighborsList()[
                                                  event_index / kSecondEventListSize]};
    // std::cout << first_jump_pair.first << ' ' << first_jump_pair.second << '\n';
    // std::cout <<"First " << first_jump_pair.first << ' ' << first_jump_pair.second << '\n';
    cfg::AtomsJump(config_, first_jump_pair);

    const auto &executed_invent = event_list_[event_index];
    // std::cout <<"Second " << executed_invent.GetJumpPair().first << ' ' << executed_invent.GetJumpPair().second << '\n';
    cfg::AtomsJump(config_, executed_invent.GetJumpPair());

    // std::pair<size_t,size_t> CheckingPair{0,0};

    if (world_rank_ == 0) {
      // update time and energy
      static std::uniform_real_distribution<double> distribution(0.0, 1.0);
      constexpr double kPrefactor = 1e14;
      auto one_step_time = log(distribution(generator_)) / (second_total_rate_ * kPrefactor);
      time_ -= one_step_time;
      one_step_change_ = executed_invent.GetEnergyChange();
      energy_ += one_step_change_;
      one_step_barrier_ = executed_invent.GetBarrier();
    }
    ++steps_;
    // world_.barrier();
  }

}

}


// namespace kmc