#include "ChainKMCSimulation.h"

#include <utility>
#include <chrono>
#include <utility>
namespace kmc {

constexpr size_t kFirstEventListSize = Al_const::kNumFirstNearestNeighbors;
constexpr size_t kSecondEventListSize = 11;

ChainKMCSimulation::ChainKMCSimulation(cfg::Config config,
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
  event_list_.resize(kFirstEventListSize);

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
ChainKMCSimulation::~ChainKMCSimulation() {
  if (MPI_GROUP_NULL != first_group_) MPI_Group_free(&first_group_);
  if (MPI_GROUP_NULL != second_group_) MPI_Group_free(&second_group_);
  if (MPI_COMM_NULL != first_comm_) MPI_Comm_free(&first_comm_);
  if (MPI_COMM_NULL != second_comm_) MPI_Comm_free(&second_comm_);
  MPI_Finalize();
}

// get the energy change and the probability from j to k, pjk by the reference.
// And return the index of j and k. Only applied to 12 sub-primary processes.
// And then pass to others, now we will have 12 different numbers for 12 different second groups.
KMCEvent ChainKMCSimulation::GetEventJK() {
  KMCEvent event_j_k;
  if (first_comm_ != MPI_COMM_NULL) {
    first_total_rate_ = 0;

    const auto k_index =
        config_.GetAtomList()[vacancy_index_].GetFirstNearestNeighborsList()[static_cast<size_t>(first_group_rank_)];

    event_j_k = KMCEvent(
        {vacancy_index_, k_index},
        lru_cache_barrier_predictor_.GetBarrierAndDiff(config_, {vacancy_index_, k_index}));

    const double first_rate = event_j_k.GetForwardRate();
    MPI_Allreduce(&first_rate, &first_total_rate_, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    event_j_k.CalculateProbability(first_total_rate_);
// MPI_Allgather(&first_probability_, 1, MPI_DOUBLE, probability_list_.data(), 1, MPI_DOUBLE, first_comm_);
  }

  MPI_Bcast(&event_j_k, sizeof(KMCEvent), MPI_BYTE, 0, second_comm_);
  return event_j_k;
}

// Return the indexed of the corresponding second neighbors
std::vector<size_t> ChainKMCSimulation::GetLIndexes() {
  std::vector<size_t> l_indexes;
  l_indexes.reserve(kSecondEventListSize);
  if (first_comm_ != MPI_COMM_NULL) {
    const auto first_neighbor_index =
        config_.GetAtomList()[vacancy_index_].
            GetFirstNearestNeighborsList()[static_cast<size_t>(first_group_rank_)];
    for (const auto second_neighbor_index
        : config_.GetAtomList()[first_neighbor_index].GetFirstNearestNeighborsList()) {
      if (second_neighbor_index != vacancy_index_) {
        l_indexes.push_back(second_neighbor_index);
      }
    }
  } else {
    l_indexes.resize(kSecondEventListSize);
  }
  MPI_Bcast(static_cast<void *>(l_indexes.data()),
            kSecondEventListSize,
            MPI_UNSIGNED_LONG,
            0,
            second_comm_);
  return l_indexes;
}
std::pair<double, double> ChainKMCSimulation::BuildEventListParallel() {
  // double first_probability, first_energy_change, first_back_rate;
  auto event_j_k = GetEventJK();
  const auto l_indexes = GetLIndexes();
  const auto probability_j_k = event_j_k.GetProbability();
  cfg::AtomsJump(config_, event_j_k.GetJumpPair());
  total_rate_k_l_ = 0.0;
  const auto l_index = l_indexes[static_cast<size_t>(second_group_rank_)];
  KMCEvent event_k_l
      ({vacancy_index_, l_index},
       lru_cache_barrier_predictor_.GetBarrierAndDiff(config_, {vacancy_index_, l_index}));
  cfg::AtomsJump(config_, event_j_k.GetJumpPair());

  // get sum r_{k to l}
  const double r_k_l = event_k_l.GetForwardRate();
  MPI_Allreduce(&r_k_l, &total_rate_k_l_, 1, MPI_DOUBLE, MPI_SUM, second_comm_);
  const auto r_k_j = event_j_k.GetBackwardRate();
  total_rate_k_l_ += r_k_j;
  const double probability_k_j = r_k_j / total_rate_k_l_;

  // gamma_{j-k}
  double gamma_jk = probability_j_k * probability_k_j;

  // collect gamma_j(i), gamma_j(0), p_{j to i}, gamma_{ij}
  double gamma_j_i = 0.0, gamma_j_0 = 0.0, probability_j_i = 0.0, gamma_ji = 0.0;
  if (first_comm_ != MPI_COMM_NULL) {
    double gamma_ji_helper = 0, probability_j_i_helper = 0;
    if (event_j_k.GetJumpPair().second == previous_i) {
      gamma_ji_helper = probability_k_j * probability_j_k;
      probability_j_i_helper = probability_j_k;
    }
    MPI_Allreduce(&gamma_jk, &gamma_j_0, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    MPI_Allreduce(&gamma_ji_helper, &gamma_ji, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    MPI_Allreduce(&probability_j_i_helper, &probability_j_i, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    gamma_j_i = (gamma_j_0 - gamma_ji) / (1 - probability_j_i);

    int flip;
    if (first_group_rank_ == 0) {
      static std::uniform_real_distribution<double> distribution(0.0, 1.0);
      const double random_number = distribution(generator_);
      if (random_number <= gamma_j_i) {
        // indirect move
        flip = int(log(random_number / gamma_j_i) / log(gamma_j_0)) + 1;
      } else {
        // direct move
        flip = 0;
      }
    }
    MPI_Bcast(&flip, 1, MPI_DOUBLE, 0, first_comm_);
    if (flip == 0) {
      double direct_possibility = (event_j_k.GetJumpPair().second == previous_i) ? 0 :
                                  (probability_j_k - gamma_jk)
                                      / (1 - probability_j_i - gamma_j_0 + gamma_ji);
      event_j_k.SetProbability(direct_possibility);
    } else {
      double indirect_possibility = (probability_j_k - gamma_jk) / (1 - gamma_j_0);
      event_j_k.SetProbability(indirect_possibility);
    }

    MPI_Allgather(&event_j_k,
                  sizeof(KMCEvent),
                  MPI_BYTE,
                  nullptr,
                  0,
                  MPI_BYTE,
                  first_comm_);
  }
  return {gamma_j_i, gamma_j_0};
}
// run this on this first process
size_t ChainKMCSimulation::SelectEvent(double gamma_j_i, double gamma_j_0) const {
  static std::uniform_real_distribution<double> distribution(0.0, 1.0);
  const double random_number = distribution(generator_);
  if (random_number <= gamma_j_i) {
    // indirect move
    const int flip = int(log(random_number / gamma_j_i) / log(gamma_j_0)) + 1;
  } else {
    // direct move
    const int flip = 0;
  }
  return 0;
}

void ChainKMCSimulation::Simulate() {
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
      auto one_step_time = log(distribution(generator_)) / (total_rate_k_l_ * kPrefactor);
      time_ -= one_step_time;
      one_step_change_ = executed_invent.GetEnergyChange();
      energy_ += one_step_change_;
      one_step_barrier_ = executed_invent.GetForwardBarrier();
    }
    ++steps_;
    // world_.barrier();
  }

}

}


// namespace kmc