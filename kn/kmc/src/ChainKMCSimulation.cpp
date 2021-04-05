#include "ChainKMCSimulation.h"

#include <utility>
#include <chrono>
#include <utility>
namespace kmc {

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
      previous_j(config_.GetAtomList()[vacancy_index_].GetFirstNearestNeighborsList()[0]),
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
void ChainKMCSimulation::Dump(std::ofstream &ofs) {
  if (steps_ % log_dump_steps_ == 0) {
    ofs << steps_ << '\t' << time_ << '\t' << energy_ << '\t' << one_step_barrier_ << '\t'
        << one_step_energy_change_ << '\t' << cfg::GetHashOfAState(config_, vacancy_index_)
        << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    cfg::Config::WriteConfig(config_, std::to_string(steps_) + ".cfg", 2);
  }
}

// get the energy change and the probability from j to k, pjk by the reference.
// And return the index of j and k. Only applied to 12 sub-primary processes.
// And then pass to others, now we will have 12 different numbers for 12 different second groups.
KMCEvent ChainKMCSimulation::GetEventI() {
  KMCEvent event_i;
  if (first_comm_ != MPI_COMM_NULL) {
    total_rate_k_ = 0;

    const auto i_index =
        config_.GetAtomList()[vacancy_index_].GetFirstNearestNeighborsList()[static_cast<size_t>(first_group_rank_)];

    event_i = KMCEvent(
        {vacancy_index_, i_index},
        lru_cache_barrier_predictor_.GetBarrierAndDiff(config_, {vacancy_index_, i_index}));

    const double first_rate = event_i.GetForwardRate();
    MPI_Allreduce(&first_rate, &total_rate_k_, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    event_i.CalculateProbability(total_rate_k_);
// MPI_Allgather(&first_probability_, 1, MPI_DOUBLE, probability_list_.data(), 1, MPI_DOUBLE, first_comm_);
  }

  MPI_Bcast(&event_i, sizeof(KMCEvent), MPI_BYTE, 0, second_comm_);
  return event_i;
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
double ChainKMCSimulation::BuildEventListParallel() {
  // double first_probability, first_energy_change, first_back_rate;
  auto event_k_i = GetEventI();
  const auto l_indexes = GetLIndexes();
  const auto probability_k_i = event_k_i.GetProbability();
  cfg::AtomsJump(config_, event_k_i.GetAtomIDJumpPair());
  total_rate_i_ = 0.0;
  const auto l_index = l_indexes[static_cast<size_t>(second_group_rank_)];
  KMCEvent event_i_l
      ({vacancy_index_, l_index},
       lru_cache_barrier_predictor_.GetBarrierAndDiff(config_, {vacancy_index_, l_index}));
  cfg::AtomsJump(config_, event_k_i.GetAtomIDJumpPair());

  // get sum r_{k to l}
  const double r_i_l = event_i_l.GetForwardRate();
  MPI_Allreduce(&r_i_l, &total_rate_i_, 1, MPI_DOUBLE, MPI_SUM, second_comm_);
  const auto r_i_k = event_k_i.GetBackwardRate();
  total_rate_i_ += r_i_k;
  const double probability_i_k = r_i_k / total_rate_i_;

  double t_2 = 0.0;
  if (first_comm_ != MPI_COMM_NULL) {
    const double beta_bar_k_i = probability_k_i * probability_i_k;
    const double beta_k_i = probability_k_i * (1 - probability_i_k);
    double beta_bar_k = 0.0, beta_k = 0.0, gamma_bar_k_j = 0.0, gamma_k_j = 0.0, beta_k_j = 0.0,
        alpha_k_j = 0.0;

    MPI_Allreduce(&beta_bar_k_i, &beta_bar_k, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    // MPI_Allreduce(&beta_k_i, &beta_k, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    beta_k = 1 - beta_bar_k;
    double gamma_bar_k_j_helper = 0.0, gamma_k_j_helper = 0.0, beta_k_j_helper = 0.0,
        alpha_k_j_helper = 0.0;
    if (event_k_i.GetAtomIDJumpPair().second == previous_j) {
      beta_k_j_helper = beta_k_i;
      alpha_k_j_helper = probability_k_i;
    } else {
      gamma_bar_k_j_helper = beta_bar_k_i;
      gamma_k_j_helper = beta_k_i;
    }
    MPI_Allreduce(&gamma_bar_k_j_helper, &gamma_bar_k_j, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    MPI_Allreduce(&gamma_k_j_helper, &gamma_k_j, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    MPI_Allreduce(&beta_k_j_helper, &beta_k_j, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    MPI_Allreduce(&alpha_k_j_helper, &alpha_k_j, 1, MPI_DOUBLE, MPI_SUM, first_comm_);

    const double one_over_one_minus_a_j = 1 / (1 - alpha_k_j);
    const double
        indirect_probability_k_j = one_over_one_minus_a_j * (gamma_bar_k_j / beta_k) * beta_k_j;
    const double
        indirect_probability_k_i = one_over_one_minus_a_j * (1 + gamma_bar_k_j / beta_k) * beta_k_i;
    if (event_k_i.GetAtomIDJumpPair().second == previous_j) {
      event_k_i.SetProbability(indirect_probability_k_j);
    } else {
      event_k_i.SetProbability(indirect_probability_k_i);
    }
    MPI_Allgather(&event_k_i,
                  sizeof(KMCEvent),
                  MPI_BYTE,
                  event_list_.data(),
                  sizeof(KMCEvent),
                  MPI_BYTE,
                  first_comm_);

    // calculate relative and cumulative probability
    double cumulative_provability = 0.0;
    for (auto &event : event_list_) {
      cumulative_provability += event.GetProbability();
      event.SetCumulativeProvability(cumulative_provability);
    }

    double t = 1 / total_rate_k_ / KMCEvent::kPrefactor;
    double t_i = 1 / total_rate_i_ / KMCEvent::kPrefactor;
    double ts_numerator = 0.0, ts_j_numerator = 0.0;
    double ts_numerator_helper = (t + t_i) * beta_bar_k_i;
    MPI_Allreduce(&ts_numerator_helper, &ts_numerator, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    double ts = ts_numerator / beta_bar_k;
    if (event_k_i.GetAtomIDJumpPair().second == previous_j) {
      ts_numerator_helper = 0;
    }
    MPI_Allreduce(&ts_numerator_helper, &ts_j_numerator, 1, MPI_DOUBLE, MPI_SUM, first_comm_);
    double ts_j = ts_j_numerator / gamma_bar_k_j;

    t_2 = one_over_one_minus_a_j
        * (gamma_k_j * t + gamma_bar_k_j * (ts_j + t + beta_bar_k / beta_k * ts));
  }
  MPI_Bcast(&t_2, 1, MPI_DOUBLE, 0, second_comm_);
  return t_2;
}
// run this on this first process
size_t ChainKMCSimulation::SelectEvent() const {
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

void ChainKMCSimulation::Simulate() {
  std::ofstream ofs("kmc_log.txt", std::ofstream::out | std::ofstream::app);
  if (world_rank_ == 0) {
    ofs << "steps\ttime\tenergy\tEa\tdE\tstate\n";
    ofs.precision(8);
  }

  while (steps_ < maximum_number_) {
    if (world_rank_ == 0) {
      Dump(ofs);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    one_step_time_change_ = BuildEventListParallel();
    KMCEvent selected_event;
    if (world_rank_ == 0) {
      selected_event = event_list_[SelectEvent()];
      static std::uniform_real_distribution<double> distribution(0.0, 1.0);
      one_step_energy_change_ *= distribution(generator_);
    }
    MPI_Bcast(&selected_event, sizeof(KMCEvent), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&one_step_energy_change_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#ifndef NDEBUG
    std::cout << "event choose " << event_index << '\n';
#endif
    atom_id_jump_pair_ = selected_event.GetAtomIDJumpPair();

    if (!CheckAndSolveEquilibrium(ofs)) {
      // update time and energy
      time_ += one_step_time_change_;
      one_step_energy_change_ = selected_event.GetEnergyChange();
      energy_ += one_step_energy_change_;
      one_step_barrier_ = selected_event.GetForwardBarrier();

      cfg::AtomsJump(config_, atom_id_jump_pair_);
      previous_j = atom_id_jump_pair_.second;
    } else {
      std::cout << world_rank_ << '\t' << cfg::GetHashOfAState(config_, vacancy_index_)
                << std::endl;
    }
    ++steps_;
  }
}
} // namespace kmc