#include "NovelKMCSimulation.h"

#include <utility>
namespace kmc {

NovelKMCSimulation::NovelKMCSimulation(const cfg::Config &config,
                                       unsigned long long int log_dump_steps,
                                       unsigned long long int config_dump_steps,
                                       unsigned long long int maximum_number,
                                       const std::set<std::string> &type_set,
                                       unsigned long long int steps,
                                       double energy,
                                       double time,
                                       const std::string &json_parameters_filename,
                                       size_t lru_size,
                                       const size_t checking_constant)
    : ChainKMCSimulation(config,
                         log_dump_steps,
                         config_dump_steps,
                         maximum_number,
                         type_set,
                         steps,
                         energy,
                         time,
                         json_parameters_filename,
                         lru_size),
      checking_constant(checking_constant) {}
NovelKMCSimulation::~NovelKMCSimulation() = default;

bool NovelKMCSimulation::GTest() const {
  for (const auto &state_count : state_count_hashmap_) {
    if (state_count.second < 10) {
      return false;
    }
  }
}
void NovelKMCSimulation::Clear() {
  state_count_hashmap_.clear();
  state_chain_.clear();
  state_hashmap_.clear();
  equilibrating_event_vector_.clear();
  state_vector_.clear();
  jump_list_.clear();
  cumulated_energy_ = 0;
  cumulated_time_ = 0;
}
size_t NovelKMCSimulation::UpdateStateVectorAndChoose() {
  double cumulated_state_rate_ = 0;
  state_vector_.clear();
  for (const auto &state:state_hashmap_) {
    cumulated_state_rate_ += state.second.state_rate_;
    state_vector_.emplace_back(state);
  }
  double cumulative_probability = 0.0;
  double total_weighted_rate = 0.0;
  for (auto &state_info : state_vector_) {
    state_info.second.state_probability_ = state_info.second.state_rate_ / cumulated_state_rate_;
    cumulative_probability += state_info.second.state_probability_;
    state_info.second.cumulated_state_probability_ = cumulative_probability;

    total_weighted_rate +=
        state_info.second.state_probability_ * state_info.second.total_absorbing_rate_;
  }
  static std::uniform_real_distribution<double> distribution(0.0, 1.0);

  solved_time_ = -log(distribution(generator_)) / (total_weighted_rate * KMCEvent::kPrefactor);

  const double random_number = distribution(generator_);
  auto it = std::lower_bound(state_vector_.cbegin(),
                             state_vector_.cend(),
                             random_number,
                             [](const auto &lhs, double value) {
                               return lhs.second.cumulated_state_probability_ < value;
                             });
  // If not find (maybe generated 1), which rarely happens, returns the last event
  if (it == state_vector_.cend()) {
    it--;
  }
  solved_energy_ = it->second.state_energy_;
  return it->first;
}

void NovelKMCSimulation::UpdateEquilibratingEventVectorAndChoose() {
  const auto state_hash = UpdateStateVectorAndChoose();
  std::cerr << "Here1";
  const auto it_state = std::find_if(state_chain_.rbegin(),
                                     state_chain_.rend(),
                                     [state_hash](const StateInfo &state_info) {
                                       return state_info.state_hash_ == state_hash;
                                     });
  for (auto it = state_chain_.rbegin(); it != (it_state - 1); ++it) {
    jump_list_.push_back(previous_j);
    // Todo check state hash same as state_hash
  }
  equilibrating_event_vector_ = it_state->quick_event_vector_;
  double total_rate = it_state->cumulated_absorbing_rate_;
  double cumulated_event_probability = 0;
  for (auto &event : equilibrating_event_vector_) {
    cumulated_event_probability += event.rate_ / total_rate;
    event.cumulated_event_probability_ = cumulated_event_probability;
  }

  static std::uniform_real_distribution<double> distribution(0.0, 1.0);
  const double random_number = distribution(generator_);
  auto it = std::lower_bound(equilibrating_event_vector_.cbegin(),
                             equilibrating_event_vector_.cend(),
                             random_number,
                             [](const auto &lhs, double value) {
                               return lhs.cumulated_event_probability_ < value;
                             });
  // If not find (maybe generated 1), which rarely happens, returns the last event
  if (it == equilibrating_event_vector_.cend()) {
    it--;
  }
  jump_list_.push_back(it->next_i);
  solved_energy_ += it->energy_change_;
  std::cerr << "Here2";
}

bool NovelKMCSimulation::CheckAndSolveEquilibrium(std::ofstream &ofs) {
  bool return_value;
  if (world_rank_ == 0) {
    cumulated_energy_ += one_step_energy_change_;
    cumulated_time_ += one_step_time_change_;
    const auto state_hash = cfg::GetHashOfAState(config_, vacancy_index_);
    const StateInfo state_info(state_hash,
                               previous_j,
                               atom_id_jump_pair_.second,
                               cumulated_energy_,
                               event_list_);
    state_count_hashmap_[state_hash]++;
    state_chain_.push_back(state_info);
    state_hashmap_[state_hash] =
        {state_info.state_energy_, state_info.state_rate_, 0.0, 0.0,
         state_info.cumulated_absorbing_rate_};
    if (state_hashmap_.size() > 20) {
      ofs << "# Stored hashmap is too large. Clear. Continuing ChainKMC" << std::endl;
      Clear();
      return_value = false;
      // Todo check if the same state hashes have the same state rate
    } else if (state_hashmap_.size() * checking_constant > state_chain_.size()) {
      return_value = false;
    } else if (!GTest()) {
      // ofs << "# G-test failed. Continuing ChainKMC" << std::endl;
      return_value = false;
    } else {
      ofs << "# G-test passed. ";
      UpdateEquilibratingEventVectorAndChoose();
      ofs << "Solved time is " << solved_time_ << std::endl;
      return_value = true;
    }
  }
  MPI_Bcast(&return_value, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
  if (return_value) {
    size_t jump_list_size = jump_list_.size();
    MPI_Bcast(&jump_list_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    jump_list_.resize(jump_list_size);
    MPI_Bcast(static_cast<void *>(jump_list_.data()),
              jump_list_size,
              MPI_UNSIGNED_LONG,
              0,
              MPI_COMM_WORLD);
    std::cerr << "Here3";

    for (auto position : jump_list_) {
      cfg::AtomsJump(config_, {vacancy_index_, position});
    }
    MPI_Bcast(&solved_time_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    time_ += solved_time_;
    MPI_Bcast(&solved_energy_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    energy_ = solved_energy_;
    previous_j = *(jump_list_.rbegin() + 1);
    Clear();
  }
  return return_value;
}

NovelKMCSimulation::StateInfo::StateInfo(size_t state_hash,
                                         size_t previous_j,
                                         size_t next_i,
                                         double state_energy,
                                         const std::vector<KMCEvent> &event_list)
    : state_hash_(state_hash),
      previous_j_(previous_j),
      next_i_(next_i),
      state_energy_(state_energy),
      state_rate_(std::exp(-state_energy * KMCEvent::kBoltzmannConstantTimesTemperatureInv)) {
  quick_event_vector_.reserve(kFirstEventListSize - 1);
  for (const auto &kmc_event : event_list) {
    auto atom_id = kmc_event.GetAtomIDJumpPair().second;
    if (atom_id == previous_j || atom_id == next_i)
      continue;
    quick_event_vector_.push_back({kmc_event.GetAtomIDJumpPair().second,
                                   kmc_event.GetEnergyChange(),
                                   kmc_event.GetForwardRate(), 0.0});
    cumulated_absorbing_rate_ += kmc_event.GetForwardRate();
  }
}

} // namespace kmc