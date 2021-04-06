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
void NovelKMCSimulation::Clear() {
  state_count_hashmap_.clear();
  state_chain_.clear();
  state_hashmap_.clear();
  position_id_hashmap_.clear();

  // equilibrating_event_vector_.clear();
  state_vector_.clear();
  jump_list_.clear();
  cumulated_energy_ = 0;
  cumulated_time_ = 0;
}
bool NovelKMCSimulation::GTest() const {
  for (const auto &state_count : state_count_hashmap_) {
    if (state_count.second < 3) {
      return false;
    }
  }
  return true;
}

void NovelKMCSimulation::ReviewAndFixRate() {
  for (auto &state : state_chain_) {
    std::erase_if(state.quick_event_vector_,
                  [&position_id_hashmap_ =
                  std::as_const(position_id_hashmap_)](const auto &quick_event) {
                    return position_id_hashmap_.contains(quick_event.site_id_);
                  });
    state.cumulated_absorbing_rate_ = std::accumulate(state.quick_event_vector_.begin(),
                                                      state.quick_event_vector_.end(),
                                                      0,
                                                      [](double cumulated_absorbing_rate,
                                                         const auto &quick_event) {
                                                        return cumulated_absorbing_rate +=
                                                                   quick_event.rate_;
                                                      });
    state_hashmap_.at(state.state_hash_).total_absorbing_rate_ = state.cumulated_absorbing_rate_;
  }
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
  const auto it_state = std::find_if(state_chain_.rbegin(),
                                     state_chain_.rend(),
                                     [state_hash](const StateInfo &state_info) {
                                       return state_info.state_hash_ == state_hash;
                                     });
  for (auto it = state_chain_.rbegin(); it != it_state; ++it) {
    jump_list_.push_back(it->previous_j_);
  }
  auto &equilibrating_event_vector = it_state->quick_event_vector_;

  if (equilibrating_event_vector.empty()) {
    std::cerr << " here" << std::endl;
    return;
  }
  double total_rate = it_state->cumulated_absorbing_rate_;
  double cumulated_event_probability = 0;
  for (auto &event : equilibrating_event_vector) {
    cumulated_event_probability += event.rate_ / total_rate;
    event.cumulated_event_probability_ = cumulated_event_probability;
  }

  static std::uniform_real_distribution<double> distribution(0.0, 1.0);
  const double random_number = distribution(generator_);
  auto it = std::lower_bound(equilibrating_event_vector.cbegin(),
                             equilibrating_event_vector.cend(),
                             random_number,
                             [](const auto &lhs, double value) {
                               return lhs.cumulated_event_probability_ < value;
                             });
  // If not find (maybe generated 1), which rarely happens, returns the last event
  if (it == equilibrating_event_vector.cend()) {
    it--;
  }
  jump_list_.push_back(it->next_i_);
  solved_energy_ += it->energy_change_;
}

bool NovelKMCSimulation::CheckAndSolveEquilibrium(std::ofstream &ofs) {
  bool return_value;
  if (world_rank_ == 0) {
    cumulated_energy_ += one_step_energy_change_;
    cumulated_time_ += one_step_time_change_;
    const auto state_hash = cfg::GetHashOfAState(config_, vacancy_index_);
    const StateInfo state_info(state_hash,
                               previous_j_,
                               atom_id_jump_pair_.second,
                               cumulated_energy_,
                               event_list_,
                               config_);
    state_count_hashmap_[state_hash]++;
    state_chain_.push_back(state_info);
    state_hashmap_[state_hash] =
        {state_info.state_energy_, state_info.state_rate_, 0.0, 0.0,
         state_info.cumulated_absorbing_rate_};
    position_id_hashmap_.insert(config_.GetSiteIdToAtomIdHashmap().at(vacancy_index_));

    if (state_hashmap_.size() > 125) {
      // ofs << "# Stored hashmap is too large. Reset. Chain size is " << state_chain_.size()
      //     << std::endl;
      Clear();
      return_value = false;
      // Todo check if the same state hashes have the same state rate
    } else if (state_hashmap_.size() * checking_constant > state_chain_.size() || !GTest()) {
      return_value = false;
    } else {
      ofs << "# G-test passed. ";
      ReviewAndFixRate();
      UpdateEquilibratingEventVectorAndChoose();
      ofs << "Solved time is " << solved_time_ << std::endl;
      return_value = true;
    }
  }
  MPI_Bcast(&return_value, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
  if (return_value) {
    size_t jump_list_size = jump_list_.size();
    MPI_Bcast(&jump_list_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    size_t jump_to_position;
    for (size_t i = 0; i < jump_list_size; ++i) {
      if (world_rank_ == 0) {
        jump_to_position = jump_list_[i];
      }
      MPI_Barrier(MPI_COMM_WORLD);

      MPI_Bcast(&jump_to_position, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
      cfg::AtomsJump(config_, {vacancy_index_, jump_to_position});
    }

    MPI_Bcast(&solved_time_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    time_ += solved_time_;
    MPI_Bcast(&solved_energy_, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    energy_ = solved_energy_;

    if (world_rank_ == 0) {
      previous_j_ = *jump_list_.rbegin();
    }
    MPI_Bcast(&previous_j_, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    Clear();
  }

  return return_value;
}

NovelKMCSimulation::StateInfo::StateInfo(size_t state_hash,
                                         size_t previous_j,
                                         size_t next_i,
                                         double state_energy,
                                         const std::vector<KMCEvent> &event_list,
                                         const cfg::Config &config)
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
    quick_event_vector_.push_back({atom_id,
                                   config.GetSiteIdToAtomIdHashmap().at(atom_id),
                                   kmc_event.GetEnergyChange(),
                                   kmc_event.GetForwardRate(), 0.0});
    cumulated_absorbing_rate_ += kmc_event.GetForwardRate();
  }
}

} // namespace kmc