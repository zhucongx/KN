#include "LSKMCSimulation.h"

#include <utility>
#include <mpi.h>
#include "KMCEvent.h"

namespace kmc {
// todo check energy
// todo debug with generated configs
LSKMCSimulation::LSKMCSimulation(const cfg::Config &config,
                                 unsigned long long int log_dump_steps,
                                 unsigned long long int config_dump_steps,
                                 unsigned long long int maximum_number,
                                 const std::set<std::string> &type_set,
                                 unsigned long long int steps,
                                 double energy,
                                 double time,
                                 const std::string &json_parameters_filename,
                                 size_t lru_size,
                                 unsigned long long int trap_steps_criteria,
                                 double one_trap_time_criteria,
                                 double energy_barrier_criteria)
    : KMCSimulation(config,
                    log_dump_steps,
                    config_dump_steps,
                    maximum_number,
                    type_set,
                    steps,
                    energy,
                    time,
                    json_parameters_filename,
                    lru_size),
      trap_steps_criteria_(trap_steps_criteria),
      one_trap_time_criteria_(one_trap_time_criteria),
      energy_barrier_criteria_(energy_barrier_criteria) {}

LSKMCSimulation::~LSKMCSimulation() = default;
void LSKMCSimulation::CheckAndSolve() {
  bool update = false;
  if (mpi_rank_ == 0) {
    update = SingleCoreReturnPathAndUpdate();
  }
  MPI_Bcast(&update, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
  if (update){
    size_t jump_list_size = path_.size();
    MPI_Bcast(&jump_list_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    size_t jump_to_position;
    for (size_t i = 0; i < jump_list_size; ++i) {
      if (mpi_rank_ == 0) {
        jump_to_position = path_[i];
      }
      MPI_Barrier(MPI_COMM_WORLD);

      MPI_Bcast(&jump_to_position, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
      cfg::AtomsJump(config_, {vacancy_index_, jump_to_position});
    }
  }
  // debug:
  // size_t state = cfg::GetHashOfAState(config_, vacancy_index_);
  // std::cerr << "rank: " << mpi_rank_ << ":\t" << state << std::endl;
}
bool LSKMCSimulation::SingleCoreReturnPathAndUpdate() {
  if (!IsTrapped()) {
    return false;
  }

  ClearAndSearch();

  if (!IsValidTrap())
    return false;

  // select event
  static std::uniform_real_distribution<double> distribution(0.0, 1.0);
  const double random_number = distribution(generator_);

  Vec_t cumulative_probability_vector;
  cumulative_probability_vector.reserve(pi_vector_.size());
  double cumulative_provability = 0.0;
  for (auto probability : pi_vector_) {
    cumulative_provability += probability;
    cumulative_probability_vector.push_back(cumulative_provability);
  }

  auto it = lower_bound(cumulative_probability_vector.begin(),
                        cumulative_probability_vector.end(),
                        random_number);

  if (it == cumulative_probability_vector.cend()) {
    it--;
  }
  auto event_index = static_cast<size_t>(std::distance(
      cumulative_probability_vector.begin(),
      it));
  size_t finial_state = matid_to_state_hashmap_[event_index];
  // std::cerr << "!!!!Debug only: finial state " << finial_state << std::endl;

  // if (state_path_hashmap_.find(finial_state) == state_path_hashmap_.end()) {
  //   std::cerr << "########################### path not found????";
  // }
  path_ = state_path_hashmap_.at(finial_state);

  if (exit_time_ > 0.0) {
    time_ += exit_time_;
  }
  std::ofstream ofs("kmc_log.txt", std::ofstream::out | std::ofstream::app);
  ofs << "# LSKMC " << steps_ << " " << time_ << " ave exit time : "
      << exit_time_ << std::endl;
  return true;
}
void LSKMCSimulation::ClearAndSearch() {
  markov_matrix_.clear();
  recurrent_matrix_.clear();
  transient_matrix_.clear();
  pi_vector_.clear();
  tau_vector_.clear();

  transient_hashset_.clear();
  absorbing_hashset_.clear();

  state_to_matid_hashmap_.clear();
  matid_to_state_hashmap_.clear();

  state_position_event_hashmap_.clear();
  state_state_event_hashmap_.clear();
  state_path_hashmap_.clear();

  Search_States_DFS();
}
KMCEvent LSKMCSimulation::CheckEventHashMapAndGet(
    const std::pair<size_t, size_t> &state_and_next_position) {
  auto it = state_position_event_hashmap_.find(state_and_next_position);

  if (it == state_position_event_hashmap_.end()) {
    const std::pair<size_t, size_t>
        jump_pair = {vacancy_index_, state_and_next_position.second};
    auto barriers_and_diff_pair = barrier_predictor_.GetBarrierAndDiff(
        config_, jump_pair);
    KMCEvent kmc_event(jump_pair,
                       barriers_and_diff_pair);
    state_position_event_hashmap_[state_and_next_position] = kmc_event;
    return kmc_event;
  } else {
    return it->second;
  }
}
void LSKMCSimulation::DFSHelper(size_t state,
                                std::vector<size_t> &path,
                                std::unordered_set<size_t> &visited) {

  if (visited.count(state)) {
    std::cerr << "I think we shouldn't get here..." << std::endl;
    return;
  }
  visited.insert(state);

  for (auto index :
      config_.GetAtomList()[vacancy_index_].GetFirstNearestNeighborsList()) {

    const auto kmc_event = CheckEventHashMapAndGet({state, index});
    const double current_barrier = kmc_event.GetForwardBarrier();
    cfg::AtomsJump(config_, {vacancy_index_, index});
    path.push_back(index);
    size_t state_after_jump = cfg::GetHashOfAState(config_, vacancy_index_);
    // std::cerr << "Jump forward" << std::endl;
    if (state_path_hashmap_.find(state_after_jump)
        == state_path_hashmap_.end()) {
      state_path_hashmap_[state_after_jump] = path;
    }
    // if (state_state_event_hashmap_.find(state)
    //     == state_state_event_hashmap_.end()) {
    //   state_state_event_hashmap_[state] = std::unordered_map<size_t, KMCEvent>{};
    // }
    state_state_event_hashmap_[state][state_after_jump] = kmc_event;

    if (visited.count(state_after_jump)) {
      // std::cerr << "state visited, jump back" << std::endl;
      cfg::AtomsJump(config_, {vacancy_index_, index});
      path.pop_back();
      continue;
    }

    // check energy
    if (current_barrier < energy_barrier_criteria_) {
      transient_hashset_.insert(state_after_jump);
      DFSHelper(state_after_jump, path, visited);
    } else {
      absorbing_hashset_.insert(state_after_jump);
    }
    cfg::AtomsJump(config_, {vacancy_index_, index});
    // std::cerr << "Done, jump back" << std::endl;
    path.pop_back();
  }
}
void LSKMCSimulation::Search_States_DFS() {
  std::unordered_set<size_t> visited;
  size_t current_state = cfg::GetHashOfAState(config_, vacancy_index_);
  std::vector<size_t> path{};
  state_path_hashmap_[current_state] = path;
  transient_hashset_.insert(current_state);
  DFSHelper(current_state, path, visited);

  // size_t end_state = cfg::GetHashOfAState(config_, vacancy_index_); //debug
  // std::cerr << "Debug only: state before " << current_state << " now "
  //           << end_state << std::endl;
  DumpBarrierStatistics(); // not necessary
}
void LSKMCSimulation::DumpBarrierStatistics() {
  // list of diffusion barriers
  std::multiset<double> barriers;
  std::transform(state_position_event_hashmap_.begin(),
                 state_position_event_hashmap_.end(),
                 std::inserter(barriers, barriers.begin()),
                 [](const auto &item) {
                   return item.second.GetForwardBarrier();
                 });
  double sum = std::accumulate(barriers.begin(), barriers.end(), 0.0);

  const int size = static_cast<int>(barriers.size());

  std::ofstream
      ofs("lskmc_barrier_stats.txt", std::ofstream::out | std::ofstream::app);
  ofs << "step " << steps_;
  ofs << " size " << size;
  ofs << " mean: " << (sum / static_cast<double>(size));

  auto it = barriers.cbegin();
  ofs << " min: " << *it;

  it = barriers.cbegin();
  std::advance(it, size / 4);
  ofs << " 25%: " << *it;

  it = barriers.cbegin();
  std::advance(it, size / 2);
  ofs << " 50%: " << *it;

  it = barriers.cbegin();
  std::advance(it, size * 3 / 4);
  ofs << " 75%: " << *it;

  ofs << " max: " << *(barriers.crbegin());
  ofs << std::endl;
}

void LSKMCSimulation::UpdateMarkovMatrix() {
  // build empty Markov matrix
  const size_t size = transient_hashset_.size() + absorbing_hashset_.size();
  // std::cerr << "Transient SIZE = " << transient_hashset_.size() << '\n';
  // std::cerr << "Absorbing SIZE = " << absorbing_hashset_.size() << '\n';

  markov_matrix_.resize(size);
  for (auto &row : markov_matrix_) {
    row.resize(size, 0);
  }
  // build maps
  state_to_matid_hashmap_.clear();
  matid_to_state_hashmap_.clear();

  size_t matid = 0;

  // absorption states MUST come ahead of transient states
  for (auto state : absorbing_hashset_) {
    state_to_matid_hashmap_[state] = matid;
    matid_to_state_hashmap_[matid] = state;
    markov_matrix_[matid][matid] = 1.0;
    ++matid;
  }

  for (auto state : transient_hashset_) {
    state_to_matid_hashmap_[state] = matid;
    matid_to_state_hashmap_[matid] = state;
    markov_matrix_[matid][matid] = 0.0;
    ++matid;
  }

  // calculate time of taking each state
  // matid_i, matid_j: matrix row and col
  // state_i, state_j: state row and col
  tau_vector_.resize(transient_hashset_.size(), 0.0);
  constexpr double kSmallRate = 1e-15;
  for (size_t matid_i = absorbing_hashset_.size(); matid_i < size;
       ++matid_i) {
    double all_rate = 0.0;
    auto state_i = matid_to_state_hashmap_[matid_i];
    const auto &all_events_from_state_i = state_state_event_hashmap_[state_i];
    // std::cerr << "Debug: all_events_from_state_i size(12?): "
    //           << all_events_from_state_i.size() <<
    //           std::endl;
    for (auto event_j : all_events_from_state_i) {
      all_rate += event_j.second.GetForwardRate();
    }
    // if (all_rate < kSmallRate)
    //   all_rate += kSmallRate;
    tau_vector_[matid_i - absorbing_hashset_.size()] = 1.0 / all_rate;
  }
  // std::cerr << "Tau: ";
  // for (auto a : tau_vector_) {
  //   std::cerr << a << '\t';
  // }
  // std::cerr << std::endl;
  // take care of i != j
  // i, j: matrix row and col
  // atom_id_i, atom_id_j: atom id
  for (auto matid_i = absorbing_hashset_.size(); matid_i < size; ++matid_i) {
    auto state_i = matid_to_state_hashmap_[matid_i];
    const auto &all_events_from_state_i = state_state_event_hashmap_[state_i];

    // const auto &atom_i_first_nearest_neighbors_list =
    //     config_.GetAtomList()[state_i].GetFirstNearestNeighborsList();
    for (size_t matid_j = 0; matid_j < size; ++matid_j) {
      if (matid_i == matid_j)
        continue;

      auto state_j = matid_to_state_hashmap_[matid_j];
      if (all_events_from_state_i.find(state_j)
          != all_events_from_state_i.end()) {
        markov_matrix_[matid_i][matid_j] =
            all_events_from_state_i.at(state_j).GetForwardRate()
                * tau_vector_[matid_i - absorbing_hashset_.size()];
      }
    }
  }

  // std::cerr << "M: \n";
  // for (const auto &a : markov_matrix_) {
  //   for (const auto &b : a) {
  //     std::cerr << b << '\t';
  //   }
  //   std::cerr << std::endl;
  // }
  // eigen_markov_matrix_ = StdVectorToVectorToEigenMatrixTranspose(markov_matrix_);
}
void LSKMCSimulation::UpdateRecurrentMatrixFromMarkovMatrix() {
  recurrent_matrix_.resize(transient_hashset_.size());
  for (size_t i = absorbing_hashset_.size(); i < markov_matrix_.size(); ++i) {
    std::vector<double> tmp(absorbing_hashset_.size(), 0.0);
    for (size_t j = 0; j < absorbing_hashset_.size(); ++j) {
      tmp[j] = markov_matrix_[i][j];
    }
    recurrent_matrix_[i - absorbing_hashset_.size()] = std::move(tmp);
  }
  // eigen_recurrent_matrix_ = StdVectorToVectorToEigenMatrixTranspose(recurrent_matrix_);
}
void LSKMCSimulation::UpdateTransientMatrixFromMarkovMatrix() {
  transient_matrix_.resize(transient_hashset_.size());
  for (size_t i = absorbing_hashset_.size(); i < markov_matrix_.size(); ++i) {
    std::vector<double> tmp(transient_hashset_.size(), 0.0);
    for (size_t j = absorbing_hashset_.size(); j < markov_matrix_.size(); ++j) {
      tmp[j - absorbing_hashset_.size()] = markov_matrix_[i][j];
    }
    transient_matrix_[i - absorbing_hashset_.size()] = std::move(tmp);
  }
  // eigen_transient_matrix_ = StdVectorToVectorToEigenMatrixTranspose(transient_matrix_);
}
void LSKMCSimulation::CalculateExitTimePi() {
  UpdateMarkovMatrix();
  UpdateRecurrentMatrixFromMarkovMatrix();
  UpdateTransientMatrixFromMarkovMatrix();

  Vec_t P0(transient_matrix_.size(), 0.0);
  P0[state_to_matid_hashmap_[vacancy_index_]] = 1.0;

  arma::vec arm_P0_T = StdVectorToArmVector(P0);
  arma::mat
      arm_transient_matrix = StdVectorToArmMatrixTranspose(transient_matrix_);

  arma::mat arm_mat =
      arma::eye(arma::size(arm_transient_matrix)) - arm_transient_matrix;
  // auto eigen_P0_T = StdVectorToEigenVector(P0);
  // Eigen::MatrixXd
  //     eigen_identity = Eigen::MatrixXd::Identity(static_cast<int>(transient_matrix_.size()),
  //                                                static_cast<int>(transient_matrix_.size()));
  // std::cerr << eigen_identity;
  // Eigen::MatrixXd
  //     eigen_mat = eigen_identity - StdVectorToVectorToEigenMatrixTranspose(transient_matrix_);
  // std::cerr << eigen_mat;

  arma::mat arm_shared_matrix = arma::solve(arm_mat, arm_P0_T);
  // auto eigen_shared_matrix = eigen_mat.colPivHouseholderQr().solve(eigen_P0_T);

  // std::cerr << arm_shared_matrix << std::endl;
  // std::cerr << StdVectorToArmVector(tau_vector_) << std::endl;
  exit_time_ = arma::dot(arm_shared_matrix, StdVectorToArmVector(tau_vector_));

  arma::mat aram_pi_matrix =
      StdVectorToArmMatrixTranspose(recurrent_matrix_) * arm_shared_matrix;

  double sum = arma::accu(aram_pi_matrix);
  aram_pi_matrix /= sum;

  pi_vector_ = ArmMatrixToStdVector(aram_pi_matrix);
}

bool LSKMCSimulation::IsValidTrap() {
  if (transient_hashset_.empty())
    return false;
  if (absorbing_hashset_.empty())
    return false;

  CalculateExitTimePi();
  if (exit_time_ > 1.0e5)
    return false;
  return true;
}

bool LSKMCSimulation::IsTrapped() {
  if (one_step_time_ < one_trap_time_criteria_) {
    ++trap_step_;
  } else {
    trap_step_ = 0;
  }
  if (trap_step_ < trap_steps_criteria_) {
    return false;
  } else {
    trap_step_ = 0;
    return true;
  }
}
} // namespace kmc
