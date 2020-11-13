#include "LSKMCSimulation.h"

#include <utility>
#include "KMCEvent.h"

namespace kmc {
LSKMCSimulation::LSKMCSimulation(
    const cfg::Config &config,
    unsigned long long int log_dump_steps,
    unsigned long long int config_dump_steps,
    unsigned long long int maximum_number,
    std::unordered_map<std::string, double> type_category_hashmap,
    unsigned long long int steps,
    double energy,
    double time,
    const unsigned long long int trap_steps_criteria,
    const double one_trap_time_criteria,
    const double energy_barrier_criteria)
    : KMCSimulation(config,
                    log_dump_steps,
                    config_dump_steps,
                    maximum_number,
                    std::move(type_category_hashmap),
                    steps,
                    energy,
                    time),
      trap_steps_criteria_(trap_steps_criteria),
      one_trap_time_criteria_(one_trap_time_criteria),
      energy_barrier_criteria_(energy_barrier_criteria) {}
LSKMCSimulation::~LSKMCSimulation() = default;

void LSKMCSimulation::CheckAndFix(double one_step_time) {
  if (!IsTrapped(one_step_time)) {
    return;
  }

  ClearAndSearch();

  if (!IsValidTrap())
    return;

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

  auto it = lower_bound(cumulative_probability_vector.begin(), cumulative_probability_vector.end(),
                        random_number);

  if (it == cumulative_probability_vector.cend()) {
    it--;
  }
  auto event_index = static_cast<size_t>(std::distance(cumulative_probability_vector.begin(), it));
  cfg::AtomsJump(config_, {vacancy_index_, mat_id_to_atom_id_hashmap_[event_index]});

  if (exit_time_ > 0.0) {
    time_ += exit_time_;
  }
  std::ofstream ofs("kmc_log.txt", std::ofstream::out | std::ofstream::app);
  ofs << "# LSKMC " << steps_ << " " << time_ << " ave exit time : "
      << exit_time_ << std::endl;
}
void LSKMCSimulation::ClearAndSearch() {
  markov_matrix_.clear();
  recurrent_matrix_.clear();
  transient_matrix_.clear();
  pi_vector_.clear();
  tau_vector_.clear();

  atom_id_to_mat_id_hashmap_.clear();
  mat_id_to_atom_id_hashmap_.clear();

  event_hashmap_.clear();

  Search_States_DFS();
}
KMCEvent LSKMCSimulation::CheckEventHashMapAndGet(const std::pair<size_t, size_t> &jump_pair) {
  auto it = event_hashmap_.find(jump_pair);

  if (it == event_hashmap_.end()) {
    auto barriers_and_diff_pair = barrier_predictor_.GetBarrierAndDiff(config_, jump_pair);
    KMCEvent event(jump_pair, barriers_and_diff_pair);
    event_hashmap_[jump_pair] = event;
    return event;
  } else {
    return it->second;
  }
}
void LSKMCSimulation::DFSHelper(size_t i,
                                std::unordered_set<int> &visited) {

  if (visited.count(i))
    return;
  visited.insert(i);

  for (auto j : config_.GetAtomList()[i].GetFirstNearestNeighborsList()) {
    if (visited.count(j))
      continue;

    double current_barrier = CheckEventHashMapAndGet({i, j}).GetBarrier();

    // check energy
    if (current_barrier < energy_barrier_criteria_) {
      trap_hashset_.insert(j);
      DFSHelper(j, visited);
    } else {
      absorb_hashset_.insert(j);
    }
  }
}
void LSKMCSimulation::Search_States_DFS() {
  std::unordered_set<int> visited;
  DFSHelper(vacancy_index_, visited);
  DumpBarrierStatistics();
}
void LSKMCSimulation::DumpBarrierStatistics() {
  // list of diffusion barriers
  std::multiset<double> barriers;
  std::transform(event_hashmap_.begin(),
                 event_hashmap_.end(),
                 std::inserter(barriers, barriers.begin()),
                 [](const auto &item) {
                   return item.second.GetBarrier();
                 });
  double sum = std::accumulate(barriers.begin(), barriers.end(), 0.0);

  const size_t size = barriers.size();

  std::ofstream ofs("barrier_stats.txt", std::ofstream::out | std::ofstream::app);
  ofs << "step " << steps_;
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

  it = barriers.cend();
  --it;
  ofs << " max: " << *(it);
  ofs << std::endl;
}

void LSKMCSimulation::UpdateMarkovMatrix() {
  // build empty Markov matrix
  const size_t size = trap_hashset_.size() + absorb_hashset_.size();
  // std::cerr << "SIZE = " << size << '\n';
  markov_matrix_.resize(size);
  for (auto &row : markov_matrix_) {
    row.resize(size, 0);
  }
  // build maps
  atom_id_to_mat_id_hashmap_.clear();
  mat_id_to_atom_id_hashmap_.clear();

  size_t mat_id = 0;

  // absorption states MUST come ahead of transient states
  for (auto atom_id : absorb_hashset_) {
    atom_id_to_mat_id_hashmap_[atom_id] = mat_id;
    mat_id_to_atom_id_hashmap_[mat_id] = atom_id;
    markov_matrix_[mat_id][mat_id] = 1.0;
    ++mat_id;
  }

  for (auto atom_id : trap_hashset_) {
    atom_id_to_mat_id_hashmap_[atom_id] = mat_id;
    mat_id_to_atom_id_hashmap_[mat_id] = atom_id;
    markov_matrix_[mat_id][mat_id] = 0.0;
    ++mat_id;
  }

  // calculate time of taking each state
  // mat_id_i, mat_id_j: matrix row and col
  // atom_id_i, atom_id_j: atom id
  tau_vector_.resize(trap_hashset_.size(), 0.0);
  constexpr double kSmallRate = 1e-15;
  for (size_t mat_id_i = absorb_hashset_.size(); mat_id_i < size; ++mat_id_i) {
    double all_rate = 0.0;
    auto atom_id_i = mat_id_to_atom_id_hashmap_[mat_id_i];
    for (auto atom_id_j : config_.GetAtomList()[atom_id_i].GetFirstNearestNeighborsList()) {
      // std::cerr << config_.GetAtomList()[atom_id_i].GetType() << " "
      //           << config_.GetAtomList()[atom_id_j].GetType() << " ";
      all_rate += CheckEventHashMapAndGet({atom_id_i, atom_id_j}).GetRate();
    }
    if (all_rate < kSmallRate)
      all_rate += kSmallRate;
    tau_vector_[mat_id_i - absorb_hashset_.size()] = 1.0 / all_rate;
  }

  // take care of i != j
  // i, j: matrix row and col
  // atom_id_i, atom_id_j: atom id
  for (auto mat_id_i = absorb_hashset_.size(); mat_id_i < size; ++mat_id_i) {
    auto atom_id_i = mat_id_to_atom_id_hashmap_[mat_id_i];
    const auto &atom_i_first_nearest_neighbors_list =
        config_.GetAtomList()[atom_id_i].GetFirstNearestNeighborsList();
    for (size_t mat_id_j = 0; mat_id_j < size; ++mat_id_j) {
      if (mat_id_i == mat_id_j)
        continue;

      auto atom_id_j = mat_id_to_atom_id_hashmap_[mat_id_j];
      if (std::find(atom_i_first_nearest_neighbors_list.begin(),
                    atom_i_first_nearest_neighbors_list.end(),
                    atom_id_j)
          != atom_i_first_nearest_neighbors_list.end()) {
        markov_matrix_[mat_id_i][mat_id_j] =
            CheckEventHashMapAndGet({atom_id_i, atom_id_j}).GetRate()
                * tau_vector_[mat_id_i - absorb_hashset_.size()];
      }
    }
  }

  // eigen_markov_matrix_ = StdVectorToVectorToEigenMatrixTranspose(markov_matrix_);
}
void LSKMCSimulation::UpdateRecurrentMatrixFromMarkovMatrix() {
  recurrent_matrix_.resize(trap_hashset_.size());
  for (size_t i = absorb_hashset_.size(); i < markov_matrix_.size(); ++i) {
    std::vector<double> tmp(absorb_hashset_.size(), 0.0);
    for (size_t j = 0; j < absorb_hashset_.size(); ++j) {
      tmp[j] = markov_matrix_[i][j];
    }
    recurrent_matrix_[i - absorb_hashset_.size()] = tmp;
  }
  // eigen_recurrent_matrix_ = StdVectorToVectorToEigenMatrixTranspose(recurrent_matrix_);
}
void LSKMCSimulation::UpdateTransientMatrixFromMarkovMatrix() {
  transient_matrix_.resize(trap_hashset_.size());
  for (size_t i = absorb_hashset_.size(); i < markov_matrix_.size(); ++i) {
    std::vector<double> tmp(trap_hashset_.size(), 0.0);
    for (size_t j = absorb_hashset_.size(); j < markov_matrix_.size(); ++j) {
      tmp[j - absorb_hashset_.size()] = markov_matrix_[i][j];
    }
    transient_matrix_[i - absorb_hashset_.size()] = tmp;
  }
  // eigen_transient_matrix_ = StdVectorToVectorToEigenMatrixTranspose(transient_matrix_);
}
void LSKMCSimulation::CalculateExitTimePi() {
  UpdateMarkovMatrix();
  UpdateRecurrentMatrixFromMarkovMatrix();
  UpdateTransientMatrixFromMarkovMatrix();

  Vec_t P0(transient_matrix_.size(), 0.0);
  P0[atom_id_to_mat_id_hashmap_[vacancy_index_]] = 1.0;

  arma::vec arm_P0_T = StdVectorToArmVector(P0);
  arma::mat arm_transient_matrix = StdVectorToArmMatrixTranspose(transient_matrix_);

  arma::mat arm_mat = arma::eye(arma::size(arm_transient_matrix)) - arm_transient_matrix;
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

  std::cerr << arm_shared_matrix << std::endl;
  std::cerr << StdVectorToArmVector(tau_vector_) << std::endl;
  exit_time_ = arma::dot(arm_shared_matrix, StdVectorToArmVector(tau_vector_));

  arma::mat aram_pi_matrix =
      StdVectorToArmMatrixTranspose(recurrent_matrix_) * arm_shared_matrix;

  double sum = arma::accu(aram_pi_matrix);
  aram_pi_matrix /= sum;

  pi_vector_ = ArmMatrixToStdVector(aram_pi_matrix);
}

bool LSKMCSimulation::IsValidTrap() {
  if (trap_hashset_.empty())
    return false;
  if (absorb_hashset_.empty())
    return false;

  CalculateExitTimePi();
  if (exit_time_ > 1.0e5)
    return false;
  return true;
}

bool LSKMCSimulation::IsTrapped(double one_step_time) {
  if (one_step_time < one_trap_time_criteria_) {
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

}

// namespace kmc
