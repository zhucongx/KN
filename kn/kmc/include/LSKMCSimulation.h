#ifndef KN_KN_KMC_INCLUDE_LSKMCSIMULATION_H_
#define KN_KN_KMC_INCLUDE_LSKMCSIMULATION_H_

/// Fichthorn, Kristen A., and Yangzheng Lin.
/// "A local superbasin kinetic Monte Carlo method."
/// The Journal of chemical physics 138.16 (2013): 164104.

// #include <unordered_map>
#include <armadillo>
// #include <eigen3/Eigen/Dense>

#include "KMCSimulation.h"

namespace kmc {
using Vec_t = std::vector<double>;
using Mat_t = std::vector<std::vector<double> >;
// inline Eigen::RowVectorXd StdVectorToEigenRowVector(const Vec_t &vec_t) {
//   Eigen::RowVectorXd eigen_row_vector =
//       Eigen::RowVectorXd::Map(vec_t.data(), vec_t.size());
//   return eigen_row_vector;
// }
// inline Eigen::VectorXd StdVectorToEigenVector(const Vec_t &vec_t) {
//   Eigen::VectorXd eigen_vector =
//       Eigen::VectorXd::Map(vec_t.data(), vec_t.size());
//   return eigen_vector;
// }
inline arma::vec StdVectorToArmVector(const Vec_t &vec_t) {
  return arma::conv_to<arma::vec>::from(vec_t);
}
inline arma::mat StdVectorToArmMatrix(const Mat_t &mat_t) {
  Vec_t continuum;
  size_t size = mat_t.size() * mat_t[0].size();
  continuum.reserve(size);
  for (size_t i = 0; i < mat_t.size(); ++i)
    for (size_t j = 0; j < mat_t[0].size(); ++j)
      continuum[j * mat_t.size() + i] = mat_t[i][j];
  return arma::mat(&continuum.front(), mat_t.size(), mat_t[0].size(), false);
}
inline arma::mat StdVectorToArmMatrixTranspose(const Mat_t &mat_t) {
  Vec_t continuum;
  size_t size = mat_t.size() * mat_t[0].size();
  continuum.reserve(size);
  for (size_t i = 0; i < mat_t.size(); ++i)
    for (size_t j = 0; j < mat_t[0].size(); ++j)
      continuum[j * mat_t.size() + i] = mat_t[i][j];
  return arma::mat(&continuum.front(),
                   mat_t.size(),
                   mat_t[0].size(),
                   false).t();
}
inline Vec_t ArmMatrixToStdVector(const arma::mat &mat) {
  Vec_t res(mat.n_rows * mat.n_cols, 0.0);
  for (size_t i = 0; i < mat.n_rows; ++i)
    for (size_t j = 0; j < mat.n_cols; ++j)
      res[i * mat.n_cols + j] = mat.at(i, j);
  return res;
}
// inline Eigen::MatrixXd StdVectorToVectorToEigenMatrix(const Mat_t &mat_t) {
//
//   Eigen::MatrixXd eigen_matrix(mat_t.size(),mat_t[0].size());
//
//   for(size_t i = 0; i < mat_t.size(); i++)
//     eigen_matrix.row(i) = Eigen::RowVectorXd::Map(mat_t[i].data(), mat_t[0].size());
//
//   return eigen_matrix;
// }
// inline Eigen::MatrixXd StdVectorToVectorToEigenMatrixTranspose(const Mat_t &mat_t) {
//
//   Eigen::MatrixXd eigen_matrix(mat_t[0].size(),mat_t.size());
//
//   for(size_t i = 0; i < mat_t[0].size(); i++)
//     eigen_matrix.col(i) = Eigen::VectorXd::Map(mat_t[i].data(), mat_t.size());
//
//   return eigen_matrix;
// }


class LSKMCSimulation : public KMCSimulation {
  public:
    LSKMCSimulation(const cfg::Config &config,
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
                    double energy_barrier_criteria);
    ~LSKMCSimulation() override;
    void CheckAndSolve() override;
  private:
    bool SingleCoreReturnPathAndUpdate();
    void ClearAndSearch();

    // get barrier from hashmap "event_hashmap_" or calculate and put to hashmap
    KMCEvent CheckEventHashMapAndGet(
        const std::pair<size_t, size_t> &state_and_next_position);

    void DFSHelper(size_t state,
                   std::vector<size_t> &path,
                   std::unordered_set<size_t> &visited);
    void Search_States_DFS();
    void DumpBarrierStatistics();
    void UpdateMarkovMatrix();
    void UpdateRecurrentMatrixFromMarkovMatrix();
    void UpdateTransientMatrixFromMarkovMatrix();
    void CalculateExitTimePi();
    bool IsValidTrap();
    bool IsTrapped();
    // simulation parameters
    const unsigned long long trap_steps_criteria_;
    const double one_trap_time_criteria_;
    const double energy_barrier_criteria_;

    // simulation statistics
    unsigned long long trap_step_{0};
    double exit_time_{0.0};

    Mat_t markov_matrix_;
    Mat_t recurrent_matrix_;
    Mat_t transient_matrix_;
    Vec_t pi_vector_;
    Vec_t tau_vector_;

    // helpful properties
    // M = [ I, 0 ; R T]
    // Eigen::MatrixXd eigen_markov_matrix_;
    // Eigen::MatrixXd eigen_recurrent_matrix_;
    // Eigen::MatrixXd eigen_transient_matrix_;
    // Eigen::MatrixXd eigen_pi_matrix_;
    // Eigen::MatrixXd eigen_tau_vector_;

    // hashset for transient states of each atom
    std::unordered_set<size_t> transient_hashset_;
    // hashset for absorbing states of each atom
    std::unordered_set<size_t> absorbing_hashset_;

    // state and final position --> event
    std::unordered_map<
        std::pair<size_t, size_t>,
        KMCEvent,
        boost::hash<std::pair<size_t, size_t> > > state_position_event_hashmap_;

    // state and next state --> event
    std::unordered_map<
        size_t,
        std::unordered_map<size_t, KMCEvent> > state_state_event_hashmap_;

    // state path hashmap: state i --> p1, p2, p3
    std::unordered_map<size_t, std::vector<size_t> > state_path_hashmap_;

    // a hashmap transient state -> matrix id
    std::unordered_map<size_t, size_t> state_to_matid_hashmap_;
    // a hashmap transient matrix id -> state
    std::unordered_map<size_t, size_t> matid_to_state_hashmap_;

    std::vector<size_t> path_{};
};
} // namespace kmc

#endif //KN_KN_KMC_INCLUDE_LSKMCSIMULATION_H_
