#ifndef KN_KN_KMC_INCLUDE_NOVELKMCSIMULATION_H_
#define KN_KN_KMC_INCLUDE_NOVELKMCSIMULATION_H_
#include "ChainKMCSimulation.h"
#include <memory>
namespace kmc {
class NovelKMCSimulation : public ChainKMCSimulation {
  public:
    NovelKMCSimulation(const cfg::Config &config,
                       unsigned long long int log_dump_steps,
                       unsigned long long int config_dump_steps,
                       unsigned long long int maximum_number,
                       const std::set<std::string> &type_set,
                       unsigned long long int steps,
                       double energy,
                       double time,
                       const std::string &json_parameters_filename,
                       size_t lru_size,
                       size_t checking_constant);
    ~NovelKMCSimulation() override;
  protected:
    struct QuickEvent {
      size_t next_i;
      double energy_change_;
      double rate_;
      double cumulated_event_probability_;
    };
    struct QuickStateInfo {
      double state_energy_;
      double state_rate_;
      double state_probability_;
      double cumulated_state_probability_;
      double total_absorbing_rate_;
    };
    struct StateInfo {
      StateInfo(size_t state_hash,
                size_t previous_j,
                size_t next_i,
                double state_energy,
                const std::vector<KMCEvent> &event_list);
      size_t state_hash_;
      size_t previous_j_;
      size_t next_i_;
      double cumulated_absorbing_rate_{0.0};
      double state_energy_;
      double state_rate_;
      // other_next_i, rate
      std::vector<QuickEvent> quick_event_vector_{};
    };

    bool GTest() const;
    void Clear();
    // Returns state hash
    size_t UpdateStateVectorAndChoose();
    void UpdateEquilibratingEventVectorAndChoose();
    bool CheckAndSolveEquilibrium(std::ofstream &ofs) override;
    const size_t checking_constant;

    std::unordered_map<size_t,size_t> state_count_hashmap_{};
    std::vector<StateInfo> state_chain_{};
    // state_hash, QuickStateInfo;
    std::unordered_map<size_t, QuickStateInfo> state_hashmap_{};
    // state_hash, next_i, probability;
    // std::vector<QuickEvent> equilibrating_event_vector_;
    // state_hash QuickStateInfo
    std::vector<std::pair<size_t, QuickStateInfo>> state_vector_;
    double cumulated_energy_{0.0};
    double cumulated_time_{0.0};
    std::vector<size_t> jump_list_{};
    double solved_time_;
    double solved_energy_;
};
} // namespace kmc
#endif //KN_KN_KMC_INCLUDE_NOVELKMCSIMULATION_H_
