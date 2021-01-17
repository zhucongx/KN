#ifndef KN_KN_KMC_INCLUDE_KMCEVENT_H_
#define KN_KN_KMC_INCLUDE_KMCEVENT_H_
#include <utility>
#include <stddef.h>
// #include <boost/serialization/serialization.hpp>
// #include <boost/serialization/utility.hpp>

namespace kmc {
constexpr double kBoltzmannConstant = 8.617333262145e-5;
constexpr double kTemperature = 300;
constexpr double kBoltzmannConstantTimesTemperatureInv = 1 / kTemperature / kBoltzmannConstant;

class KMCEvent {
    // Todo reduce KMCEvent
  public:
    // using Event_Ctor_Pair_t = std::pair<std::pair<size_t, size_t>, std::pair<double, double>>;
    /// Constructor
    KMCEvent();
    KMCEvent(std::pair<size_t, size_t> jump_pair,
             std::pair<double, double> barrier_and_diff);
    // explicit KMCEvent(const Event_Ctor_Pair_t &event_ctor_pair);
    /// Getter
    // [[nodiscard]] Event_Ctor_Pair_t GetEventCtorPair() const;
    [[nodiscard]] const std::pair<size_t, size_t> &GetJumpPair() const;
    [[nodiscard]] double GetBarrier() const;
    [[nodiscard]] double GetRate() const;
    [[nodiscard]] double GetEnergyChange() const;
    [[nodiscard]] double GetProbability() const;
    [[nodiscard]] double GetCumulativeProvability() const;
    /// Setter
    void SetJumpPair(const std::pair<size_t, size_t> &jump_pair);
    void SetBarrier(double barrier);
    void SetRate(double rate);
    void SetEnergyChange(double energy_change);
    void SetProbability(double probability);
    void SetCumulativeProvability(double cumulative_provability);
    void CalculateProbability(double total_rates);
    /// MPI serialization access
    // friend class boost::serialization::access;
    // template<class Archive>
    // // use boost serialization templates
    // void serialize(Archive &ar, [[maybe_unused]] const size_t version) {
    //   ar & jump_pair_;
    //   ar & barrier_;
    //   ar & rate_;
    //   ar & energy_change_;
    //   ar & probability_;
    //   ar & cumulative_probability_;
    // }
  protected:
    std::pair<size_t, size_t> jump_pair_{};
    double barrier_{};
    double rate_{};
    double energy_change_{};
    double probability_{};
    double cumulative_probability_{};
};
} // namespace kmc

#endif //KN_KN_KMC_INCLUDE_KMCEVENT_H_
