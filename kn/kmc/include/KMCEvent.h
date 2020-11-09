#ifndef KN_KN_KMC_INCLUDE_KMCEVENT_H_
#define KN_KN_KMC_INCLUDE_KMCEVENT_H_
#include <utility>

namespace kmc {
constexpr double kBoltzmannConstant = 8.617333262145e-5;
constexpr double kTemperature = 300;
constexpr double kBoltzmannConstantTimesTemperatureInv = 1 / kTemperature / kBoltzmannConstant;

class KMCEvent {
  public:
    /// Constructor
    explicit KMCEvent(std::pair<size_t, size_t> jump_pair);
    KMCEvent(std::pair<size_t, size_t> jump_pair,
             std::pair<double, double> barrier_and_diff);
/// Getter
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

  protected:
    std::pair<size_t, size_t> jump_pair_;
    double barrier_;
    double rate_;
    double energy_change_;
    double probability_;
    double cumulative_provability_;
};
} // namespace kmc

#endif //KN_KN_KMC_INCLUDE_KMCEVENT_H_
