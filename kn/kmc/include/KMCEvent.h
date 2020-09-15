#ifndef KN_INCLUDE_KMCEVENT_H_
#define KN_INCLUDE_KMCEVENT_H_
#include <utility>

namespace kn {
class KMCEvent {
  public:

    explicit KMCEvent(std::pair<int, int> jump_pair);

    [[nodiscard]] const std::pair<int, int> &GetJumpPair() const;
    [[nodiscard]] double GetBarrier() const;
    void SetBarrier(double barrier);
    [[nodiscard]] double GetRate() const;
    void SetRate(double rate);
    [[nodiscard]] double GetEnergyChange() const;
    void SetEnergyChange(double energy_change);
    [[nodiscard]] double GetProbability() const;
    void SetProbability(double probability);
    [[nodiscard]] double GetCumulativeProvability() const;
    void SetCumulativeProvability(double cumulative_provability);

  protected:
    std::pair<int, int> jump_pair_;
    double barrier_;
    double rate_;
    double energy_change_;
    double probability_;
    double cumulative_provability_;
};
} // namespace kn

#endif //KN_INCLUDE_KMCEVENT_H_
