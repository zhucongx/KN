#ifndef KN_INCLUDE_KMCEVENT_H_
#define KN_INCLUDE_KMCEVENT_H_
#include <utility>
namespace kn {
class KMCEvent {
 public:

  explicit KMCEvent(const std::pair<int, int> &jump_pair);

  const std::pair<int, int> &GetJumpPair() const;
  double GetBarrier() const;
  void SetBarrier(double barrier);
  double GetRate() const;
  void SetRate(double rate);
  double GetEnergyChange() const;
  void SetEnergyChange(double energy_change);
  double GetProbability() const;
  void SetProbability(double probability);
  double GetCumulativeProvability() const;
  void SetCumulativeProvability(double cumulative_provability);

 protected:
  std::pair<int, int> jump_pair_;
  double barrier_;
  double rate_;
  double energy_change_;
  double probability_;
  double cumulative_provability_;
};
}// namespace kn

#endif //KN_INCLUDE_KMCEVENT_H_
