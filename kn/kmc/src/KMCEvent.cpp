#include <cmath>
#include <utility>
#include "KMCEvent.h"

namespace kmc {
constexpr double kBoltzmannConstant = 8.617333262145e-5;
constexpr double kTemperature = 300;
constexpr double kBoltzmannConstantTimesTemperatureInv = 1 / kTemperature / kBoltzmannConstant;

KMCEvent::KMCEvent() = default;

KMCEvent::KMCEvent(std::pair<size_t, size_t> jump_pair,
                   std::pair<double, double> barrier_and_diff)
    : jump_pair_(std::move(jump_pair)),
      barrier_(barrier_and_diff.first),
      rate_(std::exp(-barrier_and_diff.first * kBoltzmannConstantTimesTemperatureInv)),
      energy_change_(barrier_and_diff.second) {}
// KMCEvent::KMCEvent(const Event_Ctor_Pair_t &event_ctor_pair)
//     : KMCEvent(event_ctor_pair.first, event_ctor_pair.second) {}
// KMCEvent::Event_Ctor_Pair_t KMCEvent::GetEventCtorPair() const {
//   return {jump_pair_, {barrier_, energy_change_}};
// }
const std::pair<size_t, size_t> &KMCEvent::GetJumpPair() const {
  return jump_pair_;
}

double KMCEvent::GetForwardBarrier() const {
  return barrier_;
}

double KMCEvent::GetForwardRate() const {
  return rate_;
}
double KMCEvent::GetBackwardRate() const {
  return std::exp(energy_change_ - barrier_ * kBoltzmannConstantTimesTemperatureInv);
}
double KMCEvent::GetEnergyChange() const {
  return energy_change_;
}

double KMCEvent::GetProbability() const {
  return probability_;
}

double KMCEvent::GetCumulativeProvability() const {
  return cumulative_probability_;
}

void KMCEvent::SetJumpPair(const std::pair<size_t, size_t> &jump_pair) {
  jump_pair_ = jump_pair;
}

void KMCEvent::SetBarrier(double barrier) {
  barrier_ = barrier;
}

void KMCEvent::SetRate(double rate) {
  rate_ = rate;
}

void KMCEvent::SetEnergyChange(double energy_change) {
  energy_change_ = energy_change;
}

void KMCEvent::SetProbability(double probability) {
  probability_ = probability;
}

void KMCEvent::SetCumulativeProvability(double cumulative_provability) {
  cumulative_probability_ = cumulative_provability;
}

void KMCEvent::CalculateProbability(double total_rates) {
  probability_ = rate_ / total_rates;
}

} // namespace kmc
