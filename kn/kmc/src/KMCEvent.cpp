#include <utility>
#include "KMCEvent.h"

namespace kn {
KMCEvent::KMCEvent(std::pair<int, int> jump_pair) : jump_pair_(std::move(jump_pair)) {
}

const std::pair<int, int> &KMCEvent::GetJumpPair() const {
  return jump_pair_;
}

double KMCEvent::GetBarrier() const {
  return barrier_;
}

void KMCEvent::SetBarrier(double barrier) {
  barrier_ = barrier;
}

double KMCEvent::GetRate() const {
  return rate_;
}

void KMCEvent::SetRate(double rate) {
  rate_ = rate;
}

double KMCEvent::GetEnergyChange() const {
  return energy_change_;
}

void KMCEvent::SetEnergyChange(double energy_change) {
  energy_change_ = energy_change;
}

double KMCEvent::GetProbability() const {
  return probability_;
}

void KMCEvent::SetProbability(double probability) {
  probability_ = probability;
}

double KMCEvent::GetCumulativeProvability() const {
  return cumulative_provability_;
}

void KMCEvent::SetCumulativeProvability(double cumulative_provability) {
  cumulative_provability_ = cumulative_provability;
}
} // namespace kn
