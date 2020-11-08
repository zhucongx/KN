#include <utility>
#include "KMCEvent.h"

namespace kmc {
KMCEvent::KMCEvent(std::pair<size_t, size_t> jump_pair) : jump_pair_(std::move(jump_pair)) {
}

const std::pair<size_t, size_t> &KMCEvent::GetJumpPair() const {
  return jump_pair_;
}

double KMCEvent::GetBarrier() const {
  return barrier_;
}

double KMCEvent::GetRate() const {
  return rate_;
}

double KMCEvent::GetEnergyChange() const {
  return energy_change_;
}

double KMCEvent::GetProbability() const {
  return probability_;
}

double KMCEvent::GetCumulativeProvability() const {
  return cumulative_provability_;
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
  cumulative_provability_ = cumulative_provability;
}
} // namespace kmc
