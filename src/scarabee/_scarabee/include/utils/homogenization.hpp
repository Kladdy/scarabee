#ifndef SCARABEE_HOMOGENIZATION_H
#define SCARABEE_HOMOGENIZATION_H
#include <data/cross_section.hpp>

#include <xtensor/containers/xtensor.hpp>

#include <functional>
#include <memory>
#include <vector>

namespace scarabee {

class HomogenizationAdaptor {
 public:
  template <typename Solver>
  HomogenizationAdaptor(const Solver& solver)
      : flux_{[&solver](std::size_t i, std::size_t g) {
          return solver.flux(i, g);
        }},
        xs_{[&solver](std::size_t i) { return solver.xs(i); }},
        volume_{[&solver](std::size_t i) { return solver.volume(i); }},
        size_{[&solver]() { return solver.size(); }},
        ngroups_{[&solver]() { return solver.ngroups(); }} {}

  template <typename Solver, typename FluxViewer>
  HomogenizationAdaptor(const Solver& solver, const FluxViewer& fv)
      : flux_{[&fv](std::size_t i, std::size_t g) { return fv(i, g); }},
        xs_{[&solver](std::size_t i) { return solver.xs(i); }},
        volume_{[&solver](std::size_t i) { return solver.volume(i); }},
        size_{[&solver]() { return solver.size(); }},
        ngroups_{[&solver]() { return solver.ngroups(); }} {}

  double flux(std::size_t i, std::size_t g) const { return flux_(i, g); }

  std::shared_ptr<CrossSection> xs(std::size_t i) const { return xs_(i); }

  double volume(std::size_t i) const { return volume_(i); }

  std::size_t size() const { return size_(); }

  std::size_t ngroups() const { return ngroups_(); }

 private:
  std::function<double(std::size_t, std::size_t)> flux_;
  std::function<std::shared_ptr<CrossSection>(std::size_t)> xs_;
  std::function<double(std::size_t)> volume_;
  std::function<std::size_t()> size_;
  std::function<std::size_t()> ngroups_;
};

std::shared_ptr<CrossSection> homogenize(
    const HomogenizationAdaptor& system,
    const std::vector<std::size_t>& regions);

xt::xtensor<double, 1> homogenize_flux_spectrum(
    const HomogenizationAdaptor& system,
    const std::vector<std::size_t>& regions);
}  // namespace scarabee

#endif
