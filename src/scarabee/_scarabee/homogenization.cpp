#include <utils/homogenization.hpp>
#include <utils/scarabee_exception.hpp>

#include <sstream>

namespace scarabee {

std::shared_ptr<CrossSection> homogenize(
    const HomogenizationAdaptor& system,
    const std::vector<std::size_t>& regions) {
  // Make sure we were actually provided with regions
  if (regions.empty()) {
    const auto mssg = "No regions were provided for homogenization.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Check all regions are valid
  if (regions.size() > system.size()) {
    const auto mssg =
        "The number of provided regions is greater than the number of regions.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (const auto m : regions) {
    if (m >= system.size()) {
      const auto mssg = "Invalid region index in homogenization list.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  // We now begin homogenization
  const std::size_t NR = regions.size();
  const std::size_t NG = system.ngroups();

  std::size_t max_l = 0;
  for (const auto m : regions) {
    const auto m_max_l = system.xs(m)->max_legendre_order();
    if (m_max_l > max_l) {
      max_l = m_max_l;
    }
  }

  xt::xtensor<double, 1> Et = xt::zeros<double>({NG});
  xt::xtensor<double, 1> Dtr = xt::zeros<double>({NG});
  xt::xtensor<double, 1> Ea = xt::zeros<double>({NG});
  xt::xtensor<double, 3> Es = xt::zeros<double>({max_l + 1, NG, NG});
  xt::xtensor<double, 1> Ef = xt::zeros<double>({NG});
  xt::xtensor<double, 1> vEf = xt::zeros<double>({NG});
  xt::xtensor<double, 1> chi = xt::zeros<double>({NG});

  // We need to calculate the total fission production in each volume for
  // generating the homogenized fission spectrum.
  std::vector<double> fiss_prod(NR, 0.);
  std::size_t j = 0;
  for (const auto i : regions) {
    const auto& mat = system.xs(i);
    const double V = system.volume(i);
    for (std::size_t g = 0; g < NG; g++) {
      fiss_prod[j] += mat->vEf(g) * system.flux(i, g) * V;
    }
    j++;
  }
  const double sum_fiss_prod =
      std::accumulate(fiss_prod.begin(), fiss_prod.end(), 0.);
  const double invs_sum_fiss_prod =
      sum_fiss_prod > 0. ? 1. / sum_fiss_prod : 1.;

  for (std::size_t g = 0; g < NG; g++) {
    // Get the sum of flux*volume for this group
    double sum_fluxV = 0.;
    for (const auto i : regions) {
      sum_fluxV += system.flux(i, g) * system.volume(i);
    }

    if (sum_fluxV == 0.) {
      std::stringstream mssg;
      mssg << "Cannot homogenize cross sections. ";
      mssg << "Sum of FSR flux*volume in group " << g << " is zero. ";
      mssg << "If you see this error while using CMFD, try skipping several "
              "MOC iterations before applying CMFD.";
      const auto err_str = mssg.str();
      spdlog::error(err_str);
      throw ScarabeeException(err_str);
    }

    const double invs_sum_fluxV = 1. / sum_fluxV;

    j = 0;
    for (const auto i : regions) {
      const auto& mat = system.xs(i);
      const double V = system.volume(i);
      const double flx = system.flux(i, g);
      const double coeff = invs_sum_fluxV * flx * V;
      Et(g) += coeff * mat->Et(g);
      Dtr(g) += coeff * mat->Dtr(g);
      Ea(g) += coeff * mat->Ea(g);
      Ef(g) += coeff * mat->Ef(g);
      vEf(g) += coeff * mat->vEf(g);

      chi(g) += invs_sum_fiss_prod * fiss_prod[j] * mat->chi(g);

      for (std::size_t l = 0; l <= max_l; l++) {
        for (std::size_t gg = 0; gg < NG; gg++) {
          Es(l, g, gg) += coeff * mat->Es(l, g, gg);
        }
      }

      j++;
    }
  }

  return std::make_shared<CrossSection>(Et, Dtr, Ea, Es, Ef, vEf, chi);
}

xt::xtensor<double, 1> homogenize_flux_spectrum(
    const HomogenizationAdaptor& system,
    const std::vector<std::size_t>& regions) {
  // Make sure we were actually provided with regions
  if (regions.empty()) {
    const auto mssg =
        "No regions were provided for homogenization of flux spectrum.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Check all regions are valid
  if (regions.size() > system.size()) {
    const auto mssg =
        "The number of provided regions is greater than the number of regions.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (const auto m : regions) {
    if (m >= system.size()) {
      const auto mssg = "Invalid region index in homogenization list.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  const std::size_t NG = system.ngroups();

  // First, calculate the sum of the volumes
  double sum_V = 0.;
  for (const auto i : regions) {
    sum_V += system.volume(i);
  }
  const double invs_sum_V = 1. / sum_V;

  xt::xtensor<double, 1> spectrum = xt::zeros<double>({NG});
  for (std::size_t g = 0; g < NG; g++) {
    for (const auto i : regions) {
      spectrum(g) += invs_sum_V * system.volume(i) * system.flux(i, g);
    }
  }

  return spectrum;
}

}  // namespace scarabee
