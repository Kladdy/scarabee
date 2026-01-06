#include <diffusion/diffusion_data.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <cereal/archives/portable_binary.hpp>

#include <filesystem>
#include <fstream>
#include <memory>

namespace scarabee {

DiffusionData::DiffusionData(std::shared_ptr<DiffusionCrossSection> xs)
    : xs_(xs), form_factors_(), adf_(), cdf_(), name_() {
  if (xs_ == nullptr) {
    auto mssg = "Diffusion cross section is None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
}

DiffusionData::DiffusionData(std::shared_ptr<DiffusionCrossSection> xs,
                             const xt::xtensor<double, 2>& form_factors)
    : xs_(xs), form_factors_(), adf_(), cdf_(), name_() {
  if (xs_ == nullptr) {
    auto mssg = "Diffusion cross section is None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  this->set_form_factors(form_factors);
}

DiffusionData::DiffusionData(std::shared_ptr<DiffusionCrossSection> xs,
                             const xt::xtensor<double, 2>& form_factors,
                             const xt::xtensor<double, 2>& adf,
                             const xt::xtensor<double, 2>& cdf)
    : xs_(xs), form_factors_(), adf_(), cdf_(), name_() {
  if (xs_ == nullptr) {
    auto mssg = "Diffusion cross section is None.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  this->set_form_factors(form_factors);
  this->set_adf(adf);
  this->set_cdf(cdf);
}

void DiffusionData::set_leakage_corrections(
    const std::optional<LeakageCorrections>& lc) {
  if (lc && lc->ngroups() != this->ngroups()) {
    auto mssg =
        "Number of groups in provided LeakageCorrections does not match number "
        "of groups in DiffusionData.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  leakage_corrections_ = lc;
}

void DiffusionData::set_form_factors(const xt::xtensor<double, 2>& ff) {
  for (std::size_t i = 0; i < ff.size(); i++) {
    if (ff.flat(i) < 0.) {
      auto mssg = "Form factors must be positive.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  form_factors_ = ff;
}

void DiffusionData::set_adf(const xt::xtensor<double, 2>& adf) {
  if (adf.shape()[0] != this->ngroups()) {
    auto mssg =
        "Number of groups in ADF array does not agree with the number of "
        "groups in the diffusion cross section.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (adf.shape()[1] != static_cast<std::size_t>(6)) {
    auto mssg = "The ADF array must have six columns.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (std::size_t i = 0; i < adf.size(); i++) {
    if (adf.flat(i) <= 0.) {
      auto mssg = "ADF values must be greater than zero.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  adf_ = adf;
}

void DiffusionData::set_cdf(const xt::xtensor<double, 2>& cdf) {
  if (cdf.shape()[0] != this->ngroups()) {
    auto mssg =
        "Number of groups in CDF array does not agree with the number of "
        "groups in the diffusion cross section.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (cdf.shape()[1] != static_cast<std::size_t>(4)) {
    auto mssg = "The CDF array must have four columns.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (std::size_t i = 0; i < cdf.size(); i++) {
    if (cdf.flat(i) <= 0.) {
      auto mssg = "CDF values must be greater than zero.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  cdf_ = cdf;
}

DiffusionData& DiffusionData::rotate_clockwise() {
  if (form_factors_.size() > 0) {
    xt::xtensor<double, 2> temp;

    // First transpose the matrix
    temp = xt::transpose(form_factors_);

    // Now we reverse each row
    form_factors_ = xt::flip(temp, 1);
  }

  // Must now swap ADF
  if (adf_.size() > 0) {
    for (std::size_t g = 0; g < this->ngroups(); g++) {
      const double temp2 = adf_(g, ADF::XN);

      adf_(g, ADF::XN) = adf_(g, ADF::YN);
      adf_(g, ADF::YN) = adf_(g, ADF::XP);
      adf_(g, ADF::XP) = adf_(g, ADF::YP);
      adf_(g, ADF::YP) = temp2;
    }
  }

  // Finally, swap CDF
  if (cdf_.size() > 0) {
    for (std::size_t g = 0; g < this->ngroups(); g++) {
      const double temp2 = cdf_(g, CDF::I);

      cdf_(g, CDF::I) = cdf_(g, CDF::II);
      cdf_(g, CDF::II) = cdf_(g, CDF::III);
      cdf_(g, CDF::III) = cdf_(g, CDF::IV);
      cdf_(g, CDF::IV) = temp2;
    }
  }

  return *this;
}

DiffusionData& DiffusionData::rotate_counterclockwise() {
  if (form_factors_.size() > 0) {
    xt::xtensor<double, 2> temp;

    // First we reverse each row
    temp = xt::flip(form_factors_, 1);

    // Now transpose the matrix
    form_factors_ = xt::transpose(temp);
  }

  // Must now swap ADF
  if (adf_.size() > 0) {
    for (std::size_t g = 0; g < this->ngroups(); g++) {
      const double temp2 = adf_(g, ADF::XN);

      adf_(g, ADF::XN) = adf_(g, ADF::YP);
      adf_(g, ADF::YP) = adf_(g, ADF::XP);
      adf_(g, ADF::XP) = adf_(g, ADF::YN);
      adf_(g, ADF::YP) = temp2;
    }
  }

  // Finally, swap CDF
  if (cdf_.size() > 0) {
    for (std::size_t g = 0; g < this->ngroups(); g++) {
      const double temp2 = cdf_(g, CDF::I);

      cdf_(g, CDF::I) = cdf_(g, CDF::IV);
      cdf_(g, CDF::IV) = cdf_(g, CDF::III);
      cdf_(g, CDF::III) = cdf_(g, CDF::II);
      cdf_(g, CDF::II) = temp2;
    }
  }

  return *this;
}

DiffusionData& DiffusionData::reflect_across_x_axis() {
  // In this case, the y-side ADFs switch, and CDFs switch along y
  if (form_factors_.size() > 0) {
    form_factors_ = xt::flip(form_factors_, 0);
  }

  // Must now swap ADF
  if (adf_.size() > 0) {
    for (std::size_t g = 0; g < this->ngroups(); g++) {
      std::swap(adf_(g, ADF::YN), adf_(g, ADF::YP));
    }
  }

  // Finally, swap CDF
  if (cdf_.size() > 0) {
    for (std::size_t g = 0; g < this->ngroups(); g++) {
      std::swap(cdf_(g, CDF::I), cdf_(g, CDF::IV));
      std::swap(cdf_(g, CDF::II), cdf_(g, CDF::III));
    }
  }

  return *this;
}

DiffusionData& DiffusionData::reflect_across_y_axis() {
  // In this case, the x-side ADFs switch, and CDFs switch along x
  if (form_factors_.size() > 0) {
    form_factors_ = xt::flip(form_factors_, 1);
  }

  // Must now swap ADF
  if (adf_.size() > 0) {
    for (std::size_t g = 0; g < this->ngroups(); g++) {
      std::swap(adf_(g, ADF::XN), adf_(g, ADF::XP));
    }
  }

  // Finally, swap CDF
  if (cdf_.size() > 0) {
    for (std::size_t g = 0; g < this->ngroups(); g++) {
      std::swap(cdf_(g, CDF::I), cdf_(g, CDF::II));
      std::swap(cdf_(g, CDF::III), cdf_(g, CDF::IV));
    }
  }

  return *this;
}

void DiffusionData::save(const std::string& fname) const {
  if (std::filesystem::exists(fname)) {
    std::filesystem::remove(fname);
  }

  std::ofstream file(fname, std::ios_base::binary);

  cereal::PortableBinaryOutputArchive arc(file);

  arc(*this);
}

std::shared_ptr<DiffusionData> DiffusionData::load(const std::string& fname) {
  if (std::filesystem::exists(fname) == false) {
    std::stringstream mssg;
    mssg << "The file \"" << fname << "\" does not exist.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  std::shared_ptr<DiffusionData> out(new DiffusionData());

  std::ifstream file(fname, std::ios_base::binary);

  cereal::PortableBinaryInputArchive arc(file);

  arc(*out);

  return out;
}

}  // namespace scarabee
