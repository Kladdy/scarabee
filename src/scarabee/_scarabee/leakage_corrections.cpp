#include <diffusion/leakage_corrections.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

#include <sstream>

namespace scarabee {

double LeakageCorrections::D(std::size_t g) const {
  if (g >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Group index of " << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return data_.at(4 * g);
}

void LeakageCorrections::set_D(std::size_t g, double val) {
  if (g >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Group index of " << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  data_.at(4 * g) = val;
}

double LeakageCorrections::Ea(std::size_t g) const {
  if (g >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Group index of " << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return data_.at(4 * g + 1);
}

void LeakageCorrections::set_Ea(std::size_t g, double val) {
  if (g >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Group index of " << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  data_.at(4 * g + 1) = val;
}

double LeakageCorrections::Ef(std::size_t g) const {
  if (g >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Group index of " << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return data_.at(4 * g + 2);
}

void LeakageCorrections::set_Ef(std::size_t g, double val) {
  if (g >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Group index of " << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  data_.at(4 * g + 2) = val;
}

double LeakageCorrections::vEf(std::size_t g) const {
  if (g >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Group index of " << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return data_.at(4 * g + 3);
}

void LeakageCorrections::set_vEf(std::size_t g, double val) {
  if (g >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Group index of " << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  data_.at(4 * g + 3) = val;
}

double LeakageCorrections::Es(std::size_t g_in, std::size_t g_out) const {
  if (g_in >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Incident group index of " << g_in << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  if (g_out >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Outgoing group index of " << g_out << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  if (g_in >= g_out) {
    std::stringstream mssg;
    mssg << "Outgoing group index must be greater than the incident group "
            "index.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  // Retrieve the linear index into the upper diagonal down-scattering matrix
  const std::size_t N = this->ngroups();
  std::size_t k = ((N * (N - 1)) / 2) - ((N - g_in) * ((N - g_in) - 1)) / 2 +
                  g_out - g_in - 1;

  return data_.at(4 * N + k);
}

void LeakageCorrections::set_Es(std::size_t g_in, std::size_t g_out,
                                double val) {
  if (g_in >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Incident group index of " << g_in << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  if (g_out >= this->ngroups()) {
    std::stringstream mssg;
    mssg << "Outgoing group index of " << g_out << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  if (g_in >= g_out) {
    std::stringstream mssg;
    mssg << "Outgoing group index must be greater than the incident group "
            "index.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  // Retrieve the linear index into the upper diagonal down-scattering matrix
  const std::size_t N = this->ngroups();
  std::size_t k = ((N * (N - 1)) / 2) - ((N - g_in) * ((N - g_in) - 1)) / 2 +
                  g_out - g_in - 1;

  data_.at(4 * N + k) = val;
}

}  // namespace scarabee
