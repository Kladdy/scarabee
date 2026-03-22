#include <diffusion/form_factors.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

#include <xtensor/views/xview.hpp>

#include <sstream>
#include <utility>

namespace scarabee {

FormFactors::FormFactors(const xt::xtensor<double, 2>& ff,
                         const xt::xtensor<double, 1>& xw,
                         const xt::xtensor<double, 1>& yw)
    : form_factors_(ff),
      x_widths_(xw),
      y_widths_(yw),
      x_width_(xt::sum(x_widths_)()),
      y_width_(xt::sum(y_widths_)()) {
  // Make sure the shapes agree !
  if (form_factors_.shape()[0] != y_widths_.size()) {
    auto mssg =
        "Shape of first dimension of form factors does not agree with the "
        "number of y widths.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (form_factors_.shape()[1] != x_widths_.size()) {
    auto mssg =
        "Shape of second dimension of form factors does not agree with the "
        "number of x widths.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure all widths are positive, as well as form factors
  for (std::size_t i = 0; i < x_widths_.size(); i++) {
    if (x_widths_(i) <= 0.) {
      std::stringstream mssg;
      mssg << "x width at index " << i << " is <= 0.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }
  }

  for (std::size_t j = 0; j < y_widths_.size(); j++) {
    if (y_widths_(j) <= 0.) {
      std::stringstream mssg;
      mssg << "y width at index " << j << " is <= 0.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }
  }

  for (std::size_t j = 0; j < y_widths_.size(); j++) {
    for (std::size_t i = 0; i < x_widths_.size(); i++) {
      if (form_factors_(j, i) < 0.) {
        auto mssg =
            "Encountered negative form factor. All form factors must be > 0.";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
    }
  }
}

FormFactors::FormFactors(const FormFactors& q1, const FormFactors& q2,
                         const FormFactors& q3, const FormFactors& q4,
                         bool half_pins) {
  // Make sure the shapes are all the same !
  if (q1.x_width() != q2.x_width() || q1.x_width() != q3.x_width() ||
      q1.x_width() != q4.x_width()) {
    auto mssg = "Width along x of all quadrants must be the same.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (q1.y_width() != q2.y_width() || q1.y_width() != q3.y_width() ||
      q1.y_width() != q4.y_width()) {
    auto mssg = "Width along y of all quadrants must be the same.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (q1.x_widths().size() != q2.x_widths().size() ||
      q1.x_widths().size() != q3.x_widths().size() ||
      q1.x_widths().size() != q3.x_widths().size() ||
      q1.x_widths().size() != q4.x_widths().size()) {
    auto mssg = "Number of tiles along x must be the same for all quadrants.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (q1.y_widths().size() != q2.y_widths().size() ||
      q1.y_widths().size() != q3.y_widths().size() ||
      q1.y_widths().size() != q3.y_widths().size() ||
      q1.y_widths().size() != q4.y_widths().size()) {
    auto mssg = "Number of tiles along y must be the same for all quadrants.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Now we can create the new arrays
  std::size_t NX = q1.x_widths().size() + q2.x_widths().size();
  std::size_t NY = q1.y_widths().size() + q4.y_widths().size();
  if (half_pins) {
    NX -= 1;
    NY -= 1;
  }
  x_widths_ = xt::zeros<double>({NX});
  y_widths_ = xt::zeros<double>({NY});
  form_factors_ = xt::zeros<double>({NY, NX});

  // Create ranges for accessing the components of the new arrays
  auto left_xrange = (NX % 2) == 0 ? xt::range(0, NX / 2)
                                   : xt::range(0, NX / 2 + 1);  // Add 1 if odd
  auto right_xrange = xt::range(NX / 2, NX);

  auto top_yrange = xt::range(NY / 2, NY);
  auto bottom_yrange = (NY % 2) == 0
                           ? xt::range(0, NY / 2)
                           : xt::range(0, NY / 2 + 1);  // Add 1 if odd
  auto ff_top_yrange = (NY % 2) == 0
                           ? xt::range(0, NY / 2)
                           : xt::range(0, NY / 2 + 1);  // Add 1 if odd
  auto ff_bottom_yrange = xt::range(NY / 2, NY);

  // Assign values. We add instead of assign so that the half pins sum to a
  // whole !
  xt::view(x_widths_, left_xrange) += q2.x_widths();
  xt::view(x_widths_, right_xrange) += q1.x_widths();
  xt::view(y_widths_, bottom_yrange) += q3.y_widths();
  xt::view(y_widths_, top_yrange) += q2.y_widths();

  // Need to make sure that we set the x and y widths
  x_width_ = xt::sum(x_widths_)();
  y_width_ = xt::sum(y_widths_)();

  // Here, we add then divide the mid-planes by 2 if we have half pins
  xt::view(form_factors_, ff_top_yrange, right_xrange) += q1.form_factors();
  xt::view(form_factors_, ff_top_yrange, left_xrange) += q2.form_factors();
  xt::view(form_factors_, ff_bottom_yrange, left_xrange) += q3.form_factors();
  xt::view(form_factors_, ff_bottom_yrange, right_xrange) += q4.form_factors();

  if (half_pins) {
    std::size_t xmid = NX / 2;
    std::size_t ymid = NY / 2;

    xt::view(form_factors_, xt::all(), xmid) *= 0.5;
    xt::view(form_factors_, ymid, xt::all()) *= 0.5;
  }
}

double FormFactors::operator()(double x, double y) const {
  // Check positions
  if (x < 0. || x > x_width_) {
    auto mssg = "x must be in [0, x_width]";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
  if (y < 0. || y > y_width_) {
    auto mssg = "y must be in [0, y_width]";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Find the x index
  std::size_t i = 0;
  double x_sum = 0.;
  for (i = 0; i < x_widths_.size(); i++) {
    x_sum += x_widths_(i);
    if (x < x_sum) break;
  }
  if (i == x_widths_.size()) i = x_widths_.size() - 1;

  // Find the y index
  std::size_t j = 0;
  double y_sum = 0.;
  for (j = 0; j < y_widths_.size(); j++) {
    y_sum += y_widths_(j);
    if (y < y_sum) break;
  }
  if (j == y_widths_.size()) j = y_widths_.size() - 1;

  // Return the form factor
  return form_factors_(form_factors_.shape()[0] - 1 - j, i);
}

double FormFactors::dx(std::size_t i) const {
  if (i >= x_widths_.size()) {
    auto mssg = "Index for x width out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  return x_widths_(i);
}

double FormFactors::dy(std::size_t j) const {
  if (j >= y_widths_.size()) {
    auto mssg = "Index for y width out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  return y_widths_(j);
}

FormFactors& FormFactors::rotate_clockwise() {
  // First transpose the matrix
  xt::xtensor<double, 2> temp = xt::transpose(form_factors_);
  // Now we reverse each row
  form_factors_ = xt::flip(temp, 1);

  // Now we must adjust the widths arrays
  std::swap(x_width_, y_width_);
  std::swap(x_widths_, y_widths_);
  y_widths_ = xt::flip(y_widths_, 0);

  return *this;
}

FormFactors& FormFactors::rotate_counterclockwise() {
  // First we reverse each row
  xt::xtensor<double, 2> temp = xt::flip(form_factors_, 1);
  // Now transpose the matrix
  form_factors_ = xt::transpose(temp);

  // Now we must adjust the widths arrays
  std::swap(x_width_, y_width_);
  std::swap(x_widths_, y_widths_);
  x_widths_ = xt::flip(x_widths_, 0);

  return *this;
}

FormFactors& FormFactors::reflect_across_x_axis() {
  form_factors_ = xt::flip(form_factors_, 0);
  y_widths_ = xt::flip(y_widths_, 0);
  return *this;
}

FormFactors& FormFactors::reflect_across_y_axis() {
  form_factors_ = xt::flip(form_factors_, 1);
  x_widths_ = xt::flip(x_widths_, 0);
  return *this;
}

FormFactors FormFactors::cut_half_pos_x() const {
  std::size_t xsize = x_widths_.size();
  std::size_t mid = xsize / 2;
  auto xrange = xt::range(mid, xsize);
  xt::xtensor<double, 1> new_x_widths = xt::view(x_widths_, xrange);
  if (xsize % 2 == 1) new_x_widths(0) *= 0.5;

  return FormFactors(xt::view(form_factors_, xt::all(), xrange), new_x_widths,
                     y_widths_);
}

FormFactors FormFactors::cut_half_neg_x() const {
  std::size_t xsize = x_widths_.size();
  std::size_t mid = xsize / 2;
  auto xrange = (xsize % 2) == 0 ? xt::range(0, mid)
                                 : xt::range(0, mid + 1);  // Add 1 if odd
  xt::xtensor<double, 1> new_x_widths = xt::view(x_widths_, xrange);
  if (xsize % 2 == 1) new_x_widths(new_x_widths.size() - 1) *= 0.5;

  return FormFactors(xt::view(form_factors_, xt::all(), xrange), new_x_widths,
                     y_widths_);
}

FormFactors FormFactors::cut_half_pos_y() const {
  std::size_t ysize = y_widths_.size();
  std::size_t mid = ysize / 2;
  auto yrange = xt::range(mid, ysize);
  xt::xtensor<double, 1> new_y_widths = xt::view(y_widths_, yrange);
  if (ysize % 2 == 1) new_y_widths(0) *= 0.5;

  auto ff_yrange = (ysize % 2) == 0 ? xt::range(0, mid)
                                    : xt::range(0, mid + 1);  // Add 1 if odd

  return FormFactors(xt::view(form_factors_, ff_yrange, xt::all()), x_widths_,
                     new_y_widths);
}

FormFactors FormFactors::cut_half_neg_y() const {
  std::size_t ysize = y_widths_.size();
  std::size_t mid = ysize / 2;
  auto yrange = (ysize % 2) == 0 ? xt::range(0, mid)
                                 : xt::range(0, mid + 1);  // Add 1 if odd
  xt::xtensor<double, 1> new_y_widths = xt::view(y_widths_, yrange);
  if (ysize % 2 == 1) new_y_widths(new_y_widths.size() - 1) *= 0.5;

  auto ff_yrange = xt::range(mid, ysize);

  return FormFactors(xt::view(form_factors_, ff_yrange, xt::all()), x_widths_,
                     new_y_widths);
}

FormFactors FormFactors::cut_quad_I() const {
  std::size_t xsize = x_widths_.size();
  std::size_t mid = xsize / 2;
  auto xrange = xt::range(mid, xsize);
  xt::xtensor<double, 1> new_x_widths = xt::view(x_widths_, xrange);
  if (xsize % 2 == 1) new_x_widths(0) *= 0.5;

  std::size_t ysize = y_widths_.size();
  mid = ysize / 2;
  auto yrange = xt::range(mid, ysize);
  xt::xtensor<double, 1> new_y_widths = xt::view(y_widths_, yrange);
  if (ysize % 2 == 1) new_y_widths(0) *= 0.5;

  auto ff_yrange = (ysize % 2) == 0 ? xt::range(0, mid)
                                    : xt::range(0, mid + 1);  // Add 1 if odd

  return FormFactors(xt::view(form_factors_, ff_yrange, xrange), new_x_widths,
                     new_y_widths);
}

FormFactors FormFactors::cut_quad_II() const {
  std::size_t xsize = x_widths_.size();
  std::size_t mid = xsize / 2;
  auto xrange = (xsize % 2) == 0 ? xt::range(0, mid)
                                 : xt::range(0, mid + 1);  // Add 1 if odd
  xt::xtensor<double, 1> new_x_widths = xt::view(x_widths_, xrange);
  if (xsize % 2 == 1) new_x_widths(new_x_widths.size() - 1) *= 0.5;

  std::size_t ysize = y_widths_.size();
  mid = ysize / 2;
  auto yrange = xt::range(mid, ysize);
  xt::xtensor<double, 1> new_y_widths = xt::view(y_widths_, yrange);
  if (ysize % 2 == 1) new_y_widths(0) *= 0.5;

  auto ff_yrange = (ysize % 2) == 0 ? xt::range(0, mid)
                                    : xt::range(0, mid + 1);  // Add 1 if odd

  return FormFactors(xt::view(form_factors_, ff_yrange, xrange), new_x_widths,
                     new_y_widths);
}

FormFactors FormFactors::cut_quad_III() const {
  std::size_t xsize = x_widths_.size();
  std::size_t mid = xsize / 2;
  auto xrange = (xsize % 2) == 0 ? xt::range(0, mid)
                                 : xt::range(0, mid + 1);  // Add 1 if odd
  xt::xtensor<double, 1> new_x_widths = xt::view(x_widths_, xrange);
  if (xsize % 2 == 1) new_x_widths(new_x_widths.size() - 1) *= 0.5;

  std::size_t ysize = y_widths_.size();
  mid = ysize / 2;
  auto yrange = (ysize % 2) == 0 ? xt::range(0, mid)
                                 : xt::range(0, mid + 1);  // Add 1 if odd
  xt::xtensor<double, 1> new_y_widths = xt::view(y_widths_, yrange);
  if (ysize % 2 == 1) new_y_widths(new_y_widths.size() - 1) *= 0.5;

  auto ff_yrange = xt::range(mid, ysize);

  return FormFactors(xt::view(form_factors_, ff_yrange, xrange), new_x_widths,
                     new_y_widths);
}

FormFactors FormFactors::cut_quad_IV() const {
  std::size_t xsize = x_widths_.size();
  std::size_t mid = xsize / 2;
  auto xrange = xt::range(mid, xsize);
  xt::xtensor<double, 1> new_x_widths = xt::view(x_widths_, xrange);
  if (xsize % 2 == 1) new_x_widths(0) *= 0.5;

  std::size_t ysize = y_widths_.size();
  mid = ysize / 2;
  auto yrange = (ysize % 2) == 0 ? xt::range(0, mid)
                                 : xt::range(0, mid + 1);  // Add 1 if odd
  xt::xtensor<double, 1> new_y_widths = xt::view(y_widths_, yrange);
  if (ysize % 2 == 1) new_y_widths(new_y_widths.size() - 1) *= 0.5;

  auto ff_yrange = xt::range(mid, ysize);

  return FormFactors(xt::view(form_factors_, ff_yrange, xrange), new_x_widths,
                     new_y_widths);
}

}  // namespace scarabee
