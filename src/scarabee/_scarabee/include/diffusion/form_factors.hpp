#ifndef SCARABEE_FORM_FACTORS_H
#define SCARABEE_FORM_FACTORS_H

#include <utils/serialization.hpp>

#include <xtensor/containers/xtensor.hpp>

#include <cereal/cereal.hpp>

namespace scarabee {

/* The form factors are indexed using (y_ind, x_ind) so that they look correct
 * on screen. Also, the y indexing goes backwards, meaning a y_ind of zero has
 * the largest y coordinate. This is also done to ensure what you see is what
 * you get on the screen. That being said, the y widths are still index from
 * low to high as normal !
 * */
class FormFactors {
 public:
  FormFactors(const xt::xtensor<double, 2>& ff,
              const xt::xtensor<double, 1>& xw,
              const xt::xtensor<double, 1>& yw);

  FormFactors(const FormFactors& q1, const FormFactors& q2,
              const FormFactors& q3, const FormFactors& q4, bool half_pins);

  const xt::xtensor<double, 2>& form_factors() const { return form_factors_; }
  const xt::xtensor<double, 1>& x_widths() const { return x_widths_; }
  const xt::xtensor<double, 1>& y_widths() const { return y_widths_; }

  double operator()(double x, double y) const;

  double dx(std::size_t i) const;
  double dy(std::size_t j) const;

  double x_width() const { return x_width_; }
  double y_width() const { return y_width_; }

  // Methods which transform the form factors
  FormFactors& rotate_clockwise();
  FormFactors& rotate_counterclockwise();
  FormFactors& reflect_across_x_axis();
  FormFactors& reflect_across_y_axis();

  // Methods which return a cut of the form factors
  FormFactors cut_half_pos_x() const;
  FormFactors cut_half_neg_x() const;
  FormFactors cut_half_pos_y() const;
  FormFactors cut_half_neg_y() const;

  FormFactors cut_quad_I() const;
  FormFactors cut_quad_II() const;
  FormFactors cut_quad_III() const;
  FormFactors cut_quad_IV() const;

 private:
  xt::xtensor<double, 2> form_factors_;
  xt::xtensor<double, 1> x_widths_;
  xt::xtensor<double, 1> y_widths_;
  double x_width_, y_width_;

  FormFactors() = default;

  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& arc) {
    arc(CEREAL_NVP(form_factors_), CEREAL_NVP(x_widths_), CEREAL_NVP(y_widths_),
        CEREAL_NVP(x_width_), CEREAL_NVP(y_width_));
  }
};

}  // namespace scarabee

#endif
