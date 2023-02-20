#if __has_include("StdAfx.h")
#	include "StdAfx.h"
#endif
#include "EdgrenNylinder.h"
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <array>

//Vilhelm Edgren och Per Nylinders "Funktioner och tabeller för bestämning av avsmalning och formkvot under bark: Tall och gran i norra och södra Sverige"
//Implemented by Jonas Hammarberg <jonas.hammarberg@skogrkaupa.com> 2023-01-26

namespace skogrkaupa {
	namespace forestry {
		namespace edgren_nylinder {
			namespace {
				constexpr size_t NUMBER_OF_FORM_CLASSES = 7;	//[0 /] 52.5% / 57.5% / 62.5% / 67.5% / 72.5% / 77.5% [/ 82.5%]

				const std::array<stem_form_constants_t, NUMBER_OF_SPECIE_TYPES * NUMBER_OF_FORM_CLASSES> stem_form_constants = { {
						//southern pine
						{  0.620, 0.8409, 15.970, 285.280, 183.440 },
						{  0.620, 0.3694, 14.948, 285.280, 458.590 },
						{  1.594, 0.4251, 14.214, 147.440, 463.070 },
						{  3.240, 1.5290, 13.646,  98.601, 171.700 },
						{  6.320, 3.9740, 13.240,  71.915,  95.286 },
						{ 13.070, 5.5100, 12.951,  54.050,  84.904 },
						{ 33.502, 6.4450, 12.755,  39.982,  83.659 },

						//northern pine
						{  0.620, 1.513, 14.233, 311.680, 123.910 },
						{  0.620, 1.228, 13.321, 311.680, 172.850 },
						{  1.594, 1.506, 12.657, 160.290, 167.680 },
						{  3.240, 2.493, 12.177, 106.670, 128.160 },
						{  6.320, 4.488, 11.880,  77.416,  94.947 },
						{ 13.056, 6.602, 11.759,  57.767,  81.725 },
						{ 32.307, 7.594, 11.753,  42.808,  80.776 },

						//southern spruce
						{  0.620, 0.892, 15.765, 287.440, 174.400 },
						{  0.620, 0.923, 14.818, 287.440, 202.650 },
						{  1.594, 1.093, 14.032, 148.970, 202.510 },
						{  3.240, 2.164, 13.479,  99.532, 132.690 },
						{  6.320, 3.324, 13.040,  72.736, 108.430 },
						{ 13.059, 4.463, 12.775,  54.618,  97.490 },
						{ 33.208, 5.586, 12.578,  40.509,  91.770 },

						//northern spruce
						{  0.620, 1.671, 16.104, 286.360, 103.060 },
						{  0.620, 1.422, 14.883, 286.360, 140.910 },
						{  1.594, 1.976, 13.784, 151.040, 127.890 },
						{  3.240, 2.906, 12.906, 102.700, 110.680 },
						{  6.320, 3.759, 12.099,  76.543, 105.150 },
						{ 13.056, 4.026, 11.321,  59.096, 112.590 },
						{ 32.012, 3.595, 10.540,  45.754, 134.760 }
					} };

				size_t form_class(const double form_quotient) {
					return std::clamp(static_cast<size_t>(trunc((form_quotient - 0.475) / 0.05)), 1U, NUMBER_OF_FORM_CLASSES) - 1U;
				}

				inline double safediv(const double t, const double d) { return (d != 0.0) ? t / d : t; }
				template<typename T> T sqr(const T lhs) { return lhs * lhs; }
				constexpr double pi = 3.14159265358979323846;
			}


			double Calculator::calculate_height_at_diameter(const double diameter_cm) const { return calculate_height_as_share_at_diameter(diameter_cm) * height_m_; }

			double Calculator::calculate_height_as_share_at_diameter(const double diameter_cm) const {
				//topp
				if(diameter_cm <= diameter_at_start_of_top_) {
					if(diameter_cm == diameter_at_start_of_top_) return START_OF_TOP_AS_HEIGHT_SHARE;
					const auto andel = 1.0 - (pow(10.0, diameter_cm * quota_ / k_.R) - 1.0) / k_.gamma;
					assert(START_OF_TOP_AS_HEIGHT_SHARE <= andel && andel <= 1.0);
					return std::clamp(andel, START_OF_TOP_AS_HEIGHT_SHARE, 1.0);
				}

				//if(diameter_cm <= rotdiameter_) {	//underdel
				if(diameter_cm <= diameter_at_end_of_root_swell_) {
					if(diameter_cm == diameter_at_end_of_root_swell_) return root_swell_height_;
					const auto andel = 1.0 - (pow(10.0, diameter_cm * quota_ / k_.Q) - 1.0) / k_.beta;
					assert(root_swell_height_ <= andel && andel <= START_OF_TOP_AS_HEIGHT_SHARE);
					return std::clamp(andel, root_swell_height_, START_OF_TOP_AS_HEIGHT_SHARE);
				}

				//nedersta delen
				const auto andel = (pow(10.0, (100.0 - diameter_cm * quota_) / k_.q) - 1.0) * INV_ALFA;
				return std::clamp(andel, 0.0, root_swell_height_);
			}


			double Calculator::calculate_diameter_at_height(const double height_m) const { return calculate_diameter_at_height_as_share(safediv(height_m, height_m_)); }

			double Calculator::calculate_diameter_at_height_as_share(const double height_as_share) const {
				if(0.0 <= height_as_share && height_as_share <= 1.0) {
					if(height_as_share >= START_OF_TOP_AS_HEIGHT_SHARE) {	//Topp
						const auto diam = (k_.R * log10(1.0 + (1.0 - height_as_share) * k_.gamma)) * inv_quota_;
						return diam;
					}

					if(height_as_share >= root_swell_height_) {	//Mellanstam
						const auto diam = (k_.Q * log10(1.0 + (1.0 - height_as_share) * k_.beta)) * inv_quota_;
						return diam;
					}

					const auto diam = (100.0 - k_.q * log10(1 + ALFA * height_as_share)) * inv_quota_;
					return diam;
				}
				return {};
			}

			double Calculator::RootSwellHeightAsShare(const specie_type_t specie, const double form_quotient) {
				struct k_t {
					double	a;
					double	b;
				};
				static const std::array<k_t, NUMBER_OF_SPECIE_TYPES> constants = { {
					{0.06873, 0.8},
					{0.05270, 0.9},
					{0.06731, 0.8},
					{0.08631, 0.5}
				} };
				const auto& k = constants[static_cast<size_t>(specie)];
				const auto n = k.a / pow((1.0 - form_quotient), k.b);
				return n;
			}

			const stem_form_constants_t& Calculator::StemFormConstants(const specie_type_t specie, const double form_quotient) {
				assert(static_cast<size_t>(specie) < NUMBER_OF_SPECIE_TYPES);
				return stem_form_constants[static_cast<size_t>(specie) * NUMBER_OF_FORM_CLASSES + form_class(form_quotient)];
			}

			double Calculator::FormQuotient(const specie_type_t tradslag, const double form_factor, const double height_m, const double diameter_ub_cm) {
				struct k_t {
					double	a;
					double	b;
					double	c;
					double	d;
				};
				static const std::array<k_t, NUMBER_OF_SPECIE_TYPES> constants = { {
					{0.372, 0.008742, 0.003263, 0.4929 },
					{0.293, 0.006690, 0.001384, 0.6348 },
					{0.209, 0.008590, 0.003157, 0.7385 },
					{0.239, 0.010460, 0.004407, 0.6532 }
				} };
				const auto& k = constants[static_cast<size_t>(tradslag)];
				const auto f = form_factor;
				const auto n = k.a + k.b * height_m - (k.c * diameter_ub_cm) + k.d * f;
				//logger::Debug(__FUNCTION__ "v {:.4}, h {:.2}, d {:.2}, a{:.3}, b{:.6}, c{:.6}, d{:.4}, f{:.4}, n {:.5}", -0.1234, height_m, diameter_ub_cm, k.a, k.b, k.c, k.d, f, n);
				return n;
			}

		}//ns edgren_nylinder
	}//ns forestry
}//ns skogrkaupa
