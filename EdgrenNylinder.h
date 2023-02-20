#pragma once
#include <math.h>

//Vilhelm Edgren och Per Nylinders "Funktioner och tabeller för bestämning av avsmalning och formkvot under bark: Tall och gran i norra och södra Sverige"
//"Functions and tables for computing taper and form quotient inside bark for pine and spruce in northern and southern Sweden"
//Implemented by Jonas Hammarberg <jonas.hammarberg@skogrkaupa.com> 2023-01-26

namespace skogrkaupa {
	namespace forestry {
		namespace edgren_nylinder {
			constexpr size_t NUMBER_OF_SPECIE_TYPES = 4;
			//pine = Pinus sylvestris
			//spruce = Picea abies
			//southern/northern -- grows south or north of latitude 60 degrees north
			enum class specie_type_t { southern_pine, northern_pine, southern_spruce, northern_spruce };
			inline int operator<=>(const specie_type_t lhs, const specie_type_t rhs) noexcept { return static_cast<int>(lhs) - static_cast<int>(rhs); }

			struct stem_form_constants_t {
				double	beta;
				double	gamma;
				double	q;
				double	Q;
				double	R;
			};

			class Calculator {
			public:
				Calculator() = delete;

				Calculator(
					const specie_type_t specie, 
					const double height_m,					//basal area weighted mean height in meter
					const double diameter_under_bark_cm,	//diameter of mean basal area in cm
					const double form_factor				//quota of cylinder volume to tree volume
				)
					: specie_{ specie }
					, height_m_{ height_m }
					, diameter_under_bark_cm_{ diameter_under_bark_cm }
					, form_factor_{ form_factor }
				{}

				double calculate_height_at_diameter(const double diameter_cm) const;

				double calculate_diameter_at_height(const double hojd_m) const;

			private:

				const specie_type_t specie_;
				const double height_m_;
				const double diameter_under_bark_cm_;
				const double form_factor_;

				const double ALFA{ 10000.0 };
				const double INV_ALFA = { 1.0 / 10000.0 };
				const double START_OF_TOP_AS_HEIGHT_SHARE = { 0.6 };

				const double form_quotient_{ FormQuotient(specie_, form_factor_, height_m_, diameter_under_bark_cm_) };
				const stem_form_constants_t& k_{ StemFormConstants(specie_, form_quotient_) };
				const double quota_{ (100.0 - k_.q * log10(1.0 + ALFA * 1.3 / height_m_)) / diameter_under_bark_cm_ };
				const double inv_quota_{ 1.0 / quota_ };

				const double root_swell_height_{ RootSwellHeightAsShare(specie_, form_quotient_) };
				const double diameter_at_end_of_root_swell_{ round(calculate_diameter_at_height_as_share(root_swell_height_)) };
				const double diameter_at_start_of_top_{ round(calculate_diameter_at_height_as_share(START_OF_TOP_AS_HEIGHT_SHARE)) };


				double calculate_height_as_share_at_diameter(const double diameter_cm) const;

				double calculate_diameter_at_height_as_share(const double height_as_share) const;

				static double RootSwellHeightAsShare(const specie_type_t tradslag, const double form_quotient);
				static double FormQuotient(const specie_type_t tradslag, const double form_factor, const double height_m, const double diameter_ub_cm);
				static const stem_form_constants_t& StemFormConstants(const specie_type_t tradslag, const double form_kvot);
			};
		}//ns edgren_nylinder
	}//ns forestry
}//ns skogrkaupa
