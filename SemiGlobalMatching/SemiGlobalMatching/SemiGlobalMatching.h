#pragma once

#include "vector.h"

#include <limits>
static constexpr float Invalid_Float = std::numeric_limits<float>::infinity();

class SemiGlobalMatching
{
public:
	SemiGlobalMatching();
	~SemiGlobalMatching();


	struct SGMOption {
		uint8_t	num_paths;
		uint16_t  min_disparity;
		uint16_t	max_disparity;

		bool	is_check_unique;
		float	uniqueness_ratio;

		bool	is_check_lr;
		float	lrcheck_thres;
		bool	is_remove_speckles;
		uint16_t		min_speckle_aera;

		// P1,P2 
		// P2 = P2_init / (Ip-Iq)
		int16_t  p1;
		int16_t  p2_init;

		SGMOption(): num_paths(8), min_disparity(0), max_disparity(64),
		             is_check_unique(true), uniqueness_ratio(0.95f),
		             is_check_lr(true), lrcheck_thres(1.0f),
		             is_remove_speckles(true), min_speckle_aera(20),
		             p1(10), p2_init(150) { }
	};

	static constexpr uint16_t MAX_IMG_WIDTH  = 450;
	static constexpr uint16_t MAX_IMG_HEIGHT = 375;
	static constexpr uint16_t MAX_DISPARITY_RANGE = 64;
	static constexpr uint8_t FILTER_WINDOW_SIZE = 3;
	static constexpr uint32_t MAX_IMG_SIZE = MAX_IMG_WIDTH * MAX_IMG_HEIGHT;
	static constexpr uint32_t MAX_DISP_IMG_SIZE = MAX_IMG_WIDTH * MAX_IMG_HEIGHT * MAX_DISPARITY_RANGE;

public:
	bool Initialize(const uint16_t& width, const uint16_t& height, const SGMOption& option);

	bool Match(const uint8_t* img_left, const uint8_t* img_right, float* disp_left);

	bool Reset(const uint16_t& width, const uint16_t& height, const SGMOption& option);

private:

	void census_transform_5x5(const uint8_t* source, uint32_t* census);
	void ComputeCost(const uint32_t* census_left, const uint32_t* census_right, uint8_t* cost_init);
	uint8_t Hamming32(const uint32_t& x, const uint32_t& y);
	void CostAggregation();
	void CostAggregate(const uint8_t* img_data, const uint8_t* cost_init, uint16_t* cost_aggr, const int8_t dx, const int8_t dy);
	void ComputeDisparity(const uint16_t* cost_aggr, float* disparity, const bool inverse);
	void LRCheck(float* disp_left, const float* disp_right);
	void MedianFilter(const float* in, float* out, const uint8_t wnd_size);
	void RemoveSpeckles(float* disparity_map, const uint8_t& diff_insame);
	
private:
	SGMOption option_;

	uint16_t width_;
	uint16_t height_;
	uint16_t disp_range_;
	bool is_initialized_;

    static uint32_t census_right_buffer[MAX_IMG_SIZE];
	static uint32_t census_left_buffer[MAX_IMG_SIZE];
    static uint8_t cost_init_buffer[MAX_DISP_IMG_SIZE];
    static uint16_t cost_aggr_buffer[MAX_DISP_IMG_SIZE];
	static float disp_left_buffer[MAX_IMG_SIZE];
	static float disp_right_buffer[MAX_IMG_SIZE];

	const uint8_t* img_left_;
	const uint8_t* img_right_;
	uint32_t* census_left_;
	uint32_t* census_right_;
	uint8_t* cost_init_;
	uint16_t* cost_aggr_;
	float* disp_left_;
	float* disp_right_;
};

