#include "SemiGlobalMatching.h"
#include <algorithm>
#include "vector.h"
#include <cassert>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdint>

uint32_t SemiGlobalMatching::census_right_buffer[MAX_IMG_SIZE];
uint32_t SemiGlobalMatching::census_left_buffer[MAX_IMG_SIZE];
uint8_t SemiGlobalMatching::cost_init_buffer[MAX_DISP_IMG_SIZE];
uint16_t SemiGlobalMatching::cost_aggr_buffer[MAX_DISP_IMG_SIZE];
float SemiGlobalMatching::disp_left_buffer[MAX_IMG_SIZE];
float SemiGlobalMatching::disp_right_buffer[MAX_IMG_SIZE];

SemiGlobalMatching::SemiGlobalMatching(): width_(0), height_(0), img_left_(nullptr), img_right_(nullptr),
                                          census_left_(nullptr), census_right_(nullptr),
                                          cost_init_(nullptr), cost_aggr_(nullptr),
                                          disp_left_(nullptr), disp_right_(nullptr),
                                          is_initialized_(false)
{
}


SemiGlobalMatching::~SemiGlobalMatching()
{
    is_initialized_ = false;
}
//static uint8_t cost_init_buffer[450 * 375 * 64];
//uint8_t* cost_init_ = cost_init_buffer;
bool SemiGlobalMatching::Initialize(const uint16_t& width, const uint16_t& height, const SGMOption& option)
{
    width_ = width;
    height_ = height;
    option_ = option;

    if(width == 0 || height == 0) {
        return false;
    }
    if(option_.max_disparity <= option_.min_disparity) {
        return false;
    }
    disp_range_ = option_.max_disparity - option_.min_disparity;

    census_right_ = census_right_buffer;
    census_left_ = census_left_buffer;

    cost_init_ = cost_init_buffer;
    cost_aggr_ = cost_aggr_buffer;
    for(uint32_t i = 0; i < MAX_DISP_IMG_SIZE; i++) {
        cost_aggr_[i] = 0;
    }

    disp_left_ = disp_left_buffer;
    disp_right_ = disp_right_buffer;

    is_initialized_ = census_left_ && census_right_ && cost_init_ && cost_aggr_ && disp_left_ && disp_right_;
    is_initialized_ = true;
    return is_initialized_;
}

bool SemiGlobalMatching::Match(const uint8_t* img_left, const uint8_t* img_right, float* disp_left)
{
    if(!is_initialized_) {
        return false;
    }
    if (img_left == nullptr || img_right == nullptr) {
        return false;
    }

    img_left_ = img_left;
    img_right_ = img_right;

    // Used: img_left_, img_right_
    census_transform_5x5(img_left_, census_left_);
    census_transform_5x5(img_right_, census_right_);
    // Write: census_left_, census_right_

    // Used: census_left_, census_right_
    ComputeCost(census_left_, census_right_, cost_init_);
    // Write: cost_init_

    // Used: cost_init_, img_left_
    CostAggregation();
    // Write: cost_aggr_

    // Used: cost_aggr_
    ComputeDisparity(cost_aggr_, disp_left_, false);
    // Write: disp_left_

    if (option_.is_check_lr) {
        // Used: img_right_, cost_init_
        ComputeDisparity(cost_aggr_, disp_right_, true);
        // Write: disp_right_

        // Used: disp_left_, disp_right_
        LRCheck(disp_left_, disp_right_);
        // Write: disp_left_
    }

    if (option_.is_remove_speckles) {
        // Used: disp_left_
        RemoveSpeckles(disp_left_, 1);
        // Write: disp_left_
    }


    MedianFilter(disp_left_, disp_left_, FILTER_WINDOW_SIZE);

    memcpy(disp_left, disp_left_, height_ * width_ * sizeof(float));

	return true;
}


bool SemiGlobalMatching::Reset(const uint16_t& width, const uint16_t& height, const SGMOption& option)
{
    is_initialized_ = false;

    return Initialize(width, height, option);
}

void SemiGlobalMatching::census_transform_5x5(const uint8_t* source, uint32_t* census)
{
	if (source == nullptr || census == nullptr || width_ <= 5 || height_ <= 5) {
		return;
	}

	for (uint16_t i = 2; i < height_ - 2; i++) {
		for (uint16_t j = 2; j < width_ - 2; j++) {
			
			const uint8_t gray_center = source[i * width_ + j];
			
			uint32_t census_val = 0u;
			for (int8_t r = -2; r <= 2; r++) {
				for (int8_t c = -2; c <= 2; c++) {
					census_val <<= 1;
					const uint8_t gray = source[(i + r) * width_ + j + c];
					if (gray < gray_center) {
						census_val += 1;
					}
				}
			}

			census[i * width_ + j] = census_val;		
		}
	}
}

void SemiGlobalMatching::ComputeCost(const uint32_t* census_left, const uint32_t* census_right, uint8_t* cost_init)
{
    for (int32_t i = 0; i < height_; i++) {
        for (int32_t j = 0; j < width_; j++) {
        	for (int32_t d = option_.min_disparity; d < option_.max_disparity; d++) {
                auto& cost = cost_init[i * width_ * disp_range_ + j * disp_range_ + (d - option_.min_disparity)];
                if (j - d < 0 || j - d >= width_) {
                    cost = UINT8_MAX/2;
                    continue;
                }
                const auto& census_val_l = census_left[i * width_ + j];
                const auto& census_val_r = census_right[i * width_ + j - d];
                cost = Hamming32(census_val_l, census_val_r);
                
            }
        }
    }
}

uint8_t SemiGlobalMatching::Hamming32(const uint32_t& x, const uint32_t& y)
{
	uint32_t dist = 0, val = x ^ y;

	// Count the number of set bits
	while (val) {
		++dist;
		val &= val - 1;
	}

	return static_cast<uint8_t>(dist);
}

void SemiGlobalMatching::CostAggregation()
{
    // ↘ ↓ ↙   5  3  7
    // →    ←	 1    2
    // ↗ ↑ ↖   8  4  6
    // Directions: (dx,dy)
    // 1: ( 1, 0)
    // 2: (-1, 0)
    // 3: ( 0, 1)
    // 4: ( 0,-1)
    // 5: ( 1, 1)
    // 6: (-1,-1)
    // 7: ( 1,-1)
    // 8: (-1, 1)

    CostAggregate(img_left_, cost_init_, cost_aggr_, static_cast<int8_t>(1), static_cast<int8_t>(0));
    CostAggregate(img_left_, cost_init_, cost_aggr_, static_cast<int8_t>(-1), static_cast<int8_t>(0));
    CostAggregate(img_left_, cost_init_, cost_aggr_, static_cast<int8_t>(0), static_cast<int8_t>(1));
    CostAggregate(img_left_, cost_init_, cost_aggr_, static_cast<int8_t>(0), static_cast<int8_t>(-1));
    CostAggregate(img_left_, cost_init_, cost_aggr_, static_cast<int8_t>(1), static_cast<int8_t>(1));
    CostAggregate(img_left_, cost_init_, cost_aggr_, static_cast<int8_t>(-1), static_cast<int8_t>(-1));
    CostAggregate(img_left_, cost_init_, cost_aggr_, static_cast<int8_t>(1), static_cast<int8_t>(-1));
    CostAggregate(img_left_, cost_init_, cost_aggr_, static_cast<int8_t>(-1), static_cast<int8_t>(1));
}

void SemiGlobalMatching::CostAggregate(const uint8_t* img_data, const uint8_t* cost_init, uint16_t* cost_aggr, const int8_t dx, const int8_t dy)
{
	assert(dx >= -1 && dx <= 1 && dy >= -1 && dy <= 1 && (dx != 0 || dy != 0));
	const bool is_forward = (dx == 1 && dy == 0) || (dx == 0 && dy == 1) || (dx == 1 && dy == 1) || (dx == -1 && dy == 1);
	const int8_t direction = is_forward ? static_cast<int8_t>(1) : static_cast<int8_t>(-1);
	// for diagonal direction, need to track current row and column
	uint16_t current_row = 0;
	uint16_t current_col = 0;

	for (uint16_t i = 0u; i < (dy == 0 ? height_ : width_); i++) { // dy == 0 -> special case: horizontal
		
		const uint8_t* cost_init_pos;
		uint16_t* cost_aggr_pos;
		const uint8_t* img_pos;
		if (dx != 0 && dy == 0) // special case: horizontal
		{
			cost_init_pos = (is_forward) ? (cost_init + i * width_ * disp_range_) : (cost_init + i * width_ * disp_range_ + (width_ - 1) * disp_range_);
			cost_aggr_pos = (is_forward) ? (cost_aggr + i * width_ * disp_range_) : (cost_aggr + i * width_ * disp_range_ + (width_ - 1) * disp_range_);
			img_pos = (is_forward) ? (img_data + i * width_) : (img_data + i * width_ + width_ - 1);
		} 
		else
		{
			cost_init_pos = (is_forward) ? (cost_init + i * disp_range_) : (cost_init + (height_ - 1) * width_ * disp_range_ + i * disp_range_);
			cost_aggr_pos = (is_forward) ? (cost_aggr + i * disp_range_) : (cost_aggr + (height_ - 1) * width_ * disp_range_ + i * disp_range_);
			img_pos = (is_forward) ? (img_data + i) : (img_data + (height_ - 1) * width_ + i);
		} 
		
		uint8_t gray = *img_pos;
		uint8_t gray_last = *img_pos;
		
		uint8_t cost_last_path[MAX_DISPARITY_RANGE + 2];
		for (size_t f = 0; f < MAX_DISPARITY_RANGE + 2; f++) {
			cost_last_path[f] = UINT8_MAX;
		}
		// memcpy(cost_aggr_pos, cost_init_pos, disp_range * sizeof(uint8));
		for(uint16_t f = 0; f < disp_range_; f++) {
			cost_aggr_pos[f] += static_cast<uint16_t>(cost_init_pos[f]);
		}
		memcpy(&cost_last_path[1], cost_init_pos, disp_range_ * sizeof(uint8_t));
		
		uint8_t mincost_last_path = UINT8_MAX;
		for (u_int8_t f = 0; f < disp_range_; f++) {
			mincost_last_path = std::min(mincost_last_path, cost_init_pos[f]);
		}

		// necessary for diagonal direction
		current_row = is_forward ? 0 : height_ - 1;
		current_col = i;

		for (uint16_t j = 0; j < (dy == 0 ? width_ : height_) - 1; j ++)  // dy == 0 -> special case: horizontal
		{	
			if (dx != 0 && dy == 0) // horizontal
			{
				cost_init_pos += direction * disp_range_;
				cost_aggr_pos += direction * disp_range_;
				img_pos += direction;
			} 
			else if (dx == 0 && dy != 0) // vertical
			{
				cost_init_pos += direction * width_ * disp_range_;
				cost_aggr_pos += direction * width_ * disp_range_;
				img_pos += direction * width_;
			} 
			else // diagonal
			{
				if ((is_forward && current_col == width_ - 1 && current_row < height_ - 1) || (!is_forward && current_col == width_ - 1 && current_row > 0)) {
					// Special case: reach right edge
					cost_init_pos = cost_init + (current_row + direction) * width_ * disp_range_;
					cost_aggr_pos = cost_aggr + (current_row + direction) * width_ * disp_range_;
					img_pos = img_data + (current_row + direction) * width_;
					current_col = 0;
				}
				else if ((!is_forward && current_col == 0 && current_row > 0) || (is_forward && current_col == 0 && current_row < height_ - 1)) {
					// Special case: reach left edge
					cost_init_pos = cost_init + (current_row + direction) * width_ * disp_range_ + (width_ - 1) * disp_range_;
					cost_aggr_pos = cost_aggr + (current_row + direction) * width_ * disp_range_ + (width_ - 1) * disp_range_;
					img_pos = img_data + (current_row + direction) * width_ + (width_ - 1);
					current_col = width_ - 1;
				} 
				else if ((dx == 1 && dy == 1) || (dx == -1 && dy == -1)) // diagonal down-right or diagonal up-left
				{
					cost_init_pos += direction * (width_ + 1) * disp_range_;
					cost_aggr_pos += direction * (width_ + 1) * disp_range_;
					img_pos += direction * (width_ + 1);
				} 
				else if ((dx == -1 && dy == 1) || (dx == 1 && dy == -1)) // diagonal down-left or diagonal up-right
				{
					cost_init_pos += direction * (width_ - 1) * disp_range_;
					cost_aggr_pos += direction * (width_ - 1) * disp_range_;
					img_pos += direction * (width_ - 1);
				}
			}
			
			gray = *img_pos;
			uint8_t min_cost = UINT8_MAX;
			uint8_t last_cost = UINT8_MAX;
			uint8_t current_cost[MAX_DISPARITY_RANGE];
			for (int32_t d = 0; d < disp_range_; d++) {
				// Lr(p,d) = C(p,d) + min( Lr(p-r,d), Lr(p-r,d-1) + P1, Lr(p-r,d+1) + P1, min(Lr(p-r))+P2 ) - min(Lr(p-r))
				const uint8_t  cost = cost_init_pos[d];
				const uint16_t l1 = cost_last_path[d + 1];
				const uint16_t l2 = cost_last_path[d] + option_.p1;
				const uint16_t l3 = cost_last_path[d + 2] + option_.p1;
				const uint16_t l4 = static_cast<uint16_t>(mincost_last_path + std::max(static_cast<int>(option_.p1), static_cast<int>(option_.p2_init) / (abs(gray - gray_last) + 1)));

				const uint8_t cost_s = cost + static_cast<uint8_t>(std::min(std::min(l1, l2), std::min(l3, l4)) - mincost_last_path);

				cost_aggr_pos[d] += cost_s;
				current_cost[d] = cost_s;
				min_cost = std::min(min_cost, cost_s);

				cost_last_path[d] = last_cost; // shift by one right
				last_cost = cost_s;
			}

			mincost_last_path = min_cost;

			// cost_aggr_pos
			// memcpy(&cost_last_path[1], current_cost, disp_range * sizeof(uint8));// error: copy uint16 to uint8
			cost_last_path[disp_range_] = last_cost;

			current_row += direction;
			if ((dx == -1 && dy == 1) || (dx == 1 && dy == -1)) // diagonal down-left or diagonal up-right
			{
				current_col -= direction;
			}
			else 
			{
				current_col += direction;
			}
			
			gray_last = gray;
		}
	}
}

void SemiGlobalMatching::ComputeDisparity(const uint16_t* cost_aggr, float* disparity, bool inverse)
{
    uint16_t cost_local[MAX_DISPARITY_RANGE];
    
	for (uint16_t i = 0; i < height_; i++) {
        for (uint16_t j = 0; j < width_; j++) {
            uint16_t min_cost = UINT16_MAX;
            uint16_t sec_min_cost = UINT16_MAX;
            uint16_t best_disparity = 0;

            for (uint16_t d = option_.min_disparity; d < option_.max_disparity; d++) {
	            const uint16_t d_idx = d - option_.min_disparity;
                if(!inverse) {
                    const auto& cost = cost_local[d_idx] = cost_aggr[i * width_ * disp_range_ + j * disp_range_ + d_idx];
                    if(min_cost > cost) {
                        min_cost = cost;
                        best_disparity = d;
                    }
                } else {
                    // wieso das? Weil ich von der anderen Seite kommen, muss ich schauen wo vom rechten Bild es ist. 
                    // Deswegen suche ich mir quasi die richtigen Pixel zusammen.
                    const int32_t col_left = j + d; 
                    if (col_left >= 0 && col_left < width_) {
                        const auto& cost = cost_local[d_idx] = cost_aggr[i * width_ * disp_range_ + col_left * disp_range_ + d_idx];
                        if (min_cost > cost) {
                            min_cost = cost;
                            best_disparity = d;
                        }
                    }
                    else {
                        cost_local[d_idx] = UINT16_MAX;
                    }
                }
            }

            if (option_.is_check_unique) {
                for (uint16_t d = option_.min_disparity; d < option_.max_disparity; d++) {
                    if (d == best_disparity) {
                        continue;
                    }
                    const auto& cost = cost_local[d - option_.min_disparity];
                    sec_min_cost = std::min(sec_min_cost, cost);
                }

                // (min-sec)/min < min*(1-uniquness)
                if (sec_min_cost - min_cost <= static_cast<uint16_t>(min_cost * (1 - option_.uniqueness_ratio))) {
                    disparity[i * width_ + j] = Invalid_Float;
                    continue;
                }
            }

            if (best_disparity == option_.min_disparity || best_disparity == option_.max_disparity - 1) {
                disparity[i * width_ + j] = Invalid_Float;
                continue;
            }
            const uint16_t idx_1 = best_disparity - 1 - option_.min_disparity;
            const uint16_t idx_2 = best_disparity + 1 - option_.min_disparity;
            const uint16_t cost_1 = cost_local[idx_1];
            const uint16_t cost_2 = cost_local[idx_2];
  
            const uint16_t denom = std::max(1, cost_1 + cost_2 - 2 * min_cost);
            disparity[i * width_ + j] = static_cast<float>(best_disparity) + static_cast<float>(cost_1 - cost_2) / (denom * 2.0f);
        }
    }
}

void SemiGlobalMatching::LRCheck(float* disp_left, const float* disp_right)
{
    for (uint16_t i = 0; i < height_; i++) {
        for (uint16_t j = 0; j < width_; j++) {
        	auto& disp = disp_left[i * width_ + j];
			if(disp == Invalid_Float){
				continue;
			}

        	const auto col_right = static_cast<int32_t>(j - disp + 0.5);
        	if(col_right >= 0 && col_right < width_) {
                const auto& disp_r = disp_right[i * width_ + col_right];
        		if (abs(disp - disp_r) > option_.lrcheck_thres) {
					disp = Invalid_Float;
                }
            }
            else {
                disp = Invalid_Float;
            }
        }
    }
}

/*
// Sorting network for median of 9 values
inline float median9(float a, float b, float c,
                     float d, float e, float f,
                     float g, float h, float i)
{
    #define SWAP(x,y) if (x > y) { float tmp = x; x = y; y = tmp; }

    SWAP(a,b); SWAP(d,e); SWAP(g,h);
    SWAP(b,c); SWAP(e,f); SWAP(h,i);
    SWAP(a,b); SWAP(d,e); SWAP(g,h);

    SWAP(a,d); SWAP(b,e); SWAP(c,f);
    SWAP(d,g); SWAP(e,h); SWAP(f,i);

    SWAP(b,d); SWAP(c,e); SWAP(f,h);
    SWAP(e,g);

    #undef SWAP

    return e; // e is median of 9 values
}*/

void SemiGlobalMatching::MedianFilter(const float* in, float* out, const uint8_t wnd_size)
{
    assert(wnd_size == FILTER_WINDOW_SIZE);
	// const int32_t radius = wnd_size / 2;

	float wnd_data[FILTER_WINDOW_SIZE * FILTER_WINDOW_SIZE];
	for (int32_t i = 1; i < height_ - 1; i++) {
		for (int32_t j = 1; j < width_ - 1; j++) {
			
			/*for (int32_t r = -radius; r <= radius; r++) {
				for (int32_t c = -radius; c <= radius; c++) {
					const int32_t row = i + r;
					const int32_t col = j + c;
					if (row >= 0 && row < height && col >= 0 && col < width) {
						wnd_data.push_back(in[row * width + col]);
					}
				}
			}*/
			wnd_data[0] = in[(i - 1) * width_ + (j - 1)];
			wnd_data[1] = in[(i - 1) * width_ + j];
			wnd_data[2] = in[(i - 1) * width_ + (j + 1)];
			wnd_data[3] = in[i * width_ + (j - 1)];
			wnd_data[4] = in[i * width_ + j];
			wnd_data[5] = in[i * width_ + (j + 1)];
			wnd_data[6] = in[(i + 1) * width_ + (j - 1)];
			wnd_data[7] = in[(i + 1) * width_ + j];
			wnd_data[8] = in[(i + 1) * width_ + (j + 1)];

			std::sort(wnd_data, wnd_data + 9);
			out[i * width_ + j] = wnd_data[4];
		}
	}
}

void SemiGlobalMatching::RemoveSpeckles(float* disparity_map, const uint8_t& diff_insame)
{

	// std::vector<bool> visited(uint32(width*height),false);
	StaticVector<uint8_t, 450*375> visited;
	RingBuffer<uint32_t, 450*375> vec; // 50 = min_speckle_aera
	for (uint32_t i = 0; i < visited.size(); i++) {
		visited[i] = 0;
	}
	for(int32_t i=0;i<height_;i++) {
		for(int32_t j=0;j<width_;j++) {
			if (visited[i * width_ + j] == 1 || disparity_map[i*width_+j] == Invalid_Float) {
				continue;
			}
			// std::vector<std::pair<int32_t, int32_t>> vec;
			vec.clear();
			vec.push(i *  width_ + j);
			visited[i * width_ + j] = 1;
			while (vec.has_next()) {
				const auto& pixel = vec.next();
				uint16_t row = static_cast<uint16_t>(pixel / width_);
				uint16_t col = static_cast<uint16_t>(pixel % width_);
				const auto& disp_base = disparity_map[pixel];
				for(int r=-1;r<=1;r++) {
					for(int c=-1;c<=1;c++) {
						if(r==0&&c==0) {
							continue;
						}
						int16_t rowr = row + r;
						int16_t colc = col + c;
						if (rowr >= 0 && rowr < height_ && colc >= 0 && colc < width_) {
							uint32_t next_pixel = rowr * width_ + colc;
							if(visited[next_pixel] == 0 &&
								(disparity_map[next_pixel] != Invalid_Float) &&
								abs(disparity_map[next_pixel] - disp_base) <= diff_insame) {
								vec.push(next_pixel);
								visited[next_pixel] = 1;
							}
						}
					}
				}
			}
			vec.reset_index();
			if(vec.size() < option_.min_speckle_aera) {
				// in this case, all (speckles) pixels are saved in vec 
				while (vec.has_next()) {
					const auto& pixel = vec.next();
					disparity_map[pixel] = Invalid_Float;
				}
			}
		}
	}
}