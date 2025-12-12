#include "SemiGlobalMatching.h"
// #include <algorithm>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h> // for exit

/* Forward declarations of static helpers (so we can call them before definitions) */
static void census_transform_5x5(const uint8_t* source, uint32_t* census);
static void ComputeCost(const uint32_t* census_left, const uint32_t* census_right, uint8_t* cost_init);
static uint8_t Hamming32(uint32_t x, uint32_t y);
static void CostAggregation(void);
static void CostAggregate(const uint8_t* img_data, const uint8_t* cost_init, uint16_t* cost_aggr, int8_t dx, int8_t dy);
static void ComputeDisparity(const uint16_t* cost_aggr, float* disparity, int inverse);
static void LRCheck(float* disp_left, const float* disp_right);
static void MedianFilter(const float* in, float* out, uint8_t wnd_size);
static void RemoveSpeckles(float* disparity_map, uint8_t diff_insame);

/* small helpers */
static inline uint8_t u8_min(uint8_t a, uint8_t b) { return a < b ? a : b; }
static inline uint16_t u16_min(uint16_t a, uint16_t b) { return a < b ? a : b; }
static inline int u_int_max_int(int a, int b) { return a > b ? a : b; }

/* SemiGlobalMatching instance (only one for now) */
SemiGlobalMatching sgm;

/* Global buffers (must match extern in header) */
uint32_t census_right_buffer[MAX_IMG_SIZE];
uint32_t census_left_buffer[MAX_IMG_SIZE];
uint8_t cost_init_buffer[MAX_DISP_IMG_SIZE];
uint16_t cost_aggr_buffer[MAX_DISP_IMG_SIZE];
float disp_left_buffer[MAX_IMG_SIZE];
float disp_right_buffer[MAX_IMG_SIZE];

bool SGM_Initialize(const uint16_t width, const uint16_t height, const SGMOption* option)
{
    sgm.width = width;
    sgm.height = height;
    sgm.option = *option;

    if(width == 0 || height == 0) {
        return false;
    }
    if(sgm.option.max_disparity <= sgm.option.min_disparity) {
        return false;
    }
    sgm.disp_range = (uint16_t)(sgm.option.max_disparity - sgm.option.min_disparity);

    sgm.census_right = census_right_buffer;
    sgm.census_left = census_left_buffer;

    sgm.cost_init = cost_init_buffer;
    sgm.cost_aggr = cost_aggr_buffer;
	/* clear aggregated cost buffer */
    memset(sgm.cost_aggr, 0, sizeof(uint16_t) * (size_t)MAX_DISP_IMG_SIZE);
    sgm.disp_left = disp_left_buffer;
    sgm.disp_right = disp_right_buffer;
    
    sgm.is_initialized = (sgm.census_left != NULL) && (sgm.census_right != NULL) &&
                         (sgm.cost_init != NULL) && (sgm.cost_aggr != NULL) &&
                         (sgm.disp_left != NULL) && (sgm.disp_right != NULL);

    return sgm.is_initialized;
}

int sum_uint32_array(const uint32_t* arr, int length) {
    int sum = 0;
    for (int i = 0; i < length; i++) {
        sum += (int)arr[i];
    }
    return sum;
}
int sum_uint8_array(const uint8_t* arr, int length) {
    int sum = 0;
    for (int i = 0; i < length; i++) {
        sum += (int)arr[i];
    }
    return sum;
}
int sum_uint16_array(const uint16_t* arr, int length) {
    int sum = 0;
    for (int i = 0; i < length; i++) {
        sum += (int)arr[i];
    }
    return sum;
}
float sum_float_array(const float* arr, int length) {
    float sum = 0;
	uint32_t valid_count = 0;
    for (int i = 0; i < length; i++) {
		if (arr[i] != INVALID_FLOAT)
		{
			sum += (float)arr[i];
			valid_count++;
		}
    }
	printf(" valid_count: %d\n", valid_count);
    return sum;
}

bool SGM_Match(const uint8_t* img_left, const uint8_t* img_right, float* disp_left)
{
    if(!sgm.is_initialized) {
        return false;
    }
    if (img_left == NULL || img_right == NULL) {
        return false;
    }

    sgm.img_left = img_left;
    sgm.img_right = img_right;

	/* Census transform */
    // Used: img_left, img_right
	printf("img_left_sum: %d\n", sum_uint8_array(sgm.img_left, MAX_IMG_SIZE));
    census_transform_5x5(sgm.img_left, sgm.census_left);
	printf("census_left_sum: %d\n", sum_uint32_array(sgm.census_left, MAX_IMG_SIZE));
	printf("img_right_sum: %d\n", sum_uint8_array(sgm.img_right, MAX_IMG_SIZE));
    census_transform_5x5(sgm.img_right, sgm.census_right);
	printf("census_right_sum: %d\n", sum_uint32_array(sgm.census_right, MAX_IMG_SIZE));
    // Write: census_left, census_right

	/* Compute matching cost (Hamming on census) */
    // Used: census_left, census_right
    ComputeCost(sgm.census_left, sgm.census_right, sgm.cost_init);
	printf("cost_init_ sum: %d\n", sum_uint8_array(sgm.cost_init, MAX_DISP_IMG_SIZE));
    // Write: cost_init

	/* Aggregate costs along directions and add into cost_aggr (uint16 sums) */
    // Used: cost_init, img_left
    CostAggregation();
	printf("cost_aggr_ sum: %d\n", sum_uint16_array(sgm.cost_aggr, MAX_DISP_IMG_SIZE));
    // Write: cost_aggr

	/* Compute disparity from aggregated costs (left-to-right) */
    // Used: cost_aggr
    ComputeDisparity(sgm.cost_aggr, sgm.disp_left, 0);
	printf("disp_left_ sum: %f\n", sum_float_array(sgm.disp_left, MAX_IMG_SIZE));
    // Write: disp_left

    if (sgm.option.is_check_lr) {
		/* compute disparity from other direction for LR-check */
        // Used: img_right, cost_init
        ComputeDisparity(sgm.cost_aggr, sgm.disp_right, 1);
		printf("disp_right_ sum: %f\n", sum_float_array(sgm.disp_right, MAX_IMG_SIZE));
        // Write: disp_right

        // Used: disp_left, disp_right
        LRCheck(sgm.disp_left, sgm.disp_right);
		printf("disp_left_ sum: %f\n", sum_float_array(sgm.disp_left, MAX_IMG_SIZE));
        // Write: disp_left
    }

    if (sgm.option.is_remove_speckles) {
        // Used: disp_left
        RemoveSpeckles(sgm.disp_left, 1);
		printf("disp_left_ sum: %f\n", sum_float_array(sgm.disp_left, MAX_IMG_SIZE));
        // Write: disp_left
    }


    MedianFilter(sgm.disp_left, sgm.disp_left, FILTER_WINDOW_SIZE);
	printf("disp_left_ sum: %f\n", sum_float_array(sgm.disp_left, MAX_IMG_SIZE));

    memcpy(disp_left, sgm.disp_left, (size_t)sgm.height * (size_t)sgm.width * sizeof(float));

	return true;
}


bool SGM_Reset(const uint16_t width, const uint16_t height, const SGMOption* option)
{
    sgm.is_initialized = false;
    return SGM_Initialize(width, height, option);
}

static void census_transform_5x5(const uint8_t* source, uint32_t* census)
{
	if (source == NULL || census == NULL || sgm.width <= 5 || sgm.height <= 5) {
		return;
	}

	for (uint16_t i = 2; i < sgm.height - 2; i++) {
		for (uint16_t j = 2; j < sgm.width - 2; j++) {
			
			const uint8_t gray_center = source[i * sgm.width + j];
			
			uint32_t census_val = 0u;
			for (int8_t r = -2; r <= 2; r++) {
				for (int8_t c = -2; c <= 2; c++) {
					census_val <<= 1;
					const uint8_t gray = source[(size_t)(i + r) * sgm.width + j + c];
					if (gray < gray_center) {
						census_val |= 1u; // change + to |
					}
				}
			}

			census[(size_t)i * sgm.width + j] = census_val;		
		}
	}
}

static void ComputeCost(const uint32_t* census_left, const uint32_t* census_right, uint8_t* cost_init)
{
	if (!census_left || !census_right || !cost_init) return;

    for (int32_t i = 0; i < sgm.height; i++) {
        for (int32_t j = 0; j < sgm.width; j++) {
        	for (int32_t d = sgm.option.min_disparity; d < sgm.option.max_disparity; d++) {
				size_t idx = (size_t)i * sgm.width * sgm.disp_range + j * sgm.disp_range + (d - sgm.option.min_disparity);
                uint8_t cost_val;
				if (j - d < 0 || j - d >= sgm.width) {
                    cost_val = UINT8_MAX/2;
                    // continue;
                } else {
					const uint32_t census_val_l = census_left[i * sgm.width + j];
                	const uint32_t census_val_r = census_right[i * sgm.width + j - d];
                	cost_val = Hamming32(census_val_l, census_val_r);
				}
                cost_init[idx] = cost_val;
                
            }
        }
    }
}

static uint8_t Hamming32(const uint32_t x, const uint32_t y)
{
	uint32_t dist = 0, val = x ^ y;

	// Count the number of set bits
	while (val) {
		++dist;
		val &= val - 1;
	}

	return (uint8_t)dist;
}

static void CostAggregation()
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

	CostAggregate(sgm.img_left, sgm.cost_init, sgm.cost_aggr, (int8_t)1,  (int8_t)0);
    CostAggregate(sgm.img_left, sgm.cost_init, sgm.cost_aggr, (int8_t)-1, (int8_t)0);
    CostAggregate(sgm.img_left, sgm.cost_init, sgm.cost_aggr, (int8_t)0,  (int8_t)1);
    CostAggregate(sgm.img_left, sgm.cost_init, sgm.cost_aggr, (int8_t)0,  (int8_t)-1);
    CostAggregate(sgm.img_left, sgm.cost_init, sgm.cost_aggr, (int8_t)1,  (int8_t)1);
    CostAggregate(sgm.img_left, sgm.cost_init, sgm.cost_aggr, (int8_t)-1, (int8_t)-1);
    CostAggregate(sgm.img_left, sgm.cost_init, sgm.cost_aggr, (int8_t)1,  (int8_t)-1);
    CostAggregate(sgm.img_left, sgm.cost_init, sgm.cost_aggr, (int8_t)-1, (int8_t)1);
}

/* CostAggregate:
   - img_data: input image (gray)
   - cost_init: initial cost C (uint8 per disparity)
   - cost_aggr: output sum (uint16 per disparity) - we add Lr into it
   - dx,dy: direction vector
*/
static void CostAggregate(const uint8_t* img_data, const uint8_t* cost_init, uint16_t* cost_aggr, const int8_t dx, const int8_t dy)
{
	assert(dx >= -1 && dx <= 1 && dy >= -1 && dy <= 1 && (dx != 0 || dy != 0));
	const bool is_forward = (dx == 1 && dy == 0) || (dx == 0 && dy == 1) || (dx == 1 && dy == 1) || (dx == -1 && dy == 1);
	const int8_t direction = is_forward ? (int8_t)1 : (int8_t)-1;
	// for diagonal direction, need to track current row and column
	uint16_t current_row = 0;
	uint16_t current_col = 0;

	for (uint16_t i = 0u; i < (dy == 0 ? sgm.height : sgm.width); i++) { // dy == 0 -> special case: horizontal
		
		const uint8_t* cost_initpos;
		uint16_t* cost_aggrpos;
		const uint8_t* img_pos;
		if (dx != 0 && dy == 0) // special case: horizontal
		{
			cost_initpos = (is_forward) ? (cost_init + i * sgm.width * sgm.disp_range) : (cost_init + i * sgm.width * sgm.disp_range + (sgm.width - 1) * sgm.disp_range);
			cost_aggrpos = (is_forward) ? (cost_aggr + i * sgm.width * sgm.disp_range) : (cost_aggr + i * sgm.width * sgm.disp_range + (sgm.width - 1) * sgm.disp_range);
			img_pos = (is_forward) ? (img_data + i * sgm.width) : (img_data + i * sgm.width + sgm.width - 1);
		} 
		else
		{
			/* vertical or diagonal starting positions */
			cost_initpos = (is_forward) ? (cost_init + i * sgm.disp_range) : (cost_init + (sgm.height - 1) * sgm.width * sgm.disp_range + i * sgm.disp_range);
			cost_aggrpos = (is_forward) ? (cost_aggr + i * sgm.disp_range) : (cost_aggr + (sgm.height - 1) * sgm.width * sgm.disp_range + i * sgm.disp_range);
			img_pos = (is_forward) ? (img_data + i) : (img_data + (sgm.height - 1) * sgm.width + i);
		} 
		
		uint8_t gray = *img_pos;
		uint8_t gray_last = *img_pos;
		
		uint8_t cost_last_path[MAX_DISPARITY_RANGE + 2];
		for (size_t f = 0; f < MAX_DISPARITY_RANGE + 2; f++) {
			cost_last_path[f] = UINT8_MAX;
		}
		/* initialize cost_aggr with cost_init for the first pixel and copy to cost_last_path */
		// memcpy(cost_aggrpos, cost_initpos, disp_range * sizeof(uint8));
		for(uint16_t f = 0; f < sgm.disp_range; f++) {
			cost_aggrpos[f] += (uint16_t)(cost_initpos[f]);
		}
		memcpy(&cost_last_path[1], cost_initpos, sgm.disp_range * sizeof(uint8_t));
		
		uint8_t mincost_last_path = UINT8_MAX;
		for (uint8_t f = 0; f < sgm.disp_range; f++) {
			uint8_t v = cost_initpos[f];
            if (v < mincost_last_path) mincost_last_path = v;
		}

		// necessary for diagonal direction
		current_row = is_forward ? 0 : sgm.height - 1;
		current_col = i;

		for (uint16_t j = 0; j < (dy == 0 ? sgm.width : sgm.height) - 1; j ++)  // dy == 0 -> special case: horizontal
		{	
			if (dx != 0 && dy == 0) // horizontal
			{
				cost_initpos += direction * sgm.disp_range;
				cost_aggrpos += direction * sgm.disp_range;
				img_pos += direction;
			} 
			else if (dx == 0 && dy != 0) // vertical
			{
				cost_initpos += direction * sgm.width * sgm.disp_range;
				cost_aggrpos += direction * sgm.width * sgm.disp_range;
				img_pos += direction * sgm.width;
			} 
			else // diagonal
			{
				if ((is_forward && current_col == sgm.width - 1 && current_row < sgm.height - 1) || (!is_forward && current_col == sgm.width - 1 && current_row > 0)) {
					// Special case: reach right edge
					cost_initpos = cost_init + (current_row + direction) * sgm.width * sgm.disp_range;
					cost_aggrpos = cost_aggr + (current_row + direction) * sgm.width * sgm.disp_range;
					img_pos = img_data + (current_row + direction) * sgm.width;
					current_col = 0;
				}
				else if ((!is_forward && current_col == 0 && current_row > 0) || (is_forward && current_col == 0 && current_row < sgm.height - 1)) {
					// Special case: reach left edge
					cost_initpos = cost_init + (current_row + direction) * sgm.width * sgm.disp_range + (sgm.width - 1) * sgm.disp_range;
					cost_aggrpos = cost_aggr + (current_row + direction) * sgm.width * sgm.disp_range + (sgm.width - 1) * sgm.disp_range;
					img_pos = img_data + (current_row + direction) * sgm.width + (sgm.width - 1);
					current_col = sgm.width - 1;
				} 
				else if ((dx == 1 && dy == 1) || (dx == -1 && dy == -1)) // diagonal down-right or diagonal up-left
				{
					cost_initpos += direction * (sgm.width + 1) * sgm.disp_range;
					cost_aggrpos += direction * (sgm.width + 1) * sgm.disp_range;
					img_pos += direction * (sgm.width + 1);
				} 
				else if ((dx == -1 && dy == 1) || (dx == 1 && dy == -1)) // diagonal down-left or diagonal up-right
				{
					cost_initpos += direction * (sgm.width - 1) * sgm.disp_range;
					cost_aggrpos += direction * (sgm.width - 1) * sgm.disp_range;
					img_pos += direction * (sgm.width - 1);
				}
			}
			
			gray = *img_pos;
			uint8_t min_cost = UINT8_MAX;
			uint8_t last_cost = UINT8_MAX;
			uint8_t current_cost[MAX_DISPARITY_RANGE];
			for (int32_t d = 0; d < sgm.disp_range; d++) {
				// Lr(p,d) = C(p,d) + min( Lr(p-r,d), Lr(p-r,d-1) + P1, Lr(p-r,d+1) + P1, min(Lr(p-r))+P2 ) - min(Lr(p-r))
				const uint8_t  cost = cost_initpos[d];
				const uint16_t l1 = cost_last_path[d + 1];
				const uint16_t l2 = cost_last_path[d] + sgm.option.p1;
				const uint16_t l3 = cost_last_path[d + 2] + sgm.option.p1;
				const uint16_t l4 = (uint16_t)(mincost_last_path + u_int_max_int((int)(sgm.option.p1), (int)(sgm.option.p2_init) / (abs(gray - gray_last) + 1)));

				// find min of l1,l2,l3,l4
				uint16_t m = l1;
                if (l2 < m) m = l2;
                if (l3 < m) m = l3;
                if (l4 < m) m = l4;

				const uint8_t cost_s = cost + m - mincost_last_path;

				cost_aggrpos[d] = cost_aggrpos[d] + cost_s;
				current_cost[d] = cost_s;
				if (cost_s < min_cost) min_cost = cost_s;

				cost_last_path[d] = last_cost; // shift by one right
				last_cost = cost_s;
			}

			mincost_last_path = min_cost;

			// cost_aggrpos
			// memcpy(&cost_last_path[1], current_cost, disp_range * sizeof(uint8));// error: copy uint16 to uint8
			cost_last_path[sgm.disp_range] = last_cost;

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

static void ComputeDisparity(const uint16_t* cost_aggr, float* disparity, int inverse)
{
	int c0 = 0;
	int c1 = 0;
	int c2 = 0;
    uint16_t cost_local[MAX_DISPARITY_RANGE];
    
	for (uint16_t i = 0; i < sgm.height; i++) {
        for (uint16_t j = 0; j < sgm.width; j++) {
            uint16_t min_cost = UINT16_MAX;
            uint16_t sec_min_cost = UINT16_MAX;
            uint16_t best_disparity = 0;

            for (uint16_t d = sgm.option.min_disparity; d < sgm.option.max_disparity; d++) {
	            const uint16_t d_idx = d - sgm.option.min_disparity;
				uint16_t cost;
                if(!inverse) {
                    cost = cost_aggr[i * sgm.width * sgm.disp_range + j * sgm.disp_range + d_idx];
                    cost_local[d_idx] = cost;
					if(min_cost > cost) {
                        min_cost = cost;
                        best_disparity = d;
                    }
                } else {
                    // wieso das? Weil ich von der anderen Seite kommen, muss ich schauen wo vom rechten Bild es ist. 
                    // Deswegen suche ich mir quasi die richtigen Pixel zusammen.
                    const int32_t col_left = j + d; 
                    if (col_left >= 0 && col_left < sgm.width) {
                        cost = cost_aggr[i * sgm.width * sgm.disp_range + col_left * sgm.disp_range + d_idx];
                        cost_local[d_idx] = cost;
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

            if (sgm.option.is_check_unique) {
                for (uint16_t d = sgm.option.min_disparity; d < sgm.option.max_disparity; d++) {
                    if (d == best_disparity) {
                        continue;
                    }
                    const uint16_t cost = cost_local[d - sgm.option.min_disparity];
					if (cost < sec_min_cost) sec_min_cost = cost;
                }

                // (min-sec)/min < min*(1-uniquness)
                if (sec_min_cost - min_cost <= (uint16_t)(min_cost * (1 - sgm.option.uniqueness_ratio))) {
                    disparity[i * sgm.width + j] = INVALID_FLOAT;
					c0++;
                    continue;
                }
            }

            if (best_disparity == sgm.option.min_disparity || best_disparity == sgm.option.max_disparity - 1) {
                disparity[i * sgm.width + j] = INVALID_FLOAT;
				c1++;
                continue;
            }
            const uint16_t idx_1 = best_disparity - 1 - sgm.option.min_disparity;
            const uint16_t idx_2 = best_disparity + 1 - sgm.option.min_disparity;
            const int16_t cost_1 = (int16_t)cost_local[idx_1];
            const int16_t cost_2 = (int16_t)cost_local[idx_2];
  
			int16_t denom = cost_1 + cost_2 - 2 * min_cost;
			if (denom < 1) denom = 1;

            disparity[(size_t)i * sgm.width + j] = (float)(best_disparity) + (float)(cost_1 - cost_2) / (denom * 2.0f);
        }
    }
	printf("ComputeDisparity invalid count (uniqueness): %d\n", c0);
	printf("ComputeDisparity invalid count (border): %d\n", c1);
}

static void LRCheck(float* disp_left, const float* disp_right)
{
	int c = 0;
	int c2 = 0;
	int c3 = 0;
	int c4 = 0;
	float sum1 = 0;
	float sum2 = 0;
	float sum3 = 0;
    for (uint16_t i = 0; i < sgm.height; i++) {
        for (uint16_t j = 0; j < sgm.width; j++) {
        	float disp = disp_left[i * sgm.width + j];
			if(disp == INVALID_FLOAT){
				c++;
				continue;
			}

			sum1 += disp;
        	const int32_t col_right = (int32_t)(j - disp + 0.5);
        	if(col_right >= 0 && col_right < sgm.width) {
                const float disp_r = disp_right[i * sgm.width + col_right];
				if (disp_r == INVALID_FLOAT)
				{
					continue;
				}
				sum2 += disp_r;
				sum3 += fabs(disp - disp_r);
        		if (fabs(disp - disp_r) > sgm.option.lrcheck_thres) {
					disp_left[i * sgm.width + j] = INVALID_FLOAT;
					c2++;
                }
				c4++;
            }
            else {
                disp_left[i * sgm.width + j] = INVALID_FLOAT;
				c3++;
            }
        }
    }
	printf("LRCheck invalid count: %d\n", c);
	printf("LRCheck mismatch count: %d\n", c2);
	printf("LRCheck out of bound count: %d\n", c3);
	printf("LRCheck valid count: %d\n", c4);
	printf("LRCheck disp_left sum: %f\n", sum1);
	printf("LRCheck disp_right sum: %f\n", sum2);
	printf("LRCheck disp diff sum: %f\n", sum3);
	printf("check_thres: %f\n", sgm.option.lrcheck_thres);
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

/* median-of-9 using insertion-style (robust and small) */
static float median9_from9(const float in[9]) // ToDO: check
{
    float s[5];
    int filled = 0;
    for (int k = 0; k < 9; ++k) {
        float v = in[k];
        if (filled < 5) {
            /* insert into sorted s[0..filled-1] */
            int i = filled - 1;
            while (i >= 0 && s[i] > v) {
                s[i+1] = s[i];
                --i;
            }
            s[i+1] = v;
            ++filled;
        } else {
            if (v < s[4]) {
                int i = 3;
                while (i >= 0 && s[i] > v) {
                    s[i+1] = s[i];
                    --i;
                }
                s[i+1] = v;
            }
        }
    }
    return s[4];
}

static void MedianFilter(const float* in, float* out, const uint8_t wnd_size)
{
    assert(wnd_size == FILTER_WINDOW_SIZE);
	// const int32_t radius = wnd_size / 2;

	float wnd_data[FILTER_WINDOW_SIZE * FILTER_WINDOW_SIZE];
	for (int32_t i = 1; i < sgm.height - 1; i++) {
		for (int32_t j = 1; j < sgm.width - 1; j++) {
			
			/*for (int32_t r = -radius; r <= radius; r++) {
				for (int32_t c = -radius; c <= radius; c++) {
					const int32_t row = i + r;
					const int32_t col = j + c;
					if (row >= 0 && row < height && col >= 0 && col < width) {
						wnd_data.push_back(in[row * width + col]);
					}
				}
			}*/
			wnd_data[0] = in[(i - 1) * sgm.width + (j - 1)];
			wnd_data[1] = in[(i - 1) * sgm.width + j];
			wnd_data[2] = in[(i - 1) * sgm.width + (j + 1)];
			wnd_data[3] = in[i * sgm.width + (j - 1)];
			wnd_data[4] = in[i * sgm.width + j];
			wnd_data[5] = in[i * sgm.width + (j + 1)];
			wnd_data[6] = in[(i + 1) * sgm.width + (j - 1)];
			wnd_data[7] = in[(i + 1) * sgm.width + j];
			wnd_data[8] = in[(i + 1) * sgm.width + (j + 1)];

			// std::sort(wnd_data, wnd_data + 9);
			out[i * sgm.width + j] = median9_from9(wnd_data); //wnd_data[4];
		}
	}
}

/* Small list helpers */
static void list_add(uint32_t* buffer, uint32_t* size, const uint32_t value)
{
	if (*size < MAX_IMG_SIZE) {
		buffer[*size] = value;
		(*size)++;
	} else {
		// overflow
		exit(-1);
	}
}
static void list_clear(uint32_t* index, uint32_t* size)
{
	*index = 0;
	*size = 0;
}
static uint32_t list_iterate(uint32_t* buffer, const uint32_t size, uint32_t* index)
{
	if(*index >= size){
		// overflow
		exit(-1);
	}
	uint32_t value = buffer[*index];
	(*index)++;
	return value;
}
static void RemoveSpeckles(float* disparity_map, const uint8_t diff_insame)
{
	int cv = 0;
	int ci = 0;
	int c0 = 0;
	int c1 = 0;
	int c2 = 0;
	int c3 = 0;
	// std::vector<bool> visited(uint32(width*height),false);
	uint8_t visited[MAX_IMG_SIZE];
	uint32_t vec[MAX_IMG_SIZE]; // 50 = min_speckle_area
	uint32_t vec_size = 0;
	uint32_t vec_index = 0;

	memset(visited, 0, sizeof(visited));
	/*for (uint32_t i = 0; i < MAX_IMG_SIZE; i++) {
		visited[i] = 0;
	}*/
	for(uint32_t i = 0; i < sgm.height; i++) {
		for(uint32_t j = 0; j < sgm.width; j++) {
			const uint32_t p = i * sgm.width + j;
			if (disparity_map[p] == INVALID_FLOAT) {
				ci++;
				continue;
			} else {
				cv++;
			}
			if (visited[p] == 1 || disparity_map[p] == INVALID_FLOAT) {
				c0++;
				continue;
			}
			c1++;
			// std::vector<std::pair<int32_t, int32_t>> vec;
			list_clear(&vec_index, &vec_size);
			list_add(vec, &vec_size, p);
			visited[p] = 1;
			while (vec_index < vec_size) {
				c2++;
				const uint32_t pixel = list_iterate(vec, vec_size, &vec_index);
				uint16_t row = (uint16_t)(pixel / sgm.width);
				uint16_t col = (uint16_t)(pixel % sgm.width);
				const float disp_base = disparity_map[pixel];

				for(int16_t r=-1;r<=1;r++) {
					for(int16_t c=-1;c<=1;c++) {
						if(r==0&&c==0) {
							continue;
						}
						int16_t rowr = row + r;
						int16_t colc = col + c;
						if (rowr >= 0 && rowr < sgm.height && colc >= 0 && colc < sgm.width) {
							uint32_t next_pixel = rowr * sgm.width + colc;
							if(visited[next_pixel] == 0 &&
								(disparity_map[next_pixel] != INVALID_FLOAT) &&
								fabs(disparity_map[next_pixel] - disp_base) <= (float)diff_insame) {
								list_add(vec, &vec_size, next_pixel);
								visited[next_pixel] = 1;
								c3++;
							}
						}
					}
				}
			}
			vec_index = 0;
			if(vec_size < sgm.option.min_speckle_area) {
				// in this case, all (speckles) pixels are saved in vec 
				while (vec_index < vec_size) {
					const uint32_t pixel = list_iterate(vec, vec_size, &vec_index);
					disparity_map[pixel] = INVALID_FLOAT;
				}
			}
		}
	}
	printf("RemoveSpeckles skipped invalid: %d (%d) = %d\n", ci, cv, cv + ci);
	printf("RemoveSpeckles skipped count: %d\n", c0);
	printf("RemoveSpeckles regions count: %d\n", c1);
	printf("RemoveSpeckles visited pixels count: %d\n", c2);
	printf("RemoveSpeckles added pixels count: %d\n", c3);
}