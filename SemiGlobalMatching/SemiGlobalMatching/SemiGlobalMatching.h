#ifndef SEMI_GLOBAL_MATCHING_H
#define SEMI_GLOBAL_MATCHING_H

#include <stdint.h>
#include <stdbool.h>
#include <float.h>      // für INFINITY

// --------------------------------------------
// Konstante Werte
// --------------------------------------------

#define INVALID_FLOAT (INFINITY)

#define MAX_IMG_WIDTH          450
#define MAX_IMG_HEIGHT         375
#define MAX_DISPARITY_RANGE    64
#define FILTER_WINDOW_SIZE     3
#define MAX_IMG_SIZE           (MAX_IMG_WIDTH * MAX_IMG_HEIGHT)
#define MAX_DISP_IMG_SIZE      (MAX_IMG_WIDTH * MAX_IMG_HEIGHT * MAX_DISPARITY_RANGE)

// --------------------------------------------
// SGMOption Struct (C-kompatibel)
// --------------------------------------------
typedef struct {
    uint8_t   num_paths;
    uint16_t  min_disparity;
    uint16_t  max_disparity;

    bool      is_check_unique;
    float     uniqueness_ratio;

    bool      is_check_lr;
    float     lrcheck_thres;

    bool      is_remove_speckles;
    uint16_t  min_speckle_area;

    int16_t   p1;
    int16_t   p2_init;
} SGMOption;

// --------------------------------------------
// SemiGlobalMatching Struct (C-kompatibel)
// --------------------------------------------
typedef struct {
    SGMOption option;

    uint16_t width;
    uint16_t height;
    uint16_t disp_range;
    bool     is_initialized;

    const uint8_t* img_left;
    const uint8_t* img_right;
    uint32_t* census_left;
    uint32_t* census_right;
    uint8_t*  cost_init;
    uint16_t* cost_aggr;
    float*    disp_left;
    float*    disp_right;
} SemiGlobalMatching;

// --------------------------------------------
// Globale Buffer
// --------------------------------------------

extern uint32_t census_right_buffer[MAX_IMG_SIZE];
extern uint32_t census_left_buffer[MAX_IMG_SIZE];
extern uint8_t  cost_init_buffer[MAX_DISP_IMG_SIZE];
extern uint16_t cost_aggr_buffer[MAX_DISP_IMG_SIZE];
extern float    disp_left_buffer[MAX_IMG_SIZE];
extern float    disp_right_buffer[MAX_IMG_SIZE];

// --------------------------------------------
// Funktionsprototypen (C-kompatibel)
// --------------------------------------------

bool SGM_Initialize(uint16_t width, uint16_t height, const SGMOption* option);
bool SGM_Reset(uint16_t width, uint16_t height, const SGMOption* option);
bool SGM_Match(const uint8_t* img_left, const uint8_t* img_right, float* disp_left);

#endif // SEMI_GLOBAL_MATCHING_H
