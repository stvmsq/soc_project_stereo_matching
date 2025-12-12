#include <stdint.h>

#ifndef __FRAME_BUFFER_H_
#define __FRAME_BUFFER_H_

#define IMG_WIDTH 1280
#define IMG_HEIGHT 720

typedef enum {
    READY,      // 0
    PROGRESS,   // 1
    DONE,       // 2
    INVALID     // 3
} FrameStatus_t;

typedef struct __attribute__((packed)) {
    float cam0[3][3];
    float cam1[3][3];
    float doffs;
    float baseline;
	float egal;
} StereoCalibration_t;

typedef struct __attribute__((packed)) {
	FrameStatus_t status;
    uint32_t frame_id;
    uint16_t width;
    uint16_t height;
	StereoCalibration_t calib;
    uint8_t left_blue[IMG_HEIGHT*IMG_WIDTH];
	uint8_t left_green[IMG_HEIGHT*IMG_WIDTH];
	uint8_t left_red[IMG_HEIGHT*IMG_WIDTH];
	uint8_t right_blue[IMG_HEIGHT*IMG_WIDTH];
	uint8_t right_green[IMG_HEIGHT*IMG_WIDTH];
	uint8_t right_red[IMG_HEIGHT*IMG_WIDTH];
} SteroPairImg_t;

typedef struct __attribute__((packed)) {
	FrameStatus_t status;
    uint32_t frame_id;
    uint16_t width;
    uint16_t height;
    float depth[IMG_HEIGHT*IMG_WIDTH];
} DepthImg_t;
// ToDo to check _Static_assert(sizeof(StereoCalibration_t) == 80, "Struct size mismatch!");

typedef struct __attribute__((packed)) {
	uint8_t current_frame_id;
	SteroPairImg_t stereo_frame[2];
	DepthImg_t last_depth_frame;
} SharedMemory_t;

void init_stereo_calib(StereoCalibration_t *calib);
void init_shared_memory(SharedMemory_t *shared_memory);
void print_shared_memory_info(SharedMemory_t *shared_memory);
uint8_t next_frame_id(SharedMemory_t *shared_memory);

#endif