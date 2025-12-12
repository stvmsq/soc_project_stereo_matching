#include "frame_buffer.h"
#include <string.h>
#include <xil_printf.h>
#include <stdint.h>

void init_stereo_calib(StereoCalibration_t *calib) {
    // Identity Matrizen
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            calib->cam0[r][c] = (r == c) ? 1.0f : 0.0f;
            calib->cam1[r][c] = (r == c) ? 1.0f : 0.0f;
        }
    }
    calib->doffs = 0.0f;
    calib->baseline = 0.1f;
}

void init_shared_memory(SharedMemory_t *shared_memory) {
    //ToDo: reicht ersteres nicht?
    // 1. Alles auf 0 setzen (auch große Arrays)
    memset(shared_memory, 0, sizeof(SharedMemory_t));

    // 2. Stereo Calibration initialisieren
    for (int i = 0; i < 2; i++) {
        shared_memory->stereo_frame[i].status = INVALID;
        shared_memory->stereo_frame[i].frame_id = 0;
        shared_memory->stereo_frame[i].width = IMG_WIDTH;
        shared_memory->stereo_frame[i].height = IMG_HEIGHT;

        init_stereo_calib(&shared_memory->stereo_frame[i].calib);
    }

    // 3. Depth Image initialisieren
    shared_memory->last_depth_frame.status = INVALID;
    shared_memory->last_depth_frame.frame_id = 0;
    shared_memory->last_depth_frame.width = IMG_WIDTH;
    shared_memory->last_depth_frame.height = IMG_HEIGHT;

    for (int i = 0; i < IMG_HEIGHT*IMG_WIDTH; i++) {
        shared_memory->last_depth_frame.depth[i] = 0.0f;
    }

    // current_frame_id initialisieren
    shared_memory->current_frame_id = 0;
}

void print_shared_memory_info(SharedMemory_t *shared_memory) {
    for (int i=0;i<2;i++){
        xil_printf("Stereo frame %d: status=%d frame_id=%u size=%dx%d doffs=%d baseline=%f\n\r",
            i,
            shared_memory->stereo_frame[i].status,
            shared_memory->stereo_frame[i].frame_id,
            shared_memory->stereo_frame[i].width,
            shared_memory->stereo_frame[i].height,
            shared_memory->stereo_frame[i].calib.doffs,
            shared_memory->stereo_frame[i].calib.baseline
        );
    }
    xil_printf("Depth frame: status=%d frame_id=%u size=%dx%d\n\r",
        shared_memory->last_depth_frame.status,
        shared_memory->last_depth_frame.frame_id,
        shared_memory->last_depth_frame.width,
        shared_memory->last_depth_frame.height
    );
}

uint8_t next_frame_id(SharedMemory_t *shared_memory) {
    return (shared_memory->current_frame_id & 1) ^ 1;
}