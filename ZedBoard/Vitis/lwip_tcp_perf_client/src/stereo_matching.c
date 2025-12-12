#include "stereo_matching.h"
#include "frame_buffer.h"
#include <stdint.h>
#include "xil_printf.h"


#include "SemiGlobalMatching.h"

#include <string.h>



void convert_to_gray(SharedMemory_t *shared_memory){
	shared_memory->last_depth_frame.frame_id = shared_memory->stereo_frame[shared_memory->current_frame_id].frame_id;
	shared_memory->last_depth_frame.width = shared_memory->stereo_frame[shared_memory->current_frame_id].width;
	shared_memory->last_depth_frame.height = shared_memory->stereo_frame[shared_memory->current_frame_id].height;
	
	for (uint32_t i = 0; i < IMG_HEIGHT * IMG_WIDTH; i++) {
        uint32_t gray = 
              76u  * shared_memory->stereo_frame[shared_memory->current_frame_id].left_red[i]
            + 150u * shared_memory->stereo_frame[shared_memory->current_frame_id].left_green[i] 
            + 29u  * shared_memory->stereo_frame[shared_memory->current_frame_id].left_blue[i];

        // Shift zurück nach 0..255
        uint8_t gray8 = (uint8_t)(gray >> 8);
        // in float32 abspeichern (für depth[] Buffer)
        shared_memory->last_depth_frame.depth[i] = (float)gray8;
        /*if(gray8 < 254){
            xil_printf("ERROR: wrong value: %d\r\n", gray8);
        }*/
    }
}

int init_stereo_matching(SharedMemory_t *shared_memory) {
    return 0;
}

int start_stereo_matching(SharedMemory_t *shared_memory) {

       return 0;
}