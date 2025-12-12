#include "frame_buffer.h"

#ifndef __STEREO_MATCHING_H_
#define __STEREO_MATCHING_H_

void convert_to_gray(SharedMemory_t *shared_memory);
int init_stereo_matching(SharedMemory_t *shared_memory) ;
int start_stereo_matching(SharedMemory_t *shared_memory) ;

#endif