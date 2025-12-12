#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "SemiGlobalMatching.h"

#define STBI_NO_LINEAR
#define STBI_NO_HDR
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


int main(int argv, char** argc)
{

    const char* path_left = "../Data/cone/im2.png"; //argc[1];
    const char* path_right = "../Data/cone/im6.png"; //argc[2];

    int w1, h1, c1;
    int w2, h2, c2;
    // load as grayscale (force 1 channel)
    unsigned char* img_left = stbi_load(path_left, &w1, &h1, &c1, 1);
    unsigned char* img_right = stbi_load(path_right, &w2, &h2, &c2, 1);

    if (!img_left || !img_right) {
        printf("Failed to load images\n");
        if (img_left) stbi_image_free(img_left);
        if (img_right) stbi_image_free(img_right);
        return -1;
    }
    if (w1 != w2 || h1 != h2) {
        printf("Images must have same dimensions\n");
        stbi_image_free(img_left);
        stbi_image_free(img_right);
        return -1;
    }


    const uint16_t width = (uint16_t)w1;
    const uint16_t height = (uint16_t)h1;
    

    printf("Loading Views...Done!\n");

    SGMOption sgm_option;
    memset(&sgm_option, 0, sizeof(SGMOption));
    sgm_option.num_paths = 8;

    sgm_option.min_disparity = 0;
    sgm_option.max_disparity = 64;

    sgm_option.is_check_lr = true;
    sgm_option.lrcheck_thres = 1.0f;

    sgm_option.is_check_unique = true;
    sgm_option.uniqueness_ratio = 0.99;

    sgm_option.is_remove_speckles = true;
    sgm_option.min_speckle_area = 50;
 
    sgm_option.p1 = 10;
    sgm_option.p2_init = 150;
 
    printf("w = %d, h = %d, d = [%d,%d]\n\n", width, height, sgm_option.min_disparity, sgm_option.max_disparity);

  
	printf("SGM Initializing...\n");
 
    if (!SGM_Initialize(width, height, &sgm_option)) {
        printf("SGM initialization failed\n");
        stbi_image_free(img_left);
        stbi_image_free(img_right);
        return -2;
    }
 
 
	printf("SGM Matching...\n");
    static float disparity[450 * 375]; // = {0};
    memset(disparity, 0, sizeof(disparity));
    if (!SGM_Match(img_left, img_right, disparity)) {
        printf("SGM matching failed\n");
        stbi_image_free(img_left);
        stbi_image_free(img_right);
        return -2;
    }
 
    // Normalize and save disparity map
    printf("Saving Disparity Map...\n");
    float min_disp = width, max_disp = -width;
    for (uint16_t i = 0; i < height; i++) {
        for (uint16_t j = 0; j < width; j++) {
            const float disp = disparity[i * width + j];
            if (disp != INVALID_FLOAT) {
                if (disp < min_disp) min_disp = disp;
                if (disp > max_disp) max_disp = disp;
            }
        }
    }
    static unsigned char disp_mat_data[450 * 375] = {0};
    const float range = (max_disp - min_disp) != 0.0f ? (max_disp - min_disp) : 1.0f;
    for (uint16_t i = 0; i < height; i++) {
        for (uint16_t j = 0; j < width; j++) {
            const float disp = disparity[i * width + j];
            if (disp == INVALID_FLOAT) {
                disp_mat_data[i * width + j] = 0;
            }
            else {
                float v = (disp - min_disp) / range * 255.0f;
                if (v < 0) v = 0;
                if (v > 255) v = 255;
                disp_mat_data[i * width + j] = (unsigned char)v;
            }
        }
    }

    const char* disp_map_path = "../Data/cone/im2.d.png";
    stbi_write_png(disp_map_path, width, height, 1, disp_mat_data, width);

    stbi_image_free(img_left);
    stbi_image_free(img_right);

    return 0;
}

