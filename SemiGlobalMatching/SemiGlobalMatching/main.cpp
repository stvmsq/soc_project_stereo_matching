
#include "SemiGlobalMatching.h"

#include <cstring>

#include <iostream>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

/**
 * \brief
 * \param argv 3
 * \param argc argc[1]:��Ӱ��·�� argc[2]: ��Ӱ��·�� argc[3]: ��С�Ӳ�[��ѡ��Ĭ��0] argc[4]: ����Ӳ�[��ѡ��Ĭ��64]
 * \param eg. ..\Data\cone\im2.png ..\Data\cone\im6.png 0 64
 * \param eg. ..\Data\Reindeer\view1.png ..\Data\Reindeer\view5.png 0 128
 * \return
 */
int main(int argv, char** argc)
{

    std::string path_left = "../Data/cone/im2.png"; //argc[1];
    std::string path_right = "../Data/cone/im6.png"; //argc[2];

    int w1, h1, c1;
    int w2, h2, c2;
    // load as grayscale (force 1 channel)
    unsigned char* img_left = stbi_load(path_left.c_str(), &w1, &h1, &c1, 1);
    unsigned char* img_right = stbi_load(path_right.c_str(), &w2, &h2, &c2, 1);

    if (!img_left || !img_right) {
        std::cout << "Failed to load images" << std::endl;
        if (img_left) stbi_image_free(img_left);
        if (img_right) stbi_image_free(img_right);
        return -1;
    }
    if (w1 != w2 || h1 != h2) {
        std::cout << "Images must have same dimensions" << std::endl;
        stbi_image_free(img_left);
        stbi_image_free(img_right);
        return -1;
    }


    const uint16_t width = static_cast<uint16_t>(w1);
    const uint16_t height = static_cast<uint16_t>(h1);
    

    printf("Loading Views...Done!\n");

    SemiGlobalMatching::SGMOption sgm_option;
    sgm_option.num_paths = 8;

    sgm_option.min_disparity = argv < 4 ? 0 : atoi(argc[3]);
    sgm_option.max_disparity = argv < 5 ? 64 : atoi(argc[4]);

    sgm_option.is_check_lr = true;
    sgm_option.lrcheck_thres = 1.0f;

    sgm_option.is_check_unique = true;
    sgm_option.uniqueness_ratio = 0.99;

    sgm_option.is_remove_speckles = true;
    sgm_option.min_speckle_aera = 50;
 
    sgm_option.p1 = 10;
    sgm_option.p2_init = 150;
 
    printf("w = %d, h = %d, d = [%d,%d]\n\n", width, height, sgm_option.min_disparity, sgm_option.max_disparity);

    SemiGlobalMatching sgm;

  
	printf("SGM Initializing...\n");
 
    if (!sgm.Initialize(width, height, sgm_option)) {
        std::cout << "SGM initialization failed" << std::endl;
        stbi_image_free(img_left);
        stbi_image_free(img_right);
        return -2;
    }
 
 
	printf("SGM Matching...\n");
    static float disparity[450 * 375] = {0};
    if (!sgm.Match(img_left, img_right, disparity)) {
        std::cout << "SGM matching failed" << std::endl;
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
            if (disp != Invalid_Float) {
                min_disp = std::min(min_disp, disp);
                max_disp = std::max(max_disp, disp);
            }
        }
    }
    static unsigned char disp_mat_data[450 * 375] = {0};
    const float range = (max_disp - min_disp) != 0.0f ? (max_disp - min_disp) : 1.0f;
    for (uint16_t i = 0; i < height; i++) {
        for (uint16_t j = 0; j < width; j++) {
            const float disp = disparity[i * width + j];
            if (disp == Invalid_Float) {
                disp_mat_data[i * width + j] = 0;
            }
            else {
                disp_mat_data[i * width + j] = static_cast<unsigned char>((disp - min_disp) / range * 255.0f);
            }
        }
    }

    /*cv::imshow("�Ӳ�ͼ", disp_mat);
    cv::Mat disp_color;
    applyColorMap(disp_mat, disp_color, cv::COLORMAP_JET);
    cv::imshow("�Ӳ�ͼ-α��", disp_color);*/

    std::string disp_map_path = "../Data/cone/im2.d.png";
    stbi_write_png(disp_map_path.c_str(), width, height, 1, disp_mat_data, width);

    stbi_image_free(img_left);
    stbi_image_free(img_right);

    return 0;
}

