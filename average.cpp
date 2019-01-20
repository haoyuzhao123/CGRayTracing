#include <cstdio>
#include <vector>
#include <cmath>

using namespace std;

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

unsigned char imgdata[3000000];

int main() {
    // if we should mark the deleted lines
    bool flag = false;
    //read image
    int w, h, bpp;
    unsigned char * img = stbi_load("result/t1.png", &w, &h, &bpp, 3);
    printf("w: %d, h: %d\n", w, h);
    for (int i = 0; i < w * h * 3; i++) {
        imgdata[i] = img[i] / 9;
    }
    stbi_image_free(img);
    img = stbi_load("result/t2.png", &w, &h, &bpp, 3);
    for (int i = 0; i < w * h * 3; i++) {
        imgdata[i] += img[i] / 9;
    }
    stbi_image_free(img);
    img = stbi_load("result/t3.png", &w, &h, &bpp, 3);
    for (int i = 0; i < w * h * 3; i++) {
        imgdata[i] += img[i] / 9;
    }
    stbi_image_free(img);
    img = stbi_load("result/t4.png", &w, &h, &bpp, 3);
    for (int i = 0; i < w * h * 3; i++) {
        imgdata[i] += img[i] / 9;
    }
    stbi_image_free(img);
    img = stbi_load("result/t5.png", &w, &h, &bpp, 3);
    for (int i = 0; i < w * h * 3; i++) {
        imgdata[i] += img[i] / 9;
    }
    stbi_image_free(img);
    img = stbi_load("result/t6.png", &w, &h, &bpp, 3);
    for (int i = 0; i < w * h * 3; i++) {
        imgdata[i] += img[i] / 9;
    }
    stbi_image_free(img);
    img = stbi_load("result/t7.png", &w, &h, &bpp, 3);
    for (int i = 0; i < w * h * 3; i++) {
        imgdata[i] += img[i] / 9;
    }
    stbi_image_free(img);
    img = stbi_load("result/t8.png", &w, &h, &bpp, 3);
    for (int i = 0; i < w * h * 3; i++) {
        imgdata[i] += img[i] / 9;
    }
    stbi_image_free(img);
    img = stbi_load("result/t9.png", &w, &h, &bpp, 3);
    for (int i = 0; i < w * h * 3; i++) {
        imgdata[i] += img[i] / 9;
    }
    stbi_image_free(img);
    stbi_write_png("test123.png", w, h, 3, imgdata, w*3);
    return 0;
}