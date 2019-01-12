#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <ctime>
#include <omp.h>

using namespace std;

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int main() {
    int width, height, bpp;
    unsigned char * texture = stbi_load("texture/iiis.png", &width, &height, &bpp, 3);
    if (texture != NULL) {
        printf("w: %d, h: %d, bpp: %d\n", width, height, bpp);
        stbi_write_png("testtexture.png", width, height, 3, texture, width * 3);
        stbi_image_free(texture);
    }
    return 0;
}