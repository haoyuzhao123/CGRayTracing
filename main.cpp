// The main cpp file for the graphics project

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

//define constants
const int width = 1024;
const int height = 1024;

char image_data[width * height * 3];

void write_image(char * filename) {
    stbi_write_png(filename, width, height, 3, image_data, width * 3);
}

int main()
{
    memset(image_data, 0, width * height * 3 / 2);
    memset(image_data + width * height * 3 / 2, 255, width * height * 3 / 2);
    write_image("test.png");
    return 0;
}
