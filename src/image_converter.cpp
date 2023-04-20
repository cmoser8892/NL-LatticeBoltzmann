#include "image_converter.h"
/// private
void imageConverter::read() {
    std::filesystem::path file{path};
    if(std::filesystem::exists(file)) {
        std::ifstream bmp_file_input;
        bmp_file_input.open(path,std::ios_base::binary);
        // read file header
        bmp_file_input.read((char*)&bmp.file_header,sizeof(BMPFileHeader_t));
        // check the type
        if(bmp.file_header.file_type != 0x4D42) {
            throw std::runtime_error("Unrecognized bmp file format!");
        }
        // read the info header
        bmp_file_input.read((char*)&bmp.info_header,sizeof(BMPInfoHeader_t));

        // color header is for transparent images
        if(bmp.info_header.bit_count == 32) {
            throw std::runtime_error("transparent colors are not implemented yet");
        }
        // go to the pixel data blocks
        bmp_file_input.seekg(bmp.file_header.offset_data, bmp_file_input.beg);
        /// todo additional checks missing

        // resize data for reading
        bmp.data.resize(bmp.info_header.bit_count * bmp.info_header.width * bmp.info_header.height /8);
        // check for row padding
        if(bmp.info_header.width % 4 == 0) {
            // just read the data
            bmp_file_input.read((char*) bmp.data.data(),(long) bmp.data.size());
            bmp.file_header.file_size += static_cast<uint32_t> (bmp.data.size());
        }
        else {
            uint32_t row_stride = bmp.info_header.width * bmp.info_header.bit_count / 8;
            uint32_t new_stride = make_stride_aligned(4,row_stride);
            std::vector<uint8_t> padding_row(new_stride - row_stride);

            for (int y = 0; y < bmp.info_header.height; ++y) {
                bmp_file_input.read((char*)(bmp.data.data() + row_stride * y), row_stride);
                bmp_file_input.read((char*)padding_row.data(), padding_row.size());
            }
            bmp.file_header.file_size += static_cast<uint32_t>(bmp.data.size())
                                         + bmp.info_header.height
                                               * static_cast<uint32_t>(padding_row.size());
        }
    }
    else {
        throw std::invalid_argument("bmp file doesnt exist!!");
    }
}

uint32_t imageConverter::make_stride_aligned(uint32_t align_stride, uint32_t row_stride) {
    uint32_t new_stride = row_stride;
    while (new_stride % align_stride != 0) {
        new_stride++;
    }
    return new_stride;
}

/// public
imageConverter::imageConverter(std::filesystem::path p) {
    path = p;
}

void imageConverter::map_colours_to_boundaries() {

}

void imageConverter::communicate_colour_decision() {

}

void imageConverter::init() {

}