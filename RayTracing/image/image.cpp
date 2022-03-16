#include"image.h"
#define STB_IMAGE_IMPLEMENTATION
#include<stb_image/stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include<stb_image/stb_image_write.h>
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include<stb_image/stb_image_resize.h>

namespace schwi {
	Color Color::operator*(const double d) {
		Color ret=*this;
		for (int i = 0; i < comp; i++) {
			ret[i] *= d;
		}
		return ret;
	}

	bool Image::load_file(const std::string& filename,bool flip) {
		stbi_set_flip_vertically_on_load(flip);
		BYTE* data = stbi_load(filename.c_str(), &w, &h, &comp, 0);
		if (data) {
			pixels = std::vector<BYTE>(data, data + w * comp * h);
			stbi_image_free(data);
			return true;
		}
		else {
			std::cerr << "Error: fail to open image" << std::endl;
			return false;
		}
	}

	void Image::write_file(const std::string& filename,bool flip) const{
		stbi_flip_vertically_on_write(flip);
		stbi_write_png(filename.c_str(), w, h, comp, pixels.data(), 0);
	}

	void Image::setColor(const int x, const int y,const Color& color) {
		if (pixels.empty() || x < 0 || y < 0 || x >= w || y >= h) return;
		int ncmp = color.comp;
		int f = y * w * comp + x * comp;
		if (ncmp < comp) {
			for (int i = ncmp; i < comp; i++) {
				pixels[f + i] = 255;
			}
		}
		for (int i = 0; i < ncmp; i++) {
			pixels[f + i] = color[i];
		}
	}

	Color Image::getColor(const int x, const int y) const {
		if (pixels.empty() || x < 0 || y < 0 || x >= w || y >= h)
			return {};
		int flag = y * w * comp + x * comp;
		return Color(pixels.data() + flag, comp);
	}
}