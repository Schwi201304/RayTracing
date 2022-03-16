#ifndef SCHWI_IMAGE_H
#define SCHWI_IMAGE_H

#include<fstream>
#include<iostream>
#include<string>
#include<vector>

namespace schwi {
	typedef unsigned char BYTE;

	class Color {
	public:
		BYTE rgba[4] = { 0,0,0,0 };
		BYTE comp;

		Color() :rgba{ 0,0,0,0 }, comp(3){};//default
		Color(BYTE r) :rgba{ r,255,255,255 }, comp(1){}//Grayscale
		Color(BYTE r, BYTE g, BYTE b) :rgba{ r,g,b,255 }, comp(3) {}//rgb
		Color(BYTE r, BYTE g, BYTE b,BYTE a) :rgba{ r,g,b,a }, comp(4) {}//rgba
		Color(const BYTE* color, int _comp) :comp(_comp) {
			for (int i = _comp; i--; rgba[i] = color[i]);
		}

		BYTE& operator[](const int i) { return rgba[i]; }
		const BYTE& operator[](const int i)const { return rgba[i]; }
		Color operator*(const double d);
	};

	class Image {
	private:
		int w, h, comp;
		std::vector<BYTE> pixels;

	public:
		enum Format { GRAY = 1, RGB = 3, RGBA = 4 };
		Image() = default;
		Image(const std::string& filename) :w(0), h(0), comp(0) {
			load_file(filename);
		}
		Image(int _w, int _h, int _comp = 3) :w(_w), h(_h), comp(_comp) {
			pixels = std::vector<BYTE>(w * comp * h);
		}

		inline const BYTE* data() { return pixels.data(); }
		inline const int width()const { return w; }
		inline const int height()const { return h; }
		inline const int channels() const { return comp; }

		bool load_file(const std::string& filename, bool flip = false);
		void write_file(const std::string& filename, bool flip = false)const;
		void setColor(const int x, const int y,const Color& color);
		Color getColor(const int x, const int y) const;
	};
}
#endif
