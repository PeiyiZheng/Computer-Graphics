#include "image.h"
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <algorithm>
using namespace std;

////////////////////////////
// Image processing stuff //
////////////////////////////
Pixel::Pixel(const Pixel32& p)
{
}
Pixel32::Pixel32(const Pixel& p)
{
}

int Image32::AddRandomNoise(const float& noise,Image32& outputImage) const
{
	srand((unsigned)time(0));
	Image32 source;
	source = *this;

	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			float rd = (float)(rand() / (double)(RAND_MAX));
			rd = (2.0f * noise) * rd - noise;
			rd *= 255.0f;

			Pixel32& temp = source.pixel(j, i);
			temp.r = min(255.0f, max(0.0f, rd + temp.r));
			temp.g = min(255.0f, max(0.0f, rd + temp.g));
			temp.b = min(255.0f, max(0.0f, rd + temp.b));
		}
	}

	outputImage = source;

	return 1;
}

int Image32::Brighten(const float& brightness,Image32& outputImage) const
{
	Image32 source;
	source = *this;

	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			Pixel32& temp = source.pixel(j, i);
			temp.r = min(255.0f, brightness * temp.r);
			temp.g = min(255.0f, brightness * temp.g);
			temp.b = min(255.0f, brightness * temp.b);
		}
	}

	outputImage = source;

	return 1;
}

int Image32::Luminance(Image32& outputImage) const
{
	Image32 source;
	source = *this;

	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			Pixel32& temp = source.pixel(j, i);

			float lp = 0.3f * temp.r + 0.59f * temp.g + 0.11f * temp.b;
			temp.r = temp.g = temp.b = lp;
		}
	}

	outputImage = source;

	return 1;
}

int Image32::Contrast(const float& contrast,Image32& outputImage) const
{
	Image32 source;
	this->Luminance(source);

	float sum = 0.0f;
	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			Pixel32& temp = source.pixel(j, i);

			sum += temp.r;
		}
	}

	source = *this;

	sum /= (float)(source.height() * source.width());

	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			Pixel32& temp = source.pixel(j, i);

			temp.r = min(255.0f, max(0.0f, (temp.r - sum) * contrast + sum));
			temp.g = min(255.0f, max(0.0f, (temp.g - sum) * contrast + sum));
			temp.b = min(255.0f, max(0.0f, (temp.b - sum) * contrast + sum));
		}
	}

	outputImage = source;

	return 1;
}

int Image32::Saturate(const float& saturation,Image32& outputImage) const
{
	Image32 source_1;
	this->Luminance(source_1);

	Image32 source_2;
	source_2 = *this;

	for (int i = 0; i < source_1.height(); ++i) {
		for (int j = 0; j < source_1.width(); ++j) {
			Pixel32& temp_1 = source_1.pixel(j, i);
			Pixel32& temp_2 = source_2.pixel(j, i);

			temp_2.r = min(255.0f, max(0.0f, (temp_2.r - temp_1.r) * saturation + temp_1.r));
			temp_2.g = min(255.0f, max(0.0f, (temp_2.g - temp_1.r) * saturation + temp_1.r));
			temp_2.b = min(255.0f, max(0.0f, (temp_2.b - temp_1.r) * saturation + temp_1.r));
		}
	}

	outputImage = source_2;

	return 1;
}

int Image32::Quantize(const int& bits,Image32& outputImage) const
{
	if (bits <= 0) {
		return 0;
	}

	float pow_of_bits = pow(2.0f, bits);
	float interval = 255.0f / (pow_of_bits - 1.0f);

	Image32 source;
	source = *this;

	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			Pixel32& temp = source.pixel(j, i);
			temp.r = floor(((float)temp.r / 256.0f) * pow_of_bits) * interval;
			temp.g = floor(((float)temp.g / 256.0f) * pow_of_bits) * interval;
			temp.b = floor(((float)temp.b / 256.0f) * pow_of_bits) * interval;
		}
	}

	outputImage = source;

	return 1;
}

int Image32::RandomDither(const int& bits,Image32& outputImage) const
{
	if (bits < 0) {
		return 0;
	}

	Image32 source;
	source = *this;

	float pow_of_bits = pow(2.0f, bits);
	float interval = 255.0f / (pow_of_bits - 1.0f);

	srand((unsigned)time(0));

	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			float rd = (float)(rand() / (double)(RAND_MAX));
			rd = 2.0f * rd - 1.0f;

			Pixel32& temp = source.pixel(j, i);

			float r = (float)temp.r, g = (float)temp.g, b = (float)temp.b;
			r = (r / 256.0f) + (rd / pow_of_bits);
			g = (g / 256.0f) + (rd / pow_of_bits);
			b = (b / 256.0f) + (rd / pow_of_bits);

			r = min(0.999999f, max(0.0f, r));
			g = min(0.999999f, max(0.0f, g));
			b = min(0.999999f, max(0.0f, b));

			r = floor(r * pow_of_bits) * interval;
			g = floor(g * pow_of_bits) * interval;
			b = floor(b * pow_of_bits) * interval;

			temp.r = (int)r;
			temp.g = (int)g;
			temp.b = (int)b;
		}
	}

	outputImage = source;

	return 1;

}
int Image32::OrderedDither2X2(const int& bits,Image32& outputImage) const
{
	Image32 source;
	source = *this;

	float d[2][2] = {{0.2f, 0.6f}, {0.8f, 0.4f}};

	float pow_of_bits = pow(2.0f, bits);
	float interval = 255.0f / (pow_of_bits - 1.0f);

	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			Pixel32& temp = source.pixel(j, i);
			int x = j % 2, y = i % 2;

			float r = (float)temp.r, g = (float)temp.g, b = (float)temp.b, e;
			r = (r / 256.0f) * (pow_of_bits - 1.0f);
			g = (g / 256.0f) * (pow_of_bits - 1.0f);
			b = (b / 256.0f) * (pow_of_bits - 1.0f);

			e = r - floor(r);
			if (e > d[x][y]) {
				r = ceil(r) * interval;
			}
			else {
				r = floor(r) * interval;
			}

			e = g - floor(g);
			if (e > d[x][y]) {
				g = ceil(g) * interval;
			}
			else {
				g = floor(g) * interval;
			}

			e = b - floor(b);
			if (e > d[x][y]) {
				b = ceil(b) * interval;
			}
			else {
				b = floor(b) * interval;
			}

			temp.r = (int)r;
			temp.g = (int)g;
			temp.b = (int)b;
		}
	}

	outputImage = source;
	return 1;
}

int Image32::FloydSteinbergDither(const int& bits,Image32& outputImage) const
{
	Image32 source;
	source = *this;
	float pow_of_bits = pow(2.0f, bits);
	float interval = 255.0f / (pow_of_bits - 1.0f);

	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			Pixel32& temp = source.pixel(j, i);
			
			float r = (float)temp.r, g = (float)temp.g, b = (float)temp.b;
			float sum = 0.0f;


			temp.r = max(0.0f, min(255.0f, (float)(floor((r / 256.0f) * pow_of_bits) * interval)));
			temp.g = max(0.0f, min(255.0f, (float)(floor((g / 256.0f) * pow_of_bits) * interval)));
			temp.b = max(0.0f, min(255.0f, (float)(floor((b / 256.0f) * pow_of_bits) * interval)));

			r -= (float)temp.r;
			g -= (float)temp.g;
			b -= (float)temp.b;

			if (i + 1 < source.height()) {
				sum += 5.0f / 16.0f; 

				if (j - 1 >= 0) {
					sum += 3.0f / 16.0f;
				}

				if (j + 1 < source.width()) {
					sum += 1.0f / 16.0f;
				}
			}
			

			if (j + 1 < source.width()) {
				sum += 7.0f / 16.0f;
			}

			if (sum == 0.0f) {
				continue;
			}

			sum = 1.0f / sum;

			if (i + 1 < source.height()) {
				Pixel32& temp1 = source.pixel(j, i + 1);
				temp1.r = (float)temp1.r + 3.0f * sum * r / 16.0f;
				temp1.g = (float)temp1.g + 3.0f * sum * g / 16.0f;
				temp1.b = (float)temp1.b + 3.0f * sum * b / 16.0f;

				if (j - 1 >= 0) {
					Pixel32& temp2 = source.pixel(j - 1, i + 1);
					temp2.r = (float)temp2.r + 3.0f * sum * r / 16.0f;
					temp2.g = (float)temp2.g + 3.0f * sum * g / 16.0f;
					temp2.b = (float)temp2.b + 3.0f * sum * b / 16.0f;
				}

				if (j + 1 < source.width()) {
					Pixel32& temp3 = source.pixel(j + 1, i + 1);
					temp3.r = (float)temp3.r + 1.0f * sum * r / 16.0f;
					temp3.g = (float)temp3.g + 1.0f * sum * g / 16.0f;
					temp3.b = (float)temp3.b + 1.0f * sum * b / 16.0f;
				}
			}
			


			if (j + 1 < source.width()) {
				Pixel32& temp4 = source.pixel(j + 1, i);
				temp4.r = (float)temp4.r + 7.0f * sum * r / 16.0f;
				temp4.g = (float)temp4.g + 7.0f * sum * g / 16.0f;
				temp4.b = (float)temp4.b + 7.0f * sum * b / 16.0f;
			}
		}
	}

	outputImage = source;
	return 1;
}

int Image32::Blur3X3(Image32& outputImage) const
{
	Image32 source;
	source = *this;
	float mask[3][3] = {{1.0f / 16.0f, 1.0f / 8.0f, 1.0f / 16.0f}, 
	                    {1.0f / 8.0f, 1.0f / 4.0f, 1.0f / 8.0f}, 
	                    {1.0f / 16.0f, 1.0f / 8.0f, 1.0f / 16.0f}};

	Image32 blur_source;
	blur_source.setSize(this->width(), this->height());

	float sum = 0.0f, r, g, b;
	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			sum = r = g = b = 0.0f;
			Pixel32& temp = blur_source.pixel(j, i);
			for (int mi = 0, p = -1; mi < 3; ++mi, ++p) {
				for (int mj = 0, q = -1; mj < 3; ++mj, ++q) {
					if (i + p < 0 ||j + q < 0 || i + p >= this->height() || j + q >= this->width()) {
						continue;
					}

					Pixel32& temp_cov = source.pixel(j + q, i + p);
					sum += mask[mi][mj];
					r += mask[mi][mj] * temp_cov.r;
					g += mask[mi][mj] * temp_cov.g;
					b += mask[mi][mj] * temp_cov.b;
				}
			}

			sum = 1.0f / sum;
			
			temp.r = min(255.0f, r * sum);
			temp.g = min(255.0f, g * sum);
			temp.b = min(255.0f, b * sum);
		}
	}

	outputImage = blur_source;

	return 1;
}


int Image32::EdgeDetect3X3(Image32& outputImage) const
{
	Image32 source;
	source = *this;
	float mask[3][3] = {{-0.125f, -0.125f, -0.125f}, 
	                    {-0.125f, 1.0f, -0.125f}, 
	                    {-0.125f, -0.125f, -0.125f}};

	Image32 edge_source;
	edge_source.setSize(this->width(), this->height());

	float r, g, b;
	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			r = g = b = 0.0f;
			Pixel32& temp = edge_source.pixel(j, i);
			for (int mi = 0, p = -1; mi < 3; ++mi, ++p) {
				for (int mj = 0, q = -1; mj < 3; ++mj, ++q) {
					if (i + p < 0 ||j + q < 0 || i + p >= this->height() || j + q >= this->width()) {
						continue;
					}

					Pixel32& temp_cov = source.pixel(j + q, i + p);
					r += mask[mi][mj] * temp_cov.r;
					g += mask[mi][mj] * temp_cov.g;
					b += mask[mi][mj] * temp_cov.b;
				}
			}

			temp.r = r;
			temp.g = g;
			temp.b = b;
		}
	}
	
	outputImage = edge_source;

	return 1;
}
int Image32::ScaleNearest(const float& scaleFactor,Image32& outputImage) const
{
	Image32 scale_source;
	if (scaleFactor == 1.0f) {
		scale_source = *this;
		outputImage = scale_source;
		return 1;
	}

	int width = (int)(scaleFactor * (float)this->width());
	int height = (int)(scaleFactor * (float)this->height());

	scale_source.setSize(width, height);

	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			scale_source.pixel(j, i) = NearestSample((float)j / scaleFactor, (float)i / scaleFactor);
		}
	}

	outputImage = scale_source;
	return 1;
}

int Image32::ScaleBilinear(const float& scaleFactor,Image32& outputImage) const
{
	Image32 scale_source;
	if (scaleFactor == 1.0f) {
		scale_source = *this;
		outputImage = scale_source;
		return 1;
	}

	int width = (int)(scaleFactor * (float)this->width());
	int height = (int)(scaleFactor * (float)this->height());

	scale_source.setSize(width, height);

	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			scale_source.pixel(j, i) = BilinearSample((float)j / scaleFactor, (float)i / scaleFactor);
		}
	}

	outputImage = scale_source;
	return 1;
}

int Image32::ScaleGaussian(const float& scaleFactor,Image32& outputImage) const
{
	Image32 scale_source;
	if (scaleFactor == 1.0f) {
		scale_source = *this;
		outputImage = scale_source;
		return 1;
	}
	int width = (int)(scaleFactor * (float)this->width());
	int height = (int)(scaleFactor * (float)this->height());

	scale_source.setSize(width, height);

	float radius,var;
	if (scaleFactor > 1) {
		radius = 1.0f;
		var = 0.25f;
	}
	else {
		radius = 1.0f / scaleFactor;
		var = 0.25 * sqrt(radius);
	}

	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			scale_source.pixel(j, i) = GaussianSample((float)j / scaleFactor, (float)i / scaleFactor, var, radius);
		}
	}

	outputImage = scale_source;
	return 1;
}

int Image32::RotateNearest(const float& angle,Image32& outputImage) const
{
	Image32 rotate_source;
	float width = (float)this->width(), height = (float)this->height();

	float diga = floor(sqrt(width * width + height * height));

	rotate_source.setSize((int)diga, (int)diga);

	float radius = angle * PI / 180.0f;
	Pixel32 black_pixel;
	black_pixel.r = black_pixel.g = black_pixel.b = 0.0f;

	for (int i = 0; i < rotate_source.height(); ++i) {
		for (int j = 0; j < rotate_source.width(); ++j) {
			float x = j, y = i, u, v;

			x -= (float)(rotate_source.width()) * 0.5;
			y -= (float)(rotate_source.height()) * 0.5;

			u = x * cos(-radius) - y * sin(-radius);
			v = x * sin(-radius) + y * cos(-radius);

			u += width * 0.5;
			v += height * 0.5;

			if (u >= -0.5f && v >= -0.5f && u <= width - 1.0f && v <= height - 1.0f) {
				rotate_source.pixel(j, i) = NearestSample(u, v);
			}
			else {
				rotate_source.pixel(j, i) = black_pixel;
			}
		}
	}

	outputImage = rotate_source;

	return 1;
}

int Image32::RotateBilinear(const float& angle,Image32& outputImage) const
{
	Image32 rotate_source;
	float width = (float)this->width(), height = (float)this->height();

	float diga = floor(sqrt(width * width + height * height));

	rotate_source.setSize((int)diga, (int)diga);

	float radius = angle * PI / 180.0f;
	Pixel32 black_pixel;
	black_pixel.r = black_pixel.g = black_pixel.b = 0.0f;

	for (int i = 0; i < rotate_source.height(); ++i) {
		for (int j = 0; j < rotate_source.width(); ++j) {
			float x = j, y = i, u, v;

			x -= (float)(rotate_source.width()) * 0.5;
			y -= (float)(rotate_source.height()) * 0.5;

			u = x * cos(-radius) - y * sin(-radius);
			v = x * sin(-radius) + y * cos(-radius);

			u += width * 0.5;
			v += height * 0.5;

			if (u >= 0.0f && v >= 0.0f && u <= width - 1.0f && v <= height - 1.0f) {
				rotate_source.pixel(j, i) = BilinearSample(u, v);
			}
			else {
				rotate_source.pixel(j, i) = black_pixel;
			}
		}
	}

	outputImage = rotate_source;

	return 1;
}
	
int Image32::RotateGaussian(const float& angle,Image32& outputImage) const
{
	Image32 rotate_source;
	float width = (float)this->width(), height = (float)this->height();

	float diga = floor(sqrt(width * width + height * height));

	rotate_source.setSize((int)diga, (int)diga);

	float radius = angle * PI / 180.0f;
	Pixel32 black_pixel;
	black_pixel.r = black_pixel.g = black_pixel.b = 0.0f;

	for (int i = 0; i < rotate_source.height(); ++i) {
		for (int j = 0; j < rotate_source.width(); ++j) {
			float x = j, y = i, u, v;

			x -= (float)(rotate_source.width()) * 0.5;
			y -= (float)(rotate_source.height()) * 0.5;

			u = x * cos(-radius) - y * sin(-radius);
			v = x * sin(-radius) + y * cos(-radius);

			u += width * 0.5;
			v += height * 0.5;

			if (u >= 0.0f && v >= 0.0f && u <= width - 1.0f && v <= height - 1.0f) {
				rotate_source.pixel(j, i) = GaussianSample(u, v, 0.25f, 1.0f);
			}
			else {
				rotate_source.pixel(j, i) = black_pixel;
			}
		}
	}

	outputImage = rotate_source;

	return 1;
}


int Image32::SetAlpha(const Image32& matte)
{
	for (int i = 0; i < this->height(); ++i) {
		for (int j = 0; j < this->width(); ++j) {
			Pixel32& temp = this->pixel(j, i);
			temp.a = matte.pixel(j, i).r / 255.0f;
		}
	}

	return 1;
}

int Image32::Composite(const Image32& overlay,Image32& outputImage) const
{
	Image32 source;
	source = *this;

	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			Pixel32& temp1 = source.pixel(j, i);
			const Pixel32& temp2 = overlay.pixel(j, i);

			float r_result = (float)temp2.a * temp2.r + (1.0f - (float)temp2.a) * temp1.r;
			float g_result = (float)temp2.a * temp2.g + (1.0f - (float)temp2.a) * temp1.g;
			float b_result = (float)temp2.a * temp2.b + (1.0f - (float)temp2.a) * temp1.b;

			temp1.r = r_result;
			temp1.g = g_result;
			temp1.b = b_result;
		}
	}

	outputImage = source;

	return 1;
}

int Image32::CrossDissolve(const Image32& source,const Image32& destination,const float& blendWeight,Image32& outputImage)
{
	Image32 temp;
	temp.setSize(source.width(), source.height());

	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			temp.pixel(j, i).r = (1.0f - blendWeight) * (float)source.pixel(j, i).r + blendWeight * (float)destination.pixel(j, i).r;
			temp.pixel(j, i).g = (1.0f - blendWeight) * (float)source.pixel(j, i).g + blendWeight * (float)destination.pixel(j, i).g;
			temp.pixel(j, i).b = (1.0f - blendWeight) * (float)source.pixel(j, i).b + blendWeight * (float)destination.pixel(j, i).b;
		}
	}

	outputImage = temp;

	return 1;
}

int Image32::Warp(const OrientedLineSegmentPairs& olsp,Image32& outputImage) const
{
	Image32 source;
	source.setSize(this->width(), this->height());

	for (int i = 0; i < source.height(); ++i) {
		for (int j = 0; j < source.width(); ++j) {
			float src_x, src_y;
			olsp.getSourcePosition(j, i, src_x, src_y);
			src_x = max(0, min(this->width() - 1, (int)src_x));
			src_y = max(0, min(this->height() - 1, (int)src_y));
			source.pixel(j, i) = this->pixel((int)src_x, (int)src_y);
		}
	}

	outputImage = source;
	
	return 1;
}

int Image32::FunFilter(Image32& outputImage) const
{
	Image32 fun_source;
	float width = (float)this->width(), height = (float)this->height();

	float diga = floor(sqrt(width * width + height * height));

	fun_source.setSize((int)diga, (int)diga);

	Pixel32 black_pixel;
	black_pixel.r = black_pixel.g = black_pixel.b = 0.0f;

	float width_center = 0.5 * width, height_center = 0.5 * height;

	for (int i = 0; i < fun_source.height(); ++i) {
		for (int j = 0; j < fun_source.width(); ++j) {
			float u = j, v = i, x, y;

			u -= diga * 0.5;
			v -= diga * 0.5;
			float r = sqrt(u * u + v * v);
			float angle = PI * r / diga;

			x = u * cos(angle) + v * sin(angle) + width_center;
			y = -u * sin(angle) + v * cos(angle) + height_center;

			if (x >= 0.0f && y >= 0.0f && y <= width - 1.0f && x <= height - 1.0f) {
				fun_source.pixel(j, i) = NearestSample(y, x);
			}
			else {
				fun_source.pixel(j, i) = black_pixel;
			}
		}
	}

	outputImage = fun_source;

	return 1;
}
int Image32::Crop(const int& x1,const int& y1,const int& x2,const int& y2,Image32& outputImage) const
{
	if (x1 < 0 || y1 < 0 || x2 < 0 || y2 < 0) {
		return 0;
	}

	if (x1 == this->width() || y1 == this->height() || x2 == this->width() || y2 == this->height()) {
		return 0;
	}

	if (x1 > x2 || y1 > y2) {
		return 0;
	}

	Image32 source;
	source.setSize(x2 - x1 + 1, y2 - y1 + 1);

	for (int i = y1, p = 0; i <= y2; ++i, ++p) {
		for (int j = x1, q = 0; j <= x2; ++j, ++q) {
			Pixel32& temp = source.pixel(p, q);

			temp.r = (this->pixel(i, j)).r;
			temp.g = (this->pixel(i, j)).g;
			temp.b = (this->pixel(i, j)).b;
		}
	}

	outputImage = source;

	return 1;
}


Pixel32 Image32::NearestSample(const float& x,const float& y) const
{
	const int dst_x = (int)trunc(x + 0.5f), dst_y = (int)trunc(y + 0.5f);
	Pixel32 temp = this->pixel(dst_x, dst_y);
	return temp;
}

Pixel32 Image32::BilinearSample(const float& x,const float& y) const
{
	if ((int)x == this->width() - 1 || (int)y == this->height() - 1) {
		return this->NearestSample(x, y);
	}

	float x1 = floor(x), y1 = floor(y);
	float x2 = x1 + 1, y2 = y1 + 1;

	float dx = x - x1;

	const int src_x1 = (int)x1, src_x2 = (int)x2;
	const int src_y1 = (int)y1, src_y2 = (int)y2;

	float part1[3], part2[3];
	part1[0] = this->pixel(src_x1, src_y1).r * (1.0f - dx) + this->pixel(src_x2, src_y1).r * dx;
	part2[0] = this->pixel(src_x1, src_y2).r * (1.0f - dx) + this->pixel(src_x2, src_y2).r * dx;

	part1[1] = this->pixel(src_x1, src_y1).g * (1.0f - dx) + this->pixel(src_x2, src_y1).b * dx;
	part2[1] = this->pixel(src_x1, src_y2).g * (1.0f - dx) + this->pixel(src_x2, src_y2).b * dx;

	part1[2] = this->pixel(src_x1, src_y1).b * (1.0f - dx) + this->pixel(src_x2, src_y1).b * dx;
	part2[2] = this->pixel(src_x1, src_y2).b * (1.0f - dx) + this->pixel(src_x2, src_y2).b * dx;

	float dy = y - y1;

	Pixel32 temp;
	temp.r = part1[0] * (1.0f - dy) + part2[0] * dy;
	temp.g = part1[1] * (1.0f - dy) + part2[1] * dy;
	temp.b = part1[2] * (1.0f - dy) + part2[2] * dy;

	return temp;
}

Pixel32 Image32::GaussianSample(const float& x,const float& y,const float& variance,const float& radius) const
{
	float sum = 0.0f;
	float gaussian[3] = {0.0f};
	int xi = (int)x, yi = (int)y;
	int r = (int)radius;
	for (int i = -r; i <= r; ++i) {
		for (int j = -r; j <= r; ++j) {
			int xx = xi + r, yy = yi + r;
			if (xx < 0 || yy < 0 || xx >= this->width() || yy >= this->height()) {
				continue;
			}

			if (sqrt((float)(i * i + j * j)) > radius) {
				continue;
			}

			float temp = pow(2.71828f, -(i * i + j * j) / (2.0f * variance)) / (2.0f * PI * variance);
			sum += temp;

			gaussian[0] += temp * this->pixel(xx, yy).r;
			gaussian[1] += temp * this->pixel(xx, yy).g;
			gaussian[2] += temp * this->pixel(xx, yy).b;

		}
	}

	sum = 1.0f / sum;

	for (int i = 0; i < 3; ++i) {
		gaussian[i] *= sum;
	}

	Pixel32 res;
	res.r = gaussian[0];
	res.g = gaussian[1];
	res.b = gaussian[2];

	return res;
}