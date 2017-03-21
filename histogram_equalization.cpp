/* Histogram Equalization on Color Images
*  HYPJUDY 2017/03/17
*/

#include <iostream>
#include <vector>
#include <string>
#include "CImg.h"
using namespace std;
using namespace cimg_library;

/* RGB to grayscale transformation */
CImg<unsigned int> rgb2gray(CImg<unsigned int> rgb_img) {
	CImg<unsigned int> gray_img(rgb_img._width, rgb_img._height, 1, 1, 0);
	cimg_forXY(rgb_img, x, y) {
		int r = rgb_img(x, y, 0);
		int g = rgb_img(x, y, 1);
		int b = rgb_img(x, y, 2);
		gray_img(x, y) = 0.299 * r + 0.587 * g + 0.114 * b;
	}
	return gray_img;
}

/* Calculate the histogram for input grayscale image. */
CImg<unsigned int> im_histogram(CImg<unsigned int> input_img) {
	CImg<unsigned int> histogram(256, 1, 1, 1, 0);
	cimg_forXY(input_img, x, y) ++histogram[input_img(x, y)];
	return histogram;
}

/* Applies histogram equalization on a gray scale image (intensity:0-255)
*  returning a gray scale image whose histogram is approximately flat. */
CImg<unsigned int> equalize_hist(CImg<unsigned int> input_img) {
	int L = 256; // number of grey levels used
	int w = input_img._width;
	int h = input_img._height;
	int number_of_pixels = w * h;

	double cdf[256] = { 0 };
	unsigned int equalized[256] = { 0 };
	CImg<unsigned int> histogram = im_histogram(input_img);

	int count = 0;
	cimg_forX(histogram, pos) { // calculate cdf and equalized transform
		count += histogram[pos];
		cdf[pos] = 1.0 * count / number_of_pixels;
		equalized[pos] = round(cdf[pos] * (L - 1));
	}

	CImg<unsigned int> output_img(w, h, 1, 1, 0);
	cimg_forXY(output_img, x, y) // calculate histogram equalization result
		output_img(x, y, 0) = equalized[input_img(x, y)];
	return output_img;
}

/* Input a rgb image and output a image equalized
*  by the average histogram of the three channel */
CImg<unsigned int> rgb_average_histogram(CImg<unsigned int> rgbImage) {
	int L = 256; // number of grey levels used
	int w = rgbImage._width;
	int h = rgbImage._height;
	int number_of_pixels = w * h;

	// Calculate the histogram on each channel separately
	CImg<unsigned int> redChannel = rgbImage.get_channel(0);
	CImg<unsigned int> greenChannel = rgbImage.get_channel(1);
	CImg<unsigned int> blueChannel = rgbImage.get_channel(2);

	double cdf[256] = { 0 };
	unsigned int equalized[256] = { 0 };
	CImg<double> histogram(256, 1, 1, 1, 0);

	cimg_forXY(redChannel, x, y) {
		++histogram[redChannel(x, y)];
		++histogram[greenChannel(x, y)];
		++histogram[blueChannel(x, y)];
	}

	// compute an average histogram of these three histograms
	cimg_forX(histogram, pos) histogram(pos) /= 3;

	int count = 0;
	cimg_forX(histogram, pos) { // calculate cdf and equalized transform
		count += histogram[pos];
		cdf[pos] = 1.0 * count / number_of_pixels;
		equalized[pos] = round(cdf[pos] * (L - 1));
	}

	// Apply this transformation to the R, G and B channels individually,
	// and again rebuild an RGB image from the three processed channels.
	CImg<unsigned int> output_img(w, h, 1, 3, 0);
	cimg_forXY(output_img, x, y) {
		output_img(x, y, 0) = equalized[redChannel(x, y)];
		output_img(x, y, 1) = equalized[greenChannel(x, y)];
		output_img(x, y, 2) = equalized[blueChannel(x, y)];
	}
	return output_img;
}

/* Convert the input image to the HSI color space, and then
*  perform histogram equalization on the intensity channel.
*  Convert the result back to the RGB color space */
CImg<double> HSI_hist(CImg<double> rgbImage) {
	int L = 256; // number of grey levels used
	int w = rgbImage._width;
	int h = rgbImage._height;
	int number_of_pixels = w * h;

	// Calculate the histogram on each channel separately
	CImg<double> R = rgbImage.get_channel(0) / 255;
	CImg<double> G = rgbImage.get_channel(1) / 255;
	CImg<double> B = rgbImage.get_channel(2) / 255;

	/* ---------- Convert the input image to the HSI color space --------- */
	// H component
	CImg<double> numerator(w, h, 1, 1, 0);
	CImg<double> denominator(w, h, 1, 1, 0);
	CImg<double> theta(w, h, 1, 1, 0);
	cimg_forXY(R, x, y) {
		numerator(x, y) = 0.5 * ((R(x, y) - G(x, y)) + (R(x, y) - B(x, y)));
		denominator(x, y) = sqrt((R(x, y) - G(x, y)) * (R(x, y) - G(x, y)) +
			(R(x, y) - B(x, y)) * (G(x, y) - B(x, y)));
		theta(x, y) = acos(numerator(x, y) / denominator(x, y));
	}
	CImg<double> H(theta);
	cimg_forXY(R, x, y)
		if (B(x, y) > G(x, y)) H(x, y) = 2 * cimg::PI - H(x, y);

	// saturation component
	CImg<double> num(w, h, 1, 1, 0);
	CImg<double> den(w, h, 1, 1, 0);
	CImg<double> S(w, h, 1, 1, 0);
	CImg<double> I(w, h, 1, 1, 0);
	cimg_forXY(R, x, y) {
		num(x, y) = cimg::min(R(x, y), G(x, y), B(x, y));
		den(x, y) = R(x, y) + G(x, y) + B(x, y);
		if (den(x, y) == 0) den(x, y) = 1e-10f; // epsilon, avoid 0 division
		S(x, y) = 1 - 3.0 * num(x, y) / den(x, y);
		if (S(x, y) == 0) H(x, y) = 0;

		// intensity component
		I(x, y) = (R(x, y) + G(x, y) + G(x, y)) / 3.0;
	}

	/* ------ Perform histogram equalization on the intensity channel ----- */
	CImg<unsigned int> I_int = equalize_hist(I * 255); // parameter must be integer
	cimg_forXY(R, x, y) I(x, y, 0) = I_int(x, y) * 1.0 / 255.0;

	/* --------- Convert the result back to the RGB color space ----------- */
	// new rgb channel
	CImg<double> r(w, h, 1, 1, 0);
	CImg<double> g(w, h, 1, 1, 0);
	CImg<double> b(w, h, 1, 1, 0);
	// RG sector(0<=H<120)
	cimg_forXY(H, x, y) {
		if (H(x, y) >= 0 && H(x, y) < 2 * cimg::PI / 3) {
			b(x, y) = I(x, y) * (1 - S(x, y));
			r(x, y) = I(x, y) * (1 + ((S(x, y) * cos(H(x, y))) / cos(cimg::PI / 3 - H(x, y))));
			g(x, y) = 3 * I(x, y) - (r(x, y) + b(x, y));
		}
	}
	// GB sector(120<=H<240)
	cimg_forXY(H, x, y) {
		if (H(x, y) >= 2 * cimg::PI / 3 && H(x, y) < 4 * cimg::PI / 3) {
			r(x, y) = I(x, y) * (1 - S(x, y));
			g(x, y) = I(x, y) * (1 + ((S(x, y) * cos(H(x, y) - 2 * cimg::PI / 3))
				/ cos(cimg::PI / 3 - (H(x, y) - 2 * cimg::PI / 3))));
			b(x, y) = 3 * I(x, y) - (r(x, y) + g(x, y));
		}
	}
	// BR sector(240<=H<=360)
	cimg_forXY(H, x, y) {
		if (H(x, y) >= 4 * cimg::PI / 3 && H(x, y) <= 2 * cimg::PI) {
			g(x, y) = I(x, y) * (1 - S(x, y));
			b(x, y) = I(x, y) * (1 + ((S(x, y) * cos(H(x, y) - 4 * cimg::PI / 3))
				/ cos(cimg::PI / 3 - (H(x, y) - 4 * cimg::PI / 3))));
			r(x, y) = 3 * I(x, y) - (g(x, y) + b(x, y));
		}
	}

	/* -------- Rebuild an RGB image from the three processed channels.--------- */
	CImg<double> output_img(w, h, 1, 3, 0);
	cimg_forXY(output_img, x, y) { //  0.0 <= value <= 1.0
		output_img(x, y, 0) = cimg::min(cimg::max(r(x, y), 0), 1);
		output_img(x, y, 1) = cimg::min(cimg::max(g(x, y), 0), 1);
		output_img(x, y, 2) = cimg::min(cimg::max(b(x, y), 0), 1);
	}

	return output_img;
}

int main() {
	CImg<unsigned int> origin_rgb_img;
	vector<const char*> num = { "0", "1", "2", "3", "4", "5", "6", "7",
		"8", "9", "10", "11", "12",  "13", "14", "15", "16", "17", "18",
		"19", "20", "21", "22", "23", "24", "25", "26", "27", "28" , "29" 
		, "30" , "31" , "32", "33" };
	for (int i = 1; i <= 20; ++i) {
		// load source image
		char inPath[80];
		strcpy(inPath, "images/");
		strcat(inPath, num[i]);
		strcat(inPath, ".bmp");
		origin_rgb_img.load_bmp(inPath);

		unsigned int w = origin_rgb_img._width;
		unsigned int h = origin_rgb_img._height;

		/* ============== Grayscale Image Histogram Equalization ================ */
		CImg<unsigned int> gray_img = rgb2gray(origin_rgb_img);
		CImg<unsigned int> equalized_gray_img = equalize_hist(gray_img);
		char inPathGray[80];
		strcpy(inPathGray, "images/");
		strcat(inPathGray, num[i]);
		strcat(inPathGray, "_gray.bmp");
		gray_img.save(inPathGray);

		char outPathGray[80];
		strcpy(outPathGray, "images/");
		strcat(outPathGray, num[i]);
		strcat(outPathGray, "_eq.bmp");
		equalized_gray_img.save(outPathGray);

		
		/* ================= Color Image Histogram Equalization ================= */

		/* --- Method 1: Independent histogram equalization based on color channel -- */
		CImg<unsigned int> R_equalized_img = equalize_hist(origin_rgb_img.get_channel(0));
		CImg<unsigned int> G_equalized_img = equalize_hist(origin_rgb_img.get_channel(1));
		CImg<unsigned int> B_equalized_img = equalize_hist(origin_rgb_img.get_channel(2));
		CImg<unsigned int> rebuild_rgb_img(w, h, 1, 3, 0);
		cimg_forXY(rebuild_rgb_img, x, y) {
			rebuild_rgb_img(x, y, 0) = R_equalized_img(x, y);
			rebuild_rgb_img(x, y, 1) = G_equalized_img(x, y);
			rebuild_rgb_img(x, y, 2) = B_equalized_img(x, y);
		}
		char outPath1[80];
		strcpy(outPath1, "images/");
		strcat(outPath1, num[i]);
		strcat(outPath1, "_1.bmp");
		rebuild_rgb_img.save(outPath1);

		/* -- Method2: Histogram equalization based on average value of color channel-- */
		CImg<unsigned int> RGB_hist_img = rgb_average_histogram(origin_rgb_img);
		char outPath2[80];
		strcpy(outPath2, "images/");
		strcat(outPath2, num[i]);
		strcat(outPath2, "_2.bmp");
		RGB_hist_img.save(outPath2);

		/* ---- Method3: Intensity component equalization based on HSI color space ---- */
		CImg<double> HSI_hist_img = HSI_hist(origin_rgb_img);
		// HSI_hist output ranges 0 - 1, perform normalization(0,255) before save
		char outPath3[80];
		strcpy(outPath3, "images/");
		strcat(outPath3, num[i]);
		strcat(outPath3, "_3.bmp");
		HSI_hist_img.normalize(0, 255).save(outPath3);
		
	}

	return 0;
}

