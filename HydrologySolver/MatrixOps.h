#pragma once

#include "DataPointer.h"
#include "Matrix.h"

#include <filesystem>
#include <fstream>
#include <iostream>


namespace MatrixOps
{
	template<int WA, int HA, int WB, int HB>
	void Multiply(
		const Matrix<WA, HA>& A,
		const Matrix<WB, HB>& B,
		Matrix<WB, HA>& C)
	{
		static_assert(WA == HB, "Matrix dimensions do not match (A.width != B.height)");

		for (int y = 0; y < HA; ++y) {
			for (int x = 0; x < WB; ++x) {
				double sum = 0.0;

				// compile-time known loop bound WA == HB
				for (int k = 0; k < WA; ++k) {
					sum += A.Get(k, y) * B.Get(x, k);
				}
				C.Set(x, y, sum);
			}
		}
	}

	template<int W, int H>
	void Transpose(
		const Matrix<W, H>& A,
		Matrix<H, W>& AT)
	{
		for (int y = 0; y < H; ++y) {
			for (int x = 0; x < W; ++x) {
				AT.Set(y, x, A.Get(x, y));
			}
		}
	}

	// calculates C = A * B where A is a ((bx * by) x (cx * cy) matrix,
	template<int W, int H, int wa, int ha, int wb, int hb>
	void Multiply(
		const Matrix<W, H>& A,
		const DataPointer<wa, ha>& B,
		DataPointer<wb, hb>& C
	)
	{
		static_assert(W == wa * ha, "Matrix dimensions do not match (A.width != B.elems)");
		static_assert(H == wb * hb, "Matrix dimensions do not match (A.height != C.elems)");
		for (int y = 0; y < H; ++y) {
			int bx = y % wb;
			int by = y / wb;
			double sum = 0.0;
			// compile-time known loop bound W == wa * ha
			for (int x = 0; x < W; ++x) {
				int ax = x % wa;
				int ay = x / wa;
				sum += A.Get(x, y) * B.Get(ax, ay);
			}
			C.Set(bx, by, sum);
		}
	}

	template<int W, int H>
	void
	Visualize(const Matrix<W, H>& A, std::filesystem::path outPath, std::string name)
	{
		std::string myString = name + ".ppm";
		outPath.replace_filename(myString);
		std::ofstream out(outPath.c_str());

		out << "P3\n" << W << " " << H << "\n255\n";
		double min = std::numeric_limits<double>::max();
		double max = std::numeric_limits<double>::min();
		for (int y = 0; y < H; ++y) {
			for (int x = 0; x < W; ++x) {
				double val = A.Get(x, y);
				if (val < min)
					min = val;
				if (val > max)
					max = val;
			}
		}
		for (int y = 0; y < H; y++) {
			for (int x = 0; x < W; x++) {
				double val = A.Get(x, y);

				int r = 0;
				int g = 0;
				int b = 0;
				if (val < 0)
					r = val / min * 255;
				else if (val > 0)
					g = val / max * 255;
				out << r << " " << g << " " << b << " ";
			}
			out << "\n";
		}
		out.close();
	}



	template<int w, int h>
	void
		Visualize(const DataPointer<w, h>& data, std::filesystem::path outPath, std::string name)
	{
		std::string myString = name + ".ppm";
		outPath.replace_filename(myString);
		std::ofstream out(outPath.c_str());

		out << "P3\n" << w << " " << h << "\n255\n";
		double min = std::numeric_limits<double>::max();
		double max = std::numeric_limits<double>::min();
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				double val = data.Get(x, y);
				if (val < min)
					min = val;
				if (val > max)
					max = val;
			}
		}
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				double val = data.Get(x, y);

				int r = 0;
				int g = 0;
				int b = 0;
				if (val < 0)
					r = val / min * 255;
				else if (val > 0)
					g = val / max * 255;
				out << r << " " << g << " " << b << " ";
			}
			out << "\n";
		}
		out.close();
	}


}