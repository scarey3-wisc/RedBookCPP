#pragma once

#include <memory>
#include "DataPointer.h"

template<int w, int h>
class FiveStencilCSRMatrix
{
public:
    static constexpr int matrixSize = w * h;
	static constexpr int nonZeros = matrixSize * 5 - 2 * (w + h); // Each row has 5 entries, except edges and corners

    int* GetRowOffsetsPtr() { return mRowOffsets; }
    int* GetColumnIndicesPtr() { return mColumnIndices; }
    double* GetValuesPtr() { return mValues; }

    // Disable copying:
    FiveStencilCSRMatrix(const FiveStencilCSRMatrix&) = delete;
    FiveStencilCSRMatrix& operator=(const FiveStencilCSRMatrix&) = delete;

    //-----------------------------------------------------------------------------

    void LoadValues(
        DataPointer<w, h>& cc, 
        DataPointer<w, h>& xp, 
        DataPointer<w, h>& xm, 
        DataPointer<w, h>& yp,
        DataPointer<w, h>& ym)
    {
        int index = 0;
        for (int j = 0; j < h; j++)
        {
            for (int i = 0; i < w; i++)
            {
                if (j > 0)
                    mValues[index++] = ym.Get(i, j); // Up
                if (i > 0)
                    mValues[index++] = xm.Get(i, j); // Left
                mValues[index++] = cc.Get(i, j); // Center
                if (i < w - 1)
                    mValues[index++] = xp.Get(i, j); // Right
                if (j < h - 1)
                    mValues[index++] = yp.Get(i, j); // Down
            }
        }
    }

    //-----------------------------------------------------------------------------

    static FiveStencilCSRMatrix AllocateNewPointer()
    {
        return FiveStencilCSRMatrix();
    }

    //-----------------------------------------------------------------------------

    void Deallocate()
    {
        delete mValues;
        delete mColumnIndices;
        delete mRowOffsets;
    }
private:
    FiveStencilCSRMatrix()
    {
        mValues = new double[nonZeros];
		mColumnIndices = new int[nonZeros];
		mRowOffsets = new int[matrixSize + 1];

        int count = 0;
        for (int i = 0; i <= matrixSize; i++)
        {
            mRowOffsets[i] = count;
            if (i == matrixSize)
                break;
			int y = i / w;
			int x = i % w;
            count++;
            if (x > 0)
				count++;
            if (x < w - 1)
                count++;
            if (y > 0)
				count++;
			if (y < h - 1)
				count++;
        }
        int index = 0;
        for (int j = 0; j < h; j++)
        {
            for (int i = 0; i < w; i++)
            {
                int colIndex = i + j * w;
                if (j > 0)
                    mColumnIndices[index++] = colIndex - w; // Up
                if (i > 0)
                    mColumnIndices[index++] = colIndex - 1; // Left
                mColumnIndices[index++] = colIndex; // Center
                if (i < w - 1)
                    mColumnIndices[index++] = colIndex + 1; // Right
                if (j < h - 1)
					mColumnIndices[index++] = colIndex + w; // Down
            }
        }
    }
private:
    double* mValues;
    int* mColumnIndices;
    int* mRowOffsets;
};