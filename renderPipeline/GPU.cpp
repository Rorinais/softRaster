#include"./GPU.h"
#include<assert.h>
GPU* GPU::mInstance = nullptr;

GPU::GPU() {}

GPU::~GPU() {
	if (mFrameBuffer)
	{
		delete mFrameBuffer;
	}
}

GPU* GPU::getInstance() {
	if (mInstance == nullptr)
	{
		mInstance = new GPU();
	}
	return mInstance;
}

void GPU::initSurface(const uint32_t& width, const uint32_t& height, void* buffer) {
	mFrameBuffer = new FrameBuffer(width, height, buffer);
}

void GPU::clear() {
	uint32_t pixelSize = mFrameBuffer->mHeight * mFrameBuffer->mWidth;
	std::fill_n(mFrameBuffer->mColorBuffer, pixelSize, RGBA(0, 0, 0, 0));
}

void GPU::drawPoint(const uint32_t& x, const uint32_t& y, const RGBA& color) {

	const bool isXValid = (x < mFrameBuffer->mWidth);
	const bool isYValid = (y < mFrameBuffer->mHeight);

	if (isXValid && isYValid) {
		const uint32_t pixelPos = y * mFrameBuffer->mWidth + x;
		assert(pixelPos < (mFrameBuffer->mWidth * mFrameBuffer->mHeight));

		mFrameBuffer->mColorBuffer[pixelPos] = color;
	}
}