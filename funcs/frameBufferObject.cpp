#include"./frameBufferObject.h"


FrameBuffer::FrameBuffer(const uint32_t& width, const uint32_t& height, void* buffer) {
	mWidth = width;
	mHeight = height;

	if (!buffer) {
		buffer = new RGBA[mWidth * mHeight];
		mExternBuffer = false;
	}
	else
	{
		mExternBuffer = true;
	}
	mColorBuffer = (RGBA*)buffer;
}

FrameBuffer::~FrameBuffer() {
	if (!mExternBuffer && mColorBuffer)
	{
		delete[] mColorBuffer;
	}
}
