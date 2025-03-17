#pragma once
#include"../global/base.h"

class FrameBuffer {
public:
	FrameBuffer(const uint32_t& width, const uint32_t& height, void* buffer = nullptr);

	~FrameBuffer();

	FrameBuffer(const FrameBuffer&) = delete;

	uint32_t mWidth{ 0 };
	uint32_t mHeight{ 0 };

	RGBA* mColorBuffer{ nullptr };
	bool mExternBuffer{ false };
};