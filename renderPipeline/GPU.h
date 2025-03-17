#pragma once
#include"../global/base.h"
#include"../funcs/frameBufferObject.h"

class GPU {
public:
	static GPU* getInstance();

	void initSurface(const uint32_t& width, const uint32_t& height,void*buffer=nullptr);
	
	void clear();

	void drawPoint(const uint32_t& x, const uint32_t& y, const RGBA& color);

	GPU(GPU& Gpu) = delete;

	GPU();

	~GPU();

private:
	static GPU* mInstance;

	FrameBuffer* mFrameBuffer{ nullptr };
};





