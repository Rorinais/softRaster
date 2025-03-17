#pragma once
#include"../global/base.h"

class RS_Windows {
public:
	static RS_Windows* getInstance();

	RS_Windows(RS_Windows& window) = delete;

	RS_Windows();

	~RS_Windows();

	bool initRSWindow(HINSTANCE hInstance, const uint32_t& width = 800, const uint32_t& height = 600);

	void handleMessage(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam);

	bool peekMessage();

	void show();

	uint32_t getWidth() {
		return mWidth;
	}

	uint32_t getHeight() {
		return mHeight;
	}

	void* getCanvas() {
		return mCanvasBuffer;
	}

private:
	BOOL createWindow(HINSTANCE hInstance);
	ATOM registerWindowClass(HINSTANCE hInstance);

private:
	static	RS_Windows* mInstance;

	bool	mAlive{ true };
	HINSTANCE mWindowInstance;
	WCHAR	mWindowClassName[100] = L"SoftRasterWindow";
	HWND	mHwnd;

	uint32_t mWidth = 800;
	uint32_t mHeight = 600;

	HDC	mhDC;
	HDC	mCanvasDC;
	HBITMAP mhBmp;
	void* mCanvasBuffer{ nullptr };
};
