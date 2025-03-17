#include"./RS_Window.h"

RS_Windows* RS_Windows::mInstance = nullptr;

RS_Windows::RS_Windows() {}

RS_Windows::~RS_Windows() {}

RS_Windows* RS_Windows::getInstance() {
	if (mInstance == nullptr)
	{
		mInstance = new RS_Windows();
	}
	return mInstance;
}

bool RS_Windows::initRSWindow(HINSTANCE hInstance, const uint32_t& width, const uint32_t& height) {
	mWidth = width;
	mHeight = height;

	if (registerWindowClass(hInstance) == 0) {
		DWORD err = GetLastError();
		std::wcerr << L"RegisterClassExW 失败，错误码: " << err << std::endl;
		return false;
	}

	if (!createWindow(hInstance)) {
		return false;
	}

	mhDC = GetDC(mHwnd);
	mCanvasDC = CreateCompatibleDC(mhDC);

	BITMAPINFO bmpInfo{};
	bmpInfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	bmpInfo.bmiHeader.biWidth = mWidth;
	bmpInfo.bmiHeader.biHeight = mHeight;
	bmpInfo.bmiHeader.biPlanes = 1;
	bmpInfo.bmiHeader.biBitCount = 32;
	bmpInfo.bmiHeader.biCompression = BI_RGB;

	mhBmp = CreateDIBSection(mCanvasDC, &bmpInfo, DIB_RGB_COLORS, (void**)&mCanvasBuffer, 0, 0);
	SelectObject(mCanvasDC, mhBmp);
	memset(mCanvasBuffer, 0, mWidth * mHeight * 4);

	return true;
}

void RS_Windows::handleMessage(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam) {
	switch (message){
		case WM_CLOSE: {
			DestroyWindow(hwnd);
			break;
		}
		case WM_PAINT: {
			PAINTSTRUCT ps;
			HDC hdc = BeginPaint(hwnd, &ps);
			EndPaint(hwnd, &ps);
			break;
		}

		case WM_DESTROY: {
			PostQuitMessage(0);
			mAlive = false;
			break;
		}
	}
}

bool RS_Windows::peekMessage() {
	MSG msg;
	if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	return mAlive;
}


void RS_Windows::show() {
	BitBlt(mhDC, 0, 0, mWidth, mHeight, mCanvasDC, 0, 0, SRCCOPY);
}

BOOL RS_Windows::createWindow(HINSTANCE hInstance) {
	mWindowInstance = hInstance;

	// 窗口样式设置
	DWORD dwExStyle = WS_EX_APPWINDOW;
	DWORD dwStyle = WS_OVERLAPPEDWINDOW | WS_CLIPSIBLINGS | WS_CLIPCHILDREN;

	// 计算客户区大小（关键修正）
	RECT clientRect = { 0, 0, static_cast<LONG>(mWidth), static_cast<LONG>(mHeight) };
	AdjustWindowRectEx(&clientRect, dwStyle, FALSE, dwExStyle);  // 注意参数顺序

	// 计算窗口实际尺寸（包含边框、标题栏等）
	int windowWidth = clientRect.right - clientRect.left;
	int windowHeight = clientRect.bottom - clientRect.top;

	// 创建窗口（确保客户区为 mWidth x mHeight）
	mHwnd = CreateWindowW(
		mWindowClassName,
		L"SoftRasterWindow",
		dwStyle,
		900, 200,  // 窗口位置
		windowWidth, windowHeight,  // 窗口尺寸（包含非客户区）
		nullptr,
		nullptr,
		hInstance,
		nullptr
	);

	if (!mHwnd) {
		DWORD err = GetLastError();
		std::cerr << "CreateWindow failed, error: " << err << std::endl;
		return FALSE;
	}

	// 显示窗口
	ShowWindow(mHwnd, SW_SHOW);
	UpdateWindow(mHwnd);

	// 设置视口原点为左下角（关键步骤）
	HDC hdc = GetDC(mHwnd);
	SetViewportOrgEx(hdc, 0, mHeight, NULL);  // 原点在左下角
	SetMapMode(hdc, MM_TEXT);                 // 使用默认映射模式
	ReleaseDC(mHwnd, hdc);
	return TRUE;
}

LRESULT CALLBACK Wndproc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam) {
	RS_Windows::getInstance()->handleMessage(hWnd, message, wParam, lParam);
	return (DefWindowProc(hWnd, message, wParam, lParam));
}

ATOM RS_Windows::registerWindowClass(HINSTANCE hInstance) {
	WNDCLASSEXW wndClass = {};
	wndClass.cbSize = sizeof(WNDCLASSEXW);
	wndClass.style = CS_HREDRAW | CS_VREDRAW;
	wndClass.lpfnWndProc = Wndproc;
	wndClass.hInstance = hInstance;
	wndClass.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	wndClass.hCursor = LoadCursor(NULL, IDC_ARROW);
	wndClass.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
	wndClass.lpszClassName = mWindowClassName;
	wndClass.hIconSm = LoadIcon(NULL, IDI_WINLOGO);
	return RegisterClassExW(&wndClass);
}


