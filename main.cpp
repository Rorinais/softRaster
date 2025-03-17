#include"./global/base.h"
#include"./renderPipeline/GPU.h"
#include"./renderPipeline/Raster.h"
#include"./renderPipeline/RS_WIndow.h"
#include<assert.h>
#include <chrono>
#include <thread>

struct renderData
{
	float time{0};
};
void render(renderData&data);

int APIENTRY wWinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE hPrevInstance, _In_ LPWSTR IpCmdline, _In_ int nCmdShow){
	auto Window = RS_Windows::getInstance();
	if (!Window->initRSWindow(hInstance, 800, 600)) {
		return -1;
	}

	GPU::getInstance()->initSurface(Window->getWidth(), Window->getHeight(), Window->getCanvas());

	renderData times;

	bool alive = true;
	while (alive) {
		Window->show();
		
		render(times);

		std::this_thread::sleep_for(std::chrono::seconds(1));

		times.time++;

		alive = Window->peekMessage();
	}
	return 0;
}

void render(renderData& data) {
	GPU::getInstance()->clear();

	std::vector<Vertex> vertices;
	std::vector<FragmentShaderOutPut> pixels;

	// 颜色定义 (RGBA)
	const float4 frontColor{ 1.0f, 0.0f, 0.0f, 1.0f };  // 前红
	const float4 backColor{ 0.0f, 1.0f, 0.0f, 1.0f };  // 后绿
	const float4 leftColor{ 0.0f, 0.0f, 1.0f, 1.0f };  // 左蓝
	const float4 rightColor{ 1.0f, 1.0f, 0.0f, 1.0f };  // 右黄
	const float4 topColor{ 0.0f, 1.0f, 1.0f, 1.0f };  // 上青
	const float4 bottomColor{ 1.0f, 0.0f, 1.0f, 1.0f };  // 下品红

	//// ================== 左表面 (X-) ==================
	// 三角形5
	vertices.emplace_back().setPosition(-1, -1, -1).setColor(leftColor).setVertexCoord(0, 0);
	vertices.emplace_back().setPosition(-1, -1, 1).setColor(leftColor).setVertexCoord(1, 0);
	vertices.emplace_back().setPosition(-1, 1, 1).setColor(leftColor).setVertexCoord(1, 1);
	// 三角形6
	vertices.emplace_back().setPosition(-1, -1, -1).setColor(leftColor).setVertexCoord(0, 0);
	vertices.emplace_back().setPosition(-1, 1, 1).setColor(leftColor).setVertexCoord(1, 1);
	vertices.emplace_back().setPosition(-1, 1, -1).setColor(leftColor).setVertexCoord(0, 1);

	// ================== 下表面 (Y-) ==================
	// 三角形11
	vertices.emplace_back().setPosition(-1, -1, -1).setColor(bottomColor).setVertexCoord(0, 1);
	vertices.emplace_back().setPosition(1, -1, -1).setColor(bottomColor).setVertexCoord(1, 1);
	vertices.emplace_back().setPosition(1, -1, 1).setColor(bottomColor).setVertexCoord(1, 0);
	// 三角形12
	vertices.emplace_back().setPosition(-1, -1, -1).setColor(bottomColor).setVertexCoord(0, 1);
	vertices.emplace_back().setPosition(1, -1, 1).setColor(bottomColor).setVertexCoord(1, 0);
	vertices.emplace_back().setPosition(-1, -1, 1).setColor(bottomColor).setVertexCoord(0, 0);

	// ================== 后表面 (Z-) ==================
	// 三角形3
	vertices.emplace_back().setPosition(-1, -1, -1).setColor(backColor).setVertexCoord(1, 0);
	vertices.emplace_back().setPosition(1, -1, -1).setColor(backColor).setVertexCoord(0, 0);
	vertices.emplace_back().setPosition(1, 1, -1).setColor(backColor).setVertexCoord(0, 1);
	// 三角形4
	vertices.emplace_back().setPosition(-1, -1, -1).setColor(backColor).setVertexCoord(1, 0);
	vertices.emplace_back().setPosition(1, 1, -1).setColor(backColor).setVertexCoord(0, 1);
	vertices.emplace_back().setPosition(-1, 1, -1).setColor(backColor).setVertexCoord(1, 1);

	// ================== 前表面 (Z+) ==================
	// 三角形1
	vertices.emplace_back().setPosition(-1, -1, 1).setColor(frontColor).setVertexCoord(0, 0);
	vertices.emplace_back().setPosition(1, -1, 1).setColor(frontColor).setVertexCoord(1, 0);
	vertices.emplace_back().setPosition(1, 1, 1).setColor(frontColor).setVertexCoord(1, 1);
	// 三角形2
	vertices.emplace_back().setPosition(-1, -1, 1).setColor(frontColor).setVertexCoord(0, 0);
	vertices.emplace_back().setPosition(1, 1, 1).setColor(frontColor).setVertexCoord(1, 1);
	vertices.emplace_back().setPosition(-1, 1, 1).setColor(frontColor).setVertexCoord(0, 1);

	// ================== 右表面 (X+) ==================
	// 三角形7
	vertices.emplace_back().setPosition(1, -1, 1).setColor(rightColor).setVertexCoord(0, 0);
	vertices.emplace_back().setPosition(1, -1, -1).setColor(rightColor).setVertexCoord(1, 0);
	vertices.emplace_back().setPosition(1, 1, -1).setColor(rightColor).setVertexCoord(1, 1);
	// 三角形8
	vertices.emplace_back().setPosition(1, -1, 1).setColor(rightColor).setVertexCoord(0, 0);
	vertices.emplace_back().setPosition(1, 1, -1).setColor(rightColor).setVertexCoord(1, 1);
	vertices.emplace_back().setPosition(1, 1, 1).setColor(rightColor).setVertexCoord(0, 1);

	// ================== 上表面 (Y+) ==================
	// 三角形9
	vertices.emplace_back().setPosition(-1, 1, 1).setColor(topColor).setVertexCoord(0, 0);
	vertices.emplace_back().setPosition(1, 1, 1).setColor(topColor).setVertexCoord(1, 0);
	vertices.emplace_back().setPosition(1, 1, -1).setColor(topColor).setVertexCoord(1, 1);
	// 三角形10
	vertices.emplace_back().setPosition(-1, 1, 1).setColor(topColor).setVertexCoord(0, 0);
	vertices.emplace_back().setPosition(1, 1, -1).setColor(topColor).setVertexCoord(1, 1);
	vertices.emplace_back().setPosition(-1, 1, -1).setColor(topColor).setVertexCoord(0, 1);

	Vector<float, 3> modelTrans{ 0,0,0 };
	Vector<float, 3> modelRot{ 0, 0, 0 };
	Vector<float, 3> modelScale{ 2,2,2 };
	auto modelMat = CreateModelMatrix(modelTrans, modelRot, modelScale);

	Vector<float, 3> eye{ 5,5,5 };
	Vector<float, 3> target{ 0,0,0 };
	Vector<float, 3> up{ 0,1,0 };
	auto viewMat = CreateViewMatrix(eye, target, up);

	float fov = 3.1415 / 3; // 60度
	auto projMat = CreatePerspectiveProjection(fov, 16.0f / 9.0f, 0.1f, 100.0f);

	for (uint32_t i = 0; i < vertices.size(); i += 3)
	{
		Vertex v0 = vertices[i];
		Vertex v1 = vertices[i + 1];
		Vertex v2 = vertices[i + 2];

		VertexShaderOutPut out0 = VertexShader(v0, modelMat, viewMat, projMat);
		VertexShaderOutPut out1 = VertexShader(v1, modelMat, viewMat, projMat);
		VertexShaderOutPut out2 = VertexShader(v2, modelMat, viewMat, projMat);
		Raster::rasterizeTriangle(pixels, out0, out1, out2);
	}

	for (auto p : pixels) {
		GPU::getInstance()->drawPoint(p.screenPosition[0], p.screenPosition[1], RGBA(p.pixelColor[0], p.pixelColor[1], p.pixelColor[2], p.pixelColor[3]));
	}

}