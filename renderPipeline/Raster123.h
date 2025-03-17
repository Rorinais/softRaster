#pragma once
#include"../global/base.h"
#include"../math/VecMath.cpp"

template <typename T>
const T& clamp(const T& value, const T& low, const T& high) {
    if (value < low) return low;
    if (value > high) return high;
    return value;
}
VertexShaderOutPut VertexShader(const Vertex& input, const Matrixf4x4& m, const Matrixf4x4& v, const Matrixf4x4& p) {
    VertexShaderOutPut output;
    Matrixf4x4 mvp = p * v * m;
    output.clipPos = mvp * input.getPosition();


    output.vertexColor = m * input.getPosition();
    return output;
}

int4 FragmentShader(const VertexShaderToRater& input) {
    // 颜色分量转换（浮点[0,1] -> 整型[0,255]）
    auto clamp_and_convert = [](float value) -> int {
        value = clamp(value, 0.0f, 1.0f);  
        return static_cast<int>(std::round(value * 255.0f)); 
        };

    return int4{
        clamp_and_convert(input.pixelColor[0]),
        clamp_and_convert(input.pixelColor[1]),
        clamp_and_convert(input.pixelColor[2]),
        clamp_and_convert(input.pixelColor[3])
    };
}


int2 NDCToScreen(const float3& ndc, int screenWidth = 800, int screenHeight = 600) {
    float clampedX = clamp(ndc[0], -1.0f, 1.0f);
    float clampedY = clamp(ndc[1], -1.0f, 1.0f);

    // 计算标准化屏幕坐标[0,1]
    float screenNormX = clampedX * 0.5f + 0.5f;
    float screenNormY = 1.0f - (clampedY * 0.5f + 0.5f);

    // 映射到实际像素坐标
    return int2{ 
        static_cast<int>(std::round(screenNormX * (screenWidth - 1))),
        static_cast<int>(std::round(screenNormY * (screenHeight - 1)))
        };
}

class Raster{
public:
    Raster(){}
    ~Raster(){}
    static void RasterLine(std::vector<screenPoint>& result, const screenPoint& v0, const screenPoint& v1);

    static void interpolantLine(screenPoint& target,const screenPoint& v0, const screenPoint& v1);

    static void rasterizeTriangle(std::vector<FragmentShaderOutPut>& result, const VertexShaderOutPut& v0, const VertexShaderOutPut& v1, const VertexShaderOutPut& v2);

    static void interpolantTriangle(VertexShaderToRater& target,
        const VertexShaderOutPut& v0, const int2& screen0,
        const VertexShaderOutPut& v1, const int2& screen1,
        const VertexShaderOutPut& v2, const int2& screen2);
};
void Raster::RasterLine(std::vector<screenPoint>& result, const screenPoint& v0, const screenPoint& v1) {
    screenPoint start = v0;
    screenPoint end = v1;

    //1.保证x方向是从小到大的
    if (start.screenVertex[0] > end.screenVertex[0])
    {
        auto temp = start;
        start = end;
        end = temp;
    }

    result.push_back(start);

    //2.保证y方向也是从大到小，如果需要翻转，必须记录
    bool flipY = false;
    if (start.screenVertex[1] > end.screenVertex[1])
    {
        start.screenVertex[1] *= -1;
        end.screenVertex[1] *= -1;

        flipY = true;
    }

    //3.保证斜率在0-1之间，如果需要调整，必须记录
    int deltaX = static_cast<int>(end.screenVertex[0] - start.screenVertex[0]);
    int deltaY = static_cast<int>(end.screenVertex[1] - start.screenVertex[1]);

    bool swapXY = false;
    if (deltaX < deltaY) {
        std::swap(start.screenVertex[0], start.screenVertex[1]);
        std::swap(end.screenVertex[0], end.screenVertex[1]);
        std::swap(deltaX, deltaY);
        swapXY = true;
    }

    //4. brensenham
    int currentX = static_cast<int>(start.screenVertex[0]);
    int currentY = static_cast<int>(start.screenVertex[1]);

    int resultX = 0;
    int resultY = 0;

    screenPoint currentPoint;
    int p = 2 * deltaY - deltaX;

    for (int i = 0; i < deltaX; i++)
    {
        if (p >= 0)
        {
            currentY += 1;
            p -= 2 * deltaX;
        }
        currentX += 1;
        p += 2 * deltaY;

        resultX = currentX;
        resultY = currentY;
        if (swapXY) {
            std::swap(resultX, resultY);
        }

        if (flipY) {
            resultY *= -1;
        }

        currentPoint.screenVertex[0] = resultX;
        currentPoint.screenVertex[1] = resultY;
        interpolantLine(start, end, currentPoint);

        result.push_back(currentPoint);
    }

}

void Raster::interpolantLine(screenPoint& target, const screenPoint& v0, const screenPoint& v1) {
    float weight = 1.0f;

    if (v0.screenVertex[0] != v1.screenVertex[0]) {
        weight = (float)(target.screenVertex[0] - v0.screenVertex[0]) / (float)(v1.screenVertex[0] - v0.screenVertex[0]);
    }
    else if (v0.screenVertex[1] != v1.screenVertex[1]) {
        weight = (float)(target.screenVertex[1] - v0.screenVertex[1]) / (float)(v1.screenVertex[1] - v0.screenVertex[1]);
    }
    int4 result;

    result[0] = static_cast<byte>(static_cast<float>(v1.screenColor[0]) * weight + (1.0f - weight) * static_cast<float>(v0.screenColor[0]));
    result[1] = static_cast<byte>(static_cast<float>(v1.screenColor[1]) * weight + (1.0f - weight) * static_cast<float>(v0.screenColor[1]));
    result[2] = static_cast<byte>(static_cast<float>(v1.screenColor[2]) * weight + (1.0f - weight) * static_cast<float>(v0.screenColor[2]));
    result[3] = static_cast<byte>(static_cast<float>(v1.screenColor[3]) * weight + (1.0f - weight) * static_cast<float>(v0.screenColor[3]));

    std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << " " << std::endl;

    target.screenColor = result;
}
//void Raster::rasterizeTriangle(std::vector<screenPoint>& results, const screenPoint& v0, const screenPoint& v1, const screenPoint& v2, bool& IsRaster)
//{
//    if (!IsRaster) {
//        int maxX = max(v0.screenVertex[0], max(v1.screenVertex[0], v2.screenVertex[0]));
//        int minX = min(v0.screenVertex[0], min(v1.screenVertex[0], v2.screenVertex[0]));
//        int maxY = max(v0.screenVertex[1], max(v1.screenVertex[1], v2.screenVertex[1]));
//        int minY = min(v0.screenVertex[1], min(v1.screenVertex[1], v2.screenVertex[1]));
//
//        int2 pv0, pv1, pv2;
//
//        screenPoint result;
//        for (int i = minX; i <= maxX; ++i) {
//            for (int j = minY; j <= maxY; ++j) {
//                pv0 = int2{ v0.screenVertex[0] - i,v0.screenVertex[1] - j };
//                pv1 = int2{ v1.screenVertex[0] - i,v1.screenVertex[1] - j };
//                pv2 = int2{ v2.screenVertex[0] - i,v2.screenVertex[1] - j };
//
//
//                auto cross1 = cross(pv0, pv1);
//                auto cross2 = cross(pv1, pv2);
//                auto cross3 = cross(pv2, pv0);
//
//                bool negativeAll = cross1 < 0 && cross2 < 0 && cross3 < 0;
//                bool positiveAll = cross1 > 0 && cross2 > 0 && cross3 > 0;
//
//                if (negativeAll || positiveAll) {
//                    result.screenVertex[0] = i;
//                    result.screenVertex[1] = j;
//                    interpolantTriangle(result, v0, v1, v2);
//
//                    fragmentShader(result);
//
//                    results.push_back(result);
//
//
//                }
//            }
//        }
//    }
//    IsRaster = true;
//
//}

//void Raster::rasterizeTriangle(std::vector<FragmentShaderOutPut>& results, const VertexShaderOutPut& v0, const VertexShaderOutPut& v1, const VertexShaderOutPut& v2, bool& IsRaster)
//{
//    float3 ndc0 = float3{ v0.clipPos[0] / v0.clipPos[3], v0.clipPos[1] / v0.clipPos[3], v0.clipPos[2] / v0.clipPos[3] };
//    float3 ndc1 = float3{ v1.clipPos[0] / v1.clipPos[3], v1.clipPos[1] / v1.clipPos[3], v1.clipPos[2] / v1.clipPos[3] };
//    float3 ndc2 = float3{ v2.clipPos[0] / v2.clipPos[3], v2.clipPos[1] / v2.clipPos[3], v2.clipPos[2] / v2.clipPos[3] };
//
//    //std::cout << ndc0[0] << " " << ndc0[1] << " " << ndc0[2] << std::endl;
//    //std::cout << ndc1[0] << " " << ndc1[1] << " " << ndc1[2] << std::endl;
//    //std::cout << ndc2[0] << " " << ndc2[1] << " " << ndc2[2] << std::endl;
//
//    // 2. 视口变换到屏幕空间
//    int2 screen0 = NDCToScreen(ndc0);
//    int2 screen1 = NDCToScreen(ndc1);
//    int2 screen2 = NDCToScreen(ndc2);
//
//    //std::cout << screen0[0] << " " << screen0[1] << std::endl;
//    //std::cout << screen1[0] << " " << screen1[1] << std::endl;
//    //std::cout << screen2[0] << " " << screen2[1] << std::endl;
//
//    int maxX = max(screen0[0], max(screen1[0], screen2[0]));
//    int minX = min(screen0[0], min(screen1[0], screen2[0]));
//    int maxY = max(screen0[1], max(screen1[1], screen2[1]));
//    int minY = min(screen0[1], min(screen1[1], screen2[1]));
//
//    int2 pv0, pv1, pv2;
//
//    FragmentShaderOutPut result1;
//    VertexShaderToRater result2;
//    for (int i = minX; i <= maxX; ++i) {
//        for (int j = minY; j <= maxY; ++j) {
//            pv0 = int2{ screen0[0] - i,screen0[1] - j };
//            pv1 = int2{ screen1[0] - i,screen1[1] - j };
//            pv2 = int2{ screen2[0] - i,screen2[1] - j };
//
//            auto cross1 = cross(pv0, pv1);
//            auto cross2 = cross(pv1, pv2);
//            auto cross3 = cross(pv2, pv0);
//
//            bool negativeAll = cross1 < 0 && cross2 < 0 && cross3 < 0;
//            bool positiveAll = cross1 > 0 && cross2 > 0 && cross3 > 0;
//
//            if (negativeAll || positiveAll) {
//                result1.screenPosition[0] = i;
//                result1.screenPosition[1] = j;
//                result2.screenPosition[0] = i;
//                result2.screenPosition[1] = j;
//                interpolantTriangle(result2, v0, v1, v2);
//
//                result1.pixelColor = FragmentShader(result2);
//
//                results.push_back(result1);
//            }
//        }
//    }
//
//}
//void Raster::interpolantTriangle(VertexShaderToRater& target, const VertexShaderOutPut& vf0, const VertexShaderOutPut& vf1, const VertexShaderOutPut& vf2)
//{
//    VertexShaderToRater v0, v1, v2;
//    v0.screenPosition[0] = vf0.clipPos[0]*800;
//    v0.screenPosition[1] = vf0.clipPos[1]*600;
//
//    v1.screenPosition[0] = vf1.clipPos[0] * 800;
//    v1.screenPosition[1] = vf1.clipPos[1] * 600;
//
//    v2.screenPosition[0] = vf2.clipPos[0] * 800;
//    v2.screenPosition[1] = vf2.clipPos[1] * 600;
//
//    v0.pixelColor[0] = vf0.vertexColor[0]*255;
//    v0.pixelColor[1] = vf0.vertexColor[1]*255;
//    v0.pixelColor[2] = vf0.vertexColor[2]*255;
//                                         
//    v1.pixelColor[0] = vf1.vertexColor[0]*255;
//    v1.pixelColor[1] = vf1.vertexColor[1]*255;
//    v1.pixelColor[2] = vf1.vertexColor[2]*255;
//                                         
//    v2.pixelColor[0] = vf2.vertexColor[0]*255;
//    v2.pixelColor[1] = vf2.vertexColor[1]*255;
//    v2.pixelColor[2] = vf2.vertexColor[2]*255;
//
//    auto e1 = int2{ v1.screenPosition[0] - v0.screenPosition[0],v1.screenPosition[1] - v0.screenPosition[1] };
//    auto e2 = int2{ v2.screenPosition[0] - v0.screenPosition[0],v2.screenPosition[1] - v0.screenPosition[1] };
//    float sumArea = abs(cross(e1, e2));
//    if (sumArea == 0) return;
//
//    auto pv0 = int2{ v0.screenPosition[0] - target.screenPosition[0],v0.screenPosition[1] - target.screenPosition[1] };
//    auto pv1 = int2{ v1.screenPosition[0] - target.screenPosition[0],v1.screenPosition[1] - target.screenPosition[1] };
//    auto pv2 = int2{ v2.screenPosition[0] - target.screenPosition[0],v2.screenPosition[1] - target.screenPosition[1] };
//
//    float v0Area = abs(cross(pv1, pv2)); // α
//    float v1Area = abs(cross(pv0, pv2)); // β
//    float v2Area = abs(cross(pv0, pv1)); // γ
//
//    float weight0 = v0Area / sumArea;
//    float weight1 = v1Area / sumArea;
//    float weight2 = v2Area / sumArea;
//
//    auto c0 = v0.pixelColor;
//    auto c1 = v1.pixelColor;
//    auto c2 = v2.pixelColor;
//
//    int4 result;
//        result[0] = static_cast<int>(c0[0] * weight0 + c1[0] * weight1 + c2[0] * weight2);
//        result[1] = static_cast<int>(c0[1] * weight0 + c1[1] * weight1 + c2[1] * weight2);
//        result[2] = static_cast<int>(c0[2] * weight0 + c1[2] * weight1 + c2[2] * weight2);
//        result[3] = static_cast<int>(c0[3] * weight0 + c1[3] * weight1 + c2[3] * weight2);
//
//    target.pixelColor = result;
//}

//void Raster::interpolantTriangle(VertexShaderOutPut& target, const VertexShaderOutPut& v0, const VertexShaderOutPut& v1, const VertexShaderOutPut& v2)
//{
//    auto e1 = int2{ v1.screenVertex[0] - v0.screenVertex[0],v1.screenVertex[1] - v0.screenVertex[1] };
//    auto e2 = int2{ v2.screenVertex[0] - v0.screenVertex[0],v2.screenVertex[1] - v0.screenVertex[1] };
//    int sumArea = abs(cross(e1, e2));
//    if (sumArea == 0) return;
//
//    auto pv0 = int2{ v0.screenVertex[0] - target.screenVertex[0],v0.screenVertex[1] - target.screenVertex[1] };
//    auto pv1 = int2{ v1.screenVertex[0] - target.screenVertex[0],v1.screenVertex[1] - target.screenVertex[1] };
//    auto pv2 = int2{ v2.screenVertex[0] - target.screenVertex[0],v2.screenVertex[1] - target.screenVertex[1] };
//
//    int v0Area = abs(cross(pv1, pv2)); // α
//    int v1Area = abs(cross(pv0, pv2)); // β
//    int v2Area = abs(cross(pv0, pv1)); // γ
//
//    float weight0 = static_cast<float>(v0Area) / sumArea;
//    float weight1 = static_cast<float>(v1Area) / sumArea;
//    float weight2 = static_cast<float>(v2Area) / sumArea;
//
//    auto c0 = v0.screenColor;
//    auto c1 = v1.screenColor;
//    auto c2 = v2.screenColor;
//
//    int4 result;
//    result[0] = static_cast<int>(c0[0] * weight0 + c1[0] * weight1 + c2[0] * weight2);
//    result[1] = static_cast<int>(c0[1] * weight0 + c1[1] * weight1 + c2[1] * weight2);
//    result[2] = static_cast<int>(c0[2] * weight0 + c1[2] * weight1 + c2[2] * weight2);
//    result[3] = static_cast<int>(c0[3] * weight0 + c1[3] * weight1 + c2[3] * weight2);
//
//    target.screenColor = result;
//}

void Raster::rasterizeTriangle(std::vector<FragmentShaderOutPut>& results,
    const VertexShaderOutPut& v0,
    const VertexShaderOutPut& v1,
    const VertexShaderOutPut& v2) {
    // 1. 透视除法得到NDC坐标
    float3 ndc0 = float3{
        v0.clipPos[0] / v0.clipPos[3],
        v0.clipPos[1] / v0.clipPos[3],
        v0.clipPos[2] / v0.clipPos[3]
    };
    float3 ndc1 = float3{
        v1.clipPos[0] / v1.clipPos[3],  // 修正分母为v1.clipPos[3]
        v1.clipPos[1] / v1.clipPos[3],  // 修正分母
        v1.clipPos[2] / v1.clipPos[3]
    };
    float3 ndc2 = float3{
        v2.clipPos[0] / v2.clipPos[3],  // 修正分母为v2.clipPos[3]
        v2.clipPos[1] / v2.clipPos[3],  // 修正分母
        v2.clipPos[2] / v2.clipPos[3]
    };

    // 2. 视口变换到屏幕空间
    int2 screen0 = NDCToScreen(ndc0);
    int2 screen1 = NDCToScreen(ndc1);
    int2 screen2 = NDCToScreen(ndc2);

    // 计算包围盒
    int maxX = max(screen0[0], max(screen1[0], screen2[0]));
    int minX = min(screen0[0], min(screen1[0], screen2[0]));
    int maxY = max(screen0[1], max(screen1[1], screen2[1]));
    int minY = min(screen0[1], min(screen1[1], screen2[1]));

    FragmentShaderOutPut result1;
    VertexShaderToRater result2;
    for (int y = minY; y <= maxY; ++y) {
        for (int x = minX; x <= maxX; ++x) {
            int2 p = { x, y };
            int2 v0p = { screen0[0] - x, screen0[1] - y };
            int2 v1p = { screen1[0] - x, screen1[1] - y };
            int2 v2p = { screen2[0] - x, screen2[1] - y };

            int cross0 = cross(v0p, v1p);
            int cross1 = cross(v1p, v2p);
            int cross2 = cross(v2p, v0p);

            // 检查是否同侧
            if ((cross0 <= 0 && cross1 <= 0 && cross2 <= 0) ||
                (cross0 >= 0 && cross1 >= 0 && cross2 >= 0)) {
                // 设置屏幕坐标
                result1.screenPosition = { x, y };
                result2.screenPosition = { x, y };

                // 透视校正插值
                interpolantTriangle(result2, v0, screen0, v1, screen1, v2, screen2);

                // 调用片段着色器并存储结果
                result1.pixelColor = FragmentShader(result2);
                results.push_back(result1);
            }
        }
    }
}

void Raster::interpolantTriangle(VertexShaderToRater& target,
    const VertexShaderOutPut& v0, const int2& screen0,
    const VertexShaderOutPut& v1, const int2& screen1,
    const VertexShaderOutPut& v2, const int2& screen2) {
    // 将int2转换为float2以提高精度
    float2 p = { static_cast<float>(target.screenPosition[0]),
             static_cast<float>(target.screenPosition[1]) };

    // 计算三角形面积（两倍面积）
    float2 a = { static_cast<float>(screen0[0]), static_cast<float>(screen0[1]) };
    float2 b = { static_cast<float>(screen1[0]), static_cast<float>(screen1[1]) };
    float2 c = { static_cast<float>(screen2[0]), static_cast<float>(screen2[1]) };

    float area = edgeFunction(a, b, c); // 整个三角形面积的两倍
    if (area == 0) return;

    // 计算各子区域面积（两倍面积）
    float w0 = edgeFunction(b, c, p);
    float w1 = edgeFunction(c, a, p);
    float w2 = edgeFunction(a, b, p);

    // 归一化权重
    float invArea = 1.0f / area;
    w0 *= invArea;
    w1 *= invArea;
    w2 *= invArea;

    // 透视校正（使用clipPos.w的倒数）
    float wa = 1.0f / v0.clipPos[3];
    float wb = 1.0f / v1.clipPos[3];
    float wc = 1.0f / v2.clipPos[3];

    float denom = w0 * wa + w1 * wb + w2 * wc;
    if (denom == 0) return;

    w0 = (w0 * wa) / denom;
    w1 = (w1 * wb) / denom;
    w2 = (w2 * wc) / denom;

    // 插值颜色（保持浮点精度）
    target.pixelColor = v0.vertexColor * w0 + v1.vertexColor * w1 + v2.vertexColor * w2;
}
