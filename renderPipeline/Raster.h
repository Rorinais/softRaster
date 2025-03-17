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

void FragmentShader(VertexShaderToRater& input) {


}

float2 NDCToScreen(const float3& ndc, float screenWidth = 800, float screenHeight = 600) {
    float clampedX = clamp(ndc[0], -1.0f, 1.0f);
    float clampedY = clamp(ndc[1], -1.0f, 1.0f);

    // 计算标准化屏幕坐标[0,1]
    float screenNormX = clampedX * 0.5f + 0.5f;
    float screenNormY = 1.0f - (clampedY * 0.5f + 0.5f);

    // 映射到实际像素坐标
    return float2{ 
        screenNormX * (screenWidth - 1),
        screenNormY * (screenHeight - 1)};
}

class Raster{
public:
    Raster(){}
    ~Raster(){}
    static void RasterLine(std::vector<screenPoint>& result, const screenPoint& v0, const screenPoint& v1);

    static void interpolantLine(screenPoint& target,const screenPoint& v0, const screenPoint& v1);

    static void rasterizeTriangle(std::vector<FragmentShaderOutPut>& results,
        const VertexShaderOutPut& v0,
        const VertexShaderOutPut& v1,
        const VertexShaderOutPut& v2,
        float screenWidth, float screenHeight);

    static void interpolantTriangle(VertexShaderToRater& target,
        const VertexShaderOutPut& v0, const float2& screen0,
        const VertexShaderOutPut& v1, const float2& screen1,
        const VertexShaderOutPut& v2, const float2& screen2,
        float screenWidth=800, float screenHeight=600);
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

void Raster::rasterizeTriangle(std::vector<FragmentShaderOutPut>& results,
    const VertexShaderOutPut& v0,
    const VertexShaderOutPut& v1,
    const VertexShaderOutPut& v2,
    float screenWidth=800, float screenHeight=600) {
    // 1. 透视除法得到NDC坐标
    float3 ndc0 = float3{
        v0.clipPos[0] / v0.clipPos[3],
        v0.clipPos[1] / v0.clipPos[3],
        v0.clipPos[2] / v0.clipPos[3]
    };
    float3 ndc1 = float3{
        v1.clipPos[0] / v1.clipPos[3],  
        v1.clipPos[1] / v1.clipPos[3],  
        v1.clipPos[2] / v1.clipPos[3]
    };
    float3 ndc2 = float3{
        v2.clipPos[0] / v2.clipPos[3], 
        v2.clipPos[1] / v2.clipPos[3], 
        v2.clipPos[2] / v2.clipPos[3]
    };

    // 2. 视口变换到屏幕空间
    float2 screen0 = NDCToScreen(ndc0);
    float2 screen1 = NDCToScreen(ndc1);
    float2 screen2 = NDCToScreen(ndc2);

    // 计算包围盒
    float maxX = max(screen0[0], max(screen1[0], screen2[0]));
    float minX = min(screen0[0], min(screen1[0], screen2[0]));
    float maxY = max(screen0[1], max(screen1[1], screen2[1]));
    float minY = min(screen0[1], min(screen1[1], screen2[1]));

    int minX_pixel = static_cast<int>(std::floor(minX));
    int maxX_pixel = static_cast<int>(std::ceil(maxX));
    int minY_pixel = static_cast<int>(std::floor(minY));
    int maxY_pixel = static_cast<int>(std::ceil(maxY));

    FragmentShaderOutPut result1;
    VertexShaderToRater result2;

    for (int y = minY_pixel; y <= maxY_pixel; ++y) {
        for (int x = minX_pixel; x <= maxX_pixel; ++x) {
            // 使用实际像素坐标进行三角形测试
            float2 p_pixel = { static_cast<float>(x), static_cast<float>(y) };
            float2 v0p = screen0 - p_pixel;
            float2 v1p = screen1 - p_pixel;
            float2 v2p = screen2 - p_pixel;

            // 叉积计算
            float cross0 = edgeFunction(screen1, screen2, p_pixel);
            float cross1 = edgeFunction(screen2, screen0, p_pixel);
            float cross2 = edgeFunction(screen0, screen1, p_pixel);

            // 判断是否在三角形内
            if ((cross0 <= 0 && cross1 <= 0 && cross2 <= 0) ||
                (cross0 >= 0 && cross1 >= 0 && cross2 >= 0)){
                VertexShaderToRater input;
                // 转换为标准化坐标
                input.screenPosition = float2{
                    static_cast<float>(x) / (screenWidth - 1.0f),
                    1.0f - (static_cast<float>(y) / (screenHeight - 1.0f))
                };

                // 插值时传入像素坐标
                interpolantTriangle(input, v0, screen0, v1, screen1, v2, screen2);
                FragmentShader(input);

                // 转换颜色
                auto clamp_and_convert = [](float value) -> int {
                    value = clamp(value, 0.0f, 1.0f);
                    return static_cast<int>(std::round(value * 255.0f));
                    };
                FragmentShaderOutPut result;
                result.screenPosition = { x, y };
                result.pixelColor = int4{
                    clamp_and_convert(input.pixelColor[0]),
                    clamp_and_convert(input.pixelColor[1]),
                    clamp_and_convert(input.pixelColor[2]),
                    clamp_and_convert(input.pixelColor[3])
                };
                results.push_back(result);
            }
        }
    }
}

void Raster::interpolantTriangle(VertexShaderToRater& target,
    const VertexShaderOutPut& v0, const float2& screen0,
    const VertexShaderOutPut& v1, const float2& screen1,
    const VertexShaderOutPut& v2, const float2& screen2,
    float screenWidth, float screenHeight) {
    // 将标准化坐标转回像素坐标用于插值
    float2 screenPos = {
        target.screenPosition[0] * (screenWidth - 1.0f),
        (1.0f - target.screenPosition[1])* (screenHeight - 1.0f) // 反转Y轴
    };

    float area = edgeFunction(screen0, screen1, screen2);
    if (area == 0) return;

    float w0 = edgeFunction(screen1, screen2, screenPos);
    float w1 = edgeFunction(screen2, screen0, screenPos);
    float w2 = edgeFunction(screen0, screen1, screenPos);

    // 归一化权重
    float invArea = 1.0f / area;
    w0 *= invArea;
    w1 *= invArea;
    w2 *= invArea;

    // 透视校正（确保分母不为零）
    float wa = 1.0f / v0.clipPos[3];
    float wb = 1.0f / v1.clipPos[3];
    float wc = 1.0f / v2.clipPos[3];
    float denom = w0 * wa + w1 * wb + w2 * wc;
    if (denom == 0) return;

    w0 = (w0 * wa) / denom;
    w1 = (w1 * wb) / denom;
    w2 = (w2 * wc) / denom;

    // 插值颜色
    target.pixelColor = v0.vertexColor * w0 + v1.vertexColor * w1 + v2.vertexColor * w2;
}
