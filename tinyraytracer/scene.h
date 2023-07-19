//
// Created by Vice on 2023/7/17.
//

#ifndef TINYRAYTRACER_SCENE_H
#define TINYRAYTRACER_SCENE_H
#include <graphics.h>
#include <cmath>
#include <vector>
using std::vector;

typedef vector<double> Point;
#define POINT(x, y, z) vector<double>{(x), (y), (z)} // 定义点
const int WIDTH = 640*1.5; // 画布宽度
const int HEIGHT = 640; // 画布高度
const double vWIDTH = 1.5; // 视口宽度
const double vHEIGHT = 1; // 视口高度
const double d = 1; // 视点到投影平面的距离
const vector<double> O = POINT(0.0, 0.0, 0.0); // 视点坐标

class Sphere { // 球体类
public:
    vector<double> center; // 球心坐标
    double radius; // 半径
    int color; // 颜色
    int specular; // 镜面反射系数
    double reflect; // 反射系数
    explicit Sphere(vector<double> center = vector<double> {0, 0, 0}, double radius = 1,
                    int color = EGERGB(255, 255, 255), int specular = 500, double reflect = 0.2)
            : center(std::move(center)), radius(radius), color(color), specular(specular), reflect(reflect){} // 构造函数
    bool intersect(vector<double> O, vector<double> D, double &t) const {
        // 计算光线与球体的交点
        // O: 光线起点坐标
        // D: 光线方向向量
        vector<double> CO = POINT(O[0] - center[0], O[1] - center[1], O[2] - center[2]);
        double a = D[0] * D[0] + D[1] * D[1] + D[2] * D[2];
        double b = 2 * (D[0] * CO[0] + D[1] * CO[1] + D[2] * CO[2]);
        double c = CO[0] * CO[0] + CO[1] * CO[1] + CO[2] * CO[2] - radius * radius;
        double delta = b * b - 4 * a * c;
        if (delta < 0) return false; // 无交点
        double t1 = (-b - sqrt(delta)) / (2 * a);
        double t2 = (-b + sqrt(delta)) / (2 * a);
        t = (t1 < t2) ? t1 : t2; // 取较小的t
        return true; // 有交点
    }
};

class Light { // 光源类
public:
    // 枚举光源类型
    enum class LightType {
        AMBIENT, // 环境光
        POINT, // 点光源
        DIRECTIONAL // 平行光
    };
    LightType type; // 光源类型
    vector<double> position; // 光源位置
    double intensity; // 光源强度
    explicit Light(LightType type = LightType::POINT, vector<double> position = vector<double> {0, 0, 0}, double intensity = 1)
            : type(type), position(std::move(position)), intensity(intensity) {} // 构造函数
};



// 场景中的物体和光源
const vector<Sphere> spheres = {Sphere(POINT(0, -1, 3), 1, EGERGB(255, 0, 0), 500, 0.1),
                                Sphere(POINT(2, 0, 4), 1, EGERGB(0, 0, 255), 500, 0.1),
                                Sphere(POINT(-2, 0, 4), 1, EGERGB(0, 255, 0), 10, 0.1),
                                Sphere(POINT(0, -5001, 0), 5000, EGERGB(255, 255, 0), 1000, 0)};

const vector<Light> lights = {Light(Light::LightType::AMBIENT, POINT(0, 0, 0), 0.2),
                              Light(Light::LightType::POINT, POINT(2, 1, 0), 0.6),
                              Light(Light::LightType::DIRECTIONAL, POINT(1, 4, 4), 0.2)};


#endif //TINYRAYTRACER_SCENE_H
