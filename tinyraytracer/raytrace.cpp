#include "raytrace.h"

void PutPixel(int x, int y, int color) {
    putpixel(x + WIDTH/2, -y + HEIGHT/2, color); // 以画布中心为原点
}

vector<double> CanvasToViewport(int x, int y) {
    return POINT(x * vWIDTH / WIDTH, y * vHEIGHT / HEIGHT, d); // 画布坐标转视口坐标
}

// ClosestIntersection
std::tuple<double, const Sphere*> ClosestIntersection(const vector<double> &O, const vector<double> &D, const double &t_min, const double &t_max) {
    // 视线与所有球体求交，返回最近的交点
    // O: 视线起点
    // D: 视线方向
    // t_min, t_max: 视线起点到视口的距离的最小值和最大值
    double closest_t = INFINITY; // 最近交点距离
    const Sphere* closest_sphere = nullptr; // 最近交点对应的球体
    for (auto &sphere : spheres) {
        double t; // 交点距离
        if (sphere.intersect(O, D, t) && t < closest_t && t > t_min && t < t_max) {
            closest_t = t; // 更新最近交点距离
            closest_sphere = &sphere; // 更新最近交点对应的球体
        }
    }
    return std::make_tuple(closest_t, closest_sphere); // 返回最近交点
}
int TraceRay(const vector<double> &O, const vector<double> &D, const double &t_min, const double &t_max) {
    // 视线与所有球体求交，返回最近的交点颜色
    // O: 视线起点
    // D: 视线方向
    // t_min, t_max: 视线起点到视口的距离的最小值和最大值
    double closest_t; // 最近交点距离
    const Sphere* closest_sphere; // 最近交点对应的球体
    std::tie(closest_t, closest_sphere) = ClosestIntersection(O, D, t_min, t_max); // 最近交点
    if (closest_sphere == nullptr) return EGERGB(125, 125, 125); // 没有交点，返回背景色
    vector<double> P = add(O, multiple(D, closest_t)); // 交点坐标
    vector<double> N = sub(P, closest_sphere->center); // 交点法向量
    N = multiple(N, 1 / sqrt(dot(N, N))); // 归一化
    double lighting = ComputeLighting(P, N, multiple(D, -1), closest_sphere->specular); // 计算光照强度
    // multiple EGERGB and double

    return EGERGB_Mul(closest_sphere->color, lighting); // 返回交点颜色

}

int TraceRay(const vector<double> &O, const vector<double> &D, const double &t_min, const double &t_max, int recursion_depth) {
    // 视线与所有球体求交，返回最近的交点颜色，带有反射
    // O: 视线起点
    // D: 视线方向
    // t_min, t_max: 视线起点到视口的距离的最小值和最大值
    // recursion_depth: 递归深度
    double closest_t; // 最近交点距离
    const Sphere* closest_sphere; // 最近交点对应的球体
    std::tie(closest_t, closest_sphere) = ClosestIntersection(O, D, t_min, t_max); // 最近交点
    if (closest_sphere == nullptr) return EGERGB(125, 125, 125); // 没有交点，返回背景色
    vector<double> P = add(O, multiple(D, closest_t)); // 交点坐标
    vector<double> N = sub(P, closest_sphere->center); // 交点法向量
    N = multiple(N, 1 / sqrt(dot(N, N))); // 归一化
    double lighting = ComputeLighting(P, N, multiple(D, -1), closest_sphere->specular); // 计算光照强度
    // multiple EGERGB and double

    // 递归
    double reflect = closest_sphere->reflect; // 反射系数
    if (recursion_depth <= 0 or reflect <= 0) return EGERGB_Mul(closest_sphere->color, lighting); // 递归深度为0，返回交点颜色
    vector<double> R =  ReflectRay(multiple(D, -1), N); // 反射光线
    int color = EGERGB_Mul(closest_sphere->color, lighting); // 交点颜色
    int reflect_color = TraceRay(P, R, 0.05, INFINITY, recursion_depth - 1); // 反射光线颜色
    color = EGERGB_Add(EGERGB_Mul(color, 1 - reflect), EGERGB_Mul(reflect_color, reflect)); // 反射光线颜色加到交点颜色上
    return color; // 返回交点颜色
}

int EGERGB_Add(int color1, int color2) {
    // EGERGB颜色相加
    // color1, color2: EGERGB颜色
    int r1 = (color1 >> 16) & 0xFF; // 取出RGB三个分量
    int g1 = (color1 >> 8) & 0xFF;
    int b1 = color1 & 0xFF;

    int r2 = (color2 >> 16) & 0xFF;
    int g2 = (color2 >> 8) & 0xFF;
    int b2 = color2 & 0xFF;

    int r = std::min(r1 + r2, 255); // RGB三个分量相加
    int g = std::min(g1 + g2, 255);
    int b = std::min(b1 + b2, 255);

    return EGERGB(r, g, b); // 返回EGERGB颜色
}

vector<double> ReflectRay(const vector<double> &R, const vector<double> &N) {
    // 计算反射光线
    // R: 入射光线
    // N: 法向量
    return sub(multiple(N, 2 * dot(R, N)), R); // 反射光线
}

double ComputeLighting(const vector<double> &P, const vector<double> &N, const vector<double> &V, const int &s) {
    // 计算光照强度
    // P: 交点坐标
    // N: 交点法向量
    // V: 视线方向
    // s: 高光系数
    double i = 0; // 光照强度
    double t_max = INFINITY; // 最大距离
    for (auto &light : lights) {
        if (light.type == Light::LightType::AMBIENT) {
            i += light.intensity; // 环境光
        } else {
            vector<double> L; // 光源方向
            if (light.type == Light::LightType::POINT) {
                // 点光源
                L = POINT(light.position[0] - P[0], light.position[1] - P[1], light.position[2] - P[2]);
                t_max = 1;
            } else {
                // 平行光
                L = POINT(light.position[0], light.position[1], light.position[2]);
                t_max = INFINITY;
            }

            // 阴影检测
            double shadow_t; // 阴影检测最近交点距离
            const Sphere* shadow_sphere; // 阴影检测最近交点对应的球体
            std::tie(shadow_t, shadow_sphere) = ClosestIntersection(P, L, 0, t_max); // 阴影检测最近交点
            if (shadow_sphere != nullptr) continue; // 有阴影，跳过该光源


            // 漫反射光
            double n_dot_l = dot(N, L); // 光源方向与法向量的点积
            if (n_dot_l > 0) {
                // 光源方向与法向量夹角小于90度
                i += light.intensity * n_dot_l / (sqrt(dot(N, N))
                        * sqrt(dot(L, L))); // 漫反射光
            }


            // 镜面反射光
            if (s != -1) {
                // 高光系数不为-1，计算镜面反射光
                vector<double> R = sub(multiple(N, 2 * n_dot_l), L); // 镜面反射光方向
                double r_dot_v = dot(R, V); // 镜面反射光方向与视线方向的点积
                if (r_dot_v > 0) {
                    // 镜面反射光方向与视线方向夹角小于90度
                    i += light.intensity * pow(r_dot_v / (sqrt(dot(R, R))
                            * sqrt(dot(V, V))), s); // 镜面反射光
                }
            }
        }
    }
    return i; // 返回光照强度
}

double dot(const vector<double> &v1, const vector<double> &v2) {
    // 向量点积
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

vector<double> multiple(const vector<double> &v, const double &k) {
    // 向量数乘
    return POINT(v[0] * k, v[1] * k, v[2] * k);
}

vector<double> add(const vector<double> &v1, const vector<double> &v2) {
    // 向量加法
    return POINT(v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);
}

vector<double> sub(const vector<double> &v1, const vector<double> &v2) {
    // 向量减法
    return POINT(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);
}

vector<double> add(const vector<double> &v1, const double &v2) {
    // 向量加法
    return POINT(v1[0] + v2, v1[1] + v2, v1[2] + v2);
}

int EGERGB_Mul(int color, double m) {
    // EGE颜色数乘
    int r = (color >> 16) & 0xFF; // 取出RGB三个分量
    int g = (color >> 8) & 0xFF;
    int b = color & 0xFF;

    int new_r = static_cast<int>(r * m); // 数乘
    int new_g = static_cast<int>(g * m);
    int new_b = static_cast<int>(b * m);

    // 防止溢出
    if (new_r > 255) new_r = 255;
    if (new_g > 255) new_g = 255;
    if (new_b > 255) new_b = 255;

    // 防止负数
    if (new_r < 0) new_r = 0;
    if (new_g < 0) new_g = 0;
    if (new_b < 0) new_b = 0;

    return EGERGB(new_r, new_g, new_b); // 返回新颜色
}
