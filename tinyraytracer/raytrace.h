//
// Created by Vice on 2023/6/14.
//

#ifndef EXERCISE_RAYTRACE_H
#define EXERCISE_RAYTRACE_H
#include <iostream>
#include <graphics.h>
#include <utility>
#include <vector>
#include <tuple>
#include <cmath>
#include "scene.h"
using std::vector;
using std::cout;
using std::endl;

void PutPixel(int x, int y, int color); // 画点
vector<double> CanvasToViewport(int x, int y); // 画布坐标转视口坐标
std::tuple<double, const Sphere*> ClosestIntersection(const vector<double> &O, const vector<double> &D, const double &t_min, const double &t_max); // 计算最近交点
int TraceRay(const vector<double> &O, const vector<double> &D, const double &t_min, const double &t_max); // 光线追踪
int TraceRay(const vector<double> &O, const vector<double> &D, const double &t_min, const double &t_max, int recursion_depth); // 光线追踪带反射
vector<double> ReflectRay(const vector<double> &L,  const vector<double> &N); // 计算反射光线
double ComputeLighting(const vector<double> &P, const vector<double> &N, const vector<double> &V, const int &s); // 计算光照
vector<double> multiple(const vector<double> &v, const double &k); // 向量数乘
vector<double> add(const vector<double> &v1, const vector<double> &v2); // 向量加法
vector<double> sub(const vector<double> &v1, const vector<double> &v2); // 向量减法
vector<double> add(const vector<double> &v1, const double &v2); // 向量数加
double dot(const vector<double> &v1, const vector<double> &v2); // 向量点乘
int EGERGB_Mul(int color, double m); // 颜色数乘
int EGERGB_Add(int color1, int color2);


#endif //EXERCISE_RAYTRACE_H
