//
// Created by Vice on 2023/6/28.
//
#include "raytrace.h"

int main() {
    int color = EGERGB(255, 255, 255);
    vector<double> D = POINT(0, 0, d);
    initgraph(WIDTH, HEIGHT); // 初始化画布

    for (int x = -WIDTH/2; x <= WIDTH/2; ++x) {
        for (int y = -HEIGHT/2; y <= HEIGHT/2; ++y) {
            D = CanvasToViewport(x, y);
            color = TraceRay(O, D, d, INFINITY, 1);
            PutPixel(x, y, color);
        }
    }

    getch();
    return 0;
}