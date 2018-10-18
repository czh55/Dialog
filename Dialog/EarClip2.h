#ifndef EARCLIP2_H
#define EARCLIP2_H

#include <pcl/point_cloud.h>
#include <vector>
#include "holefilling.h"
#include "earclip.h"

class EarClip2
{
public:
    EarClip2();
    ~EarClip2();
    // 设置输入的多边形
    void setInputPoly(POLYGON2 *pPoly);
    // 对多边形进行三角化
    bool triangulatePoly();

private:

    EarClip m_ec;
    POLYGON m_poly;
    POLYGON2* m_poly2;
};

#endif // EARCLIP2_H
