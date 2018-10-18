#include "earclip2.h"
#include "pcl/kdtree/kdtree_flann.h"

EarClip2::EarClip2()
{
}

EarClip2::~EarClip2()
{
	m_poly.clear();
}

// 设置输入的多边形
void EarClip2::setInputPoly(POLYGON2 *pPoly)
{
    if(m_poly.getSize() > 0){
        m_poly.clear();
    }

    if(pPoly->triangles.size()>0){
        pPoly->triangles.clear();
    }

    m_poly2 = pPoly;

    Vertex2 *pv2;
    int poly_index;
    for(int i = 0; i < pPoly->getSize(); ++i){
        if(i == 0){
            pv2 = pPoly->start_point;
        }else{
            poly_index = pPoly->getPolyIndex(pv2);
            pv2 = (Vertex2*)pv2->relevant_polys[poly_index].next_point;
        }
        Vertex* v = new Vertex();
        v->point.x = pv2->pos[0];
        v->point.y = pv2->pos[1];
        v->point.z = pv2->pos[2];
        if(i == 0){
            m_poly.insertFirstVertex(v);
        }else if(i == pPoly->getSize()-1){
            m_poly.insertLastVertex(v);
        }else{
            m_poly.insertMiddleVertex(v);
        }
    }
    m_poly.coeff << pPoly->coeff[0],
                    pPoly->coeff[1],
                    pPoly->coeff[2],
                    pPoly->coeff[3];

    m_ec.setInputPoly(&m_poly);
}

// 对多边形进行三角化
bool EarClip2::triangulatePoly()
{
	if(!m_ec.triangulatePoly()){
		return false;
	}

    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    vector<int> indices;
    vector<float> dists;

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
    cloud->reserve(m_poly.getSize());
    Vertex* pv;
    for(int i = 0; i < m_poly.getSize(); ++i){
        pv = i == 0 ? m_poly.start_point : pv->next_point;
        cloud->push_back(pv->point);
    }
    kdtree.setInputCloud(cloud);

    vector<Vertex2*> vertices;
    vertices.reserve(m_poly.getSize());
    Vertex2* pv2;
    int poly_index;
    for(int i = 0; i < m_poly2->getSize(); ++i){
        if(i == 0){
            pv2 = m_poly2->start_point;
        }else{
            poly_index = m_poly2->getPolyIndex(pv2);
            pv2 = (Vertex2*)pv2->relevant_polys[poly_index].next_point;
        }
        vertices.push_back(pv2);
    }

    m_poly2->triangles.reserve(m_poly.triangles.size());
    for(int i = 0; i < m_poly.triangles.size(); ++i){
		Triangle2 tri;

		for(int j = 0; j < 3; ++j){
			kdtree.nearestKSearch(m_poly.triangles[i].v[j]->point, 1, indices, dists);
			tri.v[j] = vertices[indices[0]];
		}

		m_poly2->triangles.push_back(tri);
    }

	m_poly.clear();

	return true;
}
