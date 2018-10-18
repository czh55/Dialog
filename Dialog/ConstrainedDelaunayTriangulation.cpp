#include "ConstrainedDelaunayTriangulation.h"


ConstrainedDelaunayTriangulation::ConstrainedDelaunayTriangulation(void)
{
}


ConstrainedDelaunayTriangulation::~ConstrainedDelaunayTriangulation(void)
{
}

// 输入多边形顶点序列
void ConstrainedDelaunayTriangulation::setInputPointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud_ply)
{
	m_ply.clear();

	// 计算协方差矩阵与重心    
	EIGEN_ALIGN16 Eigen::Matrix3f covariance_matrix;    
	Eigen::Vector4f xyz_centroid;    
	pcl::computeMeanAndCovarianceMatrix(*cloud_ply, covariance_matrix, xyz_centroid);
   
	// 协方差矩阵的特征值与特征向量    
	EIGEN_ALIGN16 Eigen::Vector3f eigen_values;    
	EIGEN_ALIGN16 Eigen::Matrix3f eigen_vectors;    
	pcl::eigen33(covariance_matrix, eigen_vectors, eigen_values);

	Eigen::Vector3f e1, e2;
	e1 = eigen_vectors.col(1);
	e2 = eigen_vectors.col(2);

	Eigen::MatrixXf A;
	A.resize(3,2);
	A.col(0) = e1;
	A.col(1) = e2;

	Eigen::Vector2f x;
	Eigen::Vector3f b;

	for(int i = 0; i < cloud_ply->size(); ++i){
		b << cloud_ply->points[i].x - xyz_centroid[0],
			 cloud_ply->points[i].y - xyz_centroid[1],
			 cloud_ply->points[i].z - xyz_centroid[2];
		x = A.colPivHouseholderQr().solve(b);
		Point point(x[0], x[1]);
		m_ply.push_back(point);
	}
}

void ConstrainedDelaunayTriangulation::mark_domains(CDT& ct, 
             CDT::Face_handle start, 
             int index, 
             std::list<CDT::Edge>& border )
{
  if(start->info().nesting_level != -1){
    return;
  }
  std::list<CDT::Face_handle> queue;
  queue.push_back(start);
  while(! queue.empty()){
    CDT::Face_handle fh = queue.front();
    queue.pop_front();
    if(fh->info().nesting_level == -1){
      fh->info().nesting_level = index;
      for(int i = 0; i < 3; i++){
        CDT::Edge e(fh,i);
        CDT::Face_handle n = fh->neighbor(i);
        if(n->info().nesting_level == -1){
          if(ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}
//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident 
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
void ConstrainedDelaunayTriangulation::mark_domains(CDT& cdt)
{
  for(CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
    it->info().nesting_level = -1;
  }
  std::list<CDT::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border);
  while(! border.empty()){
    CDT::Edge e = border.front();
    border.pop_front();
    CDT::Face_handle n = e.first->neighbor(e.second);
    if(n->info().nesting_level == -1){
      mark_domains(cdt, n, e.first->info().nesting_level+1, border);
    }
  }
}

// 执行CDT三角化(不插入新的顶点)
void ConstrainedDelaunayTriangulation::triangulate(std::vector<std::vector<int>> &triangles)
{
	//Insert the polygons into a constrained triangulation
	CDT cdt;
	cdt.insert_constraint(m_ply.vertices_begin(), m_ply.vertices_end(), true);

	//Mark facets that are inside the domain bounded by the polygon
	mark_domains(cdt);

	for(int i = 0; i < triangles.size(); ++i){
		triangles[i].clear();
	}
	triangles.clear();
	triangles.reserve(m_ply.size());

	//int i = 0;
	for (CDT::Finite_faces_iterator fit=cdt.finite_faces_begin();
		fit!=cdt.finite_faces_end();++fit){
			if ( fit->info().in_domain() ){
				std::vector<int> indices(3, -1);
				indices[0] = std::find(m_ply.vertices_begin(), m_ply.vertices_end(), fit->vertex(0)->point()) - m_ply.vertices_begin();
				indices[1] = std::find(m_ply.vertices_begin(), m_ply.vertices_end(), fit->vertex(1)->point()) - m_ply.vertices_begin();
				indices[2] = std::find(m_ply.vertices_begin(), m_ply.vertices_end(), fit->vertex(2)->point()) - m_ply.vertices_begin();
				triangles.push_back(indices);
			}
	}
}

// 执行CDT三角化(不插入新的顶点)
// triangles: 多边形集合(N-2个，N为多边形顶点个数)
// tri_index_with_min_angle: 具有最小内角的三角形索引
// min_angle: 三角形的最小内角
void ConstrainedDelaunayTriangulation::triangulate(std::vector<std::vector<int>> &triangles, std::vector<float> &min_angles)
{
	triangulate(triangles);
    min_angles.clear();
		
	int i0, i1, i2;
    Eigen::Vector2f e0, e1;

	for(int i = 0; i < triangles.size(); ++i){
		i0 = triangles[i][0];
		i1 = triangles[i][1];
		i2 = triangles[i][2];

        float max_cos_val = -2;
        float cos_val;
        float min_angle;

		// i0
		e0 << m_ply[i1].x() - m_ply[i0].x(),
			  m_ply[i1].y() - m_ply[i0].y();
		e1 << m_ply[i2].x() - m_ply[i0].x(),
			  m_ply[i2].y() - m_ply[i0].y();
		e0.normalize();
		e1.normalize();
		cos_val = e0.dot(e1);
        if(cos_val > max_cos_val){
            max_cos_val = cos_val;
		}

		// i1
		e0 << m_ply[i0].x() - m_ply[i1].x(),
			  m_ply[i0].y() - m_ply[i1].y();
		e1 << m_ply[i2].x() - m_ply[i1].x(),
			  m_ply[i2].y() - m_ply[i1].y();
		e0.normalize();
		e1.normalize();
		cos_val = e0.dot(e1);
        if(cos_val > max_cos_val){
            max_cos_val = cos_val;
        }

		// i2
		e0 << m_ply[i1].x() - m_ply[i2].x(),
			  m_ply[i1].y() - m_ply[i2].y();
		e1 << m_ply[i0].x() - m_ply[i2].x(),
			  m_ply[i0].y() - m_ply[i2].y();
		e0.normalize();
		e1.normalize();
		cos_val = e0.dot(e1);
        if(cos_val > max_cos_val){
            max_cos_val = cos_val;
        }

        min_angle = acos(max_cos_val);
        min_angle = min_angle/3.14159265358979 * 180.0;

        min_angles.push_back(min_angle);
	}
}

// 输入多边形顶点序列以及孤立顶点
// 所有孤立顶点默认均位于多边形的内部
void ConstrainedDelaunayTriangulation::setInputPointCloud_and_triangulate(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud, int ply_size,
                                         std::vector<std::vector<int>> &triangles, std::vector<float> &min_angles)
{
    setInputPointCloud(cloud);

    //Insert the polygons into a constrained triangulation
    CDT cdt;
    //cdt.insert_constraint(m_ply.vertices_begin(), m_ply.vertices_end(), true);
    int cloud_size = cloud->size();
    int points_size = cloud_size - ply_size;

    auto iter = m_ply.vertices_begin();
    for(int i = 0; i <= ply_size; ++i){
        if(i == 0)
            iter = m_ply.vertices_begin();
        else
            ++iter;
    }
    cdt.insert_constraint(m_ply.vertices_begin(), iter, true);  // 多边形
    auto iter_next = iter;
    ++iter_next;
    for(int i = 0; i < points_size; ++i){
        cdt.insert_constraint(iter, iter_next, true);           // 多边形内部点集
        ++iter;
        ++iter_next;
    }

    //Mark facets that are inside the domain bounded by the polygon
    mark_domains(cdt);

    for(int i = 0; i < triangles.size(); ++i){
        triangles[i].clear();
    }
    triangles.clear();
    triangles.reserve(m_ply.size());

    //int i = 0;
    for (CDT::Finite_faces_iterator fit=cdt.finite_faces_begin();
        fit!=cdt.finite_faces_end();++fit){
            if ( fit->info().in_domain() ){
                std::vector<int> indices(3, -1);
                indices[0] = std::find(m_ply.vertices_begin(), m_ply.vertices_end(), fit->vertex(0)->point()) - m_ply.vertices_begin();
                indices[1] = std::find(m_ply.vertices_begin(), m_ply.vertices_end(), fit->vertex(1)->point()) - m_ply.vertices_begin();
                indices[2] = std::find(m_ply.vertices_begin(), m_ply.vertices_end(), fit->vertex(2)->point()) - m_ply.vertices_begin();
                triangles.push_back(indices);
            }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    min_angles.clear();

    int i0, i1, i2;
    Eigen::Vector2f e0, e1;

    for(int i = 0; i < triangles.size(); ++i){
        i0 = triangles[i][0];
        i1 = triangles[i][1];
        i2 = triangles[i][2];

        float max_cos_val = 2;
        float cos_val;
        float min_angle;

        // i0
        e0 << m_ply[i1].x() - m_ply[i0].x(),
              m_ply[i1].y() - m_ply[i0].y();
        e1 << m_ply[i2].x() - m_ply[i0].x(),
              m_ply[i2].y() - m_ply[i0].y();
        e0.normalize();
        e1.normalize();
        cos_val = e0.dot(e1);
        if(cos_val > max_cos_val){
            max_cos_val = cos_val;
        }

        // i1
        e0 << m_ply[i0].x() - m_ply[i1].x(),
              m_ply[i0].y() - m_ply[i1].y();
        e1 << m_ply[i2].x() - m_ply[i1].x(),
              m_ply[i2].y() - m_ply[i1].y();
        e0.normalize();
        e1.normalize();
        cos_val = e0.dot(e1);
        if(cos_val > max_cos_val){
            max_cos_val = cos_val;
        }

        // i2
        e0 << m_ply[i1].x() - m_ply[i2].x(),
              m_ply[i1].y() - m_ply[i2].y();
        e1 << m_ply[i0].x() - m_ply[i2].x(),
              m_ply[i0].y() - m_ply[i2].y();
        e0.normalize();
        e1.normalize();
        cos_val = e0.dot(e1);
        if(cos_val > max_cos_val){
            max_cos_val = cos_val;
        }

        min_angle = acos(max_cos_val);
        min_angle = min_angle/3.14159265358979 * 180.0;

        min_angles.push_back(min_angle);
    }
}
