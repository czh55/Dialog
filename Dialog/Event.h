#pragma once
#ifndef EVENT_H
#define EVENT_H

#include "PCLViewer.h"
#include "PlaneDetect.h"

// 键盘事件处理函数
int iter_count = 0;
void keyboardEventOccurred(const pcl::visualization::KeyboardEvent &event, void* viewer_void)
{

	//前半段
	if (event.getKeySym() == "Return" && event.keyDown())
	{
		if (strcmp(state, "segment ps") == 0)
		{
			// 对参数空间
			cout << "id_sink = " << id_sink << endl;
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
			int sink_index = sinkPointsSet[id_sink].index;
			for (int i = 0; i < ps_cloud->size(); ++i)
			{
				if (param_space[i].sink == sink_index)
				{
					pcl::PointXYZRGB p = ps_cloud->points[i];
					p.r = p.g = p.b = 255;
					cloud->push_back(p);
				}
			}
			//viewer_ps.removePointCloud("sink field");
			//viewer_ps.addPointCloud(cloud, "sink field");
			cout << "ps points size = " << cloud->size() << endl;

			// 对测量空间
			PointCloudT::Ptr plane(new PointCloudT);
			for (int i = 0; i < ps_cloud->size(); ++i)
			{
				if (param_space[i].sink == sink_index)
				{
					for (int j = 0; j < param_space[i].point_indices.size(); ++j)
					{
						plane->push_back(source_cloud->points[param_space[i].point_indices[j]]);
					}
				}
			}
			((PCLViewer*)m_pclviewer)->viewer->removePointCloud("plane");
			((PCLViewer*)m_pclviewer)->viewer->addPointCloud(plane, "plane");
			((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 255, 0, "plane");
			cout << "plane has " << plane->size() << " point.s" << endl;

			id_sink = ++id_sink % sinkPointsSet.size();
		}
	}

	if (event.getKeySym() == "F1" && event.keyDown())
	{
		//viewer_ps.removePointCloud("sink field");
		((PCLViewer*)m_pclviewer)->viewer->removePointCloud("plane");
		id_sink = id_plane = 0;
	}

	//后半段
	PCLViewer *pclviewer = (PCLViewer *)viewer_void;
	if (event.getKeySym() == "Return" && event.keyDown()) {
		if (strcmp(CMD, "see poly") == 0) {
			int poly_id = iter_count++ % pclviewer->m_sv.polygons.size();
			cout << "poly id = " << poly_id << endl;
			cout << "vertices size = " << pclviewer->m_sv.polygons[poly_id].getSize() << endl;
			POLYGON *pPoly = &pclviewer->m_sv.polygons[poly_id];

			// 首先删除上一步生成的多边形
			string line_name;
			char buf[64];
			int size = 1024;
			for (int i = 0; i < size; ++i) {
				line_name = "_";
				itoa(i, buf, 16);
				line_name += buf;
				pclviewer->viewer->removeShape(line_name);
			}

			// 画出当前多边形
			Vertex *pv;
			for (int i = 0; i < pPoly->getSize(); ++i) {
				pv = i == 0 ? pPoly->start_point : pv->next_point;
				pcl::PointXYZ p1, p2;
				p1 = pv->point;
				p2 = pv->next_point->point;

				line_name = "_";
				itoa(i, buf, 16);
				line_name += buf;

				pclviewer->viewer->addLine(p1, p2, 1.0, 0.5, 0.5, line_name);
				pclviewer->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);
			}
		}
	}
}



// 点选取事件处理函数
void pointPickingEventOccurred(const pcl::visualization::PointPickingEvent &event, void* viewer_void)
{
	//前半段
	if (strcmp(state, "estimate normal") == 0 || strcmp(state, "regulate normal") == 0)
	{
		selected_point_index = event.getPointIndex();
		cout << "selected point index = " << selected_point_index << endl;
		cout << "open config.txt, and set the value of is_norm_direction_valid" << endl;

		PointCloudT::Ptr selected_point(new PointCloudT);
		selected_point->push_back(source_cloud->points[selected_point_index]);
		pcl::PointCloud<pcl::Normal>::Ptr selected_normal(new pcl::PointCloud<pcl::Normal>);
		selected_normal->push_back(source_normal->points[selected_point_index]);

		((PCLViewer*)m_pclviewer)->viewer->removePointCloud("selected normal");
		((PCLViewer*)m_pclviewer)->viewer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>
			(selected_point, selected_normal, 1, normal_length_for_selected, "selected normal");

		point_index = event.getPointIndex();
		PointCloudT::Ptr point(new PointCloudT);
		point->push_back(source_cloud->points[point_index]);
		((PCLViewer*)m_pclviewer)->viewer->removePointCloud("point");
		((PCLViewer*)m_pclviewer)->viewer->addPointCloud(point, "point");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 255, 0, 0, "point");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "point");
	}

	if (strcmp(state, "segment planes") == 0)
	{
		cout << "point_index = " << point_index << endl;
		cout << "growth_unit_id = " << growth_unit_id << endl;

		PointCloudT::Ptr plane(new PointCloudT);
		plane = plane_clouds[growth_unit_id].points_set;

		std::vector<int> indices;
		std::vector<int> knn_indices;
		std::vector<float> knn_dist;

		kdtree_source.nearestKSearch(plane->points[0], 1, indices, knn_dist);

		Eigen::Vector3f point_normal;
		point_normal << source_normal->points[indices[0]].normal_x,
			source_normal->points[indices[0]].normal_y,
			source_normal->points[indices[0]].normal_z;
		Eigen::Vector3f growth_unit_normal, local_point_normal;

		// 首先计算种子平面的法向量
		Eigen::Vector4f plane_param;
		float curvature;
		for (int i = 0; i < plane->size(); ++i)
		{
			kdtree_source.nearestKSearch(plane->points[i], 1, knn_indices, knn_dist);
			indices.push_back(knn_indices[0]);
		}
		pcl::computePointNormal(*source_cloud, indices, plane_param, curvature);
		growth_unit_normal << plane_param[0], plane_param[1], plane_param[2];
		if (growth_unit_normal.dot(point_normal) < 0)
			growth_unit_normal *= -1.0f;
		indices.clear();
		// r_local然后计算局部法向量
		kdtree_source.radiusSearch(source_cloud->points[point_index], r_local, indices, knn_dist);
		// cout << "indices.size() = " << indices.size() << endl;
		pcl::computePointNormal(*source_cloud, indices, plane_param, curvature);
		local_point_normal << plane_param[0], plane_param[1], plane_param[2];
		Eigen::Vector3f source_point_normal;
		source_point_normal << source_normal->points[point_index].normal_x,
			source_normal->points[point_index].normal_y,
			source_normal->points[point_index].normal_z;
		if (local_point_normal.dot(source_point_normal) < 0)
			local_point_normal *= -1.0f;

		float angle = arccos(growth_unit_normal.dot(local_point_normal)) * 180.0f / PI;
		cout << "angle = " << angle << endl;

		((PCLViewer*)m_pclviewer)->viewer->removeShape("sphere");
		((PCLViewer*)m_pclviewer)->viewer->addSphere(source_cloud->points[point_index], r_local, 255, 0, 0, "sphere");

		// 显示增长单元的法向量
		pcl::PointXYZ p;
		p.x = p.y = p.z = 0;
		for (int i = 0; i < plane->size(); ++i)
		{
			p.x += plane->points[i].x;
			p.y += plane->points[i].y;
			p.z += plane->points[i].z;
		}
		p.x /= float(plane->size());
		p.y /= float(plane->size());
		p.z /= float(plane->size());
		kdtree_source.nearestKSearch(p, 1, indices, knn_dist);
		p = source_cloud->points[indices[0]];

		pcl::Normal pn;
		pn.normal_x = growth_unit_normal[0];
		pn.normal_y = growth_unit_normal[1];
		pn.normal_z = growth_unit_normal[2];

		PointCloudT::Ptr cloud1(new PointCloudT);
		cloud1->push_back(p);
		pcl::PointCloud<pcl::Normal>::Ptr normal1(new pcl::PointCloud<pcl::Normal>);
		normal1->push_back(pn);
		((PCLViewer*)m_pclviewer)->viewer->removePointCloud("growth_unit_normal");
		((PCLViewer*)m_pclviewer)->viewer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>
			(cloud1, normal1, 1, normal_length_for_selected, "growth_unit_normal");
		// 显示选中点的局部法向量
		p = source_cloud->points[point_index];
		pn.normal_x = local_point_normal[0];
		pn.normal_y = local_point_normal[1];
		pn.normal_z = local_point_normal[2];
		PointCloudT::Ptr cloud2(new PointCloudT);
		cloud2->push_back(p);
		pcl::PointCloud<pcl::Normal>::Ptr normal2(new pcl::PointCloud<pcl::Normal>);
		normal2->push_back(pn);
		((PCLViewer*)m_pclviewer)->viewer->removePointCloud("local_normal");
		((PCLViewer*)m_pclviewer)->viewer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>
			(cloud2, normal2, 1, normal_length_for_selected / 2.0, "local_normal");
	}

	if (strcmp(state, "del poly") == 0 || strcmp(state, "prune poly") == 0)
	{
		if (plane_clouds_final.size() == 0) return;

		QString cloud_id = QString::number(g_selected_poly_id, 10);
		if (!g_is_poly_del[g_selected_poly_id])
		{
			// 将之前选取的多边形颜色回复
			((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 255, 255, 255, cloud_id.toStdString());
			((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, cloud_id.toStdString());
		}

		event.getPoint(g_selected_point.x, g_selected_point.y, g_selected_point.z);
		// cout << "g_selected_point = " << g_selected_point << endl;

		// 获取被选择的多边形的id
		std::vector<int> indices;
		std::vector<float> dists;
		kdtree_poly_centroids_cloud.nearestKSearch(g_selected_point, 1, indices, dists);
		g_selected_poly_id = indices[0];

		// 显示已经选取的点
		PointCloudT::Ptr point(new PointCloudT);
		point->push_back(poly_centroids_cloud->points[g_selected_poly_id]);
		((PCLViewer*)m_pclviewer)->viewer->removePointCloud("point");
		((PCLViewer*)m_pclviewer)->viewer->addPointCloud(point, "point");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 255, 0, 0, "point");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "point");

		// 用绿色显示被选中的多边形
		cloud_id = QString::number(g_selected_poly_id, 10);
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 255, 0, cloud_id.toStdString());
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, cloud_id.toStdString());
	}

	if (strcmp(state, "prune poly select points") == 0)
	{
		event.getPoint(g_selected_point.x, g_selected_point.y, g_selected_point.z);
		// 显示已经选取的点
		PointCloudT::Ptr point(new PointCloudT);
		point->push_back(g_selected_point);
		((PCLViewer*)m_pclviewer)->viewer->removePointCloud("point");
		((PCLViewer*)m_pclviewer)->viewer->addPointCloud(point, "point");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 255, 0, 0, "point");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "point");
	}

	//后半段
	PCLViewer *pclviewer = (PCLViewer *)viewer_void;
	if (strcmp(CMD, "observe poly id") == 0) {
		pcl::PointXYZ point;
		event.getPoint(point.x, point.y, point.z);
		vector<int> indices;
		vector<float> dists;
		pclviewer->m_kdtree.nearestKSearch(point, 1, indices, dists);
		int poly_id = floor(pclviewer->m_vertices->points[indices[0]].data[3] + 0.5f);
		cout << "poly id = " << poly_id << endl;

		// 高亮显示选中的顶点
		pclviewer->viewer->removePointCloud("selected point");
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->push_back(point);
		pclviewer->viewer->addPointCloud(cloud, "selected point");
		pclviewer->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "selected point");
		pclviewer->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "selected point");
	}

	if (strcmp(CMD, "select poly") == 0) {
		pcl::PointXYZ point;
		event.getPoint(point.x, point.y, point.z);
		vector<int> indices;
		vector<float> dists;
		pclviewer->m_hf.m_cloud_for_indentify_poly_id_kdtree.nearestKSearch(point, 1, indices, dists);
		pclviewer->m_hfio_dlg.selected_poly_id = floor(pclviewer->m_hf.m_cloud_for_indentify_poly_id->points[indices[0]].data[3] + 0.5);
		pclviewer->m_hfio_dlg.updateInfo();

		// 显示选取的多边形
		POLYGON2* pPoly = pclviewer->m_hf.m_polygons[pclviewer->m_hfio_dlg.selected_poly_id];
		char buf[64];
		string name;
		Vertex2* pv2, *pv2_next;
		int poly_index;
		pcl::PointXYZ p1, p2;

		// 清空之前的多边形
		int max_poly_size = 128;
		for (int i = 0; i < max_poly_size; ++i) {
			itoa(i, buf, 16);
			name = "display_";
			name.append(buf);
			name.append("_l");
			//cout << name << endl;
			pclviewer->viewer->removeShape(name);
		}
		// cout << "//////////////////////////////////////////////////////////////////////////////////" << endl;
		for (int j = 0; j < pPoly->getSize(); ++j) {
			if (j == 0) {
				pv2 = pPoly->start_point;
			}
			else {
				poly_index = pPoly->getPolyIndex(pv2);
				pv2 = (Vertex2 *)pv2->relevant_polys[poly_index].next_point;
			}

			itoa(j, buf, 16);
			name = "display_";
			name.append(buf);
			name.append("_l");
			//cout << name << endl;
			p1.x = pv2->pos[0];
			p1.y = pv2->pos[1];
			p1.z = pv2->pos[2];

			poly_index = pPoly->getPolyIndex(pv2);
			pv2_next = (Vertex2*)pv2->relevant_polys[poly_index].next_point;
			p2.x = pv2_next->pos[0];
			p2.y = pv2_next->pos[1];
			p2.z = pv2_next->pos[2];

			pclviewer->viewer->addLine(p1, p2, 0.3, 0.8, 0.7, name);
			pclviewer->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, name);
		}
	}

	if (strcmp(CMD, "select vertex1") == 0) {
		pcl::PointXYZ point;
		event.getPoint(point.x, point.y, point.z);
		vector<int> indices;
		vector<float> dists;

		pclviewer->m_hf.m_vertices_kdtree.nearestKSearch(point, 1, indices, dists);
		pclviewer->m_hfio_dlg.m_vertex1 = pclviewer->m_hf.m_vertices[indices[0]];
		pclviewer->m_hfio_dlg.updateInfo();

		point.x = pclviewer->m_hfio_dlg.m_vertex1->pos[0];
		point.y = pclviewer->m_hfio_dlg.m_vertex1->pos[1];
		point.z = pclviewer->m_hfio_dlg.m_vertex1->pos[2];

		// 显示顶点1
		pclviewer->viewer->removePointCloud("vertex1");
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->push_back(point);
		pclviewer->viewer->addPointCloud(cloud, "vertex1");
		pclviewer->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertex1");
		pclviewer->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "vertex1");
	}

	if (strcmp(CMD, "select vertex2") == 0) {
		pcl::PointXYZ point;
		event.getPoint(point.x, point.y, point.z);
		vector<int> indices;
		vector<float> dists;

		pclviewer->m_hf.m_vertices_kdtree.nearestKSearch(point, 1, indices, dists);
		pclviewer->m_hfio_dlg.m_vertex2 = pclviewer->m_hf.m_vertices[indices[0]];
		pclviewer->m_hfio_dlg.updateInfo();

		point.x = pclviewer->m_hfio_dlg.m_vertex2->pos[0];
		point.y = pclviewer->m_hfio_dlg.m_vertex2->pos[1];
		point.z = pclviewer->m_hfio_dlg.m_vertex2->pos[2];

		// 显示顶点1
		pclviewer->viewer->removePointCloud("vertex2");
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->push_back(point);
		pclviewer->viewer->addPointCloud(cloud, "vertex2");
		pclviewer->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertex2");
		pclviewer->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "vertex2");
	}

	if (strcmp(CMD, "insert dis") == 0) {
		pcl::PointXYZ point;
		event.getPoint(point.x, point.y, point.z);
		pclviewer->m_selected_p = point;

		// 高亮显示选中的顶点
		pclviewer->viewer->removePointCloud("selected point");
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->push_back(point);
		pclviewer->viewer->addPointCloud(cloud, "selected point");
		pclviewer->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "selected point");
		pclviewer->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "selected point");
	}
}


// 区域选取事件函数
void areaPickingEventOccurred(const pcl::visualization::AreaPickingEvent &event,void* viewer_void)
{
	/*
	PCLViewer *pclviewer = (PCLViewer *)viewer_void;
	if (pclviewer->m_hf.m_vertices->size() == 0) return;

	std::vector<int> indices;
	event.getPointsIndices(indices);
	pcl::PointCloud<pcl::PointXYZ>::Ptr selected_points(new pcl::PointCloud<pcl::PointXYZ>);
	for (int i = 0; i < indices.size(); ++i) selected_points->push_back(pclviewer->m_hf.m_vertices->points[indices[i]]);

	pclviewer->viewer->removePointCloud("selected points");
	pclviewer->viewer->addPointCloud(selected_points, "selected points");
	pclviewer->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1.0, 0, "selected points");
	pclviewer->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 8, "selected points");

	cout << "selected " << indices.size() << " points." << endl;

	pcl::visualization::Camera camera;
	pclviewer->viewer->getCameraParameters(camera);
	pcl::PointXYZ camera_pos;
	camera_pos.x = camera.pos[0];
	camera_pos.y = camera.pos[1];
	camera_pos.z = camera.pos[2];
	cout << "camera_pos = " << camera_pos << endl; 
	*/

	/*
	pcl::PointXYZ origin_pos;
	origin_pos.x = origin_pos.y = origin_pos.z = 0;
	pclviewer->viewer->removeShape("line");
	pclviewer->viewer->addLine(origin_pos, camera_pos, 1, 1, 1, "line");
	*/

	/*
	if(selected_points->size() < 1) return;
	pcl::Normal pn;
	Eigen::Vector3d n;
	pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);
	Eigen::Vector3d n_vec;
	n_vec << 0, 0, 0;
	for(int i = 0; i < selected_points->size(); ++i){
	n << camera_pos.x - selected_points->points[i].x,
	camera_pos.y - selected_points->points[i].y,
	camera_pos.z - selected_points->points[i].z;
	n.normalize();
	pn.normal_x = n[0];
	pn.normal_y = n[1];
	pn.normal_z = n[2];
	n_vec += n;
	normals->push_back(pn);
	}

	n_vec = 1.0/selected_points->size() * n_vec;
	for(int i = 0; i < selected_points->size(); ++i){
	normals->points[i].normal_x = n_vec[0];
	normals->points[i].normal_y = n_vec[1];
	normals->points[i].normal_z = n_vec[2];
	}

	pclviewer->viewer->removePointCloud("normals");
	pclviewer->viewer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(selected_points, normals, 1, 1, "normals");
	*/
}

#endif // !EVENT_H

