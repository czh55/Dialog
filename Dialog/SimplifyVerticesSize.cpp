#include "SimplifyVerticesSize.h"


SimplifyVerticesSize::SimplifyVerticesSize(void)
{
}


SimplifyVerticesSize::~SimplifyVerticesSize(void)
{
}

// 构造多边形顶点法向量(默认顶点排序已经满足右手法则)
// 仅利用前后顶点以及当前顶点确定当前顶点的法向量
void SimplifyVerticesSize::computeVerticesNormals(pcl::PointCloud<pcl::PointXYZ>::Ptr &contour, 
	                        std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> &vertices_normals,
							Eigen::Vector3f &plane_normal)
{
	if(vertices_normals.size() > 0) vertices_normals.clear();
    vertices_normals.reserve(contour->size());
	int prev, next;
	Eigen::Vector3f v_prev_p, v_p_next, v1, v2, v;
	float ratio1, ratio2, sum_dist, dist1, dist2;
	for(int i = 0; i < contour->size(); ++i)
	{
		prev = i == 0 ? contour->size()-1 : i-1;
		next = i == contour->size()-1 ? 0 : i+1;
		v_prev_p << contour->points[i].x - contour->points[prev].x,
					contour->points[i].y - contour->points[prev].y,
					contour->points[i].z - contour->points[prev].z;
		v_p_next << contour->points[next].x - contour->points[i].x,
					contour->points[next].y - contour->points[i].y,
					contour->points[next].z - contour->points[i].z;
		v_prev_p.normalize(); v_p_next.normalize();
		v1 = v_prev_p.cross(plane_normal);
		v2 = v_p_next.cross(plane_normal);
		v1.normalize(); v2.normalize();
		dist1 = distP2P(contour->points[i], contour->points[prev]);
		dist2 = distP2P(contour->points[i], contour->points[next]);
		sum_dist = dist1 + dist2;
		ratio1 = 1-dist1/sum_dist;
		ratio2 = 1-dist2/sum_dist;
		v = ratio1*v1+ratio2*v2;
        v.normalize();
		vertices_normals.push_back(v);
	}
}
// 区域生长构造初始的线段集合
// 根据顶点的法向量，使用区域生长法提取直线
// 对种子点的选取原则：种子点与相邻点的法向量必须具有较好地一致性(方向上与距离上)
// 随时更新的变量：拟合直线的参数，拟合直线的法向量
// 检查变量：点到拟合直线的距离，点的法向量与拟合直线的法向量的夹角
void SimplifyVerticesSize::constructInitLineSegs( pcl::PointCloud<pcl::PointXYZ>::Ptr &contour,// 多边形的轮廓线
	                        std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> &vertices_normals,// 顶点法向量
							Eigen::Vector3f &poly_normal,	// 多边形法向量
							std::vector<std::vector<int>> &line_segs )	// 作为输出，有：（1）线段的起点与终点；（2）孤立点
{
	if(line_segs.size() > 0) line_segs.clear();

	std::vector<bool> is_processed(contour->size(), false);

	pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);
	pcl::PointIndices::Ptr seg_inliers (new pcl::PointIndices);
	pcl::SACSegmentation<pcl::PointXYZ> seg;
	seg.setModelType (pcl::SACMODEL_LINE);
	seg.setMethodType (pcl::SAC_RANSAC);
	seg.setDistanceThreshold (FLT_MAX);

	pcl::PointCloud<pcl::PointXYZ>::Ptr line_cloud (new pcl::PointCloud<pcl::PointXYZ>);

	int seed_index = -1;
	while(true)
	{
		line_cloud->clear();
		getSeedIndex(contour, seed_index, is_processed, vertices_normals);
		if(seed_index == -1) break;

		int prev = seed_index == 0 ? contour->size()-1 : seed_index-1;
		int next = seed_index == contour->size()-1 ? 0 : seed_index+1;
		is_processed[prev] = is_processed[next] = is_processed[seed_index] = true;
		std::vector<int> indices;
		indices.push_back(prev); indices.push_back(seed_index); indices.push_back(next);

		// 计算拟合直线的参数
		for(int i = 0; i < indices.size(); ++i) line_cloud->push_back(contour->points[indices[i]]);
		seg.setInputCloud(line_cloud);
		seg.segment(*seg_inliers, *coefficients);
		Eigen::Vector3f p_base, line_dir, line_normal, /*poly_normal,*/ p_normal;
		p_base << coefficients->values[0], coefficients->values[1], coefficients->values[2];
		line_dir << coefficients->values[3], coefficients->values[4], coefficients->values[5];
		/*poly_normal << poly_set[poly_index].coeff[0], poly_set[poly_index].coeff[1], poly_set[poly_index].coeff[2];*/
		line_normal/*拟合直线的法向量，需指向多边形外部*/ = line_dir.cross(poly_normal); line_normal.normalize();
		p_normal = vertices_normals[seed_index];
		if(p_normal.dot(line_normal) < 0) line_normal *= -1.0f;

		// 从种子点的前后分别开始生长
		bool is_prev_stopped = false;
		bool is_next_stopped = false;

		while(true)
		{
			// update prev & next
			if(!is_prev_stopped) prev = prev == 0 ? contour->size()-1 : prev-1;
			if(!is_next_stopped) next = next == contour->size()-1 ? 0 : next+1;
			// 遇到已经被处理过的点则停止生长
			if(is_processed[prev]) is_prev_stopped = true;
			if(is_processed[next]) is_next_stopped = true;

			if(!is_prev_stopped)
			{
				float dist = distP2Line(contour->points[prev], p_base, line_dir);
				p_normal = vertices_normals[prev];/*<< normal_cloud->points[prev].normal_x, normal_cloud->points[prev].normal_y, normal_cloud->points[prev].normal_z;*/
                float angle = acos(line_normal.dot(p_normal));
				if(angle <= T_angle_betPointAndLine && dist <= T_dist_betPointAndLine)
				{
					is_processed[prev] = true;
					indices.push_back(prev);
					line_cloud->push_back(contour->points[prev]);

					// update line params
					seg.setInputCloud(line_cloud);
					seg.segment(*seg_inliers, *coefficients);
					p_base << coefficients->values[0], coefficients->values[1], coefficients->values[2];
					line_dir << coefficients->values[3], coefficients->values[4], coefficients->values[5];
					line_normal = line_dir.cross(poly_normal); line_normal.normalize();
					p_normal = vertices_normals[prev];
					if(p_normal.dot(line_normal) < 0) line_normal *= -1.0f;
				}else{ // 遇到不满足条件的点，则停止生长
					is_prev_stopped = true;
				}
			}

			if(!is_next_stopped)
			{
				float dist = distP2Line(contour->points[next], p_base, line_dir);
				p_normal = vertices_normals[next]; /*<< normal_cloud->points[next].normal_x, normal_cloud->points[next].normal_y, normal_cloud->points[next].normal_z;*/
				float angle = acos(line_normal.dot(p_normal));
				if(angle <= T_angle_betPointAndLine && dist <= T_dist_betPointAndLine)
				{
					is_processed[next] = true;
					indices.push_back(next);
					line_cloud->push_back(contour->points[next]);

					// update line params
					seg.setInputCloud(line_cloud);
					seg.segment(*seg_inliers, *coefficients);
					p_base << coefficients->values[0], coefficients->values[1], coefficients->values[2];
					line_dir << coefficients->values[3], coefficients->values[4], coefficients->values[5];
					line_normal = line_dir.cross(poly_normal); line_normal.normalize();
					p_normal = vertices_normals[next];
					if(p_normal.dot(line_normal) < 0) line_normal *= -1.0f;
				}else{ // 遇到不满足条件的点，则停止生长
					is_next_stopped = true;
				}
			}

			if(is_prev_stopped && is_next_stopped)
			{	
				line_segs.push_back(indices);
				break;
			}
		}
	}
}

void SimplifyVerticesSize::getSeedIndex( pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud,
				   int &seed_index,
				   std::vector<bool> &is_processed,
				   std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> &vertices_normals )
{
	seed_index = -1;
	Eigen::Vector3f v_prev, v_i, v_next;
	float angle1, angle2;
	for(int i = 0; i < cloud->size(); ++i)
	{
		if(is_processed[i]) continue;	// 如果该点已被处理过，则不能作为种子点

		int prev = i == 0 ? cloud->size()-1 : i-1;
		int next = i == cloud->size()-1 ? 0 : i+1;

		if(is_processed[prev] || is_processed[next]) continue;	// 如果该点之前或之后的点被处理过，则该点不能作为种子点
		
		v_prev = vertices_normals[prev];
		v_i = vertices_normals[i];
		v_next = vertices_normals[next];

		angle1 = acos(v_i.dot(v_prev));
		angle2 = acos(v_i.dot(v_next));

		if(angle1 <= T_angle_betPointAndPoint && angle2 <= T_angle_betPointAndPoint)	// 仅从法向量的方向上进行判断
		{
			seed_index = i;
			break;
		}
	}
}

float SimplifyVerticesSize::distP2Line( pcl::PointXYZ &point,
				  Eigen::Vector3f &p_base,		// 直线模型的基点
				  Eigen::Vector3f &line_dir )	// 直线模型的方向向量
{
	Eigen::Vector3f p, delta;
	p << point.x, point.y, point.z;
	delta = p_base - p;

	float lambda = -1.0 * line_dir.dot(delta);

	Eigen::Vector3f p_proj;
	p_proj = p_base + lambda * line_dir;

	return (p_proj-p).norm();
}

float SimplifyVerticesSize::distP2P(pcl::PointXYZ &p1, pcl::PointXYZ &p2)
{
	return pow((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z), 0.5f);
}

// 继续优化多边形的顶点序列
void SimplifyVerticesSize::refineLineSegs( pcl::PointCloud<pcl::PointXYZ>::Ptr &contour,			// 原始顶点
					 pcl::PointCloud<pcl::PointXYZ>::Ptr &final_vertices,	// 最终的顶点集合（x,y,z)
					 std::vector<std::vector<int>> &retrieved_line_segs, // 区域生长法提取的线段序列
					 std::vector<int> &final_line_segs )	 // 最终的线段序列(index)
{
	// final_line_segs保存最终的线段序列（仅保存线段的起点在contour中的索引)
	// 线段L=(pi,pj)的表示：只保留线段的起点pi(在final_line_segs中)
	if(final_line_segs.size()>0) final_line_segs.clear();
	if(final_vertices==NULL){
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
		final_vertices = cloud;
	} else { final_vertices->clear(); }

	// 特殊情况：retrieved_line_segs.size()==0
	if(retrieved_line_segs.size()==0)
	{
		*final_vertices = *contour;
		return;
	}

	// step1:
	// 对线段顶点进行排序(区域生长法是在种子点两旁同时生长的，因而retrieved_line_segs中某条线段内顶点的序列并不一定是连续的）
	// 一般情况仅需做一个简单的排序即可（里面不能包含多边形序列的起点和终点）
	// 特殊情况是多边形的起点和终点均被纳入同一条线段中
	for(int i = 0; i < retrieved_line_segs.size(); ++i)
	{
		std::sort(retrieved_line_segs[i].begin(), retrieved_line_segs[i].end());
		/* 索引0就一定会被包括进去吗？
		 * 在区域生长时，每次仅检查生长单元相邻的顶点，因而生成的线段序列内部一定是连续的。
		 * 如果索引0没有被纳入生长出来的线段内，则每条线段内部的点的索引必然是有序的
		 * 如果索引0被纳入生长出来的线段内，且0之前的点也有一些被纳入进来，则升序排列后，0必然出现在第一个索引上。
		 *	   且排序后的顶点序列之间有一个断层，不是加一，而是加一个比1大很多的数，后面的这一子线段应该移动到前面来。
		*/
        int i_size = retrieved_line_segs[i].size();
        if(retrieved_line_segs[i][0] == 0 && retrieved_line_segs[i][i_size-1] != i_size-1)
		{
			for(int j = 0; j < retrieved_line_segs[i].size()-1; ++j)
			{
				int next = j+1;
				if(next != retrieved_line_segs[i][j+1])
				{
					std::vector<int> tmp;
					for(int k = next; k < retrieved_line_segs[i].size(); ++k)
						tmp.push_back(retrieved_line_segs[i][k]);
					for(int k = 0; k < next; ++k)
						tmp.push_back(retrieved_line_segs[i][k]);
					for(int k = 0; k < tmp.size(); ++k)
						retrieved_line_segs[i][k] = tmp[k];
					break;
				}
			}
		}
	}

	// step2:
	// 抽取线段顶点以及两个线段之间的孤立点,最终构成初始的线段序列（保存在final_line_segs）中
	for(int i = 0; i < retrieved_line_segs.size(); ++i)
	{	// i: 当前线段在retrieved_line_segs中的索引
		int begin = retrieved_line_segs[i][0];	// 第一个线段的起点
		int end = retrieved_line_segs[i][retrieved_line_segs[i].size()-1];	// 第一个线段的终点
		// 下一条线段的起点
		int begin_next_line = i < retrieved_line_segs.size()-1 ? retrieved_line_segs[i+1][0] : retrieved_line_segs[0][0];
		final_line_segs.push_back(begin);
		final_line_segs.push_back(end);
		
		// end->end's next(end表示第i条线段之后孤立点的索引)
		end = end == contour->size()-1 ? 0 : end+1;
		while(end != begin_next_line)
		{
			final_line_segs.push_back(end);
			end = end == contour->size()-1 ? 0 : end+1;
		}
	}

	assert(final_line_segs.size() >= 3);
	if(final_line_segs.size() < 3){
		*final_vertices = *contour;
		return;
	}

	// step3:
	// 如果某条线段的长度<<线段的平均长度，则可将这条线段删除
	// 每次仅去掉一个极小线段，当不满足继续去除的条件时便停止
	std::vector<float> line_lens;
	line_lens.reserve(final_line_segs.size());
	float total_len = 0;
	for(int begin = 0; begin < final_line_segs.size(); ++begin)
	{
		int end = begin == final_line_segs.size()-1 ? 0 : begin+1;

		int begin_index/*w.r.t. contour*/ = final_line_segs[begin];
		int end_index/*w.r.t. contour*/ = final_line_segs[end];
		float dist = distP2P(contour->points[begin_index], contour->points[end_index]);
		line_lens.push_back(dist);
		total_len += dist;
	}
	float mean_len = total_len / float(final_line_segs.size());

	std::vector<bool> is_del(final_line_segs.size(), false);

	std::vector<int> indices;
	vectorSort(line_lens, indices);

	int residual_points = final_line_segs.size();
	int cur_index = indices[0];

    // float T_lenRatio = 0.25;
	int count_del = 0;
	std::vector<int> offset(indices.size(), 0);
	for(int i = 0; i < indices.size(); ++i)
	{
		cur_index = indices[i];
		if(residual_points <= 3) break;
		if(line_lens[cur_index]/mean_len < T_lenRatio){
			is_del[cur_index] = true;
			--residual_points;
			++count_del;
			offset[i]=-1;	// -1表明该索引对应的线段会被删除
		}else{
			break;
		}
	}

	// 设置offset
	int off_set = 0;
	for(int i = 0; i < offset.size(); ++i)
	{
		if(!is_del[i]){
			offset[i] = off_set;
		}else{
			++off_set;
		}
	}

	for(int i = 0; i < indices.size(); ++i)
	{
		if(offset[i] > 0 && !is_del[i]){
			final_line_segs[i - offset[i]] = final_line_segs[i];
		}
	}
	for(int i = 0; i < count_del; ++i)
		final_line_segs.pop_back();
	

	// method 2
	//// 计算线段的长度
	//std::vector<float> line_lens, line_lens_for_sort;
	//float total_len = 0;
	//for(int begin = 0; begin < final_line_segs.size(); ++begin)
	//{
	//	int end = begin == final_line_segs.size()-1 ? 0 : begin+1;

	//	int begin_index/*w.r.t. contour*/ = final_line_segs[begin];
	//	int end_index/*w.r.t. contour*/ = final_line_segs[end];
	//	float dist = distP2P(contour->points[begin_index], contour->points[end_index]);
	//	line_lens.push_back(dist);
	//	line_lens_for_sort.push_back(dist);
	//	total_len += dist;
	//}
	//float mean_len = total_len / float(final_line_segs.size());
	//std::sort(line_lens_for_sort.begin(), line_lens_for_sort.end());
	//float median_len = line_lens[line_lens.size()/2];

	//// 确定极小线段的长度相对于线段均值的最大长度比例
	//float T_lenRatio = 0.25;
	//std::vector<bool> is_del(final_line_segs.size(), false);// 用来标记对应的线段是否需要被删除
	//std::vector<int> off_set(final_line_segs.size(), 0);	// 用来标记每个顶点向前平移的位置数量
	//int count_del = 0;	// 用来标记需要删除多少条线段
	//for(int i = 0; i < final_line_segs.size(); ++i){
	//	if( line_lens[i]/mean_len < T_lenRatio ||
	//		line_lens[i]/median_len < T_lenRatio ){
	//		is_del[i] = true;
	//		++count_del;
	//		off_set[i] = -1;
	//		continue;
	//	}
	//	off_set[i] = count_del;
	//}

	//// 执行删除
	//for(int i = 0; i < final_line_segs.size(); ++i)
	//{
	//	if(off_set[i] > 0)
	//	{
	//		final_line_segs[i - off_set[i]] = final_line_segs[i];	// 移动是向前进行的
	//	}
	//}
	//for(int i = 0; i < count_del; ++i)
	//	final_line_segs.pop_back();

	//// 重要：顶点数量化简后，若多边形的顶点数量少于3个，则必须给出错误提示：
	//if(final_line_segs.size() < 3) {
	//	std::cerr << "顶点数量少于3个！" << std::endl;
	//	// 不对顶点进行化简
	//	*final_vertices = *contour;
	//	return;
	//}

	// method1:
	//// step3:
	//// 保留长度超过阈值的线段
	//// step3(a):
	//// 计算线段的长度
	//std::vector<float> line_lens;
	//for(int begin = 0; begin < final_line_segs.size(); ++begin)
	//{
	//	int end = begin == final_line_segs.size()-1 ? 0 : begin+1;

	//	int begin_index/*w.r.t. contour*/ = final_line_segs[begin];
	//	int end_index/*w.r.t. contour*/ = final_line_segs[end];
	//	float dist = distP2P(contour->points[begin_index], contour->points[end_index]);
	//	line_lens.push_back(dist);
	//}

	//// step3(b):
	//// 若有连续几条线段的长度均小于阈值，但这几条线段的长度之和又大于阈值，则合并这几条线段
	//std::vector<int> line_segs(final_line_segs.size(), 0); // 保存当前的线段
	//for(int i = 0; i < final_line_segs.size(); ++i) line_segs[i] = final_line_segs[i];
	//
	//std::vector<bool> is_len_valid(line_segs.size(), false);
	//for(int i = 0; i < line_segs.size(); ++i){
	//	if(line_lens[i] >= T_dist_refine_lines) is_len_valid[i] = true;
	//}

	//// 线段Li,Li+1,...,Li+N-1是连续N条线段，且len(Lj)<T_dist_refine_lines, j = i,i+1,...,i+N-1
	//// 若Len(Li)+Len(Li+1)+Len(Li+j)>T_dist_refine_lines
	//// 则Li...Lj可被当成一条线段
	//std::vector<std::vector<int>> continuous_lines;
	//int first_valid_line_index = -1;
	//for(int i = 0; i < line_lens.size(); ++i)
	//{// 获取第一个长度合法的线段的索引
	//	if(is_len_valid[i])
	//	{
	//		first_valid_line_index = i;
	//		break;
	//	}
	//}
	//if(first_valid_line_index == line_lens.size()-1) {std::cerr << "error: the line segment length threshold is set to be too large!" << std::endl;}
	//int cur_index = first_valid_line_index + 1;
	//while(true)
	//{
	//	if(cur_index == first_valid_line_index){break;}
	//	if(is_len_valid[cur_index])	// 若线段长度合法则继续下去
	//	{
	//		cur_index = cur_index == line_lens.size()-1 ? 0 : cur_index+1;
	//		continue;
	//	}
	//	std::vector<int> v;	// 若线段长度不合法则保存其位置
	//	while(!is_len_valid[cur_index])
	//	{
	//		v.push_back(cur_index);
	//		cur_index = cur_index == line_lens.size()-1 ? 0 : cur_index+1;
	//	}
	//	continuous_lines.push_back(v);
	//}

	//// 对于长度小于阈值的连续线段，检测其长度之和是否大于阈值，如果大于阈值则合并这些线段，否则将这些线段删去
	//std::vector<int> condition(line_segs.size(),0);	// 0: 保留； 1：删除； 2：合并（仅对于那些连续的被标记为2的点）

	//for(int i = 0; i < continuous_lines.size(); ++i)
	//{
	//	if(continuous_lines[i].size() == 1) {condition[continuous_lines[i][0]] = 1; continue;}
	//	float sum_len = 0;
	//	for(int j = 0; j < continuous_lines[i].size(); ++j)
	//	{
	//		sum_len += line_lens[continuous_lines[i][j]];
	//	}
	//	if(sum_len <= T_dist_refine_lines)
	//	{
	//		for(int j = 0; j < continuous_lines[i].size(); ++j)
	//			condition[continuous_lines[i][j]] = 1;
	//	}else{
	//		/*for(int j = 0; j < continuous_lines[i].size(); ++j)
	//			condition[continuous_lines[i][j]] = 2;*/
	//		for(int j = 1; j < continuous_lines[i].size(); ++j)
	//			condition[continuous_lines[i][j]] = 1;
	//	}
	//}
	//// 重新设置线段序列
	//final_line_segs.clear();
	//for(int i = 0; i < line_segs.size(); ++i){
	//	if(condition[i] == 0){
	//		final_line_segs.push_back(line_segs[i]);
	//		final_vertices->push_back(contour->points[line_segs[i]]);	// 这样有一个好处，线段不大可能出现与多边形其它边相交的可能性，且节省了计算量
	//	}
	//}

	//// 重要：顶点数量化简后，若多边形的顶点数量少于3个，则必须给出错误提示：
	//if(final_vertices->size() < 3) std::cerr << "顶点数量少于3个！" << std::endl;

	// step4:
	// 如果相邻线段近似平行，则合并这两个线段
	std::vector<int> line_segs;
	line_segs.clear();
	for(int i = 0; i < final_line_segs.size(); ++i) line_segs.push_back(final_line_segs[i]);

	residual_points = line_segs.size();// 保存当前多边形顶点的数量
	final_line_segs.clear();
	std::vector<int> condition;
	condition.clear();
	for(int i = 0; i < line_segs.size(); ++i) condition.push_back(0);	// 注意0是默认值，默认保留

	for(int i = 0; i < line_segs.size(); ++i)
	{
		Eigen::Vector3f line_i_begin_point, line_i_end_point, prev_line_end_point;
		line_i_begin_point << contour->points[line_segs[i]].x, contour->points[line_segs[i]].y, contour->points[line_segs[i]].z;
		int prev = i == 0 ? line_segs.size()-1 : i-1;
		int next = i == line_segs.size()-1 ? 0 : i+1;
		line_i_end_point << contour->points[line_segs[next]].x, contour->points[line_segs[next]].y, contour->points[line_segs[next]].z;
		prev_line_end_point << contour->points[line_segs[prev]].x, contour->points[line_segs[prev]].y, contour->points[line_segs[prev]].z;

		Eigen::Vector3f line1_dir, line2_dir;
		line1_dir = line_i_begin_point - prev_line_end_point;
		line2_dir = line_i_end_point - line_i_begin_point;
		line1_dir.normalize();
		line2_dir.normalize();
		float cos_val = line1_dir.dot(line2_dir);
		if(cos_val >= cos(T_angle_betLines) && residual_points > 3)// 要保证有足够的点可被删除
		{	
			condition[i] = 1;
			--residual_points;
		}
	}
	final_vertices->clear();final_line_segs.clear();
	for(int i = 0; i < line_segs.size(); ++i){
		if(condition[i]==0){
			final_line_segs.push_back(line_segs[i]);
			final_vertices->push_back(contour->points[line_segs[i]]);
		}
	}
}

bool compair(std::pair<float, int> &pa1, std::pair<float, int> &pa2)
{
	return pa1.first < pa2.first;
}
// 排序：v里面的值不变，indices保存由小到大的值对应的索引
void SimplifyVerticesSize::vectorSort(std::vector<float> &v, std::vector<int> &indices)
{
	std::vector<std::pair<float, int>> pv;
	pv.reserve(v.size());
	std::pair<float, int> pa;
	for(int i = 0; i < v.size(); ++i)
	{
		pa.first = v[i];
		pa.second = i;
		pv.push_back(pa);
	}
	std::sort(pv.begin(), pv.end(), compair);
	indices.clear();
	indices.reserve(v.size());
	for(int i = 0; i < v.size(); ++i)
	{
		indices.push_back(pv[i].second);
	}
}

//void SimplifyVerticesSize::simplifyPolyVertices(MyPoly &poly)
//{
//	// step1: 去除同一个位置上的冗余点
//	cleanPolyVer(poly.vertices_cloud);
//	// step2: 调整顶点排序，使其符合右手法则
//	validateVerticesSeq(poly.vertices_cloud, poly.poly_normal);
//	// step3：计算顶点的法向量
//	computeVerticesNormals(poly.vertices_cloud, poly.vertices_normals, poly.poly_normal);
//	// step4: 计算初始的线段序列
//	constructInitLineSegs(poly.vertices_cloud, poly.vertices_normals, poly.poly_normal, poly.line_segs);
//	// step5: 精化线段序列
//	std::vector<int> final_line_segs;
//	refineLineSegs(poly.vertices_cloud, poly.final_vertices/*这是产出物*/, poly.line_segs, final_line_segs);
//}
