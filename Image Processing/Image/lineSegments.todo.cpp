#include "lineSegments.h"
#include <math.h>

////////////////////////////
// Image processing stuff //
////////////////////////////
float OrientedLineSegment::length(void) const
{
	return sqrt(pow((float)(this->x1 - this->x2), 2) + pow((float)(this->y1 - this->y2), 2));
}
float OrientedLineSegment::distance(const int& x,const int& y) const
{
	float x1_ = (float)this->x1, x2_ = (float)this->x2, x_ = (float)x;
	float y1_ = (float)this->y1, y2_ = (float)this->y2, y_ = (float)y;
	float dx = x2_ - x1_, dy = y2_ - y1_;

	float rate = ((x_ - x1_) * dx + (y_ - y1_) * dy) / (dx * dx + dy * dy);

	if (rate >= 0.0f && rate <= 1.0f) {
		dx = x_ - (x1_ + rate * dx);
		dy = y_ - (y1_ + rate * dy);
	}
	else {
		dx = (rate < 0.0f) ? x_ - x1_ : x_ - x2_;
		dy = (rate < 0.0f) ? y_ - y1_ : y_ - y2_;
	}

	return sqrt(dx * dx + dy * dy);
}
void OrientedLineSegment::getPerpendicular(float& x,float &y) const
{
	x = y =0;
	float dx = (float)this->x1 - (float)this->x2;
	float dy = (float)this->y1 - (float)this->y2;

	x = -dy / sqrt(dx * dx + dy * dy);
	y = dx / sqrt(dx * dx + dy * dy);
}

void  OrientedLineSegment::GetSourcePosition(const OrientedLineSegment& source,const OrientedLineSegment& destination,
											 const int& targetX,const int& targetY,
											 float& sourceX,float& sourceY)
{
	sourceX = sourceY = 0;
	int tar_p_x = destination.x2, tar_p_y = destination.y2;
	int tar_q_x = destination.x1, tar_q_y = destination.y1;

	int q_p_x = tar_q_x - tar_p_x, q_p_y = tar_q_y - tar_p_y;
	int t_p_x = targetX - tar_p_x, t_p_y = targetY - tar_p_y;
	float u = (float)(t_p_x * q_p_x + t_p_y * q_p_y) / (float)(q_p_x * q_p_x + q_p_y * q_p_y);
	
	float p_qp_x, p_qp_y;
	destination.getPerpendicular(p_qp_x, p_qp_y);
	float v = (float)t_p_x * p_qp_x + (float)t_p_y * p_qp_y;

	float src_q_p_x = source.x1 - source.x2, src_q_p_y = source.y1 - source.y2;
	src_q_p_x *= u;
	src_q_p_y *= u;

	source.getPerpendicular(p_qp_x, p_qp_y);
	src_q_p_x += v * p_qp_x;
	src_q_p_y += v * p_qp_y;

	sourceX = source.x2 + src_q_p_x;
	sourceY = source.y2 + src_q_p_y; 
}
