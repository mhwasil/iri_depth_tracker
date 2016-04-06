/*
 *
 *  Created on: May 2, 2013
 *      Author: shusain
 */
#include "ros/ros.h"
#include "std_msgs/String.h"
#include <pcl_ros/point_cloud.h>
#include <sensor_msgs/PointCloud.h>
#include <pcl/point_types.h>
#include <boost/foreach.hpp>
#include <boost/progress.hpp>
#include <geometry_msgs/PoseStamped.h>

#include <image_transport/image_transport.h>
#include <cv_bridge/cv_bridge.h>
#include <sensor_msgs/image_encodings.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>
#include "Optimal_affine_tracking_3d16_fast_realtime.h"
#include "Optimal_affine_tracking_3d16_fast_realtime_initialize.h"
#include "Optimal_affine_tracking_3d16_fast_realtime_terminate.h"

#include "rt_nonfinite.h"
#include "compute_prob1.h"
//#include "LoadKinectMesh_realtime.h"
#include "init_variables.h"

#include <vector>
#include <algorithm>
//#include <omp.h>
#define PI 3.141592653589793238462643383

ros::Publisher g_pose_publisher;

real_T Ixyz[504063];
real_T corner_p[12];
static real_T mean_img[3888];
static real_T point_matrix[3888];
static real_T AR_velocity[2250000];
static real_T X_par[2250000];
static real_T dw_dp[15552];
static real_T tracked_images[38880000];
static real_T X_par_pred[2250000];
real_T centroid[3];
int32_T t;
real_T Aff_matrix[9];
int32_T i0;
int32_T i1;
int32_T i2;
int32_T row;
int32_T column;
int i;
typedef pcl::PointCloud<pcl::PointXYZRGB> PointCloud;
real_T p[12];
cv::Rect region_of_interest = cv::Rect(100, 50, 441, 381);
cv::Mat color_img(480, 640, CV_8UC3);
cv::Mat xyz_img(480, 640, CV_32FC3);
cv::Mat color_img_roi;
cv::Mat xyz_img_roi;
// for video recording
cv::Mat gray_img_roi;
// end for video recording

static float c_x=0,c_y=0;
static float v_x=0,v_y=0;
int mouse_loc_x[4];
int mouse_loc_y[4];
int mosue_click_time=0;
real_T center_y;
real_T center_x;
double xx=0,yy=0,zz=0,dxy=0,xx_C=0,yy_C=0,zz_C=0,dxy2=0,xx2=0,yy2=0,zz2=0;

void onMouse(int event, int x, int y, int flags, void* data)
{
	if (event != cv::EVENT_LBUTTONDOWN)
	return;
	mouse_loc_x[mosue_click_time] = x;
	mouse_loc_y[mosue_click_time] = y;

	//printf("%d %d\n",x,y);
	printf("%d %d\n",mouse_loc_x[mosue_click_time],mouse_loc_y[mosue_click_time]);
	++mosue_click_time;
}


void processPointCloud(const PointCloud::ConstPtr& msg) {

	boost::progress_timer timer;
	//static real_T mesh[921600];
	float x,y,z,ptrgb,yaw,pitch,roll;//,fr,fg,fb;
	uint8_t r,g,b;
	uint32_t rgb;
	i=0;

//  	cvNamedWindow("color image",CV_WINDOW_NORMAL);
//  	cvMoveWindow("color image",100,100);
  	cvNamedWindow("range image");
//  	cvMoveWindow("range image",545,100);
/*	cvNamedWindow("color image 2");
	cvMoveWindow("color image 2",100,100+415);
	cvNamedWindow("range image 2");
	cvMoveWindow("range image 2",545,100+415);
*/
	BOOST_FOREACH (const pcl::PointXYZRGB& pt, msg->points)
		{
			ptrgb = pt.rgb;
			rgb = *reinterpret_cast<int*>(&ptrgb);
			r = (rgb >> 16) & 0x0000ff;
			g = (rgb >> 8)  & 0x0000ff;
			b = (rgb)       & 0x0000ff;
			x = (float)pt.x+2;
			y = (float)pt.y+2;
			z = (float)pt.z;

			if (x != x || y != y || z != z || z == 0)
			{	x=0; y=0; z=0; r=0; g=0; b=0;
			}
			column = (int)(i / msg->width);
			row    = i % (msg->width);

/*			color_img.at<cv::Vec3b>(column,row)[0] = b;
			color_img.at<cv::Vec3b>(column,row)[1] = g;
			color_img.at<cv::Vec3b>(column,row)[2] = r;*/
			xyz_img.at<cv::Vec3f>(column,row)[0] = x;
			xyz_img.at<cv::Vec3f>(column,row)[1] = y;
			xyz_img.at<cv::Vec3f>(column,row)[2] = z;
			++i;
		}

//	color_img_roi = color_img(region_of_interest);
	xyz_img_roi = xyz_img(region_of_interest);

//	cv::cvtColor( color_img_roi, gray_img_roi, CV_RGB2GRAY );
//	cv::cvtColor( gray_img_roi, gray_img_roi, CV_GRAY2RGB );

    static int callback_counter_ = 0;
    int k,i3;
    if (callback_counter_ == 0) {
    	/* % Object template initialization */
//#pragma omp for schedule (dynamic, 8)
#pragma omp parallel for
    	  for (k = 0; k < 441; k++) {
    	    for (i3 = 0; i3 < 381; i3++) {
    	      Ixyz[i3 + 381 * k] = xyz_img_roi.at<cv::Vec3f>(i3,k)[0];
    	      Ixyz[168021 + (i3 + 381 * k)] = xyz_img_roi.at<cv::Vec3f>(i3,k)[1];
    	      Ixyz[336042 + (i3 + 381 * k)] = xyz_img_roi.at<cv::Vec3f>(i3,k)[2];
    	    }
    	  }
      real_T ptx[4];
      real_T pty[4];


//  	cv::imshow( "color image", color_img_roi );
	std::vector<cv::Mat> split_xyz(3,cv::Mat(xyz_img_roi.size(),CV_64FC1));
	cv::split(xyz_img_roi, split_xyz);
	static double min;
	static double max;
	cv::minMaxIdx(split_xyz.at(2), &min, &max);
	static cv::Mat adjMap;
	cv::convertScaleAbs(split_xyz.at(2), adjMap, 255 / max);

	cv::cvtColor( adjMap, adjMap, CV_GRAY2RGB );
  	cv::imshow( "range image", adjMap );
   	cv::setMouseCallback("range image", onMouse);
	////////////////////for video recording
	//cv::imshow( "color image 2", gray_img_roi );
	//cv::imshow( "range image 2", adjMap );
	//////////////////// end for video recording
   	cv::waitKey(0);

 	ptx[0] =  (double)mouse_loc_y[0];
  	ptx[1] =  (double)mouse_loc_y[1];
  	ptx[2] =  (double)mouse_loc_y[2];
  	ptx[3] =  (double)mouse_loc_y[3];

  	pty[0] =  (double)mouse_loc_x[0];
  	pty[1] =  (double)mouse_loc_x[1];
  	pty[2] =  (double)mouse_loc_x[2];
  	pty[3] =  (double)mouse_loc_x[3];

      init_variables(Ixyz, ptx, pty, X_par_pred, tracked_images, dw_dp, X_par,
                       AR_velocity, point_matrix, mean_img, corner_p, &center_x,
                       &center_y);
      ++callback_counter_;

    }
    else {

  	  for (k = 0; k < 441; k++) {
  	    for (i3 = 0; i3 < 381; i3++) {
  	      Ixyz[i3 + 381 * k] = xyz_img_roi.at<cv::Vec3f>(i3,k)[0];
  	      Ixyz[168021 + (i3 + 381 * k)] = xyz_img_roi.at<cv::Vec3f>(i3,k)[1];
  	      Ixyz[336042 + (i3 + 381 * k)] = xyz_img_roi.at<cv::Vec3f>(i3,k)[2];
  	    }
  	  }
    	/*     %% Particle propagation and likelihood computation  */
  	compute_prob1(X_par, X_par_pred, AR_velocity, dw_dp, (real_T)t + 2.0,
  	                  center_x, center_y, Ixyz, point_matrix, mean_img,
  	                  tracked_images, Aff_matrix, centroid);
    	++t;
    	/*     %% Display the tracking results */
//#pragma omp for// schedule (dynamic, 8)
    	for (i0 = 0; i0 < 3; i0++) {
    	 for (i1 = 0; i1 < 4; i1++) {
    	   p[i0 + 3 * i1] = 0.0;
    	   for (i2 = 0; i2 < 3; i2++) {
    		 p[i0 + 3 * i1] += Aff_matrix[i0 + 3 * i2] * corner_p[i2 + 3 * i1];
    	   }
    	 }
    	}
    	++callback_counter_;
    }

	geometry_msgs::PoseStamped pose_stamped;
	//pose_stamped.header.stamp.sec = callback_counter_;
	pose_stamped.header.frame_id = "/camera_rgb_optical_frame";

	////////////////////for video recording
	//cv::imshow( "color image 2", gray_img_roi );
	//////////////////// end for video recording

	// extract the centroid of the 4 corners of the bounding box
	c_x =  (p[4]+p[7]+p[10]+p[1])/4;
	c_y =  (p[3]+p[6]+p[9]+p[0])/4;

/*	cv::line(color_img_roi, cv::Point(p[4],p[3]), cv::Point(p[1],p[0]), cv::Scalar(0, 0, 255), 2, CV_AA);
	cv::line(color_img_roi, cv::Point(p[7],p[6]), cv::Point(p[4],p[3]), cv::Scalar(0, 0, 255), 2, CV_AA);
	cv::line(color_img_roi, cv::Point(p[10],p[9]), cv::Point(p[7],p[6]), cv::Scalar(0, 0, 255), 2, CV_AA);
	cv::line(color_img_roi, cv::Point(p[1],p[0]), cv::Point(p[10],p[9]), cv::Scalar(0, 0, 255), 2, CV_AA);
	cv::line(color_img_roi, cv::Point(c_x,c_y), cv::Point((p[7]+p[4])/2,(p[6]+p[3])/2), cv::Scalar(0, 0, 255), 2, CV_AA);*/

	v_x = c_x-((p[7]+p[4])/2);
	v_y = c_y-((p[6]+p[3])/2);

	//
//	cv::imshow( "color image", color_img_roi );
	// show depth image
	std::vector<cv::Mat> split_xyz(3,cv::Mat(xyz_img_roi.size(),CV_64FC1));
	cv::split(xyz_img_roi, split_xyz);
	static double min;
	static double max;
	cv::minMaxIdx(split_xyz.at(2), &min, &max);
	static cv::Mat adjMap;
	cv::convertScaleAbs(split_xyz.at(2), adjMap, 255 / max);
	cv::cvtColor( adjMap, adjMap, CV_GRAY2RGB );
	////////////////////for video recording
	//cv::imshow( "range image 2", adjMap );
	//////////////////// end for video recording
	cv::line(adjMap, cv::Point(p[4],p[3]), cv::Point(p[1],p[0]), cv::Scalar(0, 0, 255), 2, CV_AA);
	cv::line(adjMap, cv::Point(p[7],p[6]), cv::Point(p[4],p[3]), cv::Scalar(0, 0, 255), 2, CV_AA);
	cv::line(adjMap, cv::Point(p[10],p[9]), cv::Point(p[7],p[6]), cv::Scalar(0, 0, 255), 2, CV_AA);
	cv::line(adjMap, cv::Point(p[1],p[0]), cv::Point(p[10],p[9]), cv::Scalar(0, 0, 255), 2, CV_AA);
	cv::line(adjMap, cv::Point(c_x,c_y), cv::Point((p[7]+p[4])/2,(p[6]+p[3])/2), cv::Scalar(0, 0, 255), 2, CV_AA);
	cv::line(adjMap, cv::Point(c_x,c_y), cv::Point((p[1]+p[4])/2,(p[0]+p[3])/2), cv::Scalar(0, 255, 0), 2, CV_AA);
	cv::imshow( "range image", adjMap );
	cv::waitKey(1);


	pose_stamped.pose.position.x = centroid[0]-2;
	pose_stamped.pose.position.y = centroid[1]-2;
	pose_stamped.pose.position.z = centroid[2];
	// compute image plane rotation
	pose_stamped.pose.orientation.y = atan(v_y/v_x);
	//yaw = atan(v_y/v_x);

	// compute out of plane rotation
	cv::Point pt_C((p[4]+p[7]+p[10]+p[1])/4,(p[3]+p[6]+p[9]+p[0])/4);
	xx_C = xyz_img_roi.at<cv::Vec3f>(pt_C)[0];
	yy_C = xyz_img_roi.at<cv::Vec3f>(pt_C)[1];
	zz_C = xyz_img_roi.at<cv::Vec3f>(pt_C)[2];
	cv::Point pt((p[7]+p[4])/2,(p[6]+p[3])/2);
	xx = xyz_img_roi.at<cv::Vec3f>(pt)[0];
	yy = xyz_img_roi.at<cv::Vec3f>(pt)[1];
	zz = xyz_img_roi.at<cv::Vec3f>(pt)[2];
	dxy = pow(xx-xx_C,2)+pow(yy-yy_C,2);
	//pitch = atan((zz-zz_C)/(sqrt(dxy)));
	pose_stamped.pose.orientation.x = atan((zz-zz_C)/(sqrt(dxy)));
	//pose_stamped.pose.orientation.y = atan(zz-zz_C)*180/PI;
	//roll = 0;
	cv::Point pt2((p[1]+p[4])/2,(p[0]+p[3])/2);
	xx2 = xyz_img_roi.at<cv::Vec3f>(pt2)[0];
	yy2 = xyz_img_roi.at<cv::Vec3f>(pt2)[1];
	zz2 = xyz_img_roi.at<cv::Vec3f>(pt2)[2];
	dxy2 = pow(xx2-xx_C,2)+pow(yy2-yy_C,2);
	pose_stamped.pose.orientation.z = atan((zz2-zz_C)/(sqrt(dxy2)));

	g_pose_publisher.publish(pose_stamped);

}

int main(int argc, char **argv)
{
	ros::init(argc, argv, "depth_tracker");
	ros::NodeHandle n;
    Optimal_affine_tracking_3d16_fast_realtime_initialize();
    g_pose_publisher = n.advertise<geometry_msgs::PoseStamped>("pose_surface", 50);
    //std::cout << "Number processors: " << omp_get_num_procs() << std::endl;
	/* % Execute tracking */
    ros::Subscriber sub = n.subscribe<PointCloud>("/camera/depth_registered/points", 5, processPointCloud);
    
    Optimal_affine_tracking_3d16_fast_realtime_terminate();
    ros::spin();
    return 0;
}
