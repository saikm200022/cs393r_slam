//========================================================================
//  This software is free: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License Version 3,
//  as published by the Free Software Foundation.
//
//  This software is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  Version 3 in the file COPYING that came with this distribution.
//  If not, see <http://www.gnu.org/licenses/>.
//========================================================================
/*!
\file    slam.cc
\brief   SLAM Starter Code
\author  Joydeep Biswas, (C) 2019
*/
//========================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "shared/math/geometry.h"
#include "shared/math/math_util.h"
#include "shared/util/timer.h"

#include "slam.h"

#include "vector_map/vector_map.h"

using namespace math_util;
using Eigen::Affine2f;
using Eigen::Rotation2Df;
using Eigen::Translation2f;
using Eigen::Vector2f;
using Eigen::Vector2i;
using std::cout;
using std::endl;
using std::string;
using std::swap;
using std::vector;
using vector_map::VectorMap;

namespace slam {

SLAM::SLAM() :
    prev_odom_loc_(0, 0),
    prev_odom_angle_(0),
    odom_initialized_(false) {}

struct Pose SLAM::MostLikelyPose()
{
  // If the first pose, the cumulative transformations are all 0 therefore best pose is dx = 0, dy = 0, dtheta = 0
   if (cumulative_transform.delta_x == 0 && cumulative_transform.delta_y == 0 && cumulative_transform.delta_theta == 0)
  {
    struct Pose first_pose = {
      0, 0, 0, 0
    };
    return first_pose;
  }

  // Query 3D cube for most likely pose
  // Convert pose to map frame

  float best_x = 0.0;
  float best_y = 0.0;
  float best_theta = 0.0;
  float most_likely = -1.0;

  // Consult cube to find displacement that results in most likely pose
  for (int pixel_x = 0; pixel_x < x_width; pixel_x++)
  {
    for (int pixel_y = 0; pixel_y < y_width; pixel_y++)
    {
      for (int pixel_theta = 0; pixel_theta < theta_width; pixel_theta++)
      {
          if (most_likely < cube[pixel_x][pixel_y][pixel_theta])
          {
            best_x = pixel_x * x_incr;
            best_y = pixel_y * y_incr;
            best_theta = pixel_theta * theta_incr;
            most_likely = cube[pixel_x][pixel_y][pixel_theta];
          }
      }
    }
  }
  struct Pose new_pose = {
    best_x, best_y, best_theta, most_likely
  };
  return new_pose;
}

void SLAM::GetPose(Eigen::Vector2f* loc, float* angle) const {
  // Return the latest pose estimate of the robot.
  *loc = Vector2f(0, 0);
  *angle = 0;
  
  if (best_pose.delta_x == 0 && best_pose.delta_y == 0 && best_pose.delta_theta == 0)
    return;
  
  // Find best pose with respect to pose 1 (using the cumulative_transform variable)
  *loc = Vector2f(cumulative_transform.delta_x, cumulative_transform.delta_y);
  *angle += cumulative_transform.delta_theta;
}

Eigen::Matrix2f SLAM::GetRotationMatrix (const float angle) {
  Eigen::Matrix2f rot;
  rot(0,0) = cos(angle);
  rot(0,1) = -sin(angle);
  rot(1,0) = sin(angle);
  rot(1,1) = cos(angle);
  return rot;
}

vector<Vector2f> SLAM::GetScanPointCloud(const vector<float>& ranges,
                        float range_min,
                        float range_max,
                        float angle_min,
                        float angle_max)
  {
    vector<Vector2f> points;

    float angle = angle_min;
    int range_index = 0;
    float angle_increment = (angle_max - angle_min) / ranges.size();

    // Convert laser scans to points
    while (angle <= angle_max) 
    {
      Vector2f point;
      float range = ranges[range_index];
      point[0] = laser_off + range * cos(angle);
      point[1] = range * sin(angle);

      points.push_back(point);
      angle += angle_increment;
      range_index += 1;
    }
    return points;
  }

// x_2 = current
// x_1 previous

// estimate x_2 by trying out all x_1 + u
// Scan of x_2 ground truth and we compare that to the transformed x_1 laser scan

// Sum of gaussians approach
void SLAM::ObserveLaser(const vector<float>& ranges,
                        float range_min,
                        float range_max,
                        float angle_min,
                        float angle_max) {
  // A new laser scan has been observed. Decide whether to add it as a pose
  // for SLAM. If decided to add, align it to the scan from the last saved pose,
  // and save both the scan and the optimized pose.
  
  // For initial pose
  if (previous_pose.delta_x == -1000 && previous_pose.delta_y == -1000 && previous_pose.delta_theta == -1000)
  {
    previous_pose.delta_x = 0;
    previous_pose.delta_y = 0;
    previous_pose.delta_theta = 0;
    vector<Vector2f> points = GetScanPointCloud(ranges, range_min, range_max, angle_min, angle_max);
    previous_scan = points;
    return;
  }

  // If conditions met (0.5 m dist or 45 degrees rotated)
  if (distance_travelled <= 0 || angle_travelled <= 0)
  { 
    float current_image[x_image_width][y_image_width];
    vector<Vector2f> points = GetScanPointCloud(ranges, range_min, range_max, angle_min, angle_max);

    // Create Image by Centering gaussians around previous laser scan points
    for (auto& point : previous_scan)
    {
      // pixel_x and pixel_y are bin indices in the image
      for (int pixel_x = 0; pixel_x < x_image_width; pixel_x++)
      {
        for (int pixel_y = 0; pixel_y < y_image_width; pixel_y++)
        {
          // Take left corner of every pixel
          float x = pixel_x * x_image_incr + point[0];
          float y = pixel_y * y_image_incr + point[1];

          float prob = 1.0;
          float std_dev = k_1 * pow(pow(x, 2) + pow(y, 2), 0.5);

          // Decoupled evaluation of multivariate gaussian where product is taken along x and y
          prob *= exp(-0.5 * pow((x - point[0])/std_dev, 2));
          prob *= exp(-0.5 * pow((y - point[1])/std_dev, 2));

          // Sum of Gaussians 
          current_image[pixel_x][pixel_y] += prob;
        }
      }
    }

    for (int pixel_x = 0; pixel_x < x_width; pixel_x++)
    {
      for (int pixel_y = 0; pixel_y < y_width; pixel_y++)
      {
        for (int pixel_theta = 0; pixel_theta < theta_width; pixel_theta++)
        {
          float dx = pixel_x * x_incr;
          float dy = pixel_y * y_incr;
          float dtheta = pixel_theta * theta_incr;

          vector<Vector2f> shifted_points;
          for (auto& point : points)
          {
            // Shift point to previous pose
            float new_x = cos(-previous_pose.delta_theta) * point[0] - sin(-previous_pose.delta_theta) * point[1] - previous_pose.delta_x;
            float new_y = sin(-previous_pose.delta_theta) * point[0] + cos(-previous_pose.delta_theta) * point[1] - previous_pose.delta_y;

            // Forward shift by dx, dy, dtheta
            float final_x = cos(dtheta) * new_x - sin(dtheta) * new_y - dx;
            float final_y = sin(dtheta) * new_x + cos(dtheta) * new_y - dy;

            shifted_points.push_back(Vector2f(final_x, final_y));
          }
          
          // Calculate score by using gaussian values from image
          float score = 1.0;
          for (auto& point : shifted_points)
          {
            int point_pixel_x = point[0] / x_incr;
            int point_pixel_y = point[1] / y_incr;
            score *= current_image[point_pixel_x][point_pixel_y];
          }

          cube[pixel_x][pixel_y][pixel_theta] *= score;
        }
      }
    }

    // Update all necessary global variables
    previous_scan = points;
    best_pose = MostLikelyPose();
    cumulative_transform.delta_x += best_pose.delta_x;
    cumulative_transform.delta_y += best_pose.delta_y;
    cumulative_transform.delta_theta = best_pose.delta_theta;
    previous_pose = best_pose;
  }

  // Transform current robot scan to previous pose
  // Create image for current scan
  // For all points:
  // Get value at that point from the image
  // Calculate product
  // Update previous image to current image
}

void SLAM::ObserveOdometry(const Vector2f& odom_loc, const float odom_angle) {
  if (!odom_initialized_) {
    prev_odom_angle_ = odom_angle;
    prev_odom_loc_ = odom_loc;
    odom_initialized_ = true;
    return;
  }
  // Keep track of odometry to estimate how far the robot has moved between 
  // poses.

  // Not sure this is correct (maybe too simple?)
  Vector2f translation_hat = odom_loc - prev_odom_loc_;
  float angle_hat = odom_angle - prev_odom_angle_;

  // Try out all possible delta x, y, theta
  for (int pixel_x = 0; pixel_x < x_width; pixel_x++)
  {
    for (int pixel_y = 0; pixel_y < y_width; pixel_y++)
    {
      for (int pixel_theta = 0; pixel_theta < theta_width; pixel_theta++)
      {
        float dx = pixel_x * x_incr;
        float dy = pixel_y * y_incr;
        float dtheta = pixel_theta * theta_incr;

        float prob = 1.0;
        float std_dev = k_1 * pow(pow(translation_hat[0], 2) + pow(translation_hat[1], 2), 0.5) + k_2 * (pow(angle_hat, 2));

        // Decoupled evaluation of multivariate gaussian where product is taken along x, y, and theta
        prob *= exp(-0.5 * pow((dx - translation_hat[0])/std_dev, 2));
        prob *= exp(-0.5 * pow((dy - translation_hat[1])/std_dev, 2));
        prob *= exp(-0.5 * pow((dtheta - angle_hat)/std_dev, 2));

        // Motion model probability
        cube[pixel_x][pixel_y][pixel_theta] *= prob;
      }
    }
  }
}

vector<Vector2f> SLAM::GetMap() {
  vector<Vector2f> map;
  // Reconstruct the map as a single aligned point cloud from all saved poses
  // and their respective scans.
  return map;
}

}  // namespace slam
