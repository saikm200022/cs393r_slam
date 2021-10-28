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

void SLAM::GetPose(Eigen::Vector2f* loc, float* angle) const {
  // Return the latest pose estimate of the robot.
  *loc = Vector2f(0, 0);
  *angle = 0;

  // Query 3D cube for most likely pose
  // Convert pose to map frame
}

Eigen::Matrix2f SLAM::GetRotationMatrix (const float angle) {
  Eigen::Matrix2f rot;
  rot(0,0) = cos(angle);
  rot(0,1) = -sin(angle);
  rot(1,0) = sin(angle);
  rot(1,1) = cos(angle);
  return rot;
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

  // If conditions met (0.5 m dist or 45 degrees rotated)
  if (distance_travelled <= 0 || angle_travelled <= 0)
  { 
    slam::Pose current_image[x_image_width][y_image_width];
    vector<Vector2f> points;

    float angle = angle_min;
    int range_index = 0;
    float angle_increment = (angle_max - angle_min) / ranges.size();

    // Convert laser scans to points
    while (angle <= angle_max) 
    {
      Vector2f point;
      float range = ranges[range_index];
      point[0] = laser_offset + range * cos(angle);
      point[1] = range * sin(angle);

      points.push_back(point);
      angle += angle_increment;
      range_index += 1;
    }

    // Create Image by Centering gaussian around laser scan points
    for (auto& point : points)
    {
      // pixel_x and pixel_y are bin indices in the image
      for (int pixel_x = 0; pixel_x < x_image_width; pixel_x++)
      {
        for (int pixel_y = 0; pixel_y < y_image_width; pixel_y++)
        {
          // Take left corner of every pixel
          float x = pixel_x * x_image_incr;
          float y = pixel_y * y_image_incr;

          float prob = 1.0;
          float std_dev = k_1 * pow(pow(x, 2) + pow(y, 2), 0.5);

          // Decoupled evaluation of multivariate gaussian where product is taken along x and y
          prob *= exp(-0.5 * Math.pow((x - point[0])/obs_likelihood_stdv, 2));
          prob *= exp(-0.5 * Math.pow((y - point[1])/obs_likelihood_stdv, 2));

          // Sum of Gaussians 
          current_image[pixel_x][pixel_y] += prob;
        }
      }
    }

    // For all perturbations along delta x, delta y, delta theta (x_1 + u)
    for (float d_x = 0; d_x <= x_max; d_x += x_incr)
    {
      for (float d_y = 0; d_y <= y_max; d_y += y_incr)
      {
        for (float d_theta = 0; d_theta <= theta_max; d_theta += theta_incr)
        {
          Eigen::Matrix2f rot = GetRotationMatrix(d_theta);
          float angle = angle_min;
          int range_index = 0;

          // Transform laser scan points based on perturbations
          vector<Vector2f> transformed_points;
          float angle_increment = (angle_max - angle_min) / num_ranges;
          while (angle <= angle_max) {
            Vector2f point;
            float range = ranges[range_index];
            point[0] = 0.2 + range * cos(angle);
            point[1] = range * sin(angle);

            // Reverse Transform to previous pose
            float new_x = cos(-d_theta) * point[0] - sin(-d_theta) * point[1] - d_x;
            float new_y = sin(-d_theta) * point[0] + cos(-d_theta) * point[1] - d_y;

            transformed_points.push_back(Vector2f(new_x, new_y));
            angle += angle_increment;
            range_index += 1;
          }

          // Compare Previous Image and Transformed Scan
          float probability = 1.0;
          for (auto& point : transformed_points)
          {
            probability *= image[point[0]][point[1]];
          }

          cube[d_x / x_incr][d_y / y_incr][d_theta / theta_incr] = probability;
        }
      }
    }
    image = current_image;
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
  
  // T_2^Odom - T_1^Odom
  Vector2f delta_odom = odom_loc - prev_odom_loc;
  // R(-Theta_1^Odom)
  Eigen::Matrix2f rot = GetRotationMatrix(-prev_odom_angle);
  // R(theta_1^Map) * delta T_base_link
  Vector2f delta_T_bl = rot * delta_odom;
  // T_2^Map = T_1^Map + ...
  Vector2f translation = GetRotationMatrix(particle.angle) * delta_T_bl;

  float delta_theta_bl = odom_angle - prev_odom_angle;

  // Try out all possible delta x, y, theta
  for (float d_x = 0; d_x <= x_max; d_x += x_incr)
  {
    for (float d_y = 0; d_y <= y_max; d_y += y_incr)
    {
      for (float d_theta = 0; d_theta <= theta_max; d_theta += theta_incr)
      {
        float prob = 1.0;
        float std_dev = k_1 * pow(pow(x, 2) + pow(y, 2), 0.5) + k_2 * (delta_T_bl.norm());

        // Decoupled evaluation of multivariate gaussian where product is taken along x, y, and theta
        prob *= exp(-0.5 * Math.pow((d_x - translation[0]])/obs_likelihood_stdv, 2));
        prob *= exp(-0.5 * Math.pow((d_y - translation[1])/obs_likelihood_stdv, 2));
        prob *= exp(-0.5 * Math.pow((d_theta - delta_theta_bl)/obs_likelihood_stdv, 2));

        // Sum of Gaussians 
        cube[d_x / x_incr][d_y / y_incr][d_theta / theta_incr] *= prob;
      }
    }
  }
}

// p(x_i | x_j, u_i, u_j) = motion model

vector<Vector2f> SLAM::GetMap() {
  vector<Vector2f> map;
  // Reconstruct the map as a single aligned point cloud from all saved poses
  // and their respective scans.
  return map;
}

}  // namespace slam
