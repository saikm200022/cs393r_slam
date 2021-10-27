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

SLAM::Matrix2f ParticleFilter::GetRotationMatrix (const float angle) {
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
    // Create Image
    slam::Pose current_image[x_image_width][y_image_width];
    vector<Vector2f> points;
    while (angle <= angle_max) {
      Vector2f point;
      float range = ranges[range_index];
      point[0] = 0.2 + range * cos(angle);
      point[1] = range * sin(angle);

      points.push_back(points);
      angle += msg.angle_increment;
      range_index += 1;
    }

    for (&auto point : points)
    {
      // Find probability for "pixels" in image
      for (float d_x = 0; d_x <= x_image_max; d_x += x_image_incr)
      {
        for (float d_y = 0; d_y <= y_image_max; d_y += y_image_incr)
        {
          // Fit Gaussian over point to compute probability
          // Assign Probability values
          float x = Math.pow(Math.pow(d_x, 2) + Math.pow(d_y, 2), 0.5);
          current_image[d_x][d_y] = (1 / (s * sqrt(2 * M_PI) )) * exp(-0.5 * pow((Math.pow(x - 0)/s, obs_likelihood_stdv)));
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
          while (angle <= angle_max) {
            Vector2f point;
            float range = ranges[range_index];
            point[0] = 0.2 + range * cos(angle);
            point[1] = range * sin(angle);

            // Reverse Transform to previous pose
            float new_x = cos(-d_theta) * point[0] - sin(-d_theta) * point[1] - d_x;
            float new_y = sin(-d_theta) * point[0] + cos(-d_theta) * point[1] - d_y;

            transformed_points.push_back(Vector2f(new_x, new_y));
            angle += msg.angle_increment;
            range_index += 1;
          }

          // TODO: Compare Previous Image and Transformed Scan
          float probability = 1.0;
          for (auto& point : transformed_points)
          {
            probability *= image[point[0]][point[1]];
          }

          cube[d_x][d_y][d_theta] = probability;
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

  // Try out all possible delta x, y, theta

}

// p(x_i | x_j, u_i, u_j) = motion model

vector<Vector2f> SLAM::GetMap() {
  vector<Vector2f> map;
  // Reconstruct the map as a single aligned point cloud from all saved poses
  // and their respective scans.
  return map;
}

}  // namespace slam
