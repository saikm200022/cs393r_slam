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
#include <stdexcept>


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

void SLAM::PrintImage(float image[x_image_width][y_image_width])
{
  printf("Printing Image\n");
  for (int pixel_x = 0; pixel_x < x_image_width; pixel_x++)
  {
    for (int pixel_y = 0; pixel_y < y_image_width; pixel_y++)
      printf("[ %f ]", image[pixel_x][pixel_y]);
    
    printf("\n");
  }
}

void SLAM::PrintCube(int num_elems)
{
  printf("Printing Cube\n");
  int printed = 0;
  for (int pixel_x = 0; pixel_x < x_width; pixel_x++)
  {
    for (int pixel_y = 0; pixel_y < y_width; pixel_y++)
    {
      for (int pixel_theta = 0; pixel_theta < theta_width; pixel_theta++)
      {
        float dx = pixel_x * x_incr;
        float dy = pixel_y * y_incr;
        float dtheta = pixel_theta * theta_incr;
        printf("dx: %f dy: %f dtheta: %f --> %f\n", dx, dy, dtheta, cube[pixel_x][pixel_y][pixel_theta]);
        printed++;

        if (printed > num_elems)
          return;
      }
    }
  }
}

void SLAM::InitializeImage(float image[x_image_width][y_image_width])
{
  for (int pixel_x = 0; pixel_x < x_image_width; pixel_x++)
    for (int pixel_y = 0; pixel_y < y_image_width; pixel_y++)
      image[pixel_x][pixel_y] = 0.0;
}

void SLAM::ReinitializeCube()
{
  for (int pixel_x = 0; pixel_x < x_width; pixel_x++)
    for (int pixel_y = 0; pixel_y < y_width; pixel_y++)
      for (int pixel_theta = 0; pixel_theta < theta_width; pixel_theta++)
          cube[pixel_x][pixel_y][pixel_theta] = 1.0;
}

struct Pose SLAM::MostLikelyPose()
{
  // If the first pose, the cumulative transformations are all 0 therefore best pose is dx = 0, dy = 0, dtheta = 0
   if (previous_pose.delta_x == -1000 && previous_pose.delta_y == -1000 && previous_pose.delta_theta == -1000)
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
          if (most_likely <= cube[pixel_x][pixel_y][pixel_theta])
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

void SLAM::GetBoundingBox(float bounds[4])
{
  float xmin = 10000000000000000;
  float ymin = 10000000000000000;
  float xmax = -1000000000000000;
  float ymax = -1000000000000000;

  for (auto& point : previous_scan)
  {
    xmin = std::min(point[0], xmin);
    ymin = std::min(point[1], ymin);
    xmax = std::max(point[0], xmax);
    ymax = std::max(point[1], ymax);
  }

  bounds[0] = xmin;
  bounds[1] = xmax;
  bounds[2] = ymin;
  bounds[3] = ymax;
  // printf("BOUNDS: %f %f %f %f", xmin, xmax, ymin, ymax);
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
    ReinitializeCube();
    return;
  }

  // If conditions met (0.5 m dist or 45 degrees rotated)
  if (distance_travelled <= 0.0 || angle_travelled <= 0)
  { 
        printf("Laser\n");

    float current_image[x_image_width][y_image_width];
    InitializeImage(current_image);
    // PrintImage(current_image);
    float bounds[4];
    GetBoundingBox(bounds);
    float xmin = bounds[0];
    float xmax = bounds[1];
    float ymin = bounds[2];
    float ymax = bounds[3];
    float temp_image[(int) (xmax - xmin + image_disp)][(int) (ymax - ymin + image_disp)];
    for (unsigned int x = 0; x < sizeof(temp_image) / sizeof(temp_image[0]); x++)
      for (unsigned int y = 0; y < sizeof(temp_image[0]) / sizeof(float); y++)
        temp_image[x][y] = 0;

    for (unsigned int x = 0; x < sizeof(temp_image) / sizeof(temp_image[0]); x++)
    {
      for (unsigned int y = 0; y < sizeof(temp_image[0]) / sizeof(float); y++)
      {
        float x_val = (xmax + image_disp) - x;
        float y_val = (ymax + image_disp) - y;

        for (auto& point : previous_scan)
        {
          float prob = 1.0;
          float x_diff = x - point[0];
          float y_diff = y - point[1];

          float std_dev = 0.02 * pow(pow(x_diff, 2) + pow(y_diff, 2), 0.5);

          // Decoupled evaluation of multivariate gaussian where product is taken along x and y
          prob *= exp(-0.5 * pow((x_val - point[0])/std_dev, 2));
          prob *= exp(-0.5 * pow((y_val - point[1])/std_dev, 2));
          // printf("PROB: %f\n", prob);
          // printf("dx: %f dy: %f PROB: %f\n", x - point[0], y- point[1], prob);
          // Sum of Gaussians 
          temp_image[x][y] += prob;
          if (temp_image[x][y] >= 1.0)
            temp_image[x][y] = 1.0;
        }
      }
    }

    // for (unsigned int x = 0; x < sizeof(temp_image) / sizeof(temp_image[0]); x++)
    // {
    //   for (unsigned int y = 0; y < sizeof(temp_image[0]) / sizeof(float); y++)
    //   {
    //     printf("%f ", temp_image[x][y]);
    //   }
    //   printf("\n");
    // }
    vector<Vector2f> points = GetScanPointCloud(ranges, range_min, range_max, angle_min, angle_max);
    float obsv[x_width][y_width][theta_width];

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
            float final_x = cos(dtheta) * new_x - sin(dtheta) * new_y + dx;
            float final_y = sin(dtheta) * new_x + cos(dtheta) * new_y + dy;

            shifted_points.push_back(Vector2f(final_x, final_y));
            float score = 1.0;
            for (auto& point : shifted_points)
            {
              int point_pixel_x = (xmax - xmin + image_disp) - (point[0] - xmin);
              int point_pixel_y = (ymax - ymin + image_disp) - (point[1] - ymin);

              if (point_pixel_x < 0 || point_pixel_x >= xmax || point_pixel_y < 0 || point_pixel_y >= ymax)
                continue;

              // Not multiplying
              score = temp_image[point_pixel_x][point_pixel_y];
            }
            // printf("SCORE: %f\n", score);
            // Not properly taking max
            obsv[pixel_x][pixel_y][pixel_theta] = std::max(score, score);
          }
        }
      }
    }

    for (int pixel_x = 0; pixel_x < x_width; pixel_x++)
      for (int pixel_y = 0; pixel_y < y_width; pixel_y++)
        for (int pixel_theta = 0; pixel_theta < theta_width; pixel_theta++)
          cube[pixel_x][pixel_y][pixel_theta] *= obsv[pixel_x][pixel_y][pixel_theta];

    // Update all necessary global variables
    previous_scan = points;
    best_pose = MostLikelyPose();
    cumulative_transform.delta_x += best_pose.delta_x;
    cumulative_transform.delta_y += best_pose.delta_y;
    cumulative_transform.delta_theta = best_pose.delta_theta;
    previous_pose = best_pose;
    // PrintCube(x_width * y_width * theta_width);
    ReinitializeCube();
    distance_travelled = distance_travelled_og;
    angle_travelled = angle_travelled_og;
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

  printf("ODOM\n");

  // Keep track of odometry to estimate how far the robot has moved between 
  // poses.

  ReinitializeCube();

  // Not sure this is correct (maybe too simple?)
  Vector2f translation_hat = odom_loc - prev_odom_loc_;
  float angle_hat = odom_angle - prev_odom_angle_;
  printf("TRANSLATION: %f %f, Angle Disp: %f\n", translation_hat[0], translation_hat[1], angle_hat);
  printf("DIST: %f\n", distance_travelled);
  distance_travelled -= abs(translation_hat[0] + translation_hat[1]);
  angle_hat -= abs(angle_hat);

  prev_odom_angle_ = odom_angle;
    prev_odom_loc_ = odom_loc;

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
        // printf("PROB: %f\n", prob);
        // Motion model probability
        cube[pixel_x][pixel_y][pixel_theta] *= prob;
      }
    }
  }
  // PrintCube(5);
}

vector<Vector2f> SLAM::GetMap() {
  vector<Vector2f> map;
  // Reconstruct the map as a single aligned point cloud from all saved poses
  // and their respective scans.

  // For all points in previous scan, apply cumulative transformation to convert points to reference frame of pose 1
  for (auto& point : previous_scan)
  {
    // Transform point to reference frame of pose 1
    float new_x = cos(-cumulative_transform.delta_theta) * point[0] - sin(-cumulative_transform.delta_theta) * point[1] - cumulative_transform.delta_x;
    float new_y = sin(-cumulative_transform.delta_theta) * point[0] + cos(-cumulative_transform.delta_theta) * point[1] - cumulative_transform.delta_y;
    estimated_map.push_back(Vector2f(new_x, new_y));
  }

  return estimated_map;
}

}  // namespace slam
