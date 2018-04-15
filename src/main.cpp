#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
  int laneWidth = 4; // Width of lanes in meters
  int lanes = 3; // Total number of lanes
  int lane = 1; // 0 is far left, 1 is next to the right, etc...
  
  double ref_vel = 0; //mph, start at 0, we will accelerate later on. :)
  double speed_limit = 49.5; //Speed limit, do not exceed.
  double target_speed = 49.5; //Target speed 
  double last_lane_change = 0.0; //Distance since last lane change

  h.onMessage([&ref_vel,&lane,&lanes,&speed_limit,&target_speed,&laneWidth,&last_lane_change,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];
          	//Cost for Lane Changes and Maintain Lane, Maintain initiated at 0 to have priority
          	double leftLCCost = 99; 
          	double rightLCCost = 99;
          	double maintainLCost = 0;
          	//Count how many cars are in each lane, used for Lane Costs later.
          	int leftCarCount = 0;
          	int rightCarCount = 0;

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;
          	          	         	
          	int prev_size = previous_path_x.size();
          	
          	if(prev_size > 0)
          		car_s = end_path_s;
          	//too_close flag initiates looking into lane changes
          	bool too_close = false;
          	//Lane blocked flags to ensure the car doesn't try to change lanes into
          	//another car or off the road.
          	bool left_blocked = false;
          	bool right_blocked = false;
          	//distance to car in front to initiate lane change checks, also used for 
          	//checking to make sure we are looking at the closest car.
          	double center_closest = 60;
          	//Take speed readings to ensure we don't hit the car in front.
          	double closest_speed = speed_limit;
          	
          	//Cycle through cars and check for ones that are in our lane in front of us
          	for( int i = 0; i < sensor_fusion.size(); i++ )
          	{
          		float d = sensor_fusion[i][6];
          		if( d < laneWidth * lane + laneWidth && d > laneWidth * lane )
          		{
          			float vx = sensor_fusion[i][3];
          			float vy = sensor_fusion[i][4];
          			double check_speed = sqrt(pow(vx, 2) + pow(vy, 2)); //check speed, m/s
          			double check_car_s = sensor_fusion[i][5]; //check cars position (meter)
          			
          			check_car_s += ( (double)prev_size * .02 * check_speed);
          			double difference = check_car_s - car_s; //check distance between our cars
          			//if car is within 60m, start looking for a lane change. Calculate speed
          			//in MPH to start slowing down
          			if (check_car_s > car_s && difference < center_closest ) 
          			{
          				too_close = true;
          				center_closest = difference;
          				maintainLCost = speed_limit - ( check_speed * 2.23694 );
          				closest_speed = check_speed * 2.23694;
          			}
          		}
          		//If no car within 40 Meters, go the speed limit, if car within 40, slow to 
          		//the speed of the car, if within 20 Meters, decrease speed to give cushion space
          		if ( center_closest < 20 )
							{
								target_speed = closest_speed * .75;
							} else if ( center_closest < 40 && closest_speed <= speed_limit )
							{
								target_speed = closest_speed;
							} else
							{
								target_speed = speed_limit;
							}
          	}
						//variables to store target speed if we perform a lane change into a lane with 
						//a car going slower than speed limit.
						double left_lane_change_speed = speed_limit;
						double right_lane_change_speed = speed_limit;
						//check for cars in the left and right lanes, starting at 120 meters, also used
						//to verify the closest car in front of our vehicle.
						double left_closest = 120.0;
						double right_closest = 120.0;
						
						//check to see if we should even bother checking the left or right lane 
						//(are we already in the far lanes?)
						if ( lane == 0 )
						{
							left_blocked = true;
						} else
						{
							left_blocked = false;
						}
						if ( lane == lanes - 1 )
						{
							right_blocked = true;
						} else
						{
							right_blocked = false;
						}
							
						//Cycle through cars, check for viable lane change opportunities
						for( int i = 0; i < sensor_fusion.size(); i++ )
						{
							float d = sensor_fusion[i][6];
							//check for all cars in the left lane.
							if( ( lane != 0 ) && ( d < laneWidth * ( lane - 1 ) + laneWidth ) && ( d > laneWidth * ( lane - 1 ) ) && !left_blocked )
							{
								float vx = sensor_fusion[i][3];
								float vy = sensor_fusion[i][4];
								double check_speed = sqrt(pow(vx, 2) + pow(vy, 2));
								double check_car_s = sensor_fusion[i][5];
								check_car_s += ( (double)prev_size * .02 * check_speed);
								double difference = check_car_s - car_s;
								//Check to make sure no car is within 20 meters in front of us or 8 meters behind us.
								//Also check to make sure no car is speeding up to our position from 30 meters behind.
								if( ( difference < 20 && difference > -8 ) || ( difference < -8 && difference > -30 && check_speed * 2.23694 > ref_vel ))
								{	
									left_blocked = true;
									leftCarCount += 1;
								}
								//Check to see if the car is ahead of us, and if its the closest, and calculate its cost.
								if ( check_car_s > car_s && difference < left_closest)
								{
									leftLCCost = 10.0 + speed_limit - ( check_speed * 2.23694 );
									leftLCCost -= .2 * difference;
									left_closest = difference;
									left_lane_change_speed = check_speed * 2.23694;
									leftCarCount += 1;
								} //Check all of the above things, but now for the right lane.
							} else if( ( lane != lanes - 1 ) && ( d < laneWidth * ( lane + 1 ) + laneWidth ) && ( d > laneWidth * ( lane + 1 ) ) && !right_blocked )
							{
								float vx = sensor_fusion[i][3];
								float vy = sensor_fusion[i][4];
								double check_speed = sqrt(pow(vx, 2) + pow(vy, 2));
								double check_car_s = sensor_fusion[i][5];
								check_car_s += ( (double)prev_size * .02 * check_speed);
								double difference = check_car_s - car_s;
								if( ( difference < 20 && difference > -8 ) || ( difference < -8 && difference > -30 && check_speed * 2.23694 > ref_vel ) || lane == ( lanes - 1 ) )
								{	
									right_blocked = true;
									rightCarCount += 1;
								} 
								if ( check_car_s > car_s && difference < right_closest)
								{
									rightLCCost = 10.0 + speed_limit - ( check_speed  * 2.23694 );
									rightLCCost -= .2 * difference;
									right_closest = difference;
									right_lane_change_speed = check_speed * 2.23694;
									rightCarCount += 1;
								}
							}
						}
						//Set costs very low if there are no cars ahead of us in the respective lane
						if (leftCarCount == 0 )
						{
							leftLCCost = 5;
						}
						if( rightCarCount == 0 )
						{
							rightLCCost = 5;
						}
						
						cout <<"Current Lane: " << lane << endl << "Last Lane Change: " << car_s - last_lane_change << " Meters ago" << endl << "MAINTAIN COST: " << int(maintainLCost) << endl << "LEFT COST:     " << int(leftLCCost) << " " << left_blocked << endl << "RIGHT COST:    " << int(rightLCCost) << " " << right_blocked << endl;
						//Compare costs, make a decision whether or not to perform a lane change. If lane
						//change occurs, match speed to car ahead if its lower than ours. If a lane change
						//has occurred in the last 100 meters, do not perform one. 
						if( too_close )
						{
							if( leftLCCost < maintainLCost && ( leftLCCost <= rightLCCost || right_blocked )  && !left_blocked && abs( car_s - last_lane_change ) > 100.0 )
							{
								lane -= 1;
								last_lane_change = car_s;
								target_speed = left_lane_change_speed;
							} else if( rightLCCost < maintainLCost && !right_blocked && abs( car_s - last_lane_change ) > 100.0 )
							{
								lane += 1;
								last_lane_change = car_s;
								target_speed = right_lane_change_speed;
							}
						}
						//Use our target_speed variable to ramp our ref_vel up or down at a reasonable rate
          	if( target_speed < ref_vel )
							ref_vel -= .224;
						else if( target_speed > ref_vel )
							ref_vel += .224;
          	
          	//Pulled most of this from the video help section on the project page.
          	vector<double> ptsx;
          	vector<double> ptsy;
          	
          	double ref_x = car_x;
          	double ref_y = car_y;
          	double ref_yaw = deg2rad(car_yaw);
          	cout << "ref_x = " << ref_x << endl << "ref_y = " << ref_y << endl << "ref_yaw = " << ref_yaw << endl;
          	
          	if(prev_size < 2)
          	{
          		double prev_car_x = car_x - cos(car_yaw);
          		double prev_car_y = car_y - sin(car_yaw);
          		
          		ptsx.push_back(prev_car_x);
          		ptsx.push_back(car_x);
          		
          		ptsy.push_back(prev_car_y);
          		ptsy.push_back(prev_car_y);
          	}
          	
          	else
          	{
          		ref_x = previous_path_x[prev_size - 1];
          		ref_y = previous_path_y[prev_size - 1];
          		
          		double ref_x_prev = previous_path_x[prev_size - 2];
          		double ref_y_prev = previous_path_y[prev_size - 2];
          		ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);
          		
          		ptsx.push_back(ref_x_prev);
          		ptsx.push_back(ref_x);
          		
          		ptsy.push_back(ref_y_prev);
          		ptsy.push_back(ref_y);
          	}
          	
          	vector<double> next_wp0 = getXY(car_s + 30, (laneWidth/2 + laneWidth * lane ), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          	vector<double> next_wp1 = getXY(car_s + 60, (laneWidth/2 + laneWidth * lane ), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          	vector<double> next_wp2 = getXY(car_s + 90, (laneWidth/2 + laneWidth * lane ), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          	
          	ptsx.push_back(next_wp0[0]);
          	ptsx.push_back(next_wp1[0]);
          	ptsx.push_back(next_wp2[0]);
          	
          	ptsy.push_back(next_wp0[1]);
          	ptsy.push_back(next_wp1[1]);
          	ptsy.push_back(next_wp2[1]);
          	
          	for( int i = 0; i < ptsx.size(); i++ )
          	{
          		double shift_x = ptsx[i] - ref_x;
          		double shift_y = ptsy[i] - ref_y;
          		
          		ptsx[i] = (shift_x * cos( 0 - ref_yaw) - shift_y * sin( 0 - ref_yaw));
          		ptsy[i] = (shift_x * sin( 0 - ref_yaw) + shift_y * cos( 0 - ref_yaw));
          	}
          	
          	tk::spline s;
          	cout << "(" << ptsx[1] << "," << ptsy[1] << ")" << endl;
          	s.set_points(ptsx, ptsy);
          	
          	vector<double> next_x_vals;
          	vector<double> next_y_vals;
          	
          	for( int i = 0; i < previous_path_x.size(); i++)
          	{
          		next_x_vals.push_back(previous_path_x[i]);
          		next_y_vals.push_back(previous_path_y[i]);
          	}
          	
          	double target_x = 30.0;
          	double target_y = s(target_x);
          	double target_dist = sqrt((target_x)*(target_x)+(target_y)*(target_y));
          	
          	double x_add_on = 0;
          	
          	for( int i = 1; i <= 50 - previous_path_x.size(); i++ )
          	{
          		double N = (target_dist / ( .02 * ref_vel / 2.24 ) );
          		double x_point = x_add_on + (target_x) / N;
          		double y_point = s(x_point);
          		
          		x_add_on = x_point;
          		
          		double x_ref = x_point;
          		double y_ref = y_point;
          		
          		x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
          		y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));
          		
          		x_point += ref_x;
          		y_point += ref_y;
          		
          		next_x_vals.push_back(x_point);
							next_y_vals.push_back(y_point);
						}
						msgJson["next_x"] = next_x_vals;
						msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}