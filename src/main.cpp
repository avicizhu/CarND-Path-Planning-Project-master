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
#include <ctime>

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
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
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

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
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
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
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
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
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
enum State{
  START,
  KEEP_LANE,
  PRE_LANE_CHANGE,
  LANE_CHANE
};

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int laneForD(double d) {
  if (d < 4) return 0;
  if (d >= 4.0 & d < 8.0) return 1;
  return 2;
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
  //string map_file_ = "../data/highway_map.csv";
  string map_file_ = "../data/highway_map_bosch1.csv";
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

  const double speed_limit = 49.5;
  const float ref_dist = 20.0;
  int lane = 1;//starts in lane 1
  double ref_vel = 49.5;//speed limit in mph
  time_t prev_lane_change_time = time(0);
  State state = START;

  h.onMessage([&prev_lane_change_time, &ref_dist, &state, &speed_limit, &ref_vel, &map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy, &lane](
      uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
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

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;
          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            int prev_size = previous_path_x.size();

            if (prev_size > 0) {
              car_s = end_path_s;
            }
            //change state here

            bool ok_ahead = true;
            const int target_score = 120;
            vector<bool> ok_lanes = {true, true, true};
            vector<int> car_nums = {0, 0, 0};
            vector<double> lane_speed = {0, 0, 0}, speed_ahead = {0, 0, 0},
              dist_ahead = {200, 200, 200}, lane_score = {0, 0, 0};

            for (int i = 0; i < sensor_fusion.size(); i++) {
              double check_car_d = sensor_fusion[i][6];
              int check_car_lane = laneForD((float)check_car_d);

              double vx = sensor_fusion[i][3];
              double vy = sensor_fusion[i][4];
              double delta_theta = 1.57;
              double check_speed_v_s = sqrt(vx*vx + vy*vy);//vy * cos(delta_theta) + vx * sin(delta_theta);
              double check_speed_v_d = -vy * sin(delta_theta) + vx * cos(delta_theta);
              //cout << "v_d: "<< check_speed_v_d<<endl;
              //double check_speed = sqrt(vx * vx + vy * vy);
              double check_car_s = sensor_fusion[i][5];
              check_car_s += (double) prev_size * 0.02 * check_speed_v_s;
              check_car_d += (double) prev_size * 0.02 * check_speed_v_d;
              if(abs(check_car_s - car_s) < 5){
              //cout << "speed: " << check_car_s<<endl;
              //cout<< "car_d check_car_d: "<<car_d <<" | "<< check_car_d <<" | "<<prev_size * 0.02 * check_speed_v_d<<endl;
            }

              if(check_car_s > car_s){// get info of car in front for lane scoring
                car_nums[check_car_lane]++;
                if(check_car_s - car_s < dist_ahead[check_car_lane]){
                  dist_ahead[check_car_lane] = check_car_s - car_s;
                  speed_ahead[check_car_lane] = check_speed_v_s;
                }
                lane_speed[check_car_lane] += check_speed_v_s;
              }
              // check which lane is safe to drive
              if (check_car_lane == lane && check_car_s > car_s && check_car_s - car_s < ref_dist) {
                ok_lanes[check_car_lane] = false;
                ok_ahead = false;
              }else if (check_car_lane != lane){
                if(car_s > check_car_s && car_s - check_car_s < 0.3 * ref_dist)
                  ok_lanes[check_car_lane] = false;
                if(check_car_s > car_s && check_car_s - car_s < ref_dist)
                  ok_lanes[check_car_lane] = false;
              }
              if(abs(check_car_s - car_s) < 5 && abs(car_d - check_car_d) < 0.5 ){
                cout << "too close";
                ok_lanes[lane] = false;
                ok_ahead = false;
              }
            }
            int prefer_lane = 1;
            for(size_t i = 0 ; i < lane_score.size(); i++){
              if(car_nums[i] == 0) lane_speed[i] = speed_limit;
              else lane_speed[i] /= car_nums[i];
              lane_score[i] = lane_speed[i] + dist_ahead[i];
            }

            double best_score = lane_score[1];
            if(lane_score[2] > best_score){ prefer_lane = 2; best_score = lane_score[2]; }
            if(lane_score[0] > best_score){ prefer_lane = 0; best_score = lane_score[0]; }
            //cout<< "lane score: " << lane_score[0] <<" "<<  lane_score[1]<< " " << lane_score[2]<< endl;
            //cout<< "prefer_lane: "<<prefer_lane<<endl;
            vector<string> state_str = {"Start", "KEEP_LANE", "PRE_LANE_CHANGE", "LANE_CHANE"};
            //cout << "now state: " << state_str[state] << endl;
            int target_lane = lane + sgn(prefer_lane - lane);//make sure change lane one by one
            //choose the car to follow --- the slower one
            double min_dist_ahead = min(dist_ahead[target_lane], dist_ahead[lane]);
            double min_speed_ahead = min(speed_ahead[target_lane], speed_ahead[lane]);

            double Kp = 100;
            switch (state) {
              case START:
                lane = 1;
                ref_vel = .224;
                state = KEEP_LANE;
                break;
              case KEEP_LANE:
                if(ok_ahead && ref_vel < speed_limit)
                  ref_vel += .224;
                else if(!ok_ahead ||
                  (best_score > target_score && prefer_lane != lane && lane_score[lane] < target_score))
                  state = PRE_LANE_CHANGE;
                break;
              case PRE_LANE_CHANGE:
                // before changing lane follow the traffic first
                //cout << ((ref_vel - speed_ahead)/Kp) <<endl;
                if( min_dist_ahead > 0.7 * ref_dist && ref_vel < speed_limit)
                  ref_vel += min(.224, (ref_vel - min_speed_ahead)/Kp);
                else
                  ref_vel -= min(.224, abs(ref_vel - min_speed_ahead)/Kp);
                if(time(0) - prev_lane_change_time < 2){
                  //cout << "don't change lane too frequently"<<endl;
                  break;
                }

                if(lane != target_lane && ok_lanes[target_lane]){
                  lane = target_lane;
                  state = LANE_CHANE;
                }
                break;
              case LANE_CHANE:
              if( min_dist_ahead > 0.5 * ref_dist && ref_vel < speed_limit)
                ref_vel += min(.224, (ref_vel - min_speed_ahead)/Kp);
              else
                ref_vel -= min(.224, abs(ref_vel - min_speed_ahead)/Kp);
              if(laneForD(car_d) == lane){
                prev_lane_change_time = time(0);
                state = KEEP_LANE;
              }
                break;
            }

            vector<double> ptsx;
            vector<double> ptsy;

            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);

            if (prev_size < 2) {
              double prev_car_x = car_x - cos(car_yaw);
              double prev_car_y = car_y - sin(car_yaw);

              ptsx.push_back(prev_car_x);
              ptsx.push_back(car_x);

              ptsy.push_back(prev_car_y);
              ptsy.push_back(car_y);
            } else {
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

            vector<double> next_wp0 = getXY(car_s + 30, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp1 = getXY(car_s + 60, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp2 = getXY(car_s + 90, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

            ptsx.push_back(next_wp0[0]);
            ptsx.push_back(next_wp1[0]);
            ptsx.push_back(next_wp2[0]);

            ptsy.push_back(next_wp0[1]);
            ptsy.push_back(next_wp1[1]);
            ptsy.push_back(next_wp2[1]);

            //shift car reference angle to 0
            for (int i = 0; i < ptsx.size(); i++) {
              double shift_x = ptsx[i] - ref_x;
              double shift_y = ptsy[i] - ref_y;

              ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
              ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));

            }

            tk::spline s;
            s.set_points(ptsx, ptsy);

            for (int i = 0; i < previous_path_x.size(); i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            double target_x = 30;
            double target_y = s(target_x);
            double target_dist = sqrt((target_x) * (target_x) + (target_y) * (target_y));
            double x_add_on = 0;

            for (int i = 1; i <= 50 - previous_path_x.size(); i++) {

              double N = (target_dist / (0.02 * ref_vel / 2.24));
              double x_point = x_add_on + (target_x) / N;
              double y_point = s(x_point);

              x_add_on = x_point;

              double x_ref = x_point;
              double y_ref = y_point;
              //rotate back to normal
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
