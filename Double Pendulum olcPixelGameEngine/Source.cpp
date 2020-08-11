/*
 =============================================================================
 Simulating the double pendulum using Runge–Kutta method(RK4)
 =============================================================================
Updated on Aug 11 2020
@author: Ahmed Alkharusi

I used the olcPixelGameEngine to generate the graphics.
The "olcPixelGameEngine.h" is a single header file that enables us to draw graphics.
This is created by javidx9 (OneLoneCoder). 
please download it from his repo
https://github.com/OneLoneCoder/olcPixelGameEngine
*/



#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include "Array.h"


const double pi{ 3.14159265359 };
// These should be moved inside main() and the functions should be modified accordingly 
double m1{ 1.0 }; //mass of the 1st pendulum 
double m2{ 1.0 }; //mass of the 2nd pendulum 
double g{ 10.0 }; //gravity
double r1{ 1.0 }; //length of the 1st pendulum   
double r2{ 1.0 }; //length of the 2nd pendulum   

Array deriv_a1(Array a1_arr, Array a2_arr, Array t);
Array deriv_a2(Array a2_arr, Array a1_arr, Array t);
void rk4(Array(*deriv)(Array, Array, Array), Array& func_i, Array& func_i2, double x_i, double h);

class Pendulum : public olc::PixelGameEngine
{

public:
	Pendulum()
	{
		sAppName = " ";
	}

private:
	double step_size{ 0.002 };//step size for the RK4 method
	Array a1_arr{ pi / 2 , 0 };
	Array a2_arr{ pi / 2 , 0};
	Array temp = a1_arr;
	double t{ 0.0 };
	double pendulum1_x{0};
	double pendulum1_y{ 0 };
	double pendulum2_x{ 0 };
	double pendulum2_y{ 0 };

public:

	
	bool OnUserCreate() override
	{
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		Clear(olc::DARK_BLUE);
		
		int scale = 250;
		temp = a1_arr;	
		rk4(deriv_a1, a1_arr, a2_arr,t, step_size);
		
		rk4(deriv_a2, a2_arr, temp,t, step_size);
		t += step_size;
		pendulum1_x = scale* r1*sin(a1_arr.angle);
		pendulum1_y = scale * r1*cos(a1_arr.angle);
		pendulum2_x = (pendulum1_x + r2*scale * sin(a2_arr.angle));
		pendulum2_y = (pendulum1_y - r2*scale * cos(a2_arr.angle));

		DrawLine(pendulum1_x + ScreenWidth() / 2 , pendulum1_y + ScreenHeight() / 2 , pendulum2_x + ScreenWidth() / 2 , pendulum2_y + ScreenHeight() / 2 , olc::YELLOW);
		DrawLine(pendulum1_x + ScreenWidth() / 2 , pendulum1_y + ScreenHeight() / 2 , ScreenWidth() / 2 , ScreenHeight() / 2 , olc::YELLOW);

		FillCircle(pendulum1_x + ScreenWidth() / 2 , pendulum1_y + ScreenHeight() / 2 , 15, olc::WHITE);
		FillCircle(pendulum2_x + ScreenWidth() / 2, pendulum2_y + ScreenHeight() / 2, 15, olc::WHITE);

		DrawString(ScreenWidth()-700,10,"Menu                         | Key (Hold) \n________________________________\n\nChange initial conditions    | I\n\nChange length of the pendula | L\n\nChange the masses		          | M\n\nChange gravity	   	          | G\n\nAbout                        | A\n\n\n               This will display in a\n\n               separate window\n\n\n", olc::WHITE, 2);
		if (GetKey(olc::Key::I).bHeld) {
			std::cout << " Enter the initial conditions separated by a space :\n angle_pendulum1 angular_speed_pendulum1 angle_pendulum2 angular_speed_pendulum2 \n e.g. 45     0     60     0   " << std::endl;

			std::cin >> a1_arr.angle >> a1_arr.angular_speed >> a2_arr.angle >> a2_arr.angular_speed;
			a1_arr.angle *= pi/180.0;
			a2_arr.angle *= pi/180.0;
		}
		if (GetKey(olc::Key::A).bHeld) {
			std::cout << "\n\n=============================================================================\nSimulating the double pendulum using Runge Kutta method(RK4)\n========================================================================================\nUpdated on Aug 11 2020\n@author: Ahmed Alkharusi\nI used the olcPixelGameEngine to generate the graphics.\nThe olcPixelGameEngine.h is a single header file that enables us to draw graphics.\nThis is created by javidx9(OneLoneCoder).\nplease download it from his repo\nhttps ://github.com/OneLoneCoder/olcPixelGameEngine " << std::endl;

		}
		if (GetKey(olc::Key::G).bHeld) {
			std::cout << " Enter the value of gravity  " << std::endl;
			std::cin >> g ;
		}
		if (GetKey(olc::Key::M).bHeld) {
			std::cout << " Enter the masses of the pendula separated by a space \n e.g. 1 1  " << std::endl;
			std::cin >> m1 >> m2;
		}
		if (GetKey(olc::Key::L).bHeld) {
			std::cout << " Enter the lengths of the pendula separated by a space \n e.g. 1 1  " << std::endl;
			std::cin >> r1 >> r2;
		}
		
		return true;
	}
};


int main() {

	Pendulum pendulum;
	if (pendulum.Construct(1920, 1080, 1, 1))
		pendulum.Start();

	return 0;
}

Array deriv_a1(Array a1_arr, Array a2_arr, Array t) {
	double num = -g * (2 * m1 + m2) * std::sin(a1_arr.angle) - m2 * g * std::sin(a1_arr.angle - 2 * a2_arr.angle) - 2 * m2 * std::sin(a1_arr.angle - a2_arr.angle) * (r2 * std::pow(a2_arr.angular_speed, 2) + r1 * std::pow(a1_arr.angular_speed, 2) * std::cos(a1_arr.angle - a2_arr.angle));
	double den = r1 * (2 * m1 + m2 - m2 * std::cos(2 * a1_arr.angle - 2 * a2_arr.angle));
	// (num / den) is the angular speed
	Array temp{ a1_arr.angular_speed,num / den };
	return  temp;
}
Array deriv_a2(Array a2_arr, Array a1_arr, Array t) {
	double temp = (2 * std::sin(a1_arr.angle - a2_arr.angle));
	double num = temp * (r1 * std::pow(a1_arr.angular_speed, 2) * (m1 + m2) + g * (m1 + m2) * std::cos(a1_arr.angle) + r2 * std::pow(a2_arr.angular_speed, 2) * m2 * std::cos(a1_arr.angle - a2_arr.angle));
	double den = r2 * (2 * m1 + m2 - m2 * std::cos(2 * a1_arr.angle - 2 * a2_arr.angle));
	// (num / den) is the angular speed
	Array res{ a2_arr.angular_speed, num / den };
	return  res;
}

void rk4(Array(*deriv)(Array, Array, Array) ,Array& func_i, Array& func_i2, double x_i, double h) {
	Array k1 = deriv(func_i, func_i2, (k1 * 0) + x_i);
	Array k2 = deriv(func_i + h / 2, func_i2, k1 * (h / 2));
	Array k3 = deriv(func_i + h / 2, func_i2, k2 * (h / 2));
	Array k4 = deriv(func_i + h, func_i2, k3 * h);
	func_i = func_i + ((k1 + (k2 * 2) + (k3 * 2) + k4)) * (h / 6);
	
}
/*
 =============================================================================
 Please check the answers!!!
 =============================================================================
References:

#Implementing the RK4 method in Python
https://youtu.be/mqoqAovXxWA
by Prof. Niels Walet
#The formulas for the angular acceleration
https://www.myphysicslab.com/pendulum/double-pendulum-en.html

#Animating the double pendulum (N.B. the implementation used here is different)
https://matplotlib.org/3.2.1/gallery/animation/double_pendulum_sgskip.html
*/

