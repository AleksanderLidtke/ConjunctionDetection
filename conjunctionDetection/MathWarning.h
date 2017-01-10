/*
A custom error class that will include a time stamp for the error. To be thrown when 
Mathamatics doesn't work.

@version 1.0.0
@since 08 Aug 2014 11:43:00
@author Aleksander Lidtke
@email al11g09@soton.ac.uk, alek_l@onet.eu

CHANGELOG:
*/
#pragma once

#include <time.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

class MathWarning : public std::runtime_error{
	public:
		MathWarning(void):runtime_error("Math Warning"){}
		MathWarning(std::string msg, std::string file, int line):runtime_error(msg.c_str()) {
			char buff[25]; time_t now = time (0);
			#ifdef _MSC_VER 
				struct tm timeinfo;
				localtime_s(&timeinfo, &now);
				strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", &timeinfo);
			#else
				strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", localtime (&now));
			#endif

			//std::cerr <<buff<<" ---> MATH WARNING IN: " <<file<< " AT LINE " << line << ":\n\n";
		}
};
