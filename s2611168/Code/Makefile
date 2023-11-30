CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3
 
SRCS = main.cc basics.h camera5.h color.h cylinder.h hittable.h interval.h json_read.h lights.h material2.h motion2.h ray.h sphere.h texture.h triangle.h vec3.h

OBJS = $(SRCS:.cc=.o)

TARGET = ray_tracer

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)