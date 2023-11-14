CXX = g++
CXXFLAGS = -std=c++11

Raytracer: Raytracer.cpp
	$(CXX) $(CXXFLAGS) Raytracer.cpp -o Raytracer

clean:
	rm -f Raytracer
