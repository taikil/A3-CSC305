CXX = g++
CXXFLAGS = -std=c++20

Raytracer: Raytracer.cpp
	$(CXX) $(CXXFLAGS) Raytracer.cpp -o Raytracer

clean:
	rm -f Raytracer
