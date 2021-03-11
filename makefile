all: proj4

proj4: proj4.cpp
	g++ proj4.cpp -o proj4 -lGL -lglut -lGLU -lGLEW -lX11 -lm

clean:
	rm -f proj4
