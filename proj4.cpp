#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<string.h>
#include<cstring>
#include<fstream>
#include<sstream>
#include<math.h>

#define STB_IMAGE_IMPLEMENTATION
#include"/usr/local/include/stb_image.h"

#include<GL/glew.h>

#include<glm/glm.hpp>
#include<glm/gtc/matrix_transform.hpp>
#include<glm/gtc/type_ptr.hpp>

#include<GL/gl.h>
#include<GL/glu.h>
#include<GL/glx.h>

#include<X11/Xlib.h>
#include<X11/keysym.h>

#define rnd() (double)rand() / (double)RAND_MAX
#define pi 3.141592
#define massnum 2000

GLfloat WHITE[] = {1,1,1};

typedef double vec[3];
#define VecCross(a,b,c) \
(c)[0]=(a)[1]*(b)[2]-(a)[2]*(b)[1]; \
(c)[1]=(a)[2]*(b)[0]-(a)[0]*(b)[2]; \
(c)[2]=(a)[0]*(b)[1]-(a)[1]*(b)[0]

void vecCrossProduct(vec v0, vec v1, vec dest)
{
        dest[0] = v0[1]*v1[2] - v1[1]*v0[2];
        dest[1] = v0[2]*v1[0] - v1[2]*v0[0];
        dest[2] = v0[0]*v1[1] - v1[0]*v0[1];
}

double vecDotProduct(vec v0, vec v1)
{
        return v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
}

void vecZero(vec v)
{
        v[0] = v[1] = v[2] = 0.0;
}

void vecMake(double a, double b, double c, vec v)
{
        v[0] = a;
        v[1] = b;
        v[2] = c;
}

void vecCopy(vec source, vec dest)
{
        dest[0] = source[0];
        dest[1] = source[1];
        dest[2] = source[2];
}

double vecLength(vec v)
{
        return sqrt(vecDotProduct(v, v));
}

void vecNormalize(vec v)
{
        double len = vecLength(v);
        if (len == 0.0) {
                vecMake(0,0,1,v);
                return;
        }
        len = 1.0 / len;
        v[0] *= len;
        v[1] *= len;
        v[2] *= len;
}

void vecSub(vec v0, vec v1, vec dest)
{
        dest[0] = v0[0] - v1[0];
        dest[1] = v0[1] - v1[1];
        dest[2] = v0[2] - v1[2];
}

void getTriangleNormal(vec v1, vec v2, vec v3, vec norm) {
	vec v4, v5;
	vecSub(v2, v1, v4);
	vecSub(v3, v1, v5);
	vecSub(v3, v1, v5);

	vecCrossProduct(v4, v5, norm);
	vecNormalize(norm);
}
///////////////////////////////////SKYBOX/////////////////////////////////////////
std::vector<std::string> faces = {
        "./skybox/right.jpg",
        "./skybox/left.jpg",
        "./skybox/top.jpg",
        "./skybox/bottom.jpg",
        "./skybox/front.jpg",
        "./skybox/back.jpg"
};

std::vector<std::string> faces2 = {
	"./ely_peaks/right.tga",
	"./ely_peaks/left.tga",
	"./ely_peaks/top.tga",
	"./ely_peaks/bottom.tga",
	"./ely_peaks/front.tga",
	"./ely_peaks/back.tga"
};

std::vector<std::string> faces3 = {
	"./hw_crater/right.tga",
	"./hw_crater/left.tga",
	"./hw_crater/top.tga",
	"./hw_crater/bottom.tga",
	"./hw_crater/front.tga",
	"./hw_crater/back.tga"
};

float skyboxVerts[] = {
    -9.0f,  9.0f, -9.0f,
    -9.0f, -9.0f, -9.0f,
     9.0f, -9.0f, -9.0f,
     9.0f, -9.0f, -9.0f,
     9.0f,  9.0f, -9.0f,
    -9.0f,  9.0f, -9.0f,

    -9.0f, -9.0f,  9.0f,
    -9.0f, -9.0f, -9.0f,
    -9.0f,  9.0f, -9.0f,
    -9.0f,  9.0f, -9.0f,
    -9.0f,  9.0f,  9.0f,
    -9.0f, -9.0f,  9.0f,

     9.0f, -9.0f, -9.0f,
     9.0f, -9.0f,  9.0f,
     9.0f,  9.0f,  9.0f,
     9.0f,  9.0f,  9.0f,
     9.0f,  9.0f, -9.0f,
     9.0f, -9.0f, -9.0f,

    -9.0f, -9.0f,  9.0f,
    -9.0f,  9.0f,  9.0f,
     9.0f,  9.0f,  9.0f,
     9.0f,  9.0f,  9.0f,
     9.0f, -9.0f,  9.0f,
    -9.0f, -9.0f,  9.0f,

    -9.0f,  9.0f, -9.0f,
     9.0f,  9.0f, -9.0f,
     9.0f,  9.0f,  9.0f,
     9.0f,  9.0f,  9.0f,
    -9.0f,  9.0f,  9.0f,
    -9.0f,  9.0f, -9.0f,

    -9.0f, -9.0f, -9.0f,
    -9.0f, -9.0f,  9.0f,
     9.0f, -9.0f, -9.0f,
     9.0f, -9.0f, -9.0f,
    -9.0f, -9.0f,  9.0f,
     9.0f, -9.0f,  9.0f
};

unsigned int loadCubemap(std::vector<std::string> faces) {
        unsigned int tID;
	glGenTextures(1, &tID);
        glBindTexture(GL_TEXTURE_CUBE_MAP, tID);

        int w, h, nrChannels;
        for(unsigned int i=0; i<faces.size(); i++) {
                unsigned char* data = stbi_load(faces[i].c_str(), &w, &h, &nrChannels, 0);
                if(data && i == 0) {
                        printf("CUBEMAP TEXTURE %s LOADED\n", faces[i].c_str()); fflush(stdout);
                        glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X, 0, GL_RGB, w, h, 0, GL_RGB,
                                        GL_UNSIGNED_BYTE, data);
                        stbi_image_free(data);
                }
                else if(data && i == 1) {
                        printf("CUBEMAP TEXTURE %s LOADED\n", faces[i].c_str()); fflush(stdout);
                        glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_X, 0, GL_RGB, w, h, 0, GL_RGB,
                                        GL_UNSIGNED_BYTE, data);
                        stbi_image_free(data);
                }
                else if(data && i == 2) {
                        printf("CUBEMAP TEXTURE %s LOADED\n", faces[i].c_str()); fflush(stdout);
                        glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_Y, 0, GL_RGB, w, h, 0, GL_RGB,
                                        GL_UNSIGNED_BYTE, data);
                        stbi_image_free(data);
                }
                else if(data && i == 3) {
                        printf("CUBEMAP TEXTURE %s LOADED\n", faces[i].c_str()); fflush(stdout);
                        glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, 0, GL_RGB, w, h, 0, GL_RGB,
                                        GL_UNSIGNED_BYTE, data);
                        stbi_image_free(data);
                }
                else if(data && i == 4) {
                        printf("CUBEMAP TEXTURE %s LOADED\n", faces[i].c_str()); fflush(stdout);
                        glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_Z, 0, GL_RGB, w, h, 0, GL_RGB,
                                        GL_UNSIGNED_BYTE, data);
                        stbi_image_free(data);
                }
                else if(data && i == 5) {
                        printf("CUBEMAP TEXTURE %s LOADED\n", faces[i].c_str()); fflush(stdout);
                        glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, 0, GL_RGB, w, h, 0, GL_RGB,
                                        GL_UNSIGNED_BYTE, data);
                        stbi_image_free(data);
                }
                else {
                        printf("ERROR LOADING CUBE MAP: %s\n", faces[i].c_str()); fflush(stdout);
                        stbi_image_free(data);
                }
        }

        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

        return(tID);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Mass {
	double mass, oomass;
	vec pos;
	vec vel;
	vec force;
	GLfloat color[3];
};

struct Spring {
	int mass[2];
	double length;
	double stiffness;
};

double getLength(vec v1, vec v2) {
	double xdiff = v1[0] - v2[0];
        double ydiff = v1[1] - v2[1];
        double zdiff = v1[2] - v2[2];
        return(sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff));
}

class Bubble {
	public:
		const int latitude = 50;
		const int longitude = 25;
		GLfloat verts[26][50][3];

		Mass mass[massnum];
		Spring spring[massnum*100];
		int nsprings = 0;
		int nmasses = 0;
		
		double stiffnessLevel = -0.000004 - rnd() * 0.000006;
		int pull = 0;

		vec center;

		Bubble() {
			double offsetx =0;// (double)(rand() % xbounds) - (double)(rand() % xbounds);
			double offsetz =-5.0;// (double)(rand() % zbounds) - (double)(rand() % zbounds);
			double circle[latitude][2];
			double angle = 0.0, inc = (pi * 2.0) / (float)latitude;
			for(int i=0; i<latitude; i++) {
				circle[i][0] = cos(angle);
				circle[i][1] = sin(angle);
				angle -= inc;
			}
			for(int i=0; i<=longitude; i++) {
				for(int j=0; j<latitude; j++) {
					verts[i][j][0] = circle[j][0] * circle[i][1] + offsetx;
					verts[i][j][2] = circle[j][1] * circle[i][1] + offsetz;
					verts[i][j][1] = circle[i][0];

				}
			}
			 nmasses = 0;
			 vec tmp = {0,0,0};
			for(int i=0; i<longitude; i++) {
				for(int j=0; j<latitude; j++) {
					mass[nmasses].mass = 1.0;
					mass[nmasses].oomass = 1.0 / mass[nmasses].mass;
					mass[nmasses].pos[0] = verts[i][j][0];
					tmp[0] += verts[i][j][0];
					mass[nmasses].pos[1] = verts[i][j][1];
					tmp[1] += verts[i][j][1];
					mass[nmasses].pos[2] = verts[i][j][2];
					tmp[2] += verts[i][j][2];
					mass[nmasses].vel[0] =
					mass[nmasses].vel[1] =
					mass[nmasses].vel[2] =
					mass[nmasses].force[0] =
					mass[nmasses].force[1] =
					mass[nmasses].force[2] = 0.0;
					mass[nmasses].color[0] = 1.0;
					mass[nmasses].color[1] = 1.0;
					mass[nmasses].color[2] = 1.0;
					nmasses++;
				}
			}
			center[0] = tmp[0] / (double)nmasses;
			center[1] = tmp[1] / (double)nmasses;
			center[2] = tmp[2] / (double)nmasses;
			printf("center: %f, %f, %f\n",center[0],
					center[1],center[2]);
			printf("nmasses: %i\n", nmasses);
			nsprings = 0;
			//connects 1-1 horizontally
			for(int i=0; i<longitude; i++) {
				for(int j=0; j<latitude; j++) {
					int k = i*latitude + (j+1) % latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every other horizontally
			for(int i=0; i<longitude; i++) {
				for(int j=0; j<latitude-1; j++) {
					int k = i*latitude + (j+2) % latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every third horizontally
			for(int i=0; i<longitude; i++) {
				for(int j=0; j<latitude-2; j++) {
					int k = i*latitude + (j+3) % latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
		 //connects every fourth horizontally
			for(int i=0; i<longitude; i++) {
				for(int j=0; j<latitude-3; j++) {
					int k = i*latitude + (j+4) % latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every fifth horizontally
			for(int i=0; i<longitude; i++) {
				for(int j=0; j<latitude-4; j++) {
					int k = i*latitude + (j+5) % latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every sixth horizontally
			for(int i=0; i<longitude; i++) {
				for(int j=0; j<latitude-5; j++) {
					int k = i*latitude + (j+6) % latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every tenth  
			for(int i=0; i<longitude; i++) {
				for(int j=0; j<latitude-10; j++) {
					int k = i*latitude + (j+10) % latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects 1-1 vertically
			for(int i=0; i<longitude-1; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+1)*latitude + j;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
		  //connects every other vertically
			for(int i=0; i<longitude-2; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+2)*latitude + j;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every third vertically
			for(int i=0; i<longitude-3; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+3)*latitude + j;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every fourth vertically
			for(int i=0; i<longitude-4; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+4)*latitude + j;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every fifth vertically
			for(int i=0; i<longitude-5; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+5)*latitude + j;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every sixth vertically
			for(int i=0; i<longitude-6; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+6)*latitude + j;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
		//connects every tenth vertically
			for(int i=0; i<longitude-10; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+10)*latitude + j;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every one diagonally left-right
			for(int i=0; i<longitude-1; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+1)*latitude + (j+1)%latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every other diagonally left-right
			for(int i=0; i<longitude-2; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+2)*latitude + (j+1)%latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every third diagonally left-right
			for(int i=0; i<longitude-3; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+3)*latitude + (j+1)%latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every fourth diagonally left-right
			for(int i=0; i<longitude-4; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+4)*latitude + (j+1)%latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
		//connects every fifth diagonally left-right
			for(int i=0; i<longitude-5; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+5)*latitude + (j+1)%latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every tenth diagonally left-right
			for(int i=0; i<longitude-10; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i+10)*latitude + (j+1)%latitude;
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every one diagonally right-left
			for(int i=0; i<longitude-1; i++) {
				for(int j=latitude-1; j>0; j--) {
					int k = (i+1)*latitude + (j-1);
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every other diagonally right-left
			for(int i=0; i<longitude-1; i++) {
				for(int j=latitude-1; j>2; j--) {
					int k = (i+1)*latitude + (j-2);
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every third diagonally right-left
			for(int i=0; i<longitude-1; i++) {
				for(int j=latitude-1; j>3; j--) {
					int k = (i+1)*latitude + (j-3);
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
		//connects every fourth diagonally right-left
			for(int i=0; i<longitude-1; i++) {
				for(int j=latitude-1; j>4; j--) {
					int k = (i+1)*latitude + (j-4);
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every fifth diagonally right-left
			for(int i=0; i<longitude-1; i++) {
				for(int j=latitude-1; j>5; j--) {
					int k = (i+1)*latitude + (j-5);
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects every tenth diagonally right-left
			for(int i=0; i<longitude-1; i++) {
				for(int j=latitude-1; j>10; j--) {
					int k = (i+1)*latitude + (j-10);
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			//connects top & bottom points
			for(int j=0; j<latitude; j++) {
				int k = (longitude-1)*latitude + j;
				spring[nsprings].mass[0] = j;
				spring[nsprings].mass[1] = k;
				spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
				spring[nsprings].stiffness = stiffnessLevel;
				nsprings++;
			}
			//connects points directly across
			for(int i=0; i<longitude; i++) {
				for(int j=0; j<latitude; j++) {
					int k = (i)*latitude + (j+(latitude/2));
					spring[nsprings].mass[0] = i*latitude + j;
					spring[nsprings].mass[1] = k;
					spring[nsprings].length = getLength(mass[1].pos, mass[0].pos);
					spring[nsprings].stiffness = stiffnessLevel;
					nsprings++;
				}
			}
			printf("nsprings: %i\n", nsprings);
		}
		void maintainSprings() {
			int i, m0, m1;
			double dist, oodist, factor;
			double springVec[3];
			double springForce[3];
			const double penalty = 0.5f;

			for(i=0; i<nmasses; i++) {
				mass[i].vel[0] += mass[i].force[0] * mass[i].oomass;
				mass[i].vel[1] += mass[i].force[1] * mass[i].oomass;
				mass[i].vel[2] += mass[i].force[2] * mass[i].oomass;

				mass[i].pos[0] += mass[i].vel[0];
				mass[i].pos[1] += mass[i].vel[1];
				mass[i].pos[2] += mass[i].vel[2];

				mass[i].force[0] =
				mass[i].force[1] =
				mass[i].force[2] = 0.0;
				if(!pull) {}
					//X Bounds
/*					if(mass[i].pos[0] < -xbounds) {
						mass[i].force[0] = penalty * -mass[i].pos[0];
					}
					if(mass[i].pos[0] > xbounds) {
						mass[i].force[0] = penalty * (xbounds - mass[i].pos[0]);
					}
				}
				//Y Bounds
				if(mass[i].pos[1] < -ybounds) {
					mass[i].force[1] = penalty * (-ybounds - mass[i].pos[1]);
				}
				if(mass[i].pos[1] > ybounds) {
					mass[i].force[1] = penalty * (ybounds - mass[i].pos[1]);
				}
				//Z Bounds
				if(mass[i].pos[2] < -zbounds) {
					mass[i].force[2] = penalty * (-zbounds - mass[i].pos[2]);
				}
				if(mass[i].pos[2] > zbounds) {
					mass[i].force[2] = penalty * (zbounds - mass[i].pos[2]);
				}

				//Velocity Constraint
				if(mass[i].vel[0] > 5.0)
					mass[i].vel[0] = 5.0;
				if(mass[i].vel[1] > 5.0)
					mass[i].vel[1] = 5.0;
				if(mass[i].vel[2] > 5.0)
					mass[i].vel[2] = 5.0;
*/
			}
			//Resolve all springs...
			for (i=0; i<nsprings; i++) {
				m0 = spring[i].mass[0];
				m1 = spring[i].mass[1];
				//forces are applied here
				//vector between the two masses
				springVec[0] = mass[m0].pos[0] - mass[m1].pos[0];
				springVec[1] = mass[m0].pos[1] - mass[m1].pos[1];
				springVec[2] = mass[m0].pos[2] - mass[m1].pos[2];
				//distance between the two masses
				dist = sqrt(springVec[0]*springVec[0] + springVec[1]*springVec[1]
						+ springVec[2]*springVec[2]);
				if (dist == 0.0) dist = 0.1;
				oodist = 1.0 / dist;
				springVec[0] *= oodist;
				springVec[1] *= oodist;
				springVec[2] *= oodist;
				//the spring force is added to the mass force
				factor = -(dist - spring[i].length) * spring[i].stiffness;
				springForce[0] = springVec[0] * factor;
				springForce[1] = springVec[1] * factor;
				springForce[2] = springVec[2] * factor;
				//apply force and negative force to each end of the spring...
				mass[m0].force[0] += springForce[0];
				mass[m0].force[1] += springForce[1];
				mass[m0].force[2] += springForce[2];
				mass[m1].force[0] -= springForce[0];
				mass[m1].force[1] -= springForce[1];
				mass[m1].force[2] -= springForce[2];
				//damping
				springForce[0] = (mass[m1].vel[0] - mass[m0].vel[0]) * 0.002;
				springForce[1] = (mass[m1].vel[1] - mass[m0].vel[1]) * 0.002;
				springForce[2] = (mass[m1].vel[2] - mass[m0].vel[2]) * 0.002;
				mass[m0].force[0] += springForce[0];
				mass[m0].force[1] += springForce[1];
				mass[m0].force[2] += springForce[2];
				mass[m1].force[0] -= springForce[0];
				mass[m1].force[1] -= springForce[1];
				mass[m1].force[2] -= springForce[2];
			}
		}
		
		void updateCenter() {
			vec tmp = {0,0,0};
			for(int i=0; i<nmasses; i++) {
				tmp[0] += mass[i].pos[0];
				tmp[1] += mass[i].pos[1];
				tmp[2] += mass[i].pos[2];
			}
			center[0] = tmp[0] / (double)nmasses;
			center[1] = tmp[1] / (double)nmasses;
			center[2] = tmp[2] / (double)nmasses;
		}

		void stray() {
			static float dir = (rnd() * 0.0001) - (rnd() * 0.0001);
			for(int i=0; i<nmasses; i++)
				mass[i].pos[0] += dir;
		}	
} b;

float vertices[massnum*3];
unsigned int ind[massnum];

void initVertices() {
	for(int i=0; i<b.nmasses; i++)
		ind[i] = i;
	for(int i=0; i<b.nmasses; i+=3) {
		vertices[i] = b.mass[i].pos[0];
		vertices[i+1] = b.mass[i].pos[1];
		vertices[i+2] = b.mass[i].pos[2];
	}
}

class global {
	public:
	int xres = 1440;
	int yres = 720;
	float camx, camy, camz;
	float xangle, yangle;

	GLfloat lightAmbient[4];
        GLfloat lightDiffuse[4];
        GLfloat lightSpecular[4];
        GLfloat lightPosition[4];

	global() {
		srand(time(NULL));
		initVertices();
		camx = camy = xangle = yangle = 0.0;
		camz = -2.0;
		GLfloat la[]  = {  0.0f, 0.0f, 0.0f, 1.0f };
		GLfloat ld[]  = {  1.0f, 1.0f, 1.0f, 1.0f };
		GLfloat ls[] = {  0.5f, 0.5f, 0.5f, 1.0f };
		GLfloat lp[] = { 0.0f, 10.0f, 0.0f, 1.0f };
		lp[0] = 0.0;
		lp[1] = 1.0;//rnd() * 100.0 + 20.0;
		lp[2] = 0.0;//rnd() * 300.0 - 150.0;
		memcpy(lightAmbient, la, sizeof(GLfloat)*4);
		memcpy(lightDiffuse, ld, sizeof(GLfloat)*4);
		memcpy(lightSpecular, ls, sizeof(GLfloat)*4);
		memcpy(lightPosition, lp, sizeof(GLfloat)*4);
	}
} g;
glm::vec3 cameraPos = glm::vec3(g.camx, g.camy, g.camz);
glm::vec3 cameraFront = glm::vec3(0.0, 0.0, -1.0);
glm::vec3 cameraUp = glm::vec3(0.0, 1.0, 0.0);

class X11_wrapper {
private:
        Display *dpy;
        Window win;
        GLXContext glc;
        int xres, yres;
public:
        X11_wrapper() {
                xres = g.xres;
                yres = g.yres;
                GLint att[] = { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, None };
                //GLint att[] = { GLX_RGBA, GLX_DEPTH_SIZE, 24, None };
                XSetWindowAttributes swa;
                setup_screen_res(xres, yres);
                dpy = XOpenDisplay(NULL);
                if (!dpy) {
                        printf("cannot connect to X server\n\n");
                        exit(EXIT_FAILURE);
                }
                Window root = DefaultRootWindow(dpy);
                XVisualInfo *vi = glXChooseVisual(dpy, 0, att);
                if (!vi) {
                        printf("\n\tno appropriate visual found\n\n");
                        exit(EXIT_FAILURE);
                }
                Colormap cmap = XCreateColormap(dpy, root, vi->visual, AllocNone);
                swa.colormap = cmap;
                swa.event_mask =
                        ExposureMask |
                        StructureNotifyMask |
                        SubstructureNotifyMask |
                        KeyPressMask |
                        KeyReleaseMask |
                        PointerMotionMask |
                        ButtonPressMask |
                        ButtonReleaseMask;
                win = XCreateWindow(dpy, root, 0, 0, xres, yres, 0,
                        vi->depth, InputOutput, vi->visual, CWColormap | CWEventMask, &swa);
                set_title();
                glc = glXCreateContext(dpy, vi, NULL, GL_TRUE);
                glXMakeCurrent(dpy, win, glc);
        }
        ~X11_wrapper() { XDestroyWindow(dpy, win); XCloseDisplay(dpy); }
        int getxres() { return xres; }
        int getyres() { return yres; }
        void set_title(void) {
                //Set the window title bar.
                XMapWindow(dpy, win);
                char ts[256];
                sprintf(ts, "OpenGL Shader");
                XStoreName(dpy, win, ts);
        }
        void setup_screen_res(const int w, const int h) { xres = w; yres = h; }
	void reshape_window(int width, int height) {
                //window has been resized.
                setup_screen_res(width, height);
                glViewport(0, 0, (GLint)width, (GLint)height);
        }
        void check_resize(XEvent *e) {
                //The ConfigureNotify is sent by the
                //server if the window is resized.
                if (e->type != ConfigureNotify)
                        return;
                XConfigureEvent xce = e->xconfigure;
                if (xce.width != xres || xce.height != yres) {
                        //Window size did change.
                        reshape_window(xce.width, xce.height);
                        xres = xce.width;
                        yres = xce.height;
                }
        }
        bool getXPending() { return XPending(dpy); }
        XEvent getXNextEvent() {
                XEvent e;
                XNextEvent(dpy, &e);
                return e;
        }
        void swapBuffers() { glXSwapBuffers(dpy, win); }
} x11;

class Shader {
	public:
	GLuint ID; //program id
	GLuint vbo, vao, elem;
	char* vertexshader;
	char* fragmentshader;

	Shader(const char* vertexPath, const char* fragmentPath) {
		std::string vscode, fscode;
		std::ifstream vsFile;
		std::ifstream fsFile;
		//ensure file exceptions
		vsFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		fsFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

		try {
			vsFile.open(vertexPath);
			fsFile.open(fragmentPath);
			std::stringstream vsStream, fsStream;
			//read file contents into streams
			vsStream << vsFile.rdbuf();
			fsStream << fsFile.rdbuf();
			//close files
			vsFile.close();
			fsFile.close();
			//convert stream to string
			vscode = vsStream.str();
			fscode = fsStream.str();
		} catch(std::ifstream::failure e) {
			printf("ERROR FILE READ\n"); fflush(stdout);
		}
		printf("Files read successfully\n"); fflush(stdout);
		//printf("%s\n", vscode.c_str());
		//printf("%s\n", fscode.c_str());
		const char* vShaderCode = vscode.c_str();
		const char* fShaderCode = fscode.c_str();

		GLuint vertex, fragment;
		GLint compiled;
		//NEED TO CREATE GL CONTEXT AKA. GLEWINIT()********
		glewInit();
		////////////////////////////////////////////////

		vertex = glCreateShader(GL_VERTEX_SHADER);
		if(vertex == 0) {
			printf("ERROR: Vertex shader not created\n"); fflush(stdout);
		}
		glShaderSource(vertex, 1, &vShaderCode, NULL);
		glCompileShader(vertex);

		//Checking for compilation errors
		glGetShaderiv(vertex, GL_COMPILE_STATUS, &compiled);
		if(!compiled) {
			GLint infoLen = 0;
			glGetShaderiv(vertex, GL_INFO_LOG_LENGTH, &infoLen);
			if(infoLen > 1) {
				char* infoLog = (char*)malloc(sizeof(char) * infoLen);
				glGetShaderInfoLog(vertex, infoLen, NULL, infoLog);
				printf("ERROR COMPILING VERTEX SHADER: %s\n", infoLog); fflush(stdout);
				free(infoLog);
			}
			glDeleteShader(vertex);
		} else {
		printf("Vertex shader compiled successfully\n"); fflush(stdout);
		}
		fragment = glCreateShader(GL_FRAGMENT_SHADER);
		if(fragment == 0) {
			printf("ERROR: Fragment shader not created\n"); fflush(stdout);
		}
		glShaderSource(fragment, 1, &fShaderCode, NULL);
		glCompileShader(fragment);

		//Checking for compilation errors
		glGetShaderiv(fragment, GL_COMPILE_STATUS, &compiled);
		if(!compiled) {
			GLint infoLen = 0;
			glGetShaderiv(fragment, GL_INFO_LOG_LENGTH, &infoLen);
			if(infoLen > 1) {
				char* infoLog = (char*)malloc(sizeof(char) * infoLen);
				glGetShaderInfoLog(fragment, infoLen, NULL, infoLog);
				printf("ERROR COMPILING FRAGMENT SHADER: %s\n", infoLog); fflush(stdout);
				free(infoLog);
			}
			glDeleteShader(fragment);
		}else {
		printf("Fragment shader compiled successfully\n"); fflush(stdout);
		}
		GLint linked;
		//creating program
		ID = glCreateProgram();
		if(ID == 0) {
			printf("ERROR CREATING PROGRAM\n"); fflush(stdout);
		}

		glAttachShader(ID, vertex);
		glAttachShader(ID, fragment);

		glBindAttribLocation(ID, 0, "vPosition");
		glLinkProgram(ID);

		//Checking for linking errors
		glGetProgramiv(ID, GL_LINK_STATUS, &linked);
		if(!linked) {
			GLint infoLen = 0;
			glGetProgramiv(ID, GL_INFO_LOG_LENGTH, &infoLen);
			if(infoLen > 1) {
				char* infoLog = (char*)malloc(sizeof(char) * infoLen);
				glGetProgramInfoLog(ID, infoLen, NULL, infoLog);
				printf("ERROR LINKING PROGRAM: %s\n", infoLog); fflush(stdout);
				free(infoLog);
			}
			glDeleteProgram(ID);
		} else {
		printf("Program linked successfully\n"); fflush(stdout);
		}
	}

        void changeUniformValue1f(const char* name, float val1) {
                GLuint index = glGetUniformLocation(ID, name);
                if(index != -1) {
                        glUniform1f(index, val1);
                }
                else {
                        printf("ERROR: Changing uniform value %s\n", name); fflush(stdout);
                }
        }
        void changeUniformValue2f(const char* name, float val1, float val2) {
                GLuint index = glGetUniformLocation(ID, name);
                if(index != -1) {
                        glUniform2f(index, val1, val2);
                }
                else {
                        printf("ERROR: Changing uniform value %s\n", name); fflush(stdout);
                }
        }
        void changeUniformValue3f(const char* name, float val1, float val2, float val3) {
                GLuint index = glGetUniformLocation(ID, name);
                if(index != -1) {
                        glUniform3f(index, val1, val2, val3);
                }
                else {
                        printf("ERROR: Changing uniform value %s\n", name); fflush(stdout);
                }
        }
        void changeUniformValue4f(const char* name, float val1, float val2, float val3, float val4) {
                GLuint index = glGetUniformLocation(ID, name);
                if(index != -1) {
                        glUniform4f(index, val1, val2, val3, val4);
                }
                else {
                        printf("ERROR: Changing uniform value %s\n", name); fflush(stdout);
                }
        }
        void changeUniformMat4(const char* name, glm::mat4 mat) {
                GLuint index = glGetUniformLocation(ID, name);
                if(index != -1) {
                        glUniformMatrix4fv(index, 1, GL_FALSE, &mat[0][0]);
                }
                else {
                        printf("ERROR: Changing uniform value %s\n", name); fflush(stdout);
                }
        }
	
	void use() { //glUseProgram
		glUseProgram(ID);
	}
} bubble("./bubble.vs", "./bubble.fs"), 
	skybox("./skybox.vs", "./skybox.fs");

void initSkybox() {
	glGenBuffers(1, &skybox.vbo);
	glBindBuffer(GL_ARRAY_BUFFER, skybox.vbo);
	glBufferData(GL_ARRAY_BUFFER, 3*36*sizeof(float), &skyboxVerts, GL_STATIC_DRAW);

	glGenVertexArrays(1, &skybox.vao);
	glBindVertexArray(skybox.vao);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, skybox.vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
}

void initBubble() {
	glGenBuffers(1, &bubble.vbo);
	glGenBuffers(1, &bubble.elem);

	glBindBuffer(GL_ARRAY_BUFFER, bubble.vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bubble.elem);

	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), &vertices, GL_DYNAMIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(ind), ind, GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
}

void initOpengl() {
        glClearColor(0, 0, 0, 1.0f);
        glClearDepth(1.0);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glEnable(GL_BLEND);
        //glEnable(GL_TEXTURE_2D);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        //Enable this so material colors are the same as vert colors.
	glEnable(GL_COLOR_MATERIAL);
        glEnable( GL_LIGHTING );
        glLightfv(GL_LIGHT0, GL_AMBIENT, g.lightAmbient);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, g.lightDiffuse);
        glLightfv(GL_LIGHT0, GL_SPECULAR, g.lightSpecular);
        glLightfv(GL_LIGHT0, GL_POSITION, g.lightPosition);
        glMaterialfv(GL_FRONT, GL_SPECULAR, WHITE);
        glMaterialf(GL_FRONT, GL_SHININESS, 30);
        glEnable(GL_LIGHT0);	
}

unsigned int cubemapTexture;
void check_mouse(XEvent*);
int check_keys(XEvent*);
void render();
void physics();
void init_vbo();

int main() {

	//INITIALIZE VBO's
	initBubble();
	initSkybox();
	cubemapTexture = loadCubemap(faces);
	initOpengl();
	int done = 0;
	while(!done) {
		while(x11.getXPending()) {
			XEvent e = x11.getXNextEvent();
			x11.check_resize(&e);
			check_mouse(&e);
			done = check_keys(&e);
		}
		physics();
		render();
		x11.swapBuffers();
	}

	return(0);
}

void check_mouse(XEvent *e)
{
        //no mouse in this program
        if (e->type != ButtonPressMask && e->type != ButtonReleaseMask &&
                e->type != PointerMotionMask)
                return;
}

void screenCapture() {
        static int inc = 0;
        int xres = g.xres;
        int yres = g.yres;
        //get pixels
        unsigned char *data = new unsigned char [xres * yres * 3];
        glReadPixels(0, 0, xres, yres, GL_RGB, GL_UNSIGNED_BYTE, data);
        //write ppm file...
        char ts[256];
        sprintf(ts, "img%03i.ppm", inc++);
        FILE *fpo = fopen(ts, "w");
        fprintf(fpo, "P6\n");
        fprintf(fpo, "%i %i\n", xres, yres);
        fprintf(fpo, "255\n");
        //go backwards a row at a time...
        unsigned char *p = data;
        p = p + ((yres-1) * xres * 3);
        unsigned char *start = p;
        for (int i=0; i<yres; i++) {
                for (int j=0; j<xres*3; j++) {
                        fprintf(fpo, "%c", *p);
                        ++p;
                }
                start = start - (xres*3);
                p = start;
        }
        fclose(fpo);
}

int check_keys(XEvent *e)
{
        //Was there input from the keyboard?
        if (e->type != KeyPress && e->type != KeyRelease)
                return 0;
        int key = XLookupKeysym(&e->xkey, 0);
        if (e->type == KeyPress) {
                switch (key) {
			case XK_Left:
				g.xangle -= 1.0;
				break;
			case XK_Right:
				g.xangle += 1.0;
				break;
			case XK_Up:
				g.camz += 0.1;
				cameraPos = glm::vec3(g.camx, g.camy, g.camz);
				break;
			case XK_Down:
				g.camz -= 0.1;
				cameraPos = glm::vec3(g.camx, g.camy, g.camz);
				break;
			case XK_r:
				g.camx = g.camy = 0.0;
				g.camz = -2.0;
				cameraPos = glm::vec3(g.camx, g.camy, g.camz);
				break;
			 case XK_s:
                                for(int i=0; i<120; i++) {
                                        screenCapture();
                                        int done = 0;
                                        for(int j=0; j<5; j++) {
                                                while(x11.getXPending()) {
                                                        XEvent e = x11.getXNextEvent();
                                                        x11.check_resize(&e);
                                                        check_mouse(&e);
                                                        done = check_keys(&e);
                                                }
                                                physics();
                                                render();
                                                x11.swapBuffers();
                                        }
                                }
                                if(system("convert -loop 0 img0*.ppm proj4.gif")){}
				printf("Gif created\n");
                                break;
			case XK_Escape:
				return(1);
				break;
		}
	}

	return(0);
}

void physics() {      
	static bool less = true;
	b.updateCenter();
	double avgDist = 0.0;
	for(int i=0; i<b.nmasses; i++) {
		avgDist += getLength(b.mass[i].pos, b.center);
	}
	avgDist /= b.nmasses;
	if(avgDist > 0.7 && less) {
		for(int j=0; j<b.nsprings; j++)
			b.spring[j].stiffness *= -1;//(rnd() * 0.000001) + 0.000015;
		less ^= 1;
	}
	else if(avgDist < 0.65 && !less) {
		for(int j=0; j<b.nsprings; j++)
			b.spring[j].stiffness *= -1;//(rnd() * 0.000001) - 0.000015;
		less ^= 1;
	}
	if(b.pull < 0) {
		b.mass[0].vel[0] -= 1.5;
		++b.pull;
	}
	else if(b.pull > 0) {
		b.mass[b.nmasses-1].vel[0] += 1.5;
		--b.pull;
	}
	b.stray();
	b.maintainSprings();
}

void render() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(50, (double)g.xres/(double)g.yres, 1, 1000);

	
	//Skybox
	glDepthMask(GL_FALSE);
	skybox.use();

	//view
        glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
        skybox.changeUniformMat4("V", view);
        //Projection
        glm::mat4 projection = glm::perspective(1.0f, (float)g.xres/(float)g.yres, 0.1f, 100.0f);
        projection = glm::rotate(projection, glm::radians(g.xangle), glm::vec3(0.0, 1.0, 0.0));
        projection = glm::rotate(projection, glm::radians(g.yangle), glm::vec3(1.0, 0.0, 1.0));
        skybox.changeUniformMat4("P", projection);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
	glBindVertexArray(skybox.vao);
	glDrawArrays(GL_TRIANGLES, 0, 36);
	glDepthMask(GL_TRUE);
	//bubble	
	bubble.use();
	glm::mat4 model = glm::mat4(1.0f);
	bubble.changeUniformMat4("M", model);
	bubble.changeUniformMat4("V", view);
	glm::mat4 inverseModel = inverse(model);
	glm::mat4 TIM = transpose(inverseModel);
//	bubble.changeUniformMat4("TIM", TIM);
	bubble.changeUniformMat4("P", projection);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
	glBegin(GL_TRIANGLES);
        for(int i=0; i<b.longitude-1; i++) {
		for(int j=0; j<b.latitude; j++) {
			 glVertex3dv(b.mass[i*b.latitude + j].pos);
			 glVertex3dv(b.mass[(i+1)*b.latitude + j].pos);
			 glVertex3dv(b.mass[(i+1)*b.latitude + (j+1)%b.latitude].pos);

			 glVertex3dv(b.mass[i*b.latitude + j].pos);
			 glVertex3dv(b.mass[i*b.latitude + (j+1)%b.latitude].pos);
			 glVertex3dv(b.mass[(i+1)*b.latitude + (j+1)%b.latitude].pos);
		}
	}
        glEnd();

	glUseProgram(0);

}
