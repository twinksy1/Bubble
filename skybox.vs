#version 130
in vec3 pos;
out vec3 texcoords;
uniform mat4 P, V;
void main() {
	gl_Position = P * V * vec4(pos, 1.0);
	texcoords = pos;
}

