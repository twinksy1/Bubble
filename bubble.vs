#version 130
uniform mat4 P, V, M;

varying vec3 Normal;
varying vec3 Position;
varying vec3 I, R;
uniform vec3 cameraPos;

void main() {
	Normal = gl_NormalMatrix*gl_Normal;
	Position = vec3(gl_ModelViewMatrix*gl_Vertex);
	gl_Position = P * V * M * vec4(Position, 1.0);
	I = normalize(Position - cameraPos);
}

