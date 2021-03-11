#version 130
varying vec3 Normal;
varying vec3 Position;
varying vec3 I;
uniform vec3 cameraPos;
uniform samplerCube skybox;

const float R = 1.77;
const float G = 1.84;
const float B = 1.54;
const float fresnelPower = 2.0;
const float F = ((1.0 - G) * (1.0 - G)) / ((1.0 + G) * (1.0 + G));

void main() {
	vec3 i = normalize(I);
	vec3 n = normalize(Normal);

	float ratio = F + (1.0 - F) * pow(1.0 - dot(-i, n), fresnelPower);

	vec3 RR = vec3(gl_TextureMatrix[0] * vec4(refract(i, n, R), 1.0));
	vec3 GR = vec3(gl_TextureMatrix[0] * vec4(refract(i, n, G), 1.0));
	vec3 BR = vec3(gl_TextureMatrix[0] * vec4(refract(i, n, B), 1.0));

	vec3 reflectDir = vec3(gl_TextureMatrix[0] * vec4(reflect(i, n), 1.0));

	vec4 refractColor;
	refractColor.r = textureCube(skybox, RR).r;
	refractColor.g = textureCube(skybox, GR).g;
	refractColor.b = textureCube(skybox, BR).b;

	vec4 reflectColor = textureCube(skybox, reflectDir);
	gl_FragColor = mix(refractColor, reflectColor, ratio);
}
