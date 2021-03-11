#version 130
in vec3 texcoords;
uniform samplerCube skybox;

void main(){
	gl_FragColor = texture(skybox, texcoords);
}
