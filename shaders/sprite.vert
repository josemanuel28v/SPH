#version 430  

layout (location = 1) uniform mat4 uModelViewMatrix;
layout (location = 6) uniform vec4 uPosition;

in vec4 aPosition;
in vec4 aColor;

out vec4 vColor;
out vec4 vWPos;

void main()
{
	vColor = aColor;
	vWPos = uPosition;
	
	gl_Position = uModelViewMatrix * uPosition; 
}