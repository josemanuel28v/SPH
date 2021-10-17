#version 430

in vec3 aPosition; // Posici�n del v�rtice en el S.R. local al objeto
in vec3 aNormal; // Normal del v�rtice en el S.R. local al objeto
in vec2 aTexCoord;

out vec3 vPositionEye;
out vec3 vNormalEye;
out vec3 vColor;

layout (location = 0) uniform mat4 uModelViewProjMatrix;
layout (location = 1) uniform mat4 uModelViewMatrix;
layout (location = 2) uniform mat3 uNormalMatrix;

void main()
{
	vPositionEye = vec3(uModelViewMatrix * vec4(aPosition, 1.0)); // Posici�n en el S.R. vista
	vNormalEye = normalize(uNormalMatrix * aNormal); // Normal en el S.R. vista

	gl_Position = uModelViewProjMatrix * vec4(aPosition, 1.0);
}