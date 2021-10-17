#version 430

// Tarea por hacer: indicar el tipo de primitiva de entrada y el de salida.
layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

layout (location = 3) uniform mat4 uProjMatrix;
layout (location = 5) uniform float uParticleSize;

in vec4 vColor[];
in vec4 vWPos[];

out vec4 gWPos;
out vec2 gTexCoord;
out vec4 gColor;

void main()
{
	vec4 pos = gl_in[0].gl_Position;

	float uSize2 = uParticleSize;

	//
	//     v2    v3
	//     *-----*
	//     |\    |
	//     | \   |
	//     |  \  |
	//     |   \ |
	//     *-----*
	//    v0     v1
	//
	//  v0 = pos + vec4(-uSize2, -uSize2, 0.0, 0.0), etc.
	//  Para cada vertice de la primitiva de salida, fijar los datos de salida (gl_Position, coordenadas de textura y color) y emitir sus vertices
	
	gColor = vColor[0];
	gWPos = vWPos[0];

	gTexCoord = vec2(0, 0);
	gl_Position = uProjMatrix * (pos + vec4(-uSize2, -uSize2, 0.0, 0.0));
	EmitVertex();

	gTexCoord = vec2(1, 0);
	gl_Position = uProjMatrix * (pos + vec4(uSize2, -uSize2, 0.0, 0.0));
	EmitVertex();

	gTexCoord = vec2(0, 1);
	gl_Position = uProjMatrix * (pos + vec4(-uSize2, uSize2, 0.0, 0.0));
	EmitVertex();

	gTexCoord = vec2(1, 1);
	gl_Position = uProjMatrix * (pos + vec4(uSize2, uSize2, 0.0, 0.0));
	EmitVertex();
}
