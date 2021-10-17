#version 430

in vec4 gColor;
in vec2 gTexCoord;
in vec4 gWPos;

out vec4 fFragColor;

layout (location = 4) uniform sampler2D uSpriteTex;
layout (location = 1) uniform mat4 uModelViewMatrix;
layout (location = 5) uniform float uParticleSize;
layout (location = 3) uniform mat4 uProjMatrix;
layout (location = 0) uniform mat4 uModelViewProjMatrix;

layout (location = 7) uniform vec3 uAmbient;
layout (location = 8) uniform vec3 uDiffuse;
layout (location = 9) uniform vec3 uSpecular;
layout (location = 10) uniform float uShininess;
layout (location = 11) uniform vec4 uLightPos;
layout (location = 12) uniform vec3 uLightIntensity;
layout (location = 13) uniform vec3 uK;
layout (location = 14) uniform float uParam;
layout (location = 15) uniform int uColorMode;

const float shininess = 25;
float rad = uParticleSize;
vec3 normal;

vec3 ads(in vec3 pos, in vec3 n)
{
	vec3 lightDir = vec3(uLightPos) - pos;
	float d = length(lightDir);
	lightDir = normalize(lightDir);
	vec3 viewDir = normalize( vec3(-pos) );
	vec3 reflectDir = reflect(-lightDir, n);
	vec3 halfwayDir = normalize(lightDir + viewDir);

	vec3 intensity = uLightIntensity / (uK[0] + uK[1] * d + uK[2] * d * d);

	vec3 ambient = uAmbient;
	vec3 diff = uDiffuse * max(dot(lightDir,n), 0.0);
	vec3 spec = uSpecular * pow(max(dot(n,halfwayDir),0.0), uShininess);			// Blinn phong
	//vec3 spec = uMaterial.specular * pow(max(dot(reflectDir, viewDir),0.0), uMaterial.shininess);	// Phong

	return intensity * (ambient + diff + spec);
}

vec3 adsColor(in vec3 pos, in vec3 n, in vec3 color)
{
	vec3 lightDir = vec3(uLightPos) - pos;
	float d = length(lightDir);
	lightDir = normalize(lightDir);
	vec3 viewDir = normalize( vec3(-pos) );
	vec3 reflectDir = reflect(-lightDir, n);
	vec3 halfwayDir = normalize(lightDir + viewDir);

	vec3 intensity = uLightIntensity / (uK[0] + uK[1] * d + uK[2] * d * d);

	vec3 ambient = 0.2 * color;
	vec3 diff = color * max(dot(lightDir,n), 0.0);
	vec3 spec = vec3(0.3) * pow(max(dot(n,halfwayDir),0.0), shininess);			// Blinn phong
	//vec3 spec = uMaterial.specular * pow(max(dot(reflectDir, viewDir),0.0), uMaterial.shininess);	// Phong

	return intensity * (ambient + diff + spec);
}

void makeSphere()
{
    //clamps fragments to circle shape. 
    vec2 mapping = gTexCoord * 2.0F - 1.0F;
    float d = dot(mapping, mapping);

    float z = sqrt(1.0F - d);
    normal = /*mat3(inverse(transpose(uModelViewMatrix))) * *//*normalize(*/vec3(mapping, z)/*)*/;
    vec3 cameraPos = vec3(uModelViewMatrix * gWPos) + rad * normal;

    ////Set the depth based on the new cameraPos.
    vec4 clipPos = uProjMatrix * vec4(cameraPos, 1.0);
    float ndcDepth = clipPos.z / clipPos.w;
    gl_FragDepth = ((gl_DepthRange.diff * ndcDepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;

    //calc ambient occlusion for circle
    //ambientOcclusion = sqrt(1.0F - d * 0.5F);
	normal = normalize(normal);
}

void Impostor(out vec3 cameraPos, out vec3 cameraNormal)
{
	vec2 mapping = gTexCoord * 2.0F - 1.0F;
    //float d = dot(mapping, mapping);

    vec3 cameraSpherePos = vec3(uModelViewMatrix * gWPos);

    vec3 cameraPlanePos = vec3(mapping * uParticleSize, 0.0) + cameraSpherePos;
    vec3 rayDirection = normalize(cameraPlanePos);

    float B = 2.0 * dot(rayDirection, -cameraSpherePos);
    float C = dot(cameraSpherePos, cameraSpherePos) - (uParticleSize * uParticleSize);

    float det = (B * B) - (4 * C);
    if(det < 0.0)
        discard;

    float sqrtDet = sqrt(det);
    float posT = (-B + sqrtDet)/2;
    float negT = (-B - sqrtDet)/2;

    float intersectT = min(posT, negT);
    cameraPos = rayDirection * intersectT;
    cameraNormal = normalize(cameraPos - cameraSpherePos);

	vec4 clipPos = uProjMatrix * vec4(cameraPos, 1.0);
    float ndcDepth = clipPos.z / clipPos.w;
    gl_FragDepth = ((gl_DepthRange.diff * ndcDepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;

    //cameraNormal = mat3(transpose(uModelViewMatrix)) * cameraNormal;
	//cameraNormal = mat3(inverse(transpose(uModelViewMatrix))) * cameraNormal;
	cameraNormal = normalize(cameraNormal);
}

vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);

    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void main()
{
	if (length(gTexCoord - 0.5) > 0.38) discard;

	//makeSphere();
	
	vec3 cameraPos;
	Impostor(cameraPos, normal);
	
	vec4 lightColor;
	vec3 color;

	switch(uColorMode)
	{
		case 0: // Material 
			lightColor = vec4(ads(vec3(uModelViewMatrix * gWPos), normal), 1.0);
			break;

		case 1: // Normales
			lightColor = vec4(adsColor(vec3(uModelViewMatrix * gWPos), normal, abs(normal)), 1.0);
			break;

		case 2: // Pressure
			color = vec3(0.8, 0.8 - uParam, 0.8 - uParam); // Pressures
			lightColor = vec4(adsColor(vec3(uModelViewMatrix * gWPos), normal, color), 1.0);
			break;

		case 3: // Velocidad
			vec3 blue = vec3(0.0, 0.58, 1.0);
			vec3 white = vec3(0.95, 0.95 ,0.95);
			color = mix(blue, white, uParam);
			//color = vec3(uParam, uParam, 1.0);
			//vec3 color = abs(normalize(uParamColor3D)); // Velocity orientation
			lightColor = vec4(adsColor(vec3(uModelViewMatrix * gWPos), normal, color), 1.0);
			break;

		case 4: // Velocidad // hsv(195, 100%, 100%)
			color = hsv2rgb(vec3(uParam, 0.6, 0.85)); // Velocities hsv
			lightColor = vec4(adsColor(vec3(uModelViewMatrix * gWPos), normal, color), 1.0);
			break;	

		case 5: // Sin iluminacion
			lightColor = vec4(1.0, 1.0, 1.0, 1.0);
			break;	
	}
	
	vec4 texColor = texture(uSpriteTex, gTexCoord);

	fFragColor = texColor * lightColor;

	//fFragColor.w = 0.2;
}
