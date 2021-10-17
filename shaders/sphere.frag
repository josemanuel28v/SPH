#version 430

in vec3 vPositionEye;
in vec3 vNormalEye;
in vec3 vColor;

out vec4 fFragColor;

struct LightInfo {
	vec4 lightPos; // Posicion de la luz (S.R. de la vista)
	vec3 intensity;
	vec3 k;
};

uniform int uObject;

layout (location = 7) uniform vec3 uAmbient;
layout (location = 8) uniform vec3 uDiffuse;
layout (location = 9) uniform vec3 uSpecular;
layout (location = 10) uniform float uShininess;
layout (location = 11) uniform vec4 uLightPos;
layout (location = 12) uniform vec3 uLightIntensity;
layout (location = 13) uniform vec3 uK;
layout (location = 14) uniform float uParam;
layout (location = 15) uniform int uColorMode;

const float shininess = 16;


// ADS para la luz puntual combinada con iluminacion de hemisferio
// Parametros:	pos - posicion del vertice en el S.R. de la vista
//				n - normal del vortice en el S.R. de la vista
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

vec3 adsColor(in vec3 pos, in vec3 n, vec3 color)
{
	vec3 lightDir = vec3(uLightPos) - pos;
	float d = length(lightDir);
	lightDir = normalize(lightDir);
	vec3 viewDir = normalize( vec3(-pos) );
	vec3 reflectDir = reflect(lightDir, n);
	vec3 halfwayDir = normalize(lightDir + viewDir);

	vec3 intensity = uLightIntensity / (uK[0] + uK[1] * d + uK[2] * d * d);

	vec3 ambient = 0.2 * color;
	vec3 diff = color * max(dot(lightDir,n), 0.0);
	vec3 spec = vec3(0.3) * pow(max(dot(n,halfwayDir),0.0), shininess);			// Blinn phong
	//vec3 spec = 0.5 * color * pow(max(dot(reflectDir, viewDir),0.0), 16.0);	// Phong

	return intensity * (ambient + diff + spec);
}

vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);

    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void main()
{
	vec3 norm = normalize(vNormalEye);
	vec3 color;

	if (uObject == 0) // Plano
	{
		fFragColor = clamp(vec4(ads(vPositionEye, norm), 1.0), 0.0, 1.0);
	}
	else if (uObject == 1) // Esfera
	{
		// Color de las normales
		switch (uColorMode)
		{
			case 0: // Material
			case 5:
				fFragColor = clamp(vec4(ads(vPositionEye, norm), 1.0), 0.0, 1.0);
				break;

			case 1: // Normales
				fFragColor = clamp(vec4(adsColor(vPositionEye, norm, norm), 1.0), 0.0, 1.0); 
				break;

			case 2:	// Pressures
				color = vec3(1.0, 1.0 - uParam, 1.0 - uParam); // Pressures;
				fFragColor = clamp(vec4(adsColor(vPositionEye, norm, color), 1.0), 0.0, 1.0);
				break;

			case 3:	// Velocities
				color = vec3(uParam, uParam, 1.0);
				fFragColor = clamp(vec4(adsColor(vPositionEye, norm, color), 1.0), 0.0, 1.0);
				break;

			case 4: // Velocidad
				color = hsv2rgb(vec3(uParam, 1.0, 1.0)); // Velocities hsv
				fFragColor = clamp(vec4(adsColor(vPositionEye, norm, color), 1.0), 0.0, 1.0);
				break;							
		}

		// apply gamma correction
    	//float gamma = 1.2;
    	//fFragColor.rgb = pow(fFragColor.rgb, vec3(1.0 / gamma));	

		//fFragColor.w = 0.2;	
	}	
}