
/*
#version 330

const float dx= 0;
const float dy= 0;
const float dz= 0;

void main( )
{
    // intialiser les coordonnees des 3 sommets
    vec3 positions[3]= vec3[3]( vec3(-0.5, -0.5, 0), vec3(0.5, -0.5, 0), vec3(0, 0.5, 0) );
    
    // recuperer le sommet a traiter
    vec3 p= positions[gl_VertexID];

    // calculer le resultat 
    vec4 r;
    r.x= p.x + dx;
    r.y= p.y + dy;
    r.z= p.z + dz;
    r.w= 1;
    
    // renvoyer le sommet transforme
    gl_Position= r;
}
*/

#version 330

layout(location = 0) in vec3 position;

uniform mat4 mvp;

void main( )
{
    //vec3 positions[3]= vec3[3]( vec3(-0.5, -0.5, 0), vec3(0.5, -0.5, 0), vec3(0, 0.5, 0) );
    // recuperer le sommet a traiter
    //vec3 p= positions[gl_VertexID];
    vec4 p = mvp * vec4(position, 1);
    vec4 q = vec4(10.5, 0.0, 0.0, 0.0);
    gl_Position = p + q;
}
