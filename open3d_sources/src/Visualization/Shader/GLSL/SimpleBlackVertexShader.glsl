#version 120

attribute vec3 vertex_position;
uniform mat4 MVP;

void main()
{
    gl_Position = MVP * vec4(vertex_position, 1);
}
