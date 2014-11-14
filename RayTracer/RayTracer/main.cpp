/*
 CSCI 480
 Assignment 3 Raytracer
 
 Name: <Your name here>
 */

#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include "pic.h"
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct color {
    float rgb [3] = {0.0, 0.0, 0.0};
};

struct point {
    double x;
    double y;
    double z;
    
    point (){
        
    }
    
    point (double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }
    
    point operator - (const point &p1){
        return point(x-p1.x, y - p1.y, z - p1.z);
    }
    
    point operator + (const point &p1)
    {
        return point(x+ p1.x, y + p1.y, z + p1.z);
    }
    
    point operator * (float scalar){
        return point (x * scalar, y * scalar, z * scalar);
    }
    
    static point cross_product (point v1, point v2) {
        point v3;
        v3.x = v1.y*v2.z - v1.z*v2.y;
        v3.y = v2.x*v1.z - v1.x*v2.z;
        v3.z = v1.x*v2.y - v1.y*v2.x;
        return v3;
    }
    
    static float dot_product (point v1, point v2){
        return v1.x * v2.x +  v1.y * v2.y +  v1.z * v2.z;
    }
    
    void normalize (){
        double vectorLength = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
        this->x = this->x/vectorLength;
        this->y = this->y/vectorLength;
        this->z = this->z/vectorLength;
    }
    
    static point determine_vector (point p1, point p2){
        point v;
        v.x = p2.x-p1.x;
        v.y = p2.y-p1.y;
        v.z = p2.z-p1.z;
        v.normalize();
        return v;
    }
    
    void print (){
        cout << "X = " << x << " Y = " << y << " Z = " << z << endl;
        cout.flush();
    }
    
};

typedef point vector;

struct Vertex
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double normal[3];
    double shininess;
};

typedef struct _Triangle
{
    struct Vertex v[3];
    
    float compute_area (point A, point B, point C){
        vector cross = point::cross_product(B-A, C-A);
        return 0.5 * sqrt(point::dot_product(cross, cross));
    }
} Triangle;

typedef struct _Sphere
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double shininess;
    double radius;
} Sphere;

typedef struct _Light
{
    double position[3];
    double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

struct Ray
{
    point origin;
    vector direction;
    
    Ray (point p, vector v){
        origin = p; direction = v;
    }
    
    float check_sphere_intersection (Sphere sphere1){
        //(x0 + xd t - xc)2 + (y0 + yd t - yc)2 + (z0 + zd t - zc)2 - r2 = 0
        //Simplifies to at2 + bt + c = 0
        //a = xd2 + yd2+ zd2 = 1
        //b = 2(xd(x0-xc)+yd(y0-yc)+zd(z0-zc))
        //c = (x0-xc)2 + (y0-yc)2 +(z0-zc)2 - r2
        
        float b, c;
        b = 2*((direction.x * (origin.x - sphere1.position[0])) +
               (direction.y * (origin.y - sphere1.position[1])) +
               (direction.z * (origin.z - sphere1.position[2])));
        c = pow(origin.x-sphere1.position[0], 2) + pow(origin.y-sphere1.position[1], 2) + pow(origin.z-sphere1.position[2], 2) - pow(sphere1.radius, 2);
        
        //Solve to get t0, t1
        //t0,1 = (-b +/- sqrt(b2-4c))/2
        //If t0, t1 > 0, return min (t0, t1)
        
        float t0, t1;
        t0 = (-b + sqrt(pow(b, 2) - 4*c))/2;
        t1 = (-b - sqrt(pow(b, 2) - 4*c))/2;
        if (t0 < 0 && t1 < 0)
            return -1;
        else
            return fmax(t0, t1);
    }
    
    bool check_triangle_intersection (Triangle triangle1, color &illumination){
        //Check if ray intersects triangle plane
        
        point A(triangle1.v[0].position[0], triangle1.v[0].position[1], triangle1.v[0].position[2]);
        point B(triangle1.v[1].position[0], triangle1.v[1].position[1], triangle1.v[1].position[2]);
        point C(triangle1.v[2].position[0], triangle1.v[2].position[1], triangle1.v[2].position[2]);
        
        //Normal
        vector cross = point::cross_product(B-A, C-A);
        vector normal = cross;
        normal.normalize();
        
        float d = vector::dot_product(normal, A);
        
        /*
         Plane: 0 = p . N + d
         
         t = - (P0 .  N + d) / (V . n)
         
         t = n . (v0 - p0) / n . (p1-p0)
         
         if n . d = 0, no intersection
         if t <= 0 intersection is behind ray origin
         */
        
        if (vector::dot_product(normal, direction) == 0)
            return false;
        
        float t = vector::dot_product(normal, A-origin) / vector::dot_product(normal, direction);
        
        if (t < 0)
            return false;
        
        vector intersection = origin + (direction * t);
        
        
        /*
         Now we have 3 points instead of 2
         • Define 3 barycentric coordinates α, β, γ
         • p = α p1 + β p2 + γ p3
         • p inside triangle iff 0 ≤ α, β, γ ≤ 1, α + β + γ = 1
         • How do we calculate α, β, γ?
         
         Coordinates are ratios of triangle areas
         • Areas in these formulas should be signed
         - Clockwise (-) or anti-clockwise (+) orientation of the triangle
         - Important for point-in-triangle test
         */
        
        float alpha,beta,gamma;
        
        
        alpha = triangle1.compute_area(intersection, B, C) / triangle1.compute_area(A, B, C);
        beta = triangle1.compute_area(A, intersection, C) / triangle1.compute_area(A, B, C);
        gamma = triangle1.compute_area(A, B, intersection) / triangle1.compute_area(A, B, C);
        if (0 <= alpha && alpha <= 1 && 0 <= beta && beta <= 1 && 0 <= gamma && gamma <= 1){
            illumination = light_calculation(triangle1, intersection, alpha, beta, gamma);
            return true;
        }
        else
            return false;
    }
    
    color light_calculation (Triangle &t1, point intersection, float alpha, float beta, float gamma){
        
        color illumination;
        
        float shininess = t1.v[0].shininess*alpha
        + t1.v[1].shininess*beta
        + t1.v[2].shininess*gamma;
        
        //interpolated normal
        vector normal(t1.v[0].normal[0]*alpha + t1.v[1].normal[0]*beta + t1.v[2].normal[0]*gamma,
                      t1.v[0].normal[1]*alpha + t1.v[1].normal[1]*beta + t1.v[2].normal[1]*gamma,
                      t1.v[0].normal[2]*alpha + t1.v[1].normal[2]*beta + t1.v[2].normal[2]*gamma);
        normal.normalize();
        
        //vector viewer (-direction.x, -direction.y, -direction.z);
        vector viewer = direction;
        
        for (int rgb = 0; rgb < 3; rgb++){
            
            float diffuse = t1.v[0].color_diffuse[rgb]*alpha
            + t1.v[1].color_diffuse[rgb]*beta
            + t1.v[2].color_diffuse[rgb]*gamma;
            
            float specular = t1.v[0].color_specular[rgb]*alpha
            + t1.v[1].color_specular[rgb]*beta
            + t1.v[2].color_specular[rgb]*gamma;
            
            //illumination value
            float light = 0;
            
            //for each light source
            for (int i = 0; i < num_lights; i++){
                //unit vector to light
                vector l(lights[i].position[0] - intersection.x,
                         lights[i].position[1] - intersection.y,
                         lights[i].position[2] - intersection.z);

                l.normalize();
                
                vector reflected = (normal * 2 * vector::dot_product(l, normal)) - l;
                reflected.normalize();
                
                // I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ sh)
                
                float dif_scalar = vector::dot_product(l, normal);
                    if (dif_scalar < 0)
                        dif_scalar = 0;
                float spec_scalar = vector::dot_product(reflected, viewer);
                    if (spec_scalar < 0)
                        spec_scalar = 0;
                
                light  += lights[i].color[rgb] * ((diffuse * dif_scalar) + specular * (pow(spec_scalar, shininess)));
            }
            
            //add global ambient light
            light += ambient_light[rgb];
            
            if (light > 1)
               light = 1;
                
            illumination.rgb[rgb] = light * 255;
        }
        return illumination;
    }
};

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//MODIFY THIS FUNCTION
void draw_scene()
{
    unsigned int x,y;
    //simple output
    for(x=0; x<WIDTH; x++)
    {
        glPointSize(2.0);
        glBegin(GL_POINTS);
        for(y=0;y < HEIGHT;y++)
        {
            plot_pixel(x,y,x%256,y%256,(x+y)%256);
        }
        glEnd();
        glFlush();
    }
    printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
    glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
    glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
    buffer[HEIGHT-y-1][x][0]=r;
    buffer[HEIGHT-y-1][x][1]=g;
    buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
    point camera (0, 0, 0);
    
    //Transform screen to world space
    //Aspect ratio
    float a = WIDTH/HEIGHT;
    float scaleX = a * tan(fov/2);
    float xVal = (float) x/WIDTH;
    float scaleY = tan(fov/2);
    float yVal = (float) y/HEIGHT;
    point screen(scaleX - (2 * xVal * scaleX), scaleY - (2 * yVal * scaleY), -2);
    vector direction = vector::determine_vector(camera, screen);
    Ray r1(screen, direction);
    color color_val;
    
    bool intersects = false;
    //For every object in the scene
    //Check spheres
    for (int i = 0; i < num_spheres; i++){
        //If they intersect
        if (r1.check_sphere_intersection(spheres[i]) >= 0){
            //Plot Object's Pixel
            intersects = true;
        }
    }
    //Check Triangles
    for (int i = 0; i < num_triangles; i++){
        //If they intersect
        if (r1.check_triangle_intersection(triangles[i], color_val)){
            //Plot Object's Pixel
            intersects = true;
        }
    }
    if (intersects) {
        plot_pixel_display(x,y,color_val.rgb[0],color_val.rgb[1],color_val.rgb[2]);
    }
    else
        plot_pixel_display(x,y,0.0,0.0,0.0);
    
    
    // if(mode == MODE_JPEG)
    //   plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
    Pic *in = NULL;
    
    in = pic_alloc(640, 480, 3, NULL);
    printf("Saving JPEG file: %s\n", filename);
    
    memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
    if (jpeg_write(filename, in))
        printf("File saved Successfully\n");
    else
        printf("Error in Saving\n");
    
    pic_free(in);
    
}

void parse_check(char *expected,char *found)
{
    if(strcasecmp(expected,found))
    {
        char error[100];
        printf("Expected '%s ' found '%s '\n",expected,found);
        printf("Parse error, abnormal abortion\n");
        exit(0);
    }
    
}

void parse_doubles(FILE*file, char *check, double p[3])
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check(check,str);
    fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
    printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check("rad:",str);
    fscanf(file,"%lf",r);
    printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
    char s[100];
    fscanf(file,"%s",s);
    parse_check("shi:",s);
    fscanf(file,"%lf",shi);
    printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
    argv = "/Users/joshgreenberger/Documents/Class Material/Senior/Graphics/assign3/RayTracer/screenfile.txt";
    FILE *file = fopen(argv,"r");
    if (file == NULL) {
        printf ("Can't open file.\n");
        exit(1);
    }
    int number_of_objects;
    char type[50];
    int i;
    Triangle t;
    Sphere s;
    Light l;
    fscanf(file,"%i",&number_of_objects);
    
    printf("number of objects: %i\n",number_of_objects);
    char str[200];
    
    parse_doubles(file,"amb:",ambient_light);
    
    for(i=0;i < number_of_objects;i++)
    {
        fscanf(file,"%s\n",type);
        printf("%s\n",type);
        if(strcasecmp(type,"triangle")==0)
        {
            
            printf("found triangle\n");
            int j;
            
            for(j=0;j < 3;j++)
            {
                parse_doubles(file,"pos:",t.v[j].position);
                parse_doubles(file,"nor:",t.v[j].normal);
                parse_doubles(file,"dif:",t.v[j].color_diffuse);
                parse_doubles(file,"spe:",t.v[j].color_specular);
                parse_shi(file,&t.v[j].shininess);
            }
            
            if(num_triangles == MAX_TRIANGLES)
            {
                printf("too many triangles, you should increase MAX_TRIANGLES!\n");
                exit(0);
            }
            triangles[num_triangles++] = t;
        }
        else if(strcasecmp(type,"sphere")==0)
        {
            printf("found sphere\n");
            
            parse_doubles(file,"pos:",s.position);
            parse_rad(file,&s.radius);
            parse_doubles(file,"dif:",s.color_diffuse);
            parse_doubles(file,"spe:",s.color_specular);
            parse_shi(file,&s.shininess);
            
            if(num_spheres == MAX_SPHERES)
            {
                printf("too many spheres, you should increase MAX_SPHERES!\n");
                exit(0);
            }
            spheres[num_spheres++] = s;
        }
        else if(strcasecmp(type,"light")==0)
        {
            printf("found light\n");
            parse_doubles(file,"pos:",l.position);
            parse_doubles(file,"col:",l.color);
            
            if(num_lights == MAX_LIGHTS)
            {
                printf("too many lights, you should increase MAX_LIGHTS!\n");
                exit(0);
            }
            lights[num_lights++] = l;
        }
        else
        {
            printf("unknown type in scene description:\n%s\n",type);
            exit(0);
        }
    }
    return 0;
}

void display()
{
    
}

void init()
{
    glMatrixMode(GL_PROJECTION);
    glOrtho(0,WIDTH,0,HEIGHT,1,-1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
    //hack to make it only draw once
    static int once=0;
    if(!once)
    {
        draw_scene();
        if(mode == MODE_JPEG)
            save_jpg();
    }
    once=1;
}

int main (int argc, char ** argv)
{
    if (argc<2 || argc > 3)
    {
        printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
        exit(0);
    }
    if(argc == 3)
    {
        mode = MODE_JPEG;
        filename = argv[2];
    }
    else if(argc == 2)
        mode = MODE_DISPLAY;
    
    glutInit(&argc,argv);
    loadScene(argv[2]);
    
    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(WIDTH,HEIGHT);
    int window = glutCreateWindow("Ray Tracer");
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    init();
    glutMainLoop();
}
